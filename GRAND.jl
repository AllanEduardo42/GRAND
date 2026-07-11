################################################################################
# Allan Eduardo Feitosa
# 15 Jun 2026
# Setup for GRAND algorithm decoding simulation

using Random
using LinearAlgebra
using Plots
using LoopVectorization
using SpecialFunctions
using Polynomials
using SparseArrays
using BenchmarkTools
using DelimitedFiles

const MAXC2V = 1e3                      # saturate values for C2V (Inf approx)
const MINC2V = -MAXC2V                  # -Inf approx

function qfunc(X::Float64)
    return 0.5 * erfc(X/sqrt(2))
end

include("Koopman.jl")
include("../LDPC/GF2_functions.jl")
include("make_code.jl")
include("GRAND_sim.jl")
include("ORBGRAND_sim.jl")
include("calc_syndrome.jl")
include("../LDPC/PEG.jl")
include("../LDPC/LU_encoding.jl")
include("../LDPC/auxiliary setup functions.jl")
include("/home/allan/LDPC/Simulation core functions/encode_LDPC.jl")
include("/home/allan/LDPC/simcore.jl")
include("/home/allan/LDPC/Algorithms/Flooding.jl")

SEED::Int = 0001

PAYLOAD::Int = 234
CODE_LEN::Int = 256

REDUN::Int = CODE_LEN - PAYLOAD

TEST::Bool = false
PRINT::Bool = false
MAX_ERRORS_TEST::Int = 1

PROTOCOL::String = "CRC"

# GRAND

MAX_ERRORS::Int = 5

ABANDON::Bool = true
MAX_QUERY::Int = 1_000_000_000

MAX_ERR_LOC_VEC_LEN::Int = CODE_LEN

# EbN0 = [1.0, 1.5, 2.0, 2.5, 3.0]
EbN0 = [3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]

FULL::Bool = true
MEAN::Bool = false
MAX_DEPTH::Int = 100_000
NUM_SEGMENTS::Int = 2

### Parity-check matrix

if PROTOCOL == "PEG"
    #Generate Parity-Check Matrix by the PEG algorithm

    LAMBDA = [0.21, 0.25, 0.25, 0.29, 0]
    RO = [1.0, 0, 0, 0, 0, 0]
    H_PEG, GIRTH = PEG(LAMBDA,RO,REDUN,CODE_LEN)

    HH, LL, UU, FF = LU_encoding(H_PEG,0)

    PP = zeros(Bool,REDUN,PAYLOAD)

    for k in axes(PP,2)
        PP[:,k] = gf2_solve_LU(LL,UU,HH[:,k])
    end

elseif PROTOCOL == "CRC"

    PP, HH, CRC_POLY, HD, KOOPMAN_POLY_HEX = CRC_code(CODE_LEN,PAYLOAD)    

end

sum_P = iseven.(sum(PP,dims=1))

EVEN_CODE = iszero(sum_P)

if !EVEN_CODE
    display("Not an even code!")
end

RR = PAYLOAD/CODE_LEN

# list of checks and variables nodes
NC = make_cn2vn_list(HH)
NV = make_vn2cn_list(HH)

NTHREADS = Threads.nthreads()

RGN_SEEDS = zeros(Int,NTHREADS)
for i in 1:NTHREADS
    RGN_SEEDS[i] = SEED + i - 1
end

num_ebn0 = length(EbN0)

Errors = Matrix{Int}(undef,num_ebn0,NTHREADS)
Trials = Matrix{Int}(undef,num_ebn0,NTHREADS)

FER = zeros(num_ebn0)

### testando ints

H_COLUMNS = zeros(Int,CODE_LEN)

for i in axes(HH,2)
    for j in axes(HH,1)
        if HH[j,i]
            H_COLUMNS[i] += 2^(REDUN - j)
        end
    end
end

plotlyjs()

Y_MEAN = readdlm("$CODE_LEN.txt",'\t',Float64,'\n')

if !TEST
    PRINT = false
    for k in eachindex(EbN0)

        # transform EbN0 in standard deviations
        variance = exp10.(-EbN0[k]/10) / (2*RR)
        stdev = sqrt.(variance)

        if FULL && MEAN && NUM_SEGMENTS > 1

            anchors = Vector{Int}(undef,NUM_SEGMENTS+1)
            offsets = Vector{Int}(undef,NUM_SEGMENTS)
            alphas = Vector{Float64}(undef,NUM_SEGMENTS)
            slopes = Vector{Int}(undef,NUM_SEGMENTS)
            min_slope = line_segmentation!(anchors,offsets,slopes,alphas,
                                            Y_MEAN[2:end,k],CODE_LEN,NUM_SEGMENTS)
            # p0 = plot(Y_MEAN[2:end,k]/min_slope,legend=false,xlabel="Rank order",title="ordered |signal|")
            #     for i in 1:2
            #         if i == 1
            #             j_vec = collect(anchors[i]:anchors[i+1])
            #         else
            #             j_vec = collect((anchors[i]+1):anchors[i+1])
            #         end
            #         segment = offsets[i] .+ (j_vec .- anchors[i])*slopes[i]
            #         plot!(p0,j_vec,segment, lw = 2, color = 2)
            #     end
            # display(p0)

        else
            anchors = [0]
            offsets = [0]
            slopes = [0]
        end

        stats = @timed Threads.@threads for i in 1:NTHREADS

            errors, trials = ORBGRAND_sim(
                MAX_ERRORS,
                PP,
                RGN_SEEDS[i],
                stdev,
                PRINT,
                H_COLUMNS,
                EVEN_CODE,
                MAX_QUERY,
                ABANDON,
                FULL,
                MAX_DEPTH,
                MEAN,
                NUM_SEGMENTS,
                anchors,
                offsets,
                slopes
                )
            # errors, trials = GRAND_sim(MAX_ERRORS,PP,RGN_SEEDS[i],STDEV,false,H_COLUMNS,EVEN_CODE,MAX_ERR_LOC_VEC_LEN)
            # _,_,errors,_,trials = simcore(PAYLOAD,CODE_LEN,nothing,STDEV,HH,PP,NC,NV,[0 0],"PEG",0,"Flooding","TANH",MAX_ERRORS,50,false,0,0.0,[1],0.0,RGN_SEEDS[i],false,false) 
            Trials[k,i] = trials
            Errors[k,i] = errors
            # Errors[k,i] = errors[end]

        end
        str = """Elapsed $(round(stats.time;digits=1)) seconds ($(round(stats.gctime/stats.time*100;digits=2))% gc time, $(round(stats.compile_time/stats.time*100,digits=2))% compilation time)"""
        println(str)
    end
    Total_trials = sum(Trials,dims=2)
    FER .= sum(Errors,dims=2)
    for k = 1:num_ebn0
        FER[k] /= Total_trials[k]
    end
    title = "ORBGRAND (N = $CODE_LEN, K = $PAYLOAD)"
    if FULL
        if MEAN && NUM_SEGMENTS > 1
            label = "$NUM_SEGMENTS-line ORBGRAND using mean"
        else
            label = "$NUM_SEGMENTS-line ORBGRAND"
        end
    else
        label = "basic ORBGRAND"
    end
    if !FULL
        p = plot(EbN0,log10.(FER),lw=2,label=label,title=title,size=(750,400))
    else
        plot!(p,EbN0,log10.(FER),lw=2,label=label,title=title)
    end
    
else
    VARIANCE = exp10.(-EbN0[end]/10) / (2*RR)
    STDEV = sqrt.(VARIANCE) 

    if FULL && MEAN

        ANCHORS = Vector{Int}(undef,NUM_SEGMENTS+1)
        OFFSETS = Vector{Int}(undef,NUM_SEGMENTS)
        ALPHAS = Vector{Float64}(undef,NUM_SEGMENTS)
        SLOPES = Vector{Int}(undef,NUM_SEGMENTS)
        MIN_SLOPE = line_segmentation!(ANCHORS,OFFSETS,SLOPES,ALPHAS,
                                        Y_MEAN[2:end,end],CODE_LEN,NUM_SEGMENTS)

        # p0 = plot(Y_MEAN[2:end,end]/MIN_SLOPE,legend=false,xlabel="Rank order",title="ordered mean |signal|")
        #     for i in 1:NUM_SEGMENTS
        #         if i == 1
        #             j_vec = collect(ANCHORS[i]:ANCHORS[i+1])
        #         else
        #             j_vec = collect((ANCHORS[i]+1):ANCHORS[i+1])
        #         end
        #         segment = OFFSETS[i] .+ (j_vec .- ANCHORS[i])*SLOPES[i]
        #         plot!(p0,j_vec,segment, lw = 2, color = 2)
        #     end
        # display(p0)
    else
        ANCHORS = [0]
        OFFSETS = [0]
        SLOPES = [0]        
    end

    @time errors, trials = ORBGRAND_sim(
                        MAX_ERRORS_TEST,
                        PP,
                        RGN_SEEDS[1],
                        STDEV,
                        PRINT,
                        H_COLUMNS,
                        EVEN_CODE,
                        MAX_QUERY,
                        ABANDON,
                        FULL,
                        MAX_DEPTH,
                        MEAN,
                        NUM_SEGMENTS,
                        ANCHORS,
                        OFFSETS,
                        SLOPES,
                        )

    display((trials, errors))
    # @benchmark ORBGRAND_sim(1,$PP,$(RGN_SEEDS[1]),$STDEV,$PRINT,$H_COLUMNS,$EVEN_CODE,$MAX_QUERY,$ABANDON) seconds = 30
end
