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

const MAXC2V = 1e3                      # saturate values for C2V (Inf approx)
const MINC2V = -MAXC2V                  # -Inf approx

function qfunc(x::Float64)
    return 0.5 * erfc(x/sqrt(2))
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

PROTOCOL::String = "CRC"

# GRAND

MAX_ERRORS::Int = 5

ABANDON::Bool = true
MAX_QUERY::Int = 1_000_000_000

MAX_ERR_LOC_VEC_LEN::Int = CODE_LEN

# EbN0 = [1.0, 1.5, 2.0, 2.5, 3.0]
EbN0 = [3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]

FULL = true
NEW_PLOT::Bool = false
MAX_DEPTH = 100000

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

if !TEST
    PRINT = false
    for k in eachindex(EbN0)

        # transform EbN0 in standard deviations
        variance = exp10.(-EbN0[k]/10) / (2*RR)
        stdev = sqrt.(variance)

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
                MAX_DEPTH
                )
            # errors, trials = GRAND_sim(MAX_ERRORS,PP,RGN_SEEDS[i],stdev,false,H_COLUMNS,EVEN_CODE,MAX_ERR_LOC_VEC_LEN)
            # _,_,errors,_,trials = simcore(PAYLOAD,CODE_LEN,nothing,stdev,HH,PP,NC,NV,[0 0],"PEG",0,"Flooding","TANH",MAX_ERRORS,50,false,0,0.0,[1],0.0,RGN_SEEDS[i],false,false) 
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
        label = "2-line ORBGRAND"
    else
        label = "basic ORBGRAND"
    end
    if NEW_PLOT
        plot(EbN0, log10.(FER),label=label,title=title)
    else
        plot!(EbN0, log10.(FER),label=label,title=title)
    end
    
else
    variance = exp10.(-EbN0[end]/10) / (2*RR)
    stdev = sqrt.(variance) 
    # @benchmark GRAND_sim(1,$PP,$(RGN_SEEDS[1]),$stdev,$PRINT,$HH,$EVEN_CODE) seconds = 30
    # @time errors, trials = GRAND_sim(3,PP,RGN_SEEDS[1],stdev,PRINT,H_COLUMNS,EVEN_CODE,10)
    if PRINT
        errors, trials = ORBGRAND_sim(
                            6,
                            PP,
                            RGN_SEEDS[1],
                            stdev,
                            PRINT,
                            H_COLUMNS,
                            EVEN_CODE,
                            MAX_QUERY,
                            ABANDON,
                            FULL,
                            MAX_DEPTH)
    else
        @time @profview errors, trials = ORBGRAND_sim(
                            6,
                            PP,
                            RGN_SEEDS[1],
                            stdev,
                            PRINT,
                            H_COLUMNS,
                            EVEN_CODE,
                            MAX_QUERY,
                            ABANDON,
                            FULL,
                            MAX_DEPTH)
    end

    display((trials, errors))
    # @benchmark ORBGRAND_sim(1,$PP,$(RGN_SEEDS[1]),$stdev,$PRINT,$H_COLUMNS,$EVEN_CODE,$MAX_QUERY,$ABANDON) seconds = 30
end
