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
include("calc_syndrome.jl")
include("../LDPC/PEG.jl")
include("../LDPC/LU_encoding.jl")
include("../LDPC/auxiliary setup functions.jl")
include("/home/allan/LDPC/Simulation core functions/encode_LDPC.jl")
include("/home/allan/LDPC/simcore.jl")
include("/home/allan/LDPC/Algorithms/Flooding.jl")

SEED::Int = 1234

KK::Int = 32
NN::Int = 64

MM::Int = NN - KK

TEST::Bool = false
PRINT::Bool = false

PROTOCOL::String = "CRC"

if PROTOCOL == "PEG"
    #Generate Parity-Check Matrix by the PEG algorithm

    LAMBDA = [0.21, 0.25, 0.25, 0.29, 0]
    RO = [1.0, 0, 0, 0, 0, 0]
    H_PEG, GIRTH = PEG(LAMBDA,RO,MM,NN)

    HH, LL, UU, FF = LU_encoding(H_PEG,0)

    PP = zeros(Bool,MM,KK)

    for k in axes(PP,2)
        PP[:,k] = gf2_solve_LU(LL,UU,HH[:,k])
    end

elseif PROTOCOL == "CRC"

    PP, HH, CRC_POLY, HD, KOOPMAN_POLY_HEX = CRC_code(NN,KK)

end

sum_P = iseven.(sum(PP,dims=1))

EVEN_CODE = iszero(sum_P)

if !EVEN_CODE
    display("Not an even code!")
end

# GRAND

MAX_ERRORS::Int = 4

MAX_ERR_LOC_VEC_LEN::Int = NN

EbN0 = [1.0, 1.5, 2.0, 2.5, 3.0]
# EbN0 = [1.0]

RR = KK/NN

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

H_COLUMNS = zeros(Int,NN)

for i in axes(HH,2)
    for j in axes(HH,1)
        if HH[j,i]
            H_COLUMNS[i] += 2^(MM - j)
        end
    end
end

plotlyjs()

if !TEST
    for k in eachindex(EbN0)

        # transform EbN0 in standard deviations
        variance = exp10.(-EbN0[k]/10) / (2*RR)
        stdev = sqrt.(variance)

        stats = @timed Threads.@threads for i in 1:NTHREADS


            errors, trials = GRAND_sim(MAX_ERRORS,PP,RGN_SEEDS[i],stdev,false,H_COLUMNS,EVEN_CODE,MAX_ERR_LOC_VEC_LEN)
            # _,_,errors,_,trials = simcore(KK,NN,nothing,stdev,HH,PP,NC,NV,[0 0],"PEG",0,"Flooding","TANH",MAX_ERRORS,50,false,0,0.0,[1],0.0,RGN_SEEDS[i],false,false) 
            Trials[k,i] = trials
            Errors[k,i] = errors

        end
        str = """Elapsed $(round(stats.time;digits=1)) seconds ($(round(stats.gctime/stats.time*100;digits=2))% gc time, $(round(stats.compile_time/stats.time*100,digits=2))% compilation time)"""
        println(str)
    end
    Total_trials = sum(Trials,dims=2)
    FER .= sum(Errors,dims=2)
    for k = 1:num_ebn0
        FER[k] /= Total_trials[k]
    end

    plot!(EbN0, log10.(FER),label="L=$MAX_ERR_LOC_VEC_LEN")
else
    variance = exp10.(-EbN0[1]/10) / (2*RR)
    stdev = sqrt.(variance) 
    # @benchmark GRAND_sim(1,$PP,$(RGN_SEEDS[1]),$stdev,$PRINT,$HH,$EVEN_CODE) seconds = 30
    @time errors, trials = GRAND_sim(3,PP,RGN_SEEDS[1],stdev,PRINT,H_COLUMNS,EVEN_CODE,10)
    display((errors, trials))

end
