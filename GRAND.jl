using Random
using LinearAlgebra
using Plots

include("Koopman.jl")
include("../LDPC/GF2_poly.jl")
include("../LDPC/GF2_functions.jl")
include("make_code.jl")
include("GRAND_sim.jl")
include("calc_syndrome.jl")

SEED::Int = 1234

KK = 32
NN = 60

GG, HH, CRC_POLY, HD, KOOPMAN_POLY_HEX = CRC_code(NN,KK)

even_code = iszero(mod.(sum(GG,dims=1),2))

if !even_code
    display("Not an even code!")
    throw(error())
end

# GRAND

MAX_ERRORS = 36*3

MAX_QUERY = 0

EbN0 = [1.0, 1.5, 2.0, 2.5, 3.0]
# EbN0 = [1.0]

RR = KK/NN

PRINT = false

TEST = false

NTHREADS = Threads.nthreads()

RGN_SEEDS = zeros(Int,NTHREADS)
for i in 1:NTHREADS
    RGN_SEEDS[i] = SEED + i - 1
end

num_ebn0 = length(EbN0)

Errors = Matrix{Int}(undef,num_ebn0,NTHREADS)
Trials = Matrix{Int}(undef,num_ebn0,NTHREADS)

FER = zeros(num_ebn0)

if !TEST
    for k in eachindex(EbN0)

        # transform EbN0 in standard deviations
        variance = exp10.(-EbN0[k]/10) / (2*RR)
        stdev = sqrt.(variance)

        stats = @timed Threads.@threads for i in 1:NTHREADS

            errors, trials = GRAND_sim(MAX_ERRORS ÷ NTHREADS,KK,NN,GG,RGN_SEEDS[i],MAX_QUERY,stdev,PRINT,HH,even_code)
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

    plotlyjs()
    plot(EbN0, log10.(FER))
else
    variance = exp10.(-EbN0[1]/10) / (2*RR)
    stdev = sqrt.(variance) 
    @time GRAND_sim(10,KK,NN,GG,RGN_SEEDS[1],MAX_QUERY,stdev,PRINT,HH,even_code)
end


# MSG = rand(Bool,KK)
# CWORD = [MSG; zeros(Bool,NN-KK)]
# _,R = divide_poly(CWORD,CRC_POLY)

# CWORD = [MSG;R]

# _,R = divide_poly(CWORD,CRC_POLY)

# CWORD == GG*MSG