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

# GRAND

MAX_ERRORS = 100

MAX_QUERY = 0

EbN0 = [1.0, 1.5, 2.0, 2.5, 3.0]

RR = KK/NN

PRINT = false

NTHREADS = Threads.nthreads()

RGN_SEEDS = zeros(Int,NTHREADS)
for i in 1:NTHREADS
    RGN_SEEDS[i] = SEED + i - 1
end


num_ebn0 = length(EbN0)

Errors = Matrix{Int}(undef,num_ebn0,NTHREADS)
Trials = Matrix{Int}(undef,num_ebn0,NTHREADS)

for k in eachindex(EbN0)

    stats = @timed Threads.@threads for i in 1:NTHREADS

        errors, trials = GRAND_sim(MAX_ERRORS ÷ NTHREADS,KK,NN,GG,RGN_SEEDS[i],MAX_QUERY,EbN0[k],RR,PRINT,HH,even_code) 
        
        Trials[k,i] = trials
        Errors[k,i] = errors

    end
    str = """Elapsed $(round(stats.time;digits=1)) seconds ($(round(stats.gctime/stats.time*100;digits=2))% gc time, $(round(stats.compile_time/stats.time*100,digits=2))% compilation time)"""
    println(str)
end

Total_trials = sum(Trials,dims=2)
FER = zeros(num_ebn0)
FER .= sum(Errors,dims=2)
for k = 1:num_ebn0
    FER[k] /= Total_trials[k]
end

# @time errors, trials = GRAND_sim(MAX_ERRORS,KK,NN,GG,RGN_SEEDS[1],MAX_QUERY,EbN0[1],RR,PRINT,HH,even_code)

# @profview errors, trials = GRAND_sim(MAX_ERRORS,KK,NN,GG,SEED,MAX_QUERY,EbN0[1],RR,PRINT,HH) 

# FER[i] = log10(errors/trials)

plotlyjs()
plot(EbN0, log10.(FER))

# TEST = rand(RNG, Bool, KK)

# Y = GG*TEST

# NOISE = rand(RNG,Bool,NN)
# Y_DEMOD = Vector{Bool}(undef,NN)
# @. Y_DEMOD = Y ⊻ NOISE

# ERR_VEC = zeros(Bool,NN)
# ERR_LOC_VEC = zeros(Int,NN)
# NEW_ERR_LOC_VEC = zeros(Int,NN)
# CANDIDATE = zeros(Bool,NN)
# hard_grand!(CANDIDATE, ERR_VEC, ERR_LOC_VEC, NEW_ERR_LOC_VEC, MAX_QUERY, Y_DEMOD,NC,NN)
