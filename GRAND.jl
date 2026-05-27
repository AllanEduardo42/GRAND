using Random
using LinearAlgebra
using Plots

include("Koopman.jl")
include("../LDPC/GF2_poly.jl")
include("../LDPC/GF2_functions.jl")
include("make_code.jl")
include("GRAND_sim.jl")
include("../LDPC/auxiliary setup functions.jl")
include("calc_syndrome.jl")

SEED::Int = 1234

RNG = Xoshiro(SEED)

KK = 32
NN = 52

GG, HH, CRC_POLY, HD, KOOPMAN_POLY_HEX = CRC_code(NN,KK)

# GRAND

MAX_ERRORS = 100

MAX_QUERY = 0

NC = make_cn2vn_list(HH)

LEN = length.(NC)
INDICES = sortperm(LEN)
NC = NC[INDICES]

EbN0 = [1.0, 1.5, 2.0, 2.5, 3.0]

RR = KK/NN

FER = zeros(length(EbN0))

PRINT = false

NTHREADS = Threads.nthreads()

for i in eachindex(EbN0)

    @time errors, trials = GRAND_sim(MAX_ERRORS,KK,NN,GG,SEED,MAX_QUERY,NC,EbN0[i],RR,PRINT)
    FER[i] = log10(errors/trials)

end

# @time errors, trials = GRAND_sim(MAX_ERRORS,KK,NN,GG,SEED,MAX_QUERY,NC,EbN0[1],RR,PRINT)

@profview errors, trials = GRAND_sim(MAX_ERRORS,KK,NN,GG,SEED,MAX_QUERY,NC,EbN0[1],RR,PRINT)

# FER[i] = log10(errors/trials)

plotlyjs()
plot(EbN0, FER)

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
