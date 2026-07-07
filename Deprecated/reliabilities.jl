using Plots
using Statistics
using Random

include("line_segmentation.jl")
include("integer_splitting.jl")
include("old_integer_splitting.jl")

plotlyjs()

num_segments = 3

code_len = 64
num_trials = 1

SNR = 7.0

rng = Xoshiro(30)

variance = exp10(-SNR/10)
stdev = sqrt(variance)
x = rand(rng,Bool,code_len)
y = (2*x .- 1) + stdev*randn(rng,code_len)
y_abs_sorted = sort(abs.(y))
anchors = Vector{Int}(undef,num_segments+1)
values = Vector{Int}(undef,num_segments)
alphas = Vector{Float64}(undef,num_segments)
betas = Vector{Int}(undef,num_segments)
min_beta = line_segmentation!(anchors,values,alphas,betas,y_abs_sorted,code_len,num_segments)

p = plot(legend=false)
plot!(p,y_abs_sorted/min_beta,title="SNR = $(SNR)dB",)
for i in 1:num_segments
    if i == 1
        j_vec = collect(anchors[i]:anchors[i+1])
    else
        j_vec = collect((anchors[i]+1):anchors[i+1])
    end
    segment = values[i] .+ (j_vec .- anchors[i])*betas[i]
    plot!(p,j_vec,segment, lw = 2, color = 2)
end
plot!(p,0:code_len,200*(0:code_len)/code_len,lw=2)
display(p)

##
lim_sup = code_len*(code_len+1)÷2
W = 100
L = binomial(W + num_segments-1,num_segments-1)

##
# split = zeros(Int,num_segments)
# weights = zeros(Int,num_segments,20)
# split = [24,20,0]

# success = find_next_valid_integer_splitting!(split,J,anchors,betas,W,0,0,num_segments,weights)
# display(success)
# display(split')
##
# W_splits_old = zeros(Int,L,num_segments)
# Hamming_weights_old = zeros(Int,L,num_segments,20)
# split = zeros(Int,num_segments)
# Weights = zeros(Int,num_segments,20)

# num_splits = old_integer_splitting!(W_splits_old,values,anchors,betas,W,num_segments,split,Hamming_weights_old,Weights)

# display(W_splits_old[1:num_splits,:])

W_splits_1 = zeros(Int,L,num_segments)
Hamming_weights_1 = zeros(Int,L,num_segments,20)
split = zeros(Int,num_segments-2)
Weights = zeros(Int,num_segments,20)

num_splits = integer_splitting!(W_splits_1,values,anchors,betas,W,num_segments,Hamming_weights_1,split,Weights)
display(W_splits_1[1:num_splits,:])

# display(W_splits_old[1:num_splits,:] == W_splits_1[1:num_splits,:])
##
# num_splits = 0
# for WW in 1:W    
#     global num_splits = integer_splitting!(W_splits,J,anchors,betas,WW,num_splits,num_segments,split,hamming_weights,weights)
# end

# display(W_splits[1:num_splits,:])

include("valid_W.jl")
L = lim_sup - values[end] + 1
Weights_list = zeros(Int,L,num_segments,20)
weights = zeros(Int,20)
W_list = zeros(Int,L,num_segments)
W_lens = zeros(Int,num_segments)
for i in 1:num_segments
    valid_W!(W_list,W_lens,values,anchors,betas,i,lim_sup,Weights_list,weights)
end

include("build_splits.jl")
W_splits_2 = zeros(Int,L,num_segments)
Hamming_weights_2 = zeros(Int,L,num_segments,20)
split = zeros(Int,num_segments-2)
Weights = zeros(Int,num_segments,20)

num_splits = build_splits!(W_splits_2,W,num_segments,split,W_list,W_lens,values[end]+betas[end],Weights_list,Hamming_weights_2,Weights)

display(W_splits_2[1:num_splits,:])

display(W_splits_1[1:num_splits,:] == W_splits_2[1:num_splits,:])

# display(Hamming_weights_1[1:num_splits,:,:] == Hamming_weights_old[1:num_splits,:,:])

display(Hamming_weights_1[1:num_splits,:,:] == Hamming_weights_2[1:num_splits,:,:])