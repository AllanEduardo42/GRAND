################################################################################
# Allan Eduardo Feitosa
# 07 Jul 2026
# Function to find the valid logistic weights for full ORBGRAND

# OUTPUTS:
# W_list            :list of valid logistic weights for each segment
# W_lens            :length of the lists for each segment
# Ranges            :ranges of valid partial hamming weights

# INPUTS:
# anchors           :anchors of the segments
# offsets           :line intercepts of the piece-wise linear segment
# slopes            :slopes of each segment
# lim_sup           :superior limit for the logistic weights
# num_segments      :number of piece-wise linear segments

include("find_valid_hamming_weights.jl")

function find_valid_logistic_weights!(
    W_list::Matrix{Int},
    W_lens::Vector{Int},
    Ranges::Matrix{StepRange{Int,Int}},
    anchors::Vector{Int},
    offsets::Vector{Int},
    slopes::Vector{Int},
    W_upper_bound::Int,
    num_segments::Int
)
    @inbounds @fastmath for i in 1:num_segments
        j = 0
        for W in 0:W_upper_bound
            delta = anchors[i+1]-anchors[i]
            valid, range = find_valid_hamming_weights!(W,offsets[i],delta,slopes[i])
            if valid
                j += 1
                W_list[j,i] = W
                Ranges[j,i] = range
            end
        end
        W_lens[i] = j
    end    
end