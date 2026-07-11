################################################################################
# Allan Eduardo Feitosa
# 9 Jul 2026
# My implementation of the full ORBGRAND algorithm with one piece-wise linear
# segment to approximate the statistical model of the reliabilities

# OUTPUTS:
# err_loc_vec           : vector with the error locations 

# AUXILIARIES:
# err_loc_vecs          : 3-dimensional arry with the error location vectors
# partition_vec         : partition vector for landslide algorithm
# cum_drops             : accumulated drops for landslide algorithm

# INPUTS:
# syndrome              : base 10 correspondent of the current syndrome
# sorted_H_cols         : base-10 integer correspondents of the sorted columns of H
# max_query             : maximum number of queries
# W_list                : list of valid logistic weights for each segment 
# W_len                 : length of the list of valid W
# Ranges                : ranges of valid partial hamming weights
# anchor                : anchors of the segments
# offset                : line intercept of the piece-wise linear segment
# slope                 : slope of linear segment
# max_depth             : maximum number of error location vector at each level

include("partial_basic_ORBGRAND.jl")

function one_line_ORBGRAND!(
    err_loc_vec::Vector{Int},
    err_loc_vecs::Array{Int,3},
    partition_vec::Vector{Int},
    cum_drops::Vector{Int},
    syndrome::Int,
    sorted_H_cols::Vector{Int},
    max_query::Int,
    W_list::Matrix{Int},
    W_len::Int,
    Ranges::Matrix{StepRange{Int,Int}},
    anchor::Int,
    offset::Int,
    slope::Int,
    max_depth::Int
)

    # number of noise guesses
    n_guesses = 1
    
    @inbounds @fastmath for k in 1:W_len
        W = W_list[k,1]
        weights = Ranges[k,1]
        for w in weights
            if w > 0
                # partial basic ORBGRAND for the first segment
                WL = div(W - offset*w, slope)
                # n = number of error location vectors at segment 1
                n = partial_basic_ORBGRAND!(err_loc_vecs,
                                            partition_vec,
                                            cum_drops,
                                            1,
                                            w,
                                            WL,
                                            0,
                                            anchor,
                                            max_depth
                                            )
            else
                n = 1          # if there is no error location vector,
                                # this makes the code to enter the 
                                # following loop.
            end
            for nn in 1:n
                syn = 0        # syndrome
                if w > 0
                    # test error location vectors in segment 2
                    err_loc_vec_view = view(err_loc_vecs,1:w,nn,1)
                    for i in 1:w
                        # calculate the partial syndrome
                        syn ⊻= sorted_H_cols[err_loc_vec_view[i]]
                    end
                    if syn == syndrome
                        err_loc_vec[1:w] = err_loc_vec_view
                        return true, w
                    end
                end
            end
        end
        if n_guesses ≥ max_query
            break
        end
    end

    return false, 0

end
