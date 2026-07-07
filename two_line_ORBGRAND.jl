################################################################################
# Allan Eduardo Feitosa
# 29 Jun 2026
# My implementation of the full ORBGRAND algorithm with two piece-wise linear
# segments to approximate the statistical model of the reliabilities

# OUTPUTS:
# err_loc_vec           : vector with the error locations 

# AUXILIARIES:
# err_loc_vecs_1        : segment-1 error location vectors
# err_loc_vecs_2        : segment-2 error location vectors
# partition_vec         : partition vector for landslide algorithm
# cum_drops             : accumulated drops for landslide algorithm

# INPUTS:
# syndrome              : base 10 correspondent of the current syndrome
# sorted_H_cols         : base-10 integer correspondents of the sorted columns of H
# max_query             : maximum number of queries
# W_upper_bound         : logistic weight W ≤ upper_bound = code_len*(code_len+1)÷2
# W_list                : list of valid logistic weights for each segment 
# W_lens                : length of the lists for each segment
# Ranges                : ranges of valid partial hamming weights
# anchors               : anchors of the segments
# offsets               : line intercepts of the piece-wise linear segment
# slopes                : slopes of each segment
# max_depth             : maximum number of error location vector at each level

include("partial_basic_ORBGRAND.jl")

function two_line_ORBGRAND!(
    err_loc_vec::Vector{Int},
    err_loc_vecs_1::Matrix{Int},
    err_loc_vecs_2::Matrix{Int},
    partition_vec::Vector{Int},
    cum_drops::Vector{Int},
    syndrome::Int,
    sorted_H_cols::Vector{Int},
    max_query::Int,
    W_upper_bound::Int,
    W_list::Matrix{Int},
    W_lens::Vector{Int},
    Ranges::Matrix{StepRange{Int,Int}},
    anchors::Vector{Int},
    offsets::Vector{Int},
    slopes::Vector{Int},
    max_depth::Int
)

    # number of noise guesses
    n_guesses = 1

    a2 = anchors[2]
    a3 = anchors[3]

    o1 = offsets[1]
    o2 = offsets[2]

    s1 = slopes[1]
    s2 = slopes[2]
    
    @inbounds @fastmath for W in 0:W_upper_bound

        err_loc_vec .= 0

        # level 1 (first linear segment)
        for k1 in 1:W_lens[1]  
            W1 = W_list[k1,1]
            if W1 > W
                break
            end
            W2 = W - W1         # W = W1 + W2
            k2 = 0
            found_W2 = false
            for i in 1:W_lens[2]
                # search for W2 in the list of valid logistic weights
                k2 += 1
                inlist = W_list[i,2]
                if inlist == W2
                    found_W2 = true
                    break
                elseif inlist > W2
                    break
                end
            end
            if found_W2
                weights_1 = Ranges[k1,1]
                weights_2 = Ranges[k2,2]
                for w1 in weights_1
                    if w1 > 0
                        # partial basic ORBGRAND for the first segment
                        WL1 = div(W1 - o1*w1, s1)
                        # n1 = number of error location vectors at segment 1
                        n1 = partial_basic_ORBGRAND!(err_loc_vecs_1,
                                                     partition_vec,
                                                     cum_drops,
                                                     w1,
                                                     WL1,
                                                     0,
                                                     a2,
                                                     max_depth
                                                     )
                    else
                        n1 = 1          # if there is no error location vector,
                                        # this makes the code to enter the 
                                        # following loop.
                    end
                    for nn1 in 1:n1
                        syn1 = 0        # syndrome
                        if w1 > 0
                            # test error location vectors in segment 2
                            err_loc_vec_view_1 = view(err_loc_vecs_1,1:w1,nn1)
                            for i in 1:w1
                                # calculate the partial syndrome
                                syn1 ⊻= sorted_H_cols[err_loc_vec_view_1[i]]
                            end
                        end
                        for w2 in weights_2
                            if w2 == 0
                                # if there is no error location vector in segment 2
                                n_guesses += 1
                                if syn1 == syndrome
                                    err_loc_vec[1:w1] = err_loc_vec_view_1
                                    return true, w1
                                end
                                break   # nothing to carry on
                            end
                            # partial basic ORBGRAND for the second segment
                            WL2 = div(W2 - o2*w2, s2)
                            # n2 = number of error location vectors at segment 2
                            n2 = partial_basic_ORBGRAND!(err_loc_vecs_2,
                                                         partition_vec,
                                                         cum_drops,
                                                         w2,
                                                         WL2,
                                                         a2,
                                                         a3,
                                                         max_depth
                                                         )
                            # test error location vectors in segment 2
                            for nn2 in 1:n2
                                n_guesses += 1
                                err_loc_vec_view_2 = view(err_loc_vecs_2,1:w2,nn2)
                                syn2 = syn1
                                for i in 1:w2
                                    # completes the syndrome evaluation
                                    syn2 ⊻= sorted_H_cols[err_loc_vec_view_2[i]]
                                end
                                if syn2 == syndrome
                                    if w1 > 0
                                        for i in 1:w1
                                            err_loc_vec[i] = err_loc_vec_view_1[i]
                                        end
                                    end
                                    for i in 1:w2
                                        err_loc_vec[w1+i] = err_loc_vec_view_2[i]
                                    end
                                    return true, w1 + w2
                                end
                            end
                        end
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
