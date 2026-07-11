################################################################################
# Allan Eduardo Feitosa
# 29 Jun 2026
# My implementation of the full ORBGRAND algorithm with three piece-wise linear
# segments to approximate the statistical model of the reliabilities

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
# W_upper_bound         : logistic weight W ≤ upper_bound = code_len*(code_len+1)÷2
# W_list                : list of valid logistic weights for each segment 
# W_lens                : length of the lists for each segment
# Ranges                : ranges of valid partial hamming weights
# anchors               : anchors of the segments
# offsets               : line intercepts of the piece-wise linear segment
# slopes                : slopes of each segment
# max_depth             : maximum number of error location vector at each level

include("partial_basic_ORBGRAND.jl")

function three_line_ORBGRAND!(
    err_loc_vec::Vector{Int},
    err_loc_vecs::Array{Int,3},
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

    a1 = anchors[1]
    a2 = anchors[2]
    a3 = anchors[3]
    a4 = anchors[4]

    o1 = offsets[1]
    o2 = offsets[2]
    o3 = offsets[3]

    s1 = slopes[1]
    s2 = slopes[2]
    s3 = slopes[3]

    W1_len = W_lens[1] 
    W2_len = W_lens[2]
    W3_len = W_lens[3]
    
    @inbounds @fastmath for W in 0:W_upper_bound
        # level 1 (first linear segment)
        for k1 in 1:W1_len, k2 in 1:W2_len
            W1 = W_list[k1,1]
            if W1 > W
                break
            end
            W2 = W_list[k2,2]
            if W2 > W || W1 + W2 > W
                continue
            end
            W3 = W - W1 - W2            # W = W1 + W2 + W3
            k3 = 0
            found_W3 = false
            for i in 1:W3_len
                # search for W3 in the list of valid logistic weights
                k3 += 1
                inlist = W_list[i,3]
                if inlist == W3
                    found_W3 = true
                    break
                elseif inlist > W3
                    break
                end
            end
            if found_W3
                weights_1 = Ranges[k1,1]
                weights_2 = Ranges[k2,2]
                weights_3 = Ranges[k3,3]
                for w1 in weights_1, w2 in weights_2
                    if w1 > 0
                        # partial basic ORBGRAND for the first segment
                        WL1 = div(W1 - o1*w1, s1)
                        # n1 = number of error location vectors in segment 1
                        n1 = partial_basic_ORBGRAND!(err_loc_vecs,
                                                     partition_vec,
                                                     cum_drops,
                                                     1,
                                                     w1,
                                                     WL1,
                                                     a1,
                                                     a2,
                                                     max_depth
                                                     )
                    else
                        n1 = 1          # if there is no error location vector,
                                        # this makes the code to enter the 
                                        # following loop.
                    end
                    if w2 > 0
                        # partial basic ORBGRAND for the second segment
                        WL2 = div(W2 - o2*w2, s2)
                        # n2 = number of error location vectors in segment 2
                        n2 = partial_basic_ORBGRAND!(err_loc_vecs,
                                                     partition_vec,
                                                     cum_drops,
                                                     2,
                                                     w2,
                                                     WL2,
                                                     a2,
                                                     a3,
                                                     max_depth
                                                     )
                    else
                        n2 = 1          # if there is no error location vector,
                                        # this makes the code to enter the 
                                        # following loop.
                    end

                    for nn1 in 1:n1, nn2 in 1:n2
                        syn1 = 0        # syndrome
                        if w1 > 0
                            # test error location vectors in segment 1
                            err_loc_vec_view_1 = view(err_loc_vecs,1:w1,nn1,1)
                            for i in 1:w1
                                # calculate the partial syndrome
                                syn1 ⊻= sorted_H_cols[err_loc_vec_view_1[i]]
                            end
                        end
                        syn2 = syn1        # syndrome
                        if w2 > 0
                            # test error location vectors in segment 2
                            err_loc_vec_view_2 = view(err_loc_vecs,1:w2,nn2,2)
                            for i in 1:w2
                                # calculate the partial syndrome
                                syn2 ⊻= sorted_H_cols[err_loc_vec_view_2[i]]
                            end
                        end
                        for w3 in weights_3
                            if w3 == 0
                                # if there is no error location vector in segment 3
                                n_guesses += 1
                                if syn2 == syndrome
                                    if w1 > 0
                                        err_loc_vec[1:w1] = err_loc_vec_view_1
                                        if w2 > 0
                                            err_loc_vec[w1+1:w1+w2] = err_loc_vec_view_2
                                        end
                                    elseif w2 > 0
                                        err_loc_vec[1:w2] = err_loc_vec_view_2
                                    end
                                    return true, w1 + w2
                                end
                                break   # nothing to carry on
                            end
                            # partial basic ORBGRAND for the second segment
                            WL3 = div(W3 - o3*w3, s3)
                            # n3 = number of error location vectors in segment 3
                            n3 = partial_basic_ORBGRAND!(err_loc_vecs,
                                                         partition_vec,
                                                         cum_drops,
                                                         3,
                                                         w3,
                                                         WL3,
                                                         a3,
                                                         a4,
                                                         max_depth
                                                         )
                            # test error location vectors in segment 2
                            for nn3 in 1:n3
                                n_guesses += 1
                                err_loc_vec_view_3 = view(err_loc_vecs,1:w3,nn3,3)
                                syn3 = syn2
                                for i in 1:w3
                                    # completes the syndrome evaluation
                                    syn3 ⊻= sorted_H_cols[err_loc_vec_view_3[i]]
                                end
                                if syn3 == syndrome
                                    if w1 > 0
                                        for i in 1:w1
                                            err_loc_vec[i] = err_loc_vec_view_1[i]
                                        end
                                    end
                                    if w2 > 0
                                        for i in 1:w2
                                            err_loc_vec[w1+i] = err_loc_vec_view_2[i]
                                        end
                                    end
                                    for i in 1:w3
                                        err_loc_vec[w1+w2+i] = err_loc_vec_view_3[i]
                                    end
                                    return true, w1 + w2 + w3
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
