################################################################################
# Allan Eduardo Feitosa
# 16 Jun 2026
# My implementation of the basic ORBGRAND algorithm

# INPUTS:
#
# err_loc_vec           : vector with the error locations 
# partition_vec         : partition vector for landslide algorithm
# code_len              : codeword length
# syndrome              : base 10 correspondent of the current syndrome
# sorted_H_cols         : base-10 integer correspondents of the sorted columns of H
# max_query             : maximum number of queries (if abandon == true)
# upper_bounds          : logistic weight W ≤ upper_bound = code_len*(code_len+1)÷2 
# cum_drops             : accumulated drops for landslide algorithm
# inc                   : increment size of the length of err_loc_vec (1 or 2)
# w_init                : The hamming weight w start (1 or 2)
# abandon               : if the algorithm stops when max_query is reached

function basic_ORBGRAND!(
    err_loc_vec::Vector{Int},
    partition_vec::Vector{Int},
    code_len::Int,
    syndrome::Int,
    sorted_H_cols::Vector{Int},
    max_query::Int,
    upper_bound::Int,
    cum_drops::Vector{Int},
    inc::Int,
    w_init::Int,
    abandon::Bool
)

    # number of noise guesses
    n_guesses = 1

    # Logistic starting weight
    W = 1

    @inbounds @fastmath while W ≤ upper_bound
        w = w_init
        while w ≤ code_len
            if w*(w+1)/2 > W    # maximum value of w given W
                break
            else
                if w == 1
                    # just test the corresponding columns of sorted_H_cols
                    n_guesses += 1
                    if sorted_H_cols[W] == syndrome
                        err_loc_vec[1] = W
                        return true, w
                    end
                elseif w == 2
                    # query noise at locations [1,W-1], [2,W-2], [3,W-3],...
                    lim_sup = ceil(Int,W/2)-1   # max value of a in [a,b]
                    for i in 1:lim_sup
                        n_guesses += 1
                        syn = sorted_H_cols[i] ⊻ sorted_H_cols[W-i]
                        if syn == syndrome
                            err_loc_vec[1] = i
                            err_loc_vec[2] = W-i
                            return true, w
                        end
                    end
                else
                    w_m_1 = w - 1
                    ### Landslide algorithm
                    # With W being the target logistic weight, w being the Hamming 
                    # weight and code_len being the length of the string.

                    W1 = W - (w*(w+1))÷2    # rescaling of W
                    n1 = code_len - w       # rescaling of n (code_len)

                    
                    sum_partition = 0       # used inside landslide loop
                    partition_k = 0         # used inside landslide loop
                    landslide = true
                    while landslide
                        # Each loop generates a new integer partition
                        n_guesses += 1         

                        # montain build
                        dividend = W1 - sum_partition
                        divisor = n1 - partition_k
                        q = div(dividend, divisor)
                        r = dividend - q*divisor

                        if q != 0
                            w_m_q = w - q
                            for i in (w_m_q+1):w
                                partition_vec[i] = n1
                            end
                        else
                            w_m_q = w
                        end

                        if w_m_q > 0
                            partition_vec[w_m_q] += r
                        end

                        # evaluate accumulated drops
                        u_w = partition_vec[w]
                        for i in 1:w_m_1
                            cum_drops[i] = u_w - partition_vec[i]
                        end

                        syn = 0
                        for i in 1:w
                            syn ⊻= sorted_H_cols[partition_vec[i] + i]
                        end

                        if syn == syndrome
                            for i in 1:w
                                err_loc_vec[i] = partition_vec[i] + i
                            end
                            return true, w
                        end

                        landslide = false
                        # Find the last index with an accumulated drop >=2
                        for k in w_m_1:-1:1
                            # Since cum_drops are monotonically decreasing, we 
                            # search from the end of the list. If cum_drops[1] < 2,
                            # the landslide algorihtm halts.
                            if cum_drops[k] ≥ 2

                                # Increase its index by one.
                                partition_vec[k] += 1

                                # mountain build
                                partition_k = partition_vec[k]

                                for i in k+1:w
                                    partition_vec[i] = partition_k
                                end

                                sum_partition = partition_vec[1]
                                for i in 2:w
                                    sum_partition += partition_vec[i]
                                end

                                landslide = true
                                break
                            end
                        end
                    end
                    # clear partition vector
                    for i in 1:w
                        partition_vec[i] = 0
                    end
                end
            end
            # Increment Hamming weight 
            w += inc
        end
        W += 1
        if abandon && n_guesses ≥ max_query
            break
        end
    end

    return false, 0

end

