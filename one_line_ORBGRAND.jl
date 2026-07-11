################################################################################
# Allan Eduardo Feitosa
# 16 Jun 2026
# My implementation of the basic ORBGRAND algorithm to generate noise-effect 
# sequences

# OUTPUTS:
# err_loc_vec           : vector with the error locations 

# AUXILIARIES:
# partition_vec         : partition vector for landslide algorithm
# cum_drops             : accumulated drops for landslide algorithm

# INPUTS:
# code_len              : codeword length
# syndrome              : base 10 correspondent of the current syndrome
# sorted_H_cols         : base-10 integer correspondents of the sorted columns of H
# max_query             : maximum number of queries (if abandon == true)
# upper_bounds          : logistic weight W ≤ upper_bound = code_len*(code_len+1)÷2 
# inc                   : increment size of the length of err_loc_vec (1 or 2)
# w_init                : The hamming weight w start (1 or 2)
# cte                   : the constant 2*code_len + 1
# abandon               : if the algorithm stops when max_query is reached
# offset                : line intercept of the linear segment

function one_line_ORBGRAND!(
    err_loc_vec::Vector{Int},
    partition_vec::Vector{Int},
    cum_drops::Vector{Int},
    code_len::Int,
    syndrome::Int,
    sorted_H_cols::Vector{Int},
    max_query::Int,
    W_upper_bound::Int,
    inc::Int,
    w_init::Int,
    cte::Int,
    abandon::Bool,
    offset::Int
)

    # number of noise guesses
    n_guesses = 1
    lim_inf = 1
    cte2 = code_len 
    if offset ≠ 0 
        W_upper_bound += offset*code_len
        cte += offset
        lim_inf += offset
        cte2 += offset÷2
    end

    @inbounds @fastmath for W in lim_inf:W_upper_bound
        two_W = 2W
        w = w_init
        if W > cte2
            w += inc
            while w*(cte - w) < two_W
                w += inc
            end
        end
        WL = W - offset*w
        f_w = div(w*(w + 1),2)
        while f_w ≤ WL
            if w == 1
                # just test the corresponding columns of sorted_H_cols
                n_guesses += 1
                if sorted_H_cols[WL] == syndrome
                    err_loc_vec[1] = WL
                    return true, w
                end
            elseif w == 2
                # query noise at locations [1,WL-1], [2,WL-2], [3,WL-3],...
                b = min(WL-1,code_len)
                a = WL - b
                while a < b
                    n_guesses += 1
                    syn = sorted_H_cols[a] ⊻ sorted_H_cols[b]
                    if syn == syndrome
                        err_loc_vec[1] = a
                        err_loc_vec[2] = b
                        return true, w
                    end
                    a += 1
                    b -= 1
                end
            else
                w_m_1 = w - 1
                ### Landslide algorithm
                # With WL being the target logistic weight, w being the Hamming 
                # weight and code_len being the length of the string.

                W_prime = WL - f_w # rescaling of WL
                n_prime = code_len - w  # rescaling of n (code_len)
                
                sum_partition = 0       # used inside landslide loop
                partition_k = 0         # used inside landslide loop
                landslide = true
                while landslide
                    # Each loop generates a new integer partition
                    n_guesses += 1         

                    ### montain build
                    dividend = W_prime - sum_partition
                    divisor = n_prime - partition_k

                    if dividend > divisor
                        q = div(dividend, divisor)
                        w_m_q = w - q
                        if w_m_q > 0
                            # remainder of the division
                            r = dividend - q*divisor
                            partition_vec[w_m_q] += r
                        end
                        for i in (w_m_q+1):w
                            partition_vec[i] = n_prime
                        end
                    else
                        # w_m_q = 0, q = 0 and r = dividend
                        partition_vec[w] += dividend
                    end

                    ###

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

                            sum_partition = 0
                            for i in 1:w
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
            # Increment Hamming weight 
            w += inc
            f_w = div(w*(w + 1),2)
            WL = W - offset*w
        end
        if abandon && n_guesses ≥ max_query
            break
        end
    end

    return false, 0

end

