################################################################################
# Allan Eduardo Feitosa
# 29 Jun 2026
# This is a version of the basic ORBGRAND algorithm for partial noise-effect
# sequence generation to be used in the full ORBGRAND algortihm.

# OUTPUTS:
# err_loc_vecs          : all error locations vectors

# AUXILIARIES:
# partition_vec         : partition vector for landslide algorithm
# cum_drops             : accumulated drops for landslide algorithm

# INPUTS:
# w                     : hamming weight
# WL                    : normalized logistic weight div(W - offset*w, slope)
# anchor_1              : anchor of the previous segment
# anchor_2              : anchor of the current segment
# max_depth             : maximum number of error location vector at each level


function partial_basic_ORBGRAND!(
    err_loc_vecs::Matrix{Int},
    partition_vec::Vector{Int},
    cum_drops::Vector{Int},
    w::Int,
    WL::Int,
    anchor_1::Int,
    anchor_2::Int,
    max_depth::Int
)
    @fastmath begin
        n = 0
        segment_len = anchor_2 - anchor_1
        if w == 1
            n += 1
            err_loc_vecs[1,n] = WL + anchor_1
        elseif w == 2
            # query noise at locations [anchor_1,W-1], [anchor_1+1,W-2], [anchor_1+2,W-3],...
            b = min(WL-1,anchor_2-anchor_1)
            a = WL - b
            while a < b
                n += 1
                err_loc_vecs[1,n] = a + anchor_1
                err_loc_vecs[2,n] = b + anchor_1
                a += 1
                b -= 1            
            end
        else
            w_m_1 = w - 1
            ### Landslide algorithm
            # With W being the target logistic weight, w being the Hamming 
            # weight and code_len being the length of the string.

            W_prime = WL - div(w*(w+1),2) # rescaling of W
            n_prime = segment_len - w   # rescaling of n (code_len)
            
            sum_partition = 0       # used inside landslide loop
            partition_k = 0         # used inside landslide loop
            landslide = true
            while landslide
                # Each loop generates a new integer partition
                n += 1   

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
    
                for i in 1:w
                    err_loc_vecs[i,n] = partition_vec[i] + i + anchor_1
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
                if n == max_depth
                    break
                end
            end
            # clear partition vector
            for i in 1:w
                partition_vec[i] = 0
            end
        end
    end

    return n

end