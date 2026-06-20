################################################################################
# Allan Eduardo Feitosa
# 16 Jun 2026
# My implementation of the soft decision GRAND algorithm (ORBGRAND)

function soft_grand!(
    err_loc_vec::Vector{Int},
    u::Vector{Int},
    code_len::Int,
    syndrome::Int,
    test_H_cols::Vector{Int},
    even_code::Bool,
    demod_parity::Int,
    max_query::Int,
    upper_bound::Int,
    b::Float64,
    bb::Float64,
    D::Vector{Int}
)

    n_guesses = 1

    # Logistic starting weight
    W = 1

    @inbounds @fastmath while n_guesses ≤ max_query && W ≤ upper_bound

        # Hamming weight
        w = max(1,ceil(Int,(b-sqrt(bb - 2*W))))
        # If the code is even and the parity is incorrect, increment 
        if even_code && rem(w,2) != demod_parity
            w += 1
        end
        while w ≤ code_len
            if W < w*(w+1)/2
                break
            else
                if w == 1
                    if test_H_cols[W] == syndrome
                        err_loc_vec[1] = W
                        return true, 1
                    end
                elseif w == 2
                    lim_sup = ceil(Int,W/2)-1
                    for i in 1:lim_sup
                        syn = test_H_cols[i] ⊻ test_H_cols[W-i]
                        if syn == syndrome
                            err_loc_vec[1] = i
                            err_loc_vec[2] = W-i
                            return true, w
                        end
                    end
                else
                    w_m_1 = w - 1
                    ### Landslide algorithm:
                    # With W being the target logistic weight, w being the Hamming 
                    # weight and code_len being the length of the string.

                    W1 = W - (w*(w+1))÷2
                    n1 = code_len - w 
                    k = 1
                    D_1 = 2 # Force the loop to run at least once

                    # Each loop generates a new integer partition
                    while D_1 ≥ 2

                        n_guesses = n_guesses +1           
                        
                        # mountain build
                        u_k = u[k]

                        for i in k+1:w
                            u[i] = u_k
                        end

                        sum_u = u[1]
                        for i in 2:w
                            sum_u += u[i]
                        end

                        dividend = W1 - sum_u
                        divisor = n1 - u_k
                        q = div(dividend, divisor)
                        r = dividend - q*divisor

                        w_m_q = w - q

                        if q != 0
                            for i in (w_m_q+1):w
                                u[i] = n1
                            end
                        end

                        if w_m_q > 0
                            u[w_m_q] += r
                        end

                        # evaluate acumulated drops
                        u_w = u[w]
                        for i in 1:w_m_1
                            D[i] = u_w - u[i]
                        end

                        syn = 0
                        for i in 1:w
                            idx = u[i] + i
                            err_loc_vec[i] = idx
                            syn ⊻= test_H_cols[idx]
                        end

                        if syn == syndrome
                            return true, w
                        end
                        # Find the last index with an accumulated drop >=2
                        k = 0
                        for outer k in w_m_1:-1:1
                            if D[k] ≥ 2
                                break
                            end
                        end
                        # Increase its index by one.
                        u[k] = u[k] + 1    
                        D_1 = D[1]
                    end
                    # clear vectors
                    for i in 1:w_m_1
                        u[i] = 0
                        D[i] = 0
                    end
                    u[w] = 0
                end
            end

            # Increment Hamming weight 
            w = w + 1
            # If the code is even
            if even_code && rem(w,2) != demod_parity
                w = w + 1
            end
        end
        W = W + 1
    end

    return false, 0

end