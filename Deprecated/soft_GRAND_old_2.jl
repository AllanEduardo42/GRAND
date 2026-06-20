################################################################################
# Allan Eduardo Feitosa
# 16 Jun 2026
# My implementation of the soft decision GRAND algorithm (ORBGRAND)

function soft_grand_2!(
    candidate::Vector{Bool},
    err_vec::Vector{Bool},
    err_loc_vec::Vector{Int},
    err_vec_perm::Vector{Bool},
    u::Vector{Int},
    code_len::Int,
    syndrome::Int,
    test_H_cols::Vector{Int},
    even_code::Bool,
    demod_parity::Int,
    max_query::Int,
    inv_perm::Vector{Int},
    c::Int,
    d::Vector{Int},
    D::Vector{Int}
    # idx_order::Vector{Int},
    # pG1::Float64,
    # sPG::Float64,
    # prob_demod::Vector{Float64},
    # K::Int,
)

    # Logistic starting weight
    wt = c + 1

    n_noises = 1

    cte = c*code_len + code_len*(code_len+1)/2

    n_guesses = 1

    @inbounds @fastmath while n_guesses ≤ max_query && wt ≤ cte

        # display(wt)

        # Hamming weight
        w = max(1,ceil(Int,(1+2*(code_len+c)-sqrt((1+2*(code_len+c))^2-8*wt))/2))
        # If the code is even and the parity is incorrect, increment 
        if even_code && rem(w,2) != demod_parity
            w += 1
        end
        while w ≤ code_len

            # display(" w = $w")
            # Logistic weight
            W = wt - c*w
            if W < w*(w+1)/2
                break
            else
                # display("w = $w, W = $W")
                # Make error vectors
                # Internally converts W and code_len to W' and code_len'.

                ### landslide
                W1 = W - (w*(w+1))÷2
                n1 = code_len - w 
                # Create the first integer partition
                j = 1
                # Start with empty vector and breaking at first index
                u .= 0
                k = 1

                # mountain_build!(u,k,w,W1,n1)
                for i in k+1:w
                    u[i] = u[k]
                end
                sum_u = 0
                for i in 1:w
                    sum_u += u[i]
                end
                W2 = W1 - sum_u
                if n1 == u[k]
                    throw(error("n1 cannot be equal to u[k]"))
                end
                q = floor(Int,W2/(n1-u[k]))
                # display("q = $q")
                r = W2 - q*(n1-u[k])
                if q != 0
                    for i in (w-q+1):w
                        u[i] = n1
                    end
                end
                if (w-q)> 0
                    u[w-q] += r
                end

                # Evaluate drops
                u_2 = u[1]
                for i in 1:w-1
                    u_1 = u_2
                    u_2 = u[i+1]
                    d[i] = u_2 - u_1
                end
                d[w] = 0
                # Evaluate accumulated drops
                sum = d[w]
                D[w] = sum
                for i in w-1:-1:1
                    sum += d[i]
                    D[i] = sum
                end
                # display("fora")
                # display([u d D])

                # for i in 1:w
                #     noise_locations[i,j] = u[i] + i
                # end

                n_guesses = n_guesses +1
                # err_vec_perm .= false
                # for i in 1:w
                #     err_vec_perm[u[i] + i] = true
                # end

                # syn = fast_gf2_mat_mul(test_H_cols,err_vec_perm)
                syn = 0
                for i in 1:w
                    err_loc_vec[i] = u[i] + i
                    syn ⊻= test_H_cols[err_loc_vec[i]]
                end

                if syn == syndrome
                    err_vec_perm .= false
                    for i in 1:w
                        err_vec_perm[err_loc_vec[i]] = true
                    end
                    for i in eachindex(err_vec)
                        err_vec[i] = err_vec_perm[inv_perm[i]]
                    end
                    candidate .⊻= err_vec
                    return true
                end

                # Each loop generates a new integer partition
                while D[1] ≥ 2
                    # Find the last index with an accumulated drop >=2
                    # k = findlast(x -> x ≥ 2, D)
                    k = 0
                    for outer k in w:-1:1
                        if D[k] ≥ 2
                            break
                        end
                    end
                    # Increase its index by one.
                    u[k] = u[k] + 1
                    # mountain_build!(u,k,w,W1,n1)
                    for i in k+1:w
                        u[i] = u[k]
                    end
                    sum_u = 0
                    for i in 1:w
                        sum_u += u[i]
                    end
                    W2 = W1 - sum_u
                    denominator = n1 - u[k]
                    if denominator == 0
                        throw(error("n1 cannot be equal to u[k]"))
                    end
                    # q = fld(W2,denominator)
                    q, r = divrem(W2,denominator)
                    # display("q = $q")
                    # r = W2 - q*denominator
                    if q != 0
                        for i in (w-q+1):w
                            u[i] = n1
                        end
                    end
                    if (w-q)> 0
                        u[w-q] += r
                    end
                    # Record the partition
                    j += 1
                    # Evaluate drops
                    # Evaluate drops
                    u_2 = u[1]
                    for i in 1:w-1
                        u_1 = u_2
                        u_2 = u[i+1]
                        d[i] = u_2 - u_1
                    end
                    d[w] = 0
                    # Evaluate acumuated drops
                    # D[w] = d[w]
                    sum = d[w]
                    D[w] = sum
                    for i in w-1:-1:1
                        sum += d[i]
                        D[i] = sum
                    end
                    # display("dentro")
                    # display([u d D])
                    # for i in 1:w
                    #     noise_locations[i,j] = u[i] + i
                    # end

                    # n_guesses = n_guesses +1
                    # err_vec_perm .= false
                    # for i in 1:w
                    #     err_vec_perm[u[i] + i] = true
                    # end

                    # syn = fast_gf2_mat_mul(test_H_cols,err_vec_perm)

                    syn = 0
                    for i in 1:w
                        err_loc_vec[i] = u[i] + i
                        syn ⊻= test_H_cols[err_loc_vec[i]]
                    end

                    if syn == syndrome
                        err_vec_perm .= false
                        for i in 1:w
                            err_vec_perm[err_loc_vec[i]] = true
                        end
                        for i in eachindex(err_vec)
                            err_vec[i] = err_vec_perm[inv_perm[i]]
                        end
                        candidate .⊻= err_vec
                        return true
                    end

                end

                # display(noise_locations[1:w,1:j])

                # n_noises = j
               
                # For each error vector
                # for j in 2:n_noises
                #     n_guesses = n_guesses +1
                #     err_vec_perm .= false
                #     for i in 1:w
                #         err_vec_perm[noise_locations[i,j]] = true
                #     end

                #     syn = fast_gf2_mat_mul(test_H_cols,err_vec_perm)

                #     if syn == syndrome
                #         for i in eachindex(err_vec)
                #             err_vec[i] = err_vec_perm[inv_perm[i]]
                #         end
                #         candidate .⊻= err_vec
                #         return true
                #     end
                # end
            end
            # Increment Hamming weight 
            w = w + 1
            # If the code is even
            if even_code && rem(w,2) != demod_parity
                w = w + 1
            end
        end
        wt = wt + 1
    end
    # If we max out on queries or total weight
    err_vec .= false

    return false

end