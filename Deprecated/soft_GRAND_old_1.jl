################################################################################
# Allan Eduardo Feitosa
# 16 Jun 2026
# My implementation of the soft decision GRAND algorithm (ORBGRAND)

function soft_grand!(
    candidate::Vector{Bool},
    err_vec::Vector{Bool},
    err_vec_perm::Vector{Bool},
    noise_locations::Matrix{Int},
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

    while n_guesses ≤ max_query && wt ≤ cte

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
                n_noises = landslide!(noise_locations,u,W,w,code_len,d,D) 
                # display(noise_locations[1:w,1:n_noises])                
                # For each error vector
                for j in 1:n_noises
                    n_guesses = n_guesses + 1
                    err_vec_perm .= false
                    for i in 1:w
                        err_vec_perm[noise_locations[i,j]] = true
                    end
                    # pG = pG1*prod(prob_demod[idx_order[indexes[1:w]]]./(1 .- prob_demod[idx_order[indexes[1:w]]]))
                    # sPG += pG

                    syn = fast_gf2_mat_mul(test_H_cols,err_vec_perm)

                    if syn == syndrome
                        # display("zero syndrome")
                        # display(err_vec_perm)
                        for i in eachindex(err_vec)
                            err_vec[i] = err_vec_perm[inv_perm[i]]
                        end
                        # display(err_vec)
                        candidate .⊻= err_vec
                        # display(err_vec)
                        # if even_code == 0
                        #     return pG/(pG+(1-sPG)*(2^K-1)/(2^code_len-n_guesses))
                        # else
                        #     return pG/(pG+(1-sPG)*(2^K-1)/(2^(code_len-1)-n_guesses))
                        # end
                        return true
                    end
                end
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

# With W being the target logistic weight, w being the Hamming weight and code_len
# being the length of the string, W1 = W-w(w+1)/2 and n1 = code_len-w.

function landslide!(
    noise_locations::Matrix{Int},
    u::Vector{Int},
    W::Int,
    w::Int,
    code_len::Int,
    d::Vector{Int},
    D::Vector{Int}
)

    W1 = W - (w*(w+1))÷2
    n1 = code_len - w 
    # Create the first integer partition
    j = 1
    # Start with empty vector and breaking at first index
    u .= 0
    k = 1
    mountain_build!(u,k,w,W1,n1)
    for i in 1:w
        noise_locations[i,j] = u[i]
    end
    # Evaluate drops
    # d = circshift(u,-1) - u
    for i in 1:w-1
        d[i] = u[i+1] - u[i]
    end
    d[w] = 0
    # Evaluate accumulated drops
    # D = cumsum(d[end:-1:1])
    D[w] = d[w]
    for i in w-1:-1:1
        D[i] = d[i] + D[i+1]
    end
    # display("j=$j")
    # display([u[1:w] d[1:w] D[1:w]])
    # Each loop generates a new integer partition
    while D[1] ≥ 2
        # Find the last index with an accumulated drop >=2
        k = findlast(x -> x ≥ 2, D)
        # Increase its index by one.
        u[k] = u[k] + 1
        mountain_build!(u,k,w,W1,n1)
        # Record the partition
        j += 1
        for i in 1:w
            noise_locations[i,j] = u[i]
        end
        # Evaluate drops
        for i in 1:w-1
            d[i] = u[i+1] - u[i]
        end
        d[w] = 0
        # Evaluate acumuated drops
        D[w] = d[w]
        for i in w-1:-1:1
            D[i] = d[i] + D[i+1]
        end
        # display("j=$j")
        # display([u[1:w] d[1:w] D[1:w]])
    end
    for i in 1:j
        for k in 1:w
            noise_locations[k,i] += k
        end
    end

    return j
end

function mountain_build!(
    u::Vector{Int},
    k::Int,
    w::Int,
    W1::Int,
    n1::Int
)
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
end