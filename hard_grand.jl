function hard_grand!(
    candidate::Vector{Bool},
    err_loc_vec::Vector{Int},
    err_loc_vec_len::Int,
    max_query::Int,
    y_demod::Vector{Bool}, 
    N::Int,
    syndrome::Vector{Bool},
    H::Matrix{Bool},
    inc::Int
)

    n_guesses = 0

    if max_query == 0
        while_one = true
    else
        while_one = false
    end

    zerosym = false

    @fastmath @inbounds while while_one || n_guesses < max_query

        syn_flag = true
        for j in eachindex(syndrome)
            count = 0
            H_line = view(H,j,:)
            @simd for i in 1:err_loc_vec_len
                if H_line[err_loc_vec[i]]
                    count += 1
                end               
            end
            sym = isodd(count)
            if sym != syndrome[j]
                syn_flag = false
                break
            end
        end   

        if syn_flag
            candidate .= y_demod
            for i in err_loc_vec
                if i == 0
                    break
                end
                candidate[i] ⊻= true
            end
            return true
        end
           
        err_loc_vec_len = increase_error!(err_loc_vec,err_loc_vec_len,inc,N)
        
        n_guesses += 1

    end

    return zerosym
    
end

# This function generates the next error location vector given the previous one

function increase_error!(
    err_loc_vec::Vector{Int},
    err_loc_vec_len::Int,
    inc::Int,
    N::Int
)

    inbounds = false

    @fastmath @inbounds begin
    
        for i = err_loc_vec_len:-1:1

            if err_loc_vec[i] == N
                continue
            end
            inbounds = true
            index = err_loc_vec[i]
            err_loc_vec[i] = index + 1
            if i < err_loc_vec_len
                i_p = i + 1
                offset = index - i + 1
                for j = i_p : err_loc_vec_len
                    new_index = offset + j
                    if new_index ≤ N
                        err_loc_vec[j] = new_index
                    else
                        inbounds = false
                        break
                    end
                end
            end
            if inbounds
                break
            end
        end

        if !inbounds
            if err_loc_vec_len == N
                err_loc_vec_len = 0
                return err_loc_vec_len
            end
            err_loc_vec_len += inc      # increments the length of err_loc_vec
            @simd for i = 1:err_loc_vec_len
                err_loc_vec[i] = i
            end  
        
        end
    end

    return err_loc_vec_len
end