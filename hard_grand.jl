function hard_grand!(
    candidate::Vector{Bool},
    err_loc_vec::Vector{Int},
    err_loc_vec_len::Int,
    max_query::Int,
    y_demod::Vector{Bool}, 
    N::Int,
    syndrome::Vector{Bool},
    H::Matrix{Bool},
    inc::Int,
    Sum_syndromes::Matrix{Bool}
)

    n_guesses = 0

    if max_query == 0
        while_one = true
    else
        while_one = false
    end

    start_index = 1

    @fastmath @inbounds while while_one || n_guesses < max_query

        syn_flag = true

        if err_loc_vec_len == 1  
            for i=1:N   
                n_guesses += 1
                syn_flag = true
                h_line = view(H,:,i)       
                for j in eachindex(syndrome)
                    sym = h_line[j]
                    if sym != syndrome[j]
                        syn_flag = false
                        break
                    end
                end
                if syn_flag
                    err_loc_vec[1] = i
                    break
                end
            end
        else
            if start_index == 1
                idx = err_loc_vec[1]
                h_line = view(H,:,idx) 
                @simd for j in axes(Sum_syndromes,1)
                    Sum_syndromes[j,1] = h_line[j]
                end
                start_index += 1
            end
            if err_loc_vec_len > 2
                for i in start_index:err_loc_vec_len-1
                    idx = err_loc_vec[i]
                    h_line = view(H,:,idx)  
                    @simd for j in axes(Sum_syndromes,1)
                        sym = Sum_syndromes[j,i-1]
                        sym ⊻= h_line[j]
                        Sum_syndromes[j,i] = sym
                    end
                end
            end
            last_idx = err_loc_vec[err_loc_vec_len]
            for i = last_idx:N
                n_guesses += 1
                syn_flag = true 
                h_line = view(H,:,i)   
                for j in eachindex(syndrome)
                    sym = Sum_syndromes[j,err_loc_vec_len-1] 
                    sym ⊻= h_line[j]           
                    if sym != syndrome[j]
                        syn_flag = false
                        break
                    end
                end
                if syn_flag
                    err_loc_vec[err_loc_vec_len] = i
                    break
                end
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
        else
            err_loc_vec[err_loc_vec_len] = N
            err_loc_vec_len, start_index = increase_error!(err_loc_vec,err_loc_vec_len,inc,N)
        end

        # for i in 1:err_loc_vec_len
        #     print(err_loc_vec[i])
        #     print(" ")
        # end
        # println()        

    end

    return false
    
end

# This function generates the next error location vector given the previous one

function increase_error!(
    err_loc_vec::Vector{Int},
    err_loc_vec_len::Int,
    inc::Int,
    N::Int
)

    inbounds = false

    start_index = 1

    @fastmath @inbounds begin
    
        for i = err_loc_vec_len-1:-1:1

            inbounds = true
            index = err_loc_vec[i]
            err_loc_vec[i] = index + 1
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
            if inbounds
                start_index = i
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

    return err_loc_vec_len, start_index
end