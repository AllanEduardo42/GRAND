function hard_grand!(
    candidate::Vector{Bool},
    err_vec::Vector{Bool},
    err_loc_vec::Vector{Int},
    max_query::Int,
    y_demod::Vector{Bool},
    Nc::Vector{Vector{Int}},    
    N::Int
)

    err_loc_vec_len = 1

    inc = 1

    n_guesses = 0

    if max_query == 0
        while_one = true
    else
        while_one = false
    end

    zerosym = false

    @fastmath @inbounds while while_one || n_guesses < max_query

        @simd for i in eachindex(candidate)
            candidate[i] = y_demod[i] ⊻ err_vec[i]
        end
        zerosym = iszerosyndrome(candidate,Nc)   
        if zerosym
            return zerosym
        end
           
        err_loc_vec_len = gen_next_err!(err_vec, err_loc_vec, err_loc_vec_len, inc, N)
        
        n_guesses += 1

    end

    return zerosym
    
end

# This function generates the new error location vector, given a previous error
# location vector
#
# Inputs:
#   N           - Code length
#   err_loc_vec - Previous error location vector
#
#Outputs:
#   err_loc_vec - New error location vector. [] if a new error location vector cannot be generated
#   err_vec     - Binary error vector that corresponds to err_loc_vec. Zero vector if no error could be generated

function gen_next_err!(
    err_vec::Vector{Bool},
    err_loc_vec::Vector{Int},
    err_loc_vec_len::Int,
    inc::Int,
    N::Int
)
    @fastmath @inbounds begin

        err_vec .= false

        err_loc_vec_len = increase_error!(err_loc_vec,err_loc_vec_len,inc,N)

        for i in err_loc_vec
            if i == 0
                break
            end
            err_vec[i] = true
        end         
    end

    return err_loc_vec_len

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