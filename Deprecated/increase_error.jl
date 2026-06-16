################################################################################
# Allan Eduardo Feitosa
# 15 Jun 2026
# This function generates the next error location vector given the previous one

function increase_error!(
    err_loc_vec::Vector{Int},
    err_loc_vec_len::Int,
    max_err_loc_vec_len::Int,
    inc::Int,
    N::Int
)

    inbounds = false

    start_index = 1

    @inbounds @fastmath begin
    
        for i = err_loc_vec_len-1:(-1):1
            new_loc = (err_loc_vec[i] + 1) + (err_loc_vec_len - i)
            if new_loc > N
                continue
            else
                inbounds = true
                for j = err_loc_vec_len:(-1):i
                    err_loc_vec[j] = new_loc
                    new_loc -= 1
                end
                start_index = i
                break
            end
        end

        if !inbounds
            err_loc_vec_len += inc      # increments the length of err_loc_vec
            if err_loc_vec_len > max_err_loc_vec_len
                return err_loc_vec_len, start_index, false
            end
            for i = 1:err_loc_vec_len
                err_loc_vec[i] = i
            end        
        end
    end

    return err_loc_vec_len, start_index, true
    
end