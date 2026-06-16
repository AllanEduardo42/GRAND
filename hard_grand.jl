################################################################################
# Allan Eduardo Feitosa
# 15 Jun 2026
# My implementation of the hard decision GRAND algorithm.

# INPUTS:
#
# candidate             : the candidate transmitted codeword
# err_loc_vec           : vector with the error locations 
# err_loc_vec_len       : number of error locations in err_loc_vec
# max_err_loc_vec_len   : maximum length of err_loc_vec
# code_len              : codeword length
# syndrome              : the base 10 correspondent of the current syndrome
# H_cols                : the base-10 integer correspondents of the columns of H
# inc                   : increment size of the length of err_loc_vec (1 or 2)
# Sum_H_cols            : sum of the columns of H

# This function searches for the binary  vector noise 'w' such that 
# H*w = syndrome.
#
# The indices were w are equal to 1 (true) are the error locations represented in
# err_loc_vec.
#
# For performance, the columns of the parity-check matrix H are transformed to 
# their corresponding base-10 integers and stored into the vector H_cols.
# 
# Since H*w = the sum of the columns of H corresponding to the indices where w
# is equal 1, we use Sum_H_cols to store the sums of the columns of H. Each
# index 'idx' of Sum_H_cols corresponds to the sum of the first 'idx' error
# locations currently contained in err_loc_vec. 
# Example, if err_loc_vec = [1,3,7,10,27], then
# Sum_H_cols[1] = H_cols[1]
# Sum_H_cols[2] = H_cols[1] ⊻ H_cols[3]              (⊻: xor operation)
# Sum_H_cols[3] = H_cols[1] ⊻ H_cols[3] ⊻ H_cols[7]
# Sum_H_cols[4] = H_cols[1] ⊻ H_cols[3] ⊻ H_cols[7] ⊻ H_cols[10]
# Therefore, Sum_H_cols[idx] = Sum_H_cols[idx-1] ⊻ H_cols[err_loc_vec[idx]]
#
# Note that H_cols[27] is not added (see lines 74 and following below).  

function hard_grand!(
    candidate::Vector{Bool},
    err_loc_vec::Vector{Int},
    err_loc_vec_len::Int,
    max_err_loc_vec_len::Int,
    code_len::Int,
    syndrome::Int,
    H_cols::Vector{Int},
    inc::Int,
    Sum_H_cols::Vector{Int}
)

    # Unless the code is even and we query only an even number of errors,
    # the next step is to search "one-error" noises.
    
    if err_loc_vec_len == 1
        # In the case where there is only one error, just look at the 
        # column values of H.
        for i in 1:code_len
            if H_cols[i] == syndrome
                err_loc_vec[1] = i
                candidate[i] ⊻= true
                return true
            end
        end
        # Jump to the "two-errors" noise guesses
        err_loc_vec[1] = 1
        err_loc_vec[2] = 2
        err_loc_vec_len = 2
        if inc == 2
            # if the code is even and we query only an odd number of errors,
            # skip to the "three-errors" noises.
            err_loc_vec[3] = 3
            err_loc_vec_len = 3
        end                       
    end    

    start_idx = 1 # the smaller index of err_loc_vec that has changed between updates
    
    loop = true

    @inbounds @fastmath while loop

        syn_flag = false          # indicates if the right noise guess was found

        # We must add the column values of H according to the error locations.
        if start_idx == 1
            # i.e., if all indices of err_loc_vec have changed, then we must 
            # change the value stored at Sum_H_cols[1]
            Sum_H_cols[1] = H_cols[err_loc_vec[1]]
            start_idx = 2   # the next index of Sum_H_cols to update       
        end
        # add the columns of H and stores in Sum_H_cols accordingly
        for idx in start_idx:(err_loc_vec_len-1)
            Sum_H_cols[idx] = Sum_H_cols[idx-1] ⊻ H_cols[err_loc_vec[idx]]
        end

        # Finally, we add the column of H corresponding to the last error
        # location, and compares the result with the syndrome
        last_err_loc = err_loc_vec[err_loc_vec_len]
        sum_H_cols = Sum_H_cols[err_loc_vec_len-1]
        for i in last_err_loc:code_len
            # we test all error locations from 'last_err_loc' until 'code_len'
            syn = sum_H_cols ⊻ H_cols[i]  
            if syn == syndrome
                err_loc_vec[err_loc_vec_len] = i
                syn_flag = true
                break
            end
        end

        if syn_flag 
            # if the right noise guess was found, sum 1 at each error location
            for idx in 1:err_loc_vec_len
                candidate[err_loc_vec[idx]] ⊻= true
            end
            return true
        else   
            # updates err_loc_vec
            
            success = false             # if err_loc_vec was successfully updated
            start_idx = 1
            
            # At this point, the last error location guess is always at 'code_len'.
            # We must find the index of err_loc_vec from which we start tp update it.
            # Example: if err_loc_vec = [1,3,5,62,63,64], and code_len = 64, then
            # 1) from index = 5, the updated would result [1,3,5,62,64,65], which is not allowed
            # 2) from index = 4, we have [1,3,5,63,64,65], again not allowed.
            # 3) from index = 3, we have [1,3,6,7,8,9], which is the valid update.
            # The routine below looks for this index, testing if the resulting 
            # last error location is larger than code_len.

            for idx = (err_loc_vec_len-1):(-1):1
                new_loc = (err_loc_vec[idx] + 1) + (err_loc_vec_len - idx)
                if new_loc > code_len
                    continue
                else
                    success = true      # we found the index from which to start 
                    # Next, we update err_loc_ vec from the last to the starting index
                    for j = err_loc_vec_len:(-1):idx    
                        err_loc_vec[j] = new_loc
                        new_loc -= 1
                    end
                    start_idx = idx     # we need the starting index outside the for loop
                    break
                end
            end

            if !success # if was not possible to update err_loc_vec with its current length               
                err_loc_vec_len += inc   
                if err_loc_vec_len > max_err_loc_vec_len
                    # then we must terminate simulation
                    loop = false 
                else  
                    # then we increment the length of err_loc_vec and initiates it
                    # accordinly
                    for idx = 1:err_loc_vec_len
                        err_loc_vec[idx] = idx
                    end   
                end   
            end
        end
    end

    return false
    
end