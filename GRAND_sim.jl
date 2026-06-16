################################################################################
# Allan Eduardo Feitosa
# 15 Jun 2026
# Core function to simulate the performance of the GRAND algorithm

include("hard_grand.jl")

function GRAND_sim(
    max_errors::Int,
    P::Matrix{Bool},
    rgn_seed::Int,
    stdev::Float64,
    print::Bool,
    H_cols::Vector{Int},
    even_code::Bool,
    max_err_loc_vec_len::Int
)

    M,K = size(P)

    code_len = M + K

    # Set the random seeds
    rng = Xoshiro(rgn_seed)

    # payload
    msg = Vector{Bool}(undef,K)

    # codeword
    cword = Vector{Bool}(undef,code_len)

    biterror = Vector{Bool}(undef,code_len)

    # received signal
    signal = Vector{Float64}(undef,code_len)

    # demodulated signal
    y_demod = Vector{Bool}(undef,code_len)

    # true error vector
    true_err_vec = Vector{Bool}(undef,code_len)

    # found hard decision candidate
    candidate = Vector{Bool}(undef,code_len)

    # vector with the error locations
    err_loc_vec = zeros(Int,max_err_loc_vec_len)

    # parity bits
    w = Vector{Bool}(undef,M)

    errors = 0
    trials = 0

    Sum_H_cols = zeros(Int,max_err_loc_vec_len-1)

    if print
        err_vec = zeros(Bool,code_len)         # error vector
    end

    @fastmath @inbounds while min(errors,trials - errors) < max_errors

        decoded = false
        trials += 1

        ### 1) generate the random message
        rand!(rng,msg,Bool)

        ### 2) generate the codeword
        _gf2_mat_mult!(w,P,msg)         # w = P*msg
        cword[1:K] .= msg
        cword[K+1:end] .= w

        # test encoding
        if fast_gf2_mat_mul(H_cols,cword) != 0
            throw(error(lazy"encoding error"))
        end

        ### 3) sum the noise to the modulated cword to produce the received signal
        randn!(rng,signal)              # put the noise in the vector 'signal'
        lmul!(stdev,signal)             # multiply by the standard deviation

        @simd for i in eachindex(cword)
            u = 1 - 2*cword[i]  
            signal[i] += u                  # sum the modulated signal to the noise              
            y_demod[i] = signbit(signal[i]) # demodulated signal
        end

        if print
            println("""
________________________________________________________________________________

                                    TRIAL #$trials
________________________________________________________________________________
""")
            display("trial = $trials, errors = $errors")
            print_test("Message",msg)
            print_test("Codeword",cword)
            print_test("Demodulated",y_demod)
             # True error vector
            @simd for i in eachindex(true_err_vec)
                true_err_vec[i] = cword[i] ⊻ y_demod[i]
            end
            print_test("True error Vector", true_err_vec)
        end

        # Only if the code is even, then check the parity of 'y_demod' to decide
        # whether to query even or odd numbers of error locations; i,e., whether
        # length(err_loc_vec) is odd or even.

        len_one = true # if or not the "one-error" noise guesses are performed
        if even_code
            inc = 2     # Hamming weight increment
            if mod(sum(y_demod),2) == 0 
                # this means that err_loc_vec_len ∈ {2,4,6,8,...}
                len_one = false
            end
        else
            inc = 1
        end  

        ### 4) noise guessing

        # The first noise guess is the "all-zeros" noise
        candidate .= y_demod            
        err_loc_vec .= 0   
        syndrome = fast_gf2_mat_mul(H_cols,candidate)
        zerosyn = iszero(syndrome)        
        
        if !zerosyn
            if len_one
                # In the case where there is only one error, just look at the 
                # column values of H.
                for i in 1:code_len
                    if H_cols[i] == syndrome
                        candidate[i] ⊻= true
                        zerosyn = true
                        break
                    end
                end                
            end
        end                    

        if !zerosyn
            # The "hard_grand" function starts with err_loc_vec = [1,2]
            err_loc_vec[1] = 1
            err_loc_vec[2] = 2
            zerosyn = hard_grand!(
                candidate,
                err_loc_vec,
                max_err_loc_vec_len,
                code_len,
                syndrome,
                H_cols,
                inc,
                Sum_H_cols
            )
        end
            
        ### 5) Calculate bit error and verify decoding
        @turbo for i in eachindex(biterror)
            biterror[i] = candidate[i] ⊻ cword[i]
        end

        if zerosyn
            if iszero(biterror)
                decoded = true
            end
        end

        if !decoded
            errors += 1
        end  

        if print
            err_vec .= false
            for i in eachindex(err_loc_vec)
                if err_loc_vec[i] == 0
                    break
                end
                err_vec[err_loc_vec[i]] = true
            end
            print_test("Estimated Error Vector", err_vec)      
            print_test("Bit Error",biterror)
        end      

    end

    return errors, trials

end

function print_test(
    text::String,
    array::Vector{Bool}

)    
    println()
    print("$text (L = $(length(array))):")
    for i in eachindex(array)
        if i%80 == 1
            println()
        end
        print(Int(array[i]))
    end
    println()
end 

function fast_gf2_mat_mul(
    int_vec::Vector{Int},
    bool_vec::Vector{Bool}
)

    if length(int_vec) != length(bool_vec)
        throw(error(lazy"dimensions must be equal"))
    end

    result = 0
    for i in eachindex(int_vec)
        if bool_vec[i]
            result ⊻= int_vec[i]
        end
    end

    return result

end