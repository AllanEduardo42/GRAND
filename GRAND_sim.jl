################################################################################
# Allan Eduardo Feitosa
# 15 Jun 2026
# Core function to simulate the performance of the GRAND algorithm

include("hard_grand.jl")
include("auxiliary functions.jl")

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
    err_loc_vec = zeros(Int,max(2,max_err_loc_vec_len))

    # parity bits
    w = Vector{Bool}(undef,M)

    errors = 0
    trials = 0

    Sum_H_cols = zeros(Int,max_err_loc_vec_len-1)

    if print
        err_vec = zeros(Bool,code_len)         # error vector
    end

    @inbounds @fastmath while min(errors,trials - errors) < max_errors

        err_loc_vec .= 0

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

        for i in eachindex(cword)
            u = 1 - 2*cword[i]  
            signal[i] += u                  # sum the modulated signal to the noise              
            y_demod[i] = signbit(signal[i]) # hard demodulated signal
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
            for i in eachindex(true_err_vec)
                true_err_vec[i] = cword[i] ⊻ y_demod[i]
            end
            print_test("True error Vector", true_err_vec)
        end

        # Only if the code is even, then check the parity of 'y_demod' to decide
        # whether to query even or odd numbers of error locations; i,e., whether
        # length(err_loc_vec) is odd or even.

        even_errors = false # if or not to query "one-error" noises
        if even_code
            inc = 2     # Hamming weight increment
            if mod(sum(y_demod),2) == 0 
                # this means that err_loc_vec_len ∈ {2,4,6,8,...}
                even_errors = true
            end
        else
            inc = 1
        end  

        ### 4) noise guessing

        # The first noise guess is the "all-zeros" noise
        candidate .= y_demod               
        syndrome = fast_gf2_mat_mul(H_cols,candidate)
        zerosyn = iszero(syndrome)                 

        if !zerosyn

            # The next noise guesses are the "one-error" noises
            err_loc_vec[1] = 1          
            err_loc_vec_len = 1

            # If the code is even and we query only an even number of errors,
            # skip to "two-errors" noises.
            if even_code && even_errors
                err_loc_vec[2] = 2      # err_loc_vec = [1,2]
                err_loc_vec_len = 2
            end

            if err_loc_vec_len ≤ max_err_loc_vec_len

                zerosyn = hard_grand!(
                    candidate,
                    err_loc_vec,
                    err_loc_vec_len,
                    max_err_loc_vec_len,
                    code_len,
                    syndrome,
                    H_cols,
                    inc,
                    Sum_H_cols
                )
            end
        end
            
        ### 5) Calculate bit error and verify decoding
        for i in eachindex(biterror)
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
            if zerosyn
                err_vec .= false
                for i in eachindex(err_loc_vec)
                    if err_loc_vec[i] == 0
                        break
                    end
                    err_vec[err_loc_vec[i]] = true
                end
                print_test("Estimated Error Vector", err_vec)      
                print_test("Bit Error",biterror)
            else
                println()
                println("No noise candidate have been found!")
            end
        end      

    end

    return errors, trials

end