include("hard_grand.jl")

function GRAND_sim(
    max_errors::Int,
    P::Matrix{Bool},
    rgn_seed::Int,
    max_query::Int,
    stdev::Float64,
    print::Bool,
    H::Matrix{Bool},
    even_code::Bool
)

    M,N = size(H)

    K = N - M

    # Set the random seeds
    rng = Xoshiro(rgn_seed)

    # payload
    msg = Vector{Bool}(undef,K)

    # codeword
    cword = Vector{Bool}(undef,N)

    biterror = Vector{Bool}(undef,N)

    # received signal
    signal = Vector{Float64}(undef,N)

    y_demod = Vector{Bool}(undef,N)
    gabarito = Vector{Bool}(undef,N)
    candidate = Vector{Bool}(undef,N)
    err_loc_vec = zeros(Int,N)

    w = Vector{Bool}(undef,M)

    errors = 0
    trials = 0

    syndrome = zeros(Bool,N-K)

    Sum_syndromes = zeros(Bool,M,N-1)

    if print
        err_vec = zeros(Bool,N)         # error vector
    end

    @fastmath @inbounds while min(errors,trials - errors) < max_errors

        decoded = false
        trials += 1

        ### 1) generate the random message
        rand!(rng,msg,Bool)

        ### 2) generate the codeword
        cword[1:K] .= msg
        _gf2_mat_mult!(w,P,msg)
        cword[K+1:end] .= w

        # if print
            if !iszero(H*cword)
                throw(error(lazy"encoding error"))
            end
        # end

        ### 3) sum the noise to the modulated cword to produce the received signal
        randn!(rng,signal)              # put the noise in the vector 'signal'
        lmul!(stdev,signal)             # multiply by the standard deviation

        @simd for i in eachindex(cword)
            u = 1 - 2*cword[i]  
            signal[i] += u              # sum the modulated signal
            y_demod[i] = signbit(signal[i])
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
            @simd for i in eachindex(gabarito)
                gabarito[i] = cword[i] ⊻ y_demod[i]
            end
            print_test("True error Vector", gabarito)
        end

        _gf2_mat_mult!(syndrome,H,y_demod)           # syndrome

        # If the code is even, then check the parity of the demod to decide
        # whether to query even or odd

        test_zero_noise = true
        init_type_1 = true
        if even_code
            y_demod_parity = mod(sum(y_demod),2)
            if y_demod_parity == 1
                test_zero_noise = false
            else
                init_type_1 = false
            end
            inc = 2                     # Hamming weight increment
        else
            inc = 1
        end

        zerosym = false

        if test_zero_noise
            zerosym = iszero(syndrome)
        end

        err_loc_vec .= 0                
        
        if zerosym 
            @simd for i in eachindex(biterror)
                biterror[i] = y_demod[i] ⊻ cword[i]
            end
        else      
            err_loc_vec[1] = 1
            if init_type_1
                err_loc_vec_len = 1
            else
                err_loc_vec[2] = 2
                err_loc_vec_len = 2
            end

            zerosym = hard_grand!(
                candidate,
                err_loc_vec,
                err_loc_vec_len,
                max_query,
                y_demod,
                N,
                syndrome,
                H,
                inc,
                Sum_syndromes
            )
            
            # Calculate bit error
            @simd for i in eachindex(biterror)
                biterror[i] = candidate[i] ⊻ cword[i]
            end

        end   

        if zerosym
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