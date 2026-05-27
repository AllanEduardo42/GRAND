include("hard_grand.jl")

function GRAND_sim(
    max_errors::Int,
    K::Int,
    N::Int,
    G::Matrix{Bool},
    rgn_seed::Int,
    max_query::Int,
    ebn0::Float64,
    R::Float64,
    print::Bool,
    H::Matrix{Bool},
    even_code::Bool
)

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

    errors = 0
    trials = 0

    # transform EbN0 in standard deviations
    variance = exp10.(-ebn0/10) / (2*R)
    stdev = sqrt.(variance)

    syndrome = zeros(Bool,N-K)

    @fastmath @inbounds while min(errors,trials - errors) < max_errors

        decoded = false
        trials += 1

        ### 1) generate the random message
        rand!(rng,msg,Bool)

        ### 2) generate the codeword
        cword = G*msg

        ### 3) sum the noise to the modulated cword to produce the received signal
        randn!(rng,signal)              # put the noise in the vector 'signal'
        lmul!(stdev,signal)             # multiply by the standard deviation
        # 
        @simd for i in eachindex(cword)
            u = 1 - 2*cword[i]  
            signal[i] += u              # sum the modulated signal
        end

        @simd for i in eachindex(y_demod)
            y_demod[i] = signbit(signal[i])
        end

        # If the code is even, then check the parity of the demod to decide
        # whether to query even or odd

        # if even_code
        #     y_demod_parity = mod(sum(y_demod),2)
        #     if y_demod_parity == 1
        #         test_zero_noise = false
        #     end
        #     inc = 2 # Hamming weight increment
        # else
        #     inc = 1
        # end

        syndrome .= H*y_demod

        zerosym = iszero(syndrome)
        
        if zerosym 
            @simd for i in eachindex(biterror)
                biterror[i] = y_demod[i] ⊻ cword[i]
            end
        else          
            err_loc_vec .= 0      
            err_loc_vec[1] = 1
            err_loc_vec_len = 1
            inc = 1
            zerosym = hard_grand!(candidate, err_loc_vec,err_loc_vec_len,max_query,y_demod,N,syndrome,H,inc)
            
            @simd for i in eachindex(biterror)
                biterror[i] = candidate[i] ⊻ cword[i]
            end

        end

        if print
            @simd for i in eachindex(gabarito)
                gabarito[i] = cword[i] ⊻ y_demod[i]
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
            println("""
________________________________________________________________________________

                                    TRIAL #$trials
________________________________________________________________________________
""")
            display("trial = $trials, errors = $errors")
            print_test("Message",msg)
            print_test("Codeword",cword)
            print_test("Demodulated",y_demod)
            print_test("True error Vector", gabarito)
            # print_test("Estimated Error Vector", err_vec)
            display("zerosym = $zerosym")       
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