################################################################################
# Allan Eduardo Feitosa
# 15 Jun 2026
# Core function to simulate the performance of the ORBGRAND algorithm

include("basic_ORBGRAND.jl")
include("two_line_ORBGRAND.jl")
include("line_segmentation.jl")
include("find_valid_logistic_weights.jl")
include("auxiliary functions.jl")

function ORBGRAND_sim(
    max_errors::Int,
    P::Matrix{Bool},
    rgn_seed::Int,
    stdev::Float64,
    printtest::Bool,
    H_cols::Vector{Int},
    even_code::Bool,
    max_query::Int,
    abandon::Bool,
    full::Bool,
    max_depth::Int,
    mean::Bool,
    anchors_mean::Vector{Int},
    offsets_mean::Vector{Int},
    slopes_mean::Vector{Int}
)

    M,K = size(P)

    code_len = M + K

    # ORBGRAND constants
    W_upper_bound = code_len*(code_len+1)÷2   # logistic weight W ≤ upper_bound
    cte = 2*code_len + 1

    # Set the random seeds
    rng = Xoshiro(rgn_seed)

    # payload
    msg = Vector{Bool}(undef,K)

    # codeword
    cword = Vector{Bool}(undef,code_len)

    biterror = Vector{Bool}(undef,code_len)

    # received signal
    signal = Vector{Float64}(undef,code_len)

    # absolute offsets of the signal
    abs_signal = Vector{Float64}(undef,code_len)

    # demodulated signal
    y_demod = Vector{Bool}(undef,code_len)

    # true error vector
    true_err_vec = Vector{Bool}(undef,code_len)

    # found hard decision candidate
    candidate = Vector{Bool}(undef,code_len)

    # parity bits
    parity_bits = Vector{Bool}(undef,M)

    errors = 0
    trials = 0

    err_vec = zeros(Bool,code_len)              # error vector

    err_loc_vec = Vector{Int}(undef,code_len)   # error location vector

    partition_vec = Vector{Int}(undef,code_len) # partition vector

    cum_drops = Vector{Int}(undef,code_len)     # accumulated drops

    sorted_ind = Vector{Int}(undef,code_len)    # sorted indices

    sorted_H_cols = Vector{Int}(undef,code_len)

    sorted_abs_signal = Vector{Float64}(undef,code_len)

    if full  
        err_loc_vecs_1 = zeros(Int,code_len,max_depth)
        err_loc_vecs_2 = zeros(Int,code_len,max_depth)
        W_lens = zeros(Int,2)
        Ranges = Matrix{StepRange{Int,Int}}(undef,W_upper_bound,2)
        W_list = Matrix{Int}(undef,W_upper_bound,2)
        if !mean      
            anchors = Vector{Int}(undef,3)
            offsets = Vector{Int}(undef,2)
            alphas = Vector{Float64}(undef,2)
            slopes = Vector{Int}(undef,2)
        else
            anchors = copy(anchors_mean)
            offsets = copy(offsets_mean)
            slopes = copy(slopes_mean)
            W_list .= 0
            find_valid_logistic_weights!(W_list,W_lens,Ranges,anchors,offsets,
                                                        slopes,W_upper_bound,2)
        end        
    end

    @inbounds @fastmath while min(errors,trials - errors) < max_errors

        err_loc_vec .= false

        decoded = false
        trials += 1

        ### 1) generate the random message
        rand!(rng,msg,Bool)

        ### 2) generate the codeword
        _gf2_mat_mult!(parity_bits,P,msg)         # w = P*msg
        cword[1:K] .= msg
        cword[K+1:end] .= parity_bits

        # test encoding
        if fast_gf2_mat_mul(H_cols,cword) != 0
            throw(error(lazy"encoding error"))
        end

        ### 3) sum the noise to the modulated cword to produce the received signal
        randn!(rng,signal)              # put the noise in the vector 'signal'
        lmul!(stdev,signal)             # multiply by the standard deviation

        for i in eachindex(cword)
            bpsk = 1 - 2*cword[i]  
            signal[i] += bpsk                  # sum the modulated signal to the noise     
            abs_signal[i] = abs(signal[i]) 
            y_demod[i] = signbit(signal[i]) # hard demodulated signal
        end

        if printtest
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

        if even_code
            inc = 2
            if reduce(xor,y_demod)
                w_init = 1
            else
                w_init = 2
            end
        else
            inc = 1
            w_init = 1
        end  

        ### 4) noise guessing

        # The first noise guess is the "all-zeros" noise
        err_vec .= false
        candidate .= y_demod               
        syndrome = fast_gf2_mat_mul(H_cols,candidate)
        zerosyn = iszero(syndrome)   

        if !zerosyn && (!abandon || max_query > 1)

            sortperm!(sorted_ind,abs_signal,alg=QuickSort)

            # This is the H columns reordered to put in ML order
            for i in eachindex(sorted_ind)
                sorted_abs_signal[i] = abs_signal[sorted_ind[i]]
                sorted_H_cols[i] = H_cols[sorted_ind[i]]
            end

            if printtest
                sorted_true_err_vec = true_err_vec[sorted_ind]
                print_test("Sorted true error vector", sorted_true_err_vec)

                sorted_true_err_loc_vec = zeros(Int,code_len)
                
                k = 0
                for i in eachindex(sorted_true_err_vec)
                    if sorted_true_err_vec[i]
                        k+= 1
                        sorted_true_err_loc_vec[k] = i
                    end 
                end
                println()
                display("Sorted true error location vector")
                display(sorted_true_err_loc_vec[1:k]')
            end

            partition_vec .= 0
            cum_drops .= 0

            if full
                if !mean
                    min_slope = line_segmentation!(anchors,offsets,slopes,alphas,
                                                sorted_abs_signal,code_len,2)
                    W_list .= 0
                    find_valid_logistic_weights!(W_list,W_lens,Ranges,anchors,offsets,
                                                            slopes,W_upper_bound,2)
                end                

                # display(W_list)

                # if sorted_true_err_loc_vec[1:k] == [4,5,9,10,11,120]
                #     display("anchors: ")
                #     display(anchors')
                #     display("offsets: ")
                #     display(offsets')
                #     display("slopes: ")
                #     display(slopes')
                # end

                zerosyn, w = two_line_ORBGRAND!(
                    err_loc_vec,
                    err_loc_vecs_1,
                    err_loc_vecs_2,
                    partition_vec,
                    cum_drops,
                    syndrome,
                    sorted_H_cols,
                    max_query,
                    W_upper_bound,
                    W_list,
                    W_lens,
                    Ranges,
                    anchors,
                    offsets,
                    slopes,
                    max_depth
                )
            else
                zerosyn, w = basic_ORBGRAND!(
                        err_loc_vec,
                        partition_vec,
                        cum_drops,
                        code_len,
                        syndrome,
                        sorted_H_cols,
                        max_query,
                        W_upper_bound,
                        inc,
                        w_init,
                        cte,
                        abandon
                    )
            end

            if zerosyn
                for i in 1:w
                    err_vec[sorted_ind[err_loc_vec[i]]] = true
                end
                candidate .⊻= err_vec
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

        if printtest
            if zerosyn
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