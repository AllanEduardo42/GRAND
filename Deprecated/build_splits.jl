function build_splits!(
    W_splits::Matrix{Int},
    W::Int,
    num_splits::Int,
    sum_W::Int,
    i::Int,
    num_segments::Int,
    split::Vector{Int},
    W_list::Matrix{Int},
    W_lens::Vector{Int},
    lim_inf::Int,
    Weights_list::Array{Int,3},
    Hamming_weights::Array{Int,3},
    Weights::Matrix{Int}
)
    i += 1  # i-th segment
    for k in 1:W_lens[i]  
        Wi = W_list[k,i]
        if Wi > W - sum_W
            return num_splits
        end
        if i == num_segments - 1
            j = i + 1
            Wi_p1 = W - (Wi + sum_W)
            if Wi_p1 == 0 || Wi_p1 ≥ lim_inf
                num_splits += 1
                W_splits[num_splits,i] = Wi
                Hamming_weights[num_splits,i,:] = Weights_list[k,i,:]
                W_splits[num_splits,j] = Wi_p1
                if Wi_p1 == 0
                    kk = 1
                else
                    kk = Wi_p1 - lim_inf + 2
                end
                Hamming_weights[num_splits,j,:] = Weights_list[kk,j,:]
                for ii in 1:i-1
                    W_splits[num_splits,ii] = split[ii]
                    Hamming_weights[num_splits,ii,:] = Weights[ii,:]
                end
            end
        else
            split[i] = Wi
            Weights[i,:] = Weights_list[k,i,:]
            num_splits = build_splits!(W_splits,W,num_splits,sum_W + Wi,i,num_segments,split,W_list,W_lens,lim_inf,Weights_list,Hamming_weights,Weights)
        end
    end

    return num_splits
end

function build_splits!(
    W_splits::Matrix{Int},
    W::Int,
    num_segments::Int,
    split::Vector{Int},
    W_list::Matrix{Int},
    W_lens::Vector{Int},
    lim_inf::Int,
    Weights_list::Array{Int,3},
    Hamming_weights::Array{Int,3},
    Weights::Matrix{Int}
)
    return build_splits!(W_splits,W,0,0,0,num_segments,split,W_list,W_lens,lim_inf,Weights_list,Hamming_weights,Weights)
end
