include("valid_partial_hamming_weights.jl")

function integer_splitting!(
    W_splits::Matrix{Int},
    values::Vector{Int},
    anchors::Vector{Int},
    betas::Vector{Int},
    W::Int,
    num_splits::Int,
    sum_W::Int,
    i::Int,
    num_segments::Int,
    Hamming_weights::Array{Int,3},
    split::Vector{Int},
    Weights::Matrix{Int}
)
    i += 1  # i-th segment
    for Wi in 0:(W-sum_W)
        if Wi == 0 || valid_partial_hamming_weights!(Wi,values[i],anchors[i],anchors[i+1],betas[i],Weights,i)
            if i == num_segments - 1
                j = i + 1
                Wi_p1 = W - (Wi + sum_W)
                if Wi_p1 == 0 || valid_partial_hamming_weights!(Wi_p1,values[j],anchors[j],anchors[j+1],betas[j],Weights,j)
                    num_splits += 1
                    W_splits[num_splits,i] = Wi
                    Hamming_weights[num_splits,i,:] = Weights[i,:]
                    W_splits[num_splits,j] = Wi_p1
                    Hamming_weights[num_splits,j,:] = Weights[j,:]
                    for ii in 1:i-1
                        W_splits[num_splits,ii] = split[ii]
                        Hamming_weights[num_splits,ii,:] = Weights[ii,:]
                    end
                end
                Weights[num_segments,:] .= 0
            else
                split[i] = Wi
                num_splits = integer_splitting!(W_splits,values,anchors,betas,W,num_splits,sum_W + Wi,i,num_segments,Hamming_weights,split,Weights)
            end
            Weights[i,:] .= 0
        end
    end

    return num_splits
end

function integer_splitting!(
    W_splits::Matrix{Int},
    values::Vector{Int},
    anchors::Vector{Int},
    betas::Vector{Int},
    W::Int,
    num_segments::Int,
    Hamming_weights::Array{Int,3},
    split::Vector{Int},
    Weights::Matrix{Int}
)

    return integer_splitting!(W_splits,values,anchors,betas,W,0,0,0,num_segments,Hamming_weights,split,Weights)

end