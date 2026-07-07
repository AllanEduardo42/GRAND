function old_valid_partial_hamming_weights!(
    Wi::Int,
    Ji::Int,
    Ii::Int,
    Ii_p1::Int,
    beta_i::Int,
    Weights::Matrix{Int},
    i::Int
)
    num_w = 0
    lim_sup = floor(Int,(sqrt(1 + 8*Wi)-1)/2)
    for w in 1:lim_sup
        a = w*(w+1)/2
        b = Wi - w*Ji
        c = (Ii_p1 - Ii + 1)*w - a
        if b ≥ a
            if b ≤ c 
                if rem(b,beta_i) == 0
                    num_w += 1
                    # if i == 1
                    #     display("i = $i, Wi = $Wi, w = $w")
                    # elseif i == 2
                    #     display("    i = $i, Wi = $Wi, w = $w")
                    # elseif i == 3
                    #     display("        i = $i, Wi = $Wi, w = $w")
                    # end
                    Weights[i,num_w] = w
                end
            end
        end
    end
    if num_w == 0
        # if i == 1
        #     display("i = $i, Wi = $Wi failed")
        # elseif i == 2
        #     display("    i = $i, Wi = $Wi failed")
        # elseif i == 3
        #     display("        i = $i, Wi = $Wi failed")
        # end
        return false # fail
    else
        return true
    end
end

function old_integer_splitting!(
    W_splits::Matrix{Int},
    J::Vector{Int},
    anchors::Vector{Int},
    betas::Vector{Int},
    W::Int,
    sum_W::Int,
    num_splits::Int,
    i::Int,
    num_segments::Int,
    split::Vector{Int},
    Hamming_weights::Array{Int,3},
    Weights::Matrix{Int}
)
    i += 1
    for Wi in 0:(W-sum_W)
        if Wi == 0 || old_valid_partial_hamming_weights!(Wi,J[i],anchors[i],anchors[i+1],betas[i],Weights,i)
            if i == num_segments - 1
                j = i + 1
                Wi_p1 = W - (Wi + sum_W)
                if Wi_p1 == 0 || old_valid_partial_hamming_weights!(Wi_p1,J[j],anchors[j],anchors[j+1],betas[j],Weights,j)
                    num_splits += 1
                    W_splits[num_splits,i] = Wi
                    Hamming_weights[num_splits,i,:] .= Weights[i,:]
                    W_splits[num_splits,j] = Wi_p1
                    Hamming_weights[num_splits,j,:] .= Weights[j,:]
                    for j in 1:i-1
                        W_splits[num_splits,j] = split[j]
                        Hamming_weights[num_splits,j,:] .= Weights[j,:]
                    end
                end
                Weights[num_segments,:] .= 0
            else
                split[i] = Wi
                num_splits = old_integer_splitting!(W_splits,J,anchors,betas,W,sum_W + Wi,num_splits,i,num_segments,split,Hamming_weights,Weights)
            end
            Weights[i,:] .= 0
        end
    end

    return num_splits
end

function old_integer_splitting!(
    W_splits::Matrix{Int},
    J::Vector{Int},
    anchors::Vector{Int},
    betas::Vector{Int},
    W::Int,
    num_segments::Int,
    split::Vector{Int},
    Hamming_weights::Array{Int,3},
    Weights::Matrix{Int}
)

    return old_integer_splitting!(W_splits,J,anchors,betas,W,0,0,0,num_segments,split,Hamming_weights,Weights)

end