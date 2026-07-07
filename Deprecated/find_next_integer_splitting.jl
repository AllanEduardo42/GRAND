include("integer_splitting.jl")

function find_next_valid_integer_splitting!(
    split::Vector{Int},
    J::Vector{Int},
    anchors::Vector{Int},
    betas::Vector{Int},
    W::Int,
    sum_W::Int,
    i::Int,
    num_segments::Int,
    weights::Matrix{Int}
)   

    i += 1
    W_init = split[i]
    for Wi in W_init:(W-sum_W)
        if Wi == 0 || valid_partial_hamming_weights!(Wi,J[i],anchors[i],anchors[i+1],betas[i],weights,i)
            if i == num_segments - 1
                j = i + 1
                Wi_p1 = W - (Wi + sum_W)
                if Wi_p1 == 0 || valid_partial_hamming_weights!(Wi_p1,J[j],anchors[j],anchors[j+1],betas[j],weights,j)
                    split[i] = Wi
                    split[j] = Wi_p1
                    return true
                end
                weights[num_segments,:] .= 0
            else
                split[i] = Wi
                success = find_next_valid_integer_splitting!(split,J,anchors,betas,W,sum_W + Wi,i,num_segments,weights)
                if success
                    return true
                end
            end
            weights[i,:] .= 0
        end
    end
    split[i] = 0

    return false

end