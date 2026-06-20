N = 5
# Sizes = zeros(Int,N)

W = N*(N+1)÷2
u = zeros(Int,N)
d = zeros(Int,N)
D = zeros(Int,N)

noise_locations = zeros(Int,N,W)

for wt = 1:W
    w = max(1,ceil(Int,(1+2*(N)-sqrt((1+2*(N))^2-8*wt))/2))
    while w ≤ N
        if wt < w*(w+1)/2
            break
        end
        display("wt = $wt, w = $w")
        len = landslide!(noise_locations,u,wt,w,N,d,D) 
        display(noise_locations[1:w,1:len])
        # Sizes[w] = len
        w += 1
    end
end

function landslide_2!(
    u::Vector{Int},
    W::Int,
    w::Int,
    code_len::Int,
    d::Vector{Int},
    D::Vector{Int}
)

    W1 = W - (w*(w+1))÷2
    n1 = code_len - w 
    # Create the first integer partition
    j = 1
    # Start with empty vector and breaking at first index
    u .= 0
    k = 1
    mountain_build!(u,k,w,W1,n1)
    # Evaluate drops
    # d = circshift(u,-1) - u
    for i in 1:w-1
        d[i] = u[i+1] - u[i]
    end
    d[w] = 0
    # Evaluate accumulated drops
    # D = cumsum(d[end:-1:1])
    D[w] = d[w]
    for i in w-1:-1:1
        D[i] = d[i] + D[i+1]
    end
    # display("j=$j")
    # display([u[1:w] d[1:w] D[1:w]])
    # Each loop generates a new integer partition
    while D[1] ≥ 2
        # Find the last index with an accumulated drop >=2
        k = findlast(x -> x ≥ 2, D)
        # Increase its index by one.
        u[k] = u[k] + 1
        mountain_build!(u,k,w,W1,n1)
        # Record the partition
        j += 1
        # Evaluate drops
        for i in 1:w-1
            d[i] = u[i+1] - u[i]
        end
        d[w] = 0
        # Evaluate acumuated drops
        D[w] = d[w]
        for i in w-1:-1:1
            D[i] = d[i] + D[i+1]
        end
        # display("j=$j")
        # display([u[1:w] d[1:w] D[1:w]])
    end

    return j
end