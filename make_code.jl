

function hex_to_bitvector(koopman_poly_hex::Integer)
    
    # Determine number of bits (4 bits per hex digit)
    num_bits = length(string(koopman_poly_hex,base=16)) * 4
    
    # Extract bits from most significant to least significant
    bitvector = zeros(Bool,num_bits)
    for i in 1:num_bits
        bitvector[i] = (koopman_poly_hex >> (num_bits - i)) & 1
    end

    index = 1
    for i in eachindex(bitvector)
        if bitvector[i]
            break
        end
        index += 1
    end
    
    return bitvector[index:end]
end


function koopman_poly(
    N::Int, 
    K::Int
)

    hd, koopman_poly_hex = koopman_CRC_hex(N,K)
    

    return hd, koopman_poly_hex, [hex_to_bitvector(koopman_poly_hex); true]

end

function vizualize_poly(bitvector::Vector{Bool})

    L = length(bitvector)

    i = 0
    while i < L
        i += 1
        if bitvector[i]
            print("x^$(L-i)")
            break
        end
    end
    while i < L-2
        i += 1
        if bitvector[i]
            print(" + x^$(L-i)")
        end
    end
    if bitvector[L-1]
        print(" + x")
    end
    if bitvector[L]
        print(" + 1")
    end

end

function vizualize_poly(args::Tuple{Int64,UInt16,Vector{Bool}})

    hd, koopman_poly_hex, bitvector = args

    println()
    print("Hamming distance = $hd, Koopman poly hexadecimal = ")
    display(koopman_poly_hex)
    println()
    vizualize_poly(bitvector)
    println()

end

function CRC_code(
    N::Int,
    K::Int
)

    hd, koopman_poly_hex, crc_poly = koopman_poly(N,K)

    G = zeros(Bool,N,K)
    u = [zeros(Bool,K-1); true]
    uu = zeros(Bool,N)
    for i =1:K
        uu[1:K] = circshift(u,i)
        _,r = divide_poly(uu,crc_poly) 
        G[:,i] = [uu[1:K];r]
    end
    H = [G[K+1:end,:] I(N-K)]

    return G, H, crc_poly, hd, koopman_poly_hex

end