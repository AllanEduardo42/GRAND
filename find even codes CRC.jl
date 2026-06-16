even_codes = []

for m in 3:32
    for k in 1:(64 - m)
        n = m + k
        PP, HH, CRC_POLY, HD, KOOPMAN_POLY_HEX = CRC_code(n,k)
        if iszero(iseven.(sum(PP,dims=1)))
            push!(even_codes,(k,n))
        end
    end
end

