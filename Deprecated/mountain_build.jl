function mountain_build(
    partition_vec::Vector{Int},
    k::Int,
    w::Int,
    W_prime::Int,
    n_prime::Int
)

    partition_k = partition_vec[k]
    for i in k+1:w
        partition_vec[i] = partition_k
    end
    sum_partition = 0
    for i in 1:w
        sum_partition += partition_vec[i]
    end
    dividend = W_prime - sum_partition
    divisor = n_prime - partition_k
    if dividend > divisor
        q = div(dividend, divisor)
        w_m_q = w - q
        if w_m_q > 0
            # remainder of the division
            r = dividend - q*divisor
            partition_vec[w_m_q] += r
        end
        for i in (w_m_q+1):w
            partition_vec[i] = n_prime
        end
    else
        # w_m_q = 0, q = 0 and r = dividend
        partition_vec[w] += dividend
    end
end