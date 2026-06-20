################################################################################
# Allan Eduardo Feitosa
# 15 Jun 2026
# Auxiliary functions

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

function fast_gf2_mat_mul(
    int_vec::Vector{Int},
    bool_vec::Vector{Bool}
)

    if length(int_vec) != length(bool_vec)
        throw(error(lazy"dimensions must be equal"))
    end

    result = 0
    for i in eachindex(int_vec)
        if bool_vec[i]
            result ⊻= int_vec[i]
        end
    end

    return result

end

function vector_diff!(
    d::Vector{Int},
    u::Vector{Int}
)

    b = u[1]
    for i in 1:w-1
        u_1 = u_2
        u_2 = u[i+1]
        d[i] = u_2 - u_1
    end

end