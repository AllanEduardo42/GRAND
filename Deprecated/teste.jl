# N = 10
# err_loc_vec = zeros(Int,N)
# err_loc_vec_2 = zeros(Int,N)
# err_loc_vec_len = 1
# err_loc_vec_len_2 = 1
# inc = 1

# include("hard_grand.jl")

# err_loc_vec[1] = 1
# err_loc_vec_2[1] = 1

# iters = 1022

# @time for i = 1:iters

#     global err_loc_vec_len = increase_error!(err_loc_vec,err_loc_vec_len,inc,N)

# end

# @time for i = 1:iters

#     global err_loc_vec_len_2 = increase_error_2!(err_loc_vec_2,err_loc_vec_len_2,inc,N)

# end

# display([err_loc_vec err_loc_vec_2])

# display(err_loc_vec == err_loc_vec_2)

include("hard_grand.jl")

N = 10
nthreads = 12
err_vec = zeros(Bool,N,nthreads)
err_loc_vec = zeros(Int,N,nthreads)
err_loc_vec_len = zeros(Int,nthreads)
inc = 1

iters = 5

for thread = 1:nthreads

    for i = 1:iters

        global err_loc_vec_len[thread] = gen_next_err_par!(err_vec,err_loc_vec,err_loc_vec_len[thread],inc,N,thread,nthreads)

    end

end

display(err_loc_vec)