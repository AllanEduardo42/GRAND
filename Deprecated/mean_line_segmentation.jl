function mean_line_segmentation(
    stdev::Float64,
    num_segments::Int,
    code_len::Int,
    num_trials::Int
)

    x = zeros(code_len)
    y = zeros(code_len)
    y_abs_sorted = zeros(code_len,num_trials)
    anchors = zeros(Int,num_segments+1)
    J = zeros(Int,num_segments)
    alphas = zeros(num_segments)
    betas = zeros(Int,num_segments)

    for j in 1:num_trials
        x .= rand(Bool,code_len)
        y .= (2*x .- 1) + stdev*randn(code_len)
        y_abs_sorted[:,j] .= sort(abs.(y))
    end

    y_mean = mean(y_abs_sorted,dims=2)[:]

    end_point = code_len÷2
    anchors[1] = 1
    anchors[end] = end_point

    L1 = y_mean[1]

    for i in 2:num_segments
        m = (y_mean[end_point]-L1)/(end_point - 1)
        max = 0
        idx = 1
        for j in 1:end_point
            dif = y_mean[j] - (L1 + (j - 1)*m)
            if dif > max
                max = dif
                idx = j
            end
        end
        anchors[i] = idx
        end_point = idx
    end
    sort!(anchors)

    min_beta = Inf
    for i = 1:num_segments
        alphas[i] = (y_mean[anchors[i+1]] - y_mean[anchors[i]])/(anchors[i+1]-anchors[i])
        if alphas[i] < min_beta
            min_beta = alphas[i]
        end
    end

    for i in 1:num_segments
        betas[i] = round(Int,alphas[i]/min_beta)
        J[i] = round(Int,y_mean[anchors[i]]/min_beta)
    end

    return anchors, J, betas, y_mean, min_beta
    
end