################################################################################
# Allan Eduardo Feitosa
# 07 Jul 2026
# Function to find the anchors, offsets and slopes of the piece-wise linear
# statistical model of the reliabilities for full ORBGRAND algorithm

# OUTPUTS:
# anchors               :Integer anchors of the line segments
# offsets               :Integer offsets of the line segments
# slopes                :Integer slopes of the line segments

# AUXILIARIES:
# alphas                :intermediary Float slopes

# INPUTS:
# sorted_abs_signal     :sorted absolute values of the received signal
# code_len              :codeword length (i.e., N)
# num_segments          :number of linear segments



function line_segmentation!(
    anchors::Vector{Int},
    offsets::Vector{Int},    
    slopes::Vector{Int},
    alphas::Vector{Float64},
    sorted_abs_signal::Vector{Float64},
    code_len::Int,
    num_segments::Int
)
    
    @inbounds @fastmath begin

        # 1) Find anchor points

        end_point = code_len÷2
        anchors[1] = 1
        anchors[end] = end_point

        L = sorted_abs_signal[1]

        for i in 2:num_segments
            m = (sorted_abs_signal[end_point]-L)/(end_point - 1)
            max_dif = 0
            idx = 1
            for j in 1:end_point
                dif = abs(sorted_abs_signal[j] - (L + (j - 1)*m))
                if dif > max_dif
                    max_dif = dif
                    idx = j
                end
            end
            anchors[i] = idx
            end_point = idx
        end
        sort!(anchors)

        # 2) Find minimum slope
        min_slope = sorted_abs_signal[end]
        for i = 1:num_segments
            x2 = anchors[i+1]
            x1 = anchors[i]
            alphas[i] = (sorted_abs_signal[x2] - sorted_abs_signal[x1])/(x2-x1)
            if alphas[i] < min_slope
                min_slope = alphas[i]
            end
        end

        # 3) Quantized (integer) slopes and offsets
        for i in 1:num_segments
            slopes[i] = round(Int,alphas[i]/min_slope)
            offsets[i] = round(Int,sorted_abs_signal[anchors[i]]/min_slope)
            if i == 1
                offsets[i] -= slopes[i]
            end
        end

        anchors[1] = 0
        # extention of the last segment
        anchors[end] = code_len
    end

    return min_slope

end