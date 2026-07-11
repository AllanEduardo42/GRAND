################################################################################
# Allan Eduardo Feitosa
# 07 Jul 2026
# Function to find the valid partial hamming weights for full ORBGRAND

# OUTPUTS:
# valid             :if the given logistic weight W is valid
# range             :range of the valid partial hamming weights

# INPUTS:
# W                 :value of the logistic weight
# offset            :line intercept of the piece-wise linear statistical model
# size              :size of the linear segment
# slope             :slope of the piece-wise linear statistical model

# Conditions for valid logistic weight W:
#   i) W - offset*w ≥ slope*w(w+1)/2   
#  ii) W - offset*w ≤ slope*(delta + 1)*w - slope*w(w-1)/2
# iii) rem(W - offset*w,slope) = 0
# Obs: W = 0 is always valid

function find_valid_hamming_weights!(
    W::Int,
    offset::Int,
    size::Int,
    slope::Int
)

    @inbounds @fastmath begin
        # default outputs
        valid = false
        range = 0:1:0
        if offset == 0 && rem(W,slope) != 0
            # if the line intercept is zero, condition
            return valid, range
        end
        # The range of valid partial hamming weights is given by
        #
        #       w_min ≤ w ≤ w_max
        #s
        # where w_min = max(0,ceil(Int,A)), where
        #              _________________
        #        b₂ - √(b₂)² - 8*slope*W
        # A =    ----------------------- , 
        #                2*slope
        #
        # where b₂ = b₁ + 2*slope*size, b₁ = 2*offset + slope, and
        #
        #       w_max = min(floor(Int,B),floor(Int,C)),
        #
        # where
        #              _________________
        #        b₂ + √(b₂)² - 8*slope*W
        # B =    ----------------------- , and
        #                2*slope
        #               _________________
        #        -b₁ + √(b₁)² + 8*slope*W
        # C =    ------------------------ ,
        #                2*slope
        
        two_slopes = 2*slope
        b1 = 2*offset + slope        
        b2 = b1 + two_slopes*size
        c = 8*slope*W
        d = b2^2 - c
        if d ≥ 0    # there is no range of valid hamming weights otherwise
            sq_d = sqrt(d)
            if isinteger(sq_d)
                sq_d = Int(sq_d)        # to avoid round errors (0.99999999... ≠ 1)
            end
            x = b2 - sq_d
            if x > 0
                w_min = unsafe_trunc(Int,ceil(x/two_slopes))
            else
                w_min = 0
            end
            sq_e = sqrt(b1^2 + c)
            if isinteger(sq_e)
                sq_e = Int(sq_e)        # to avoid round errors (0.99999999... ≠ 1)
            end
            w_max =  unsafe_trunc(Int,min(floor((b2 + sq_d)/two_slopes),
                                          floor((sq_e - b1)/two_slopes)))
            if w_min ≤ w_max    # verify if the range is valid
                # next, find the mininum w valid inside the range w_min:w_max
                w = w_min               # start at w_min
                wL = W - w*offset
                r = rem(wL,slope)       # condition iii)
                success = true          # flag indicating that the minimum had 
                                        # been found
                # if r == 0, w = w_min is a valid hamming weight
                while r ≠ 0
                    if w == w_max
                        # i.e., no valid w was found in the range w_min:w_max
                        success = false
                        break
                    end
                    # search for w_min+1,w_min+2,...,w_max, until r != 0
                    w += 1
                    wL = W - w*offset
                    r = rem(wL,slope)
                end              
                if success
                    # if we found a valid w inside w_min:w_max

                    # Next, we find the step delta_w of valid hamming weights,
                    # i.e., w:delta_w:w_max

                    # delta_w must obey the relation (if offset ≠ 0)
                    # delta_w      slope
                    # ------- = - ------- ,
                    # delta_k      offset
                    # where delta_k = k₁ - k₂, such that
                    #
                    # W - w*offset = k₁*slope, 
                    #
                    # and
                    #
                    # W - (w + delta_w)*offset = k₂*slope, 
                    #
                    # with k₁,k₂ ∈ Z.

                    if offset == 0
                        delta_w = 1     # all w in the range are valid
                    elseif offset == 1
                        delta_w = slope
                    else
                        g = gcd(slope,offset)
                        if g == 1
                            delta_w = slope
                        else
                            delta_w = div(slope,g)
                        end
                    end
                    # updated output
                    range = w:delta_w:w_max
                    valid = true
                    return valid, range
                end
            end
        end 
    end

    if W > 0
        return valid, range
    else
        valid = true
        return valid, range
    end

end