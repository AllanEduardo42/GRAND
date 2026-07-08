using DelimitedFiles

EbN0 = 3.0:0.5:8.0

L = 100_000

x = zeros(CODE_LEN)
y = zeros(CODE_LEN)
y_abs_sorted = zeros(CODE_LEN,L)
y_mean = zeros(CODE_LEN,length(EbN0))


p = plot(
    legend_title="EbN0 (dB)",
    legendtitlefontsize = 8,
    title="Ordered Statistical model for |Yᵢ|",
    xlabel = "rank order",
    legend = :outertopright
)

k = 0 
for ebn0 in EbN0
    global k += 1

    # transform EbN0 in standard deviations
    variance = exp10.(-ebn0/10) / (2*RR)
    stdev = sqrt.(variance)    

    for j in 1:L
        x .= rand(Bool,CODE_LEN)
        y .= (2*X .- 1) + stdev*randn(CODE_LEN)
        y_abs_sorted[:,j] .= sort(abs.(y))
    end

    y_mean[:,k] = mean(y_abs_sorted,dims=2)
    plot!(p,1:CODE_LEN,y_mean[:,k],label="$ebn0")
end

display(p)

open("$CODE_LEN.txt","w") do io
    writedlm(io,[collect(EbN0)'; y_mean])
end