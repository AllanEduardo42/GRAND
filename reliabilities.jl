using Plots
using Statistics

plotlyjs()

N = 10000

p = plot(legend = :topleft)
for snr=0:2:16
    variance = exp10(-snr/10)
    x = rand(Bool,N)
    y = (2*x .- 1) + sqrt(variance)*randn(N)
    plot!(p,sort(abs.(y)),label="SN = $snr")
end
display(p)
plot!(1:N,2*(1:N)/N, ls=:dash, lw = 2, color = :black)