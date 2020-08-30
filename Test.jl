include("BesselYzero.jl")

using Plots

ν_p = 0
xlim = 12
xs = range(0.00001,stop = xlim,length=500);
ys = bessely.(ν_p,xs)

ms = 1:4
roots_p = initial_guess_yν.(ν_p, ms)

plot(xs,ys,label ="")
plot!([0,xlim],[0,0],label="")
scatter!(roots_p,zeros(length(roots_p)),label="Guess")
xlabel!("x")
ylabel!("Y$ν_p")



indices = rand(100,100) * (1000-0) .+ 0;
indices_ra = Int.(floor.(rand(100,100) * (10-0) .+ 0));

pistas = [initial_guess_yν(indices[i,j],indices_ra[i,j]) for i in 1:100, j in 1:100]

limsup = [try 
    bessely(indices[i,j],real(pistas[i,j]+0.001pistas[i,j])) 
catch 
    0 
end 
for i in 1:100, j in 1:100]

    limsin = [try 
    bessely(indices[i,j],real(pistas[i,j]-0.001pistas[i,j])) 
catch 
    0 
end 
for i in 1:100, j in 1:100]

signos=    sign.(limsin .* limsup)

#using PyPlot

mes = heatmap(signos)