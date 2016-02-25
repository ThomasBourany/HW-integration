
module HW_int

	# question 1 b) 
	# here are the packages I used

	using FastGaussQuadrature
	using Roots
	using Sobol
	using PyPlot
	using Distributions

	# here are some functions I defined for useage 
	# in several sub questions

	# demand function

p = linspace(0.2, 5, 1000)
q(p)=2p.^(-0.5)
q1 = linspace(q(1), q(1), 1000)
q2 = linspace(q(4), q(4), 1000)

plot(p, q(p))
plot(p, q1)
plot(p, q2)
title("Inverse demand function")
xlabel("Prices")
ylabel("Quantities")


	# gauss-legendre adjustment factors for map change

	# eqm condition for question 2
	# this is the equilibrium condition: total demand = supply, 


	# weighted sum for integration from the slides.



	function question_1b(n)

using FastGaussQuadrature 
q(p)=2p^(-0.5)
nodes, weights = gausslegendre(n)
x = []
nods = []
points = []
for i = 1:n
push!(x, weights[i]q(1.5nodes[i] + 2.5))
push!(nods, 1.5nodes[i]+2.5)
push!(points, q(1.5nodes[i]+2.5))
end
approx_legendre = 1.5sum(x)
println("Result of the Gauss-Legendre approximation is $approx_legendre")

using PyPlot
p = linspace(0.5, 5, 1000)
quant(p)=2p.^(-0.5)
plot(p, quant(p))
scatter(nods, points, s=10, color = "green")
title("Inverse demand function, with the nodes of the Gauss-Legendre Approximation")
xlabel("Prices")
ylabel("Quantities")
	end


	function question_1c(n)


using FastGaussQuadrature 
using PyPlot
q(p)=2p^(-0.5)
x = []
for i = 1:n
    nodes = 3*rand() + 1
    push!(x, (1/n)q(nodes))
end
approx_montecarlo = 3sum(x)
println("Result of the Monte-Carlo approximation is $approx_montecarlo")

p = linspace(0.5, 5, 1000)
quant(p)=2p.^(-0.5)
nodes = 3*rand(n)+1
points = quant(nodes)
plot(p, quant(p))
scatter(nodes, points, s=10, color = "purple")
title("Inverse demand function, with the nodes of the Monte-Carlo Approximation")
xlabel("Prices")
ylabel("Quantities")

	end


	function question_1d(n)

using PyPlot
using Sobol
s = SobolSeq(1, 1, 4)
q(p)=2p^(-0.5)
points = []
nods = []
for i = 1:n 
	x = next(s)[1]
	nods = push!(nods, x)
	push!(points, q(x))
end 
approx_quasimontecarlo = 3(1/n)sum(points)
println("Result of the Quasi-Monte-Carlo approximation is $approx_quasimontecarlo")


p = linspace(0.5, 5, 1000)
quant(p)=2p.^(-0.5)
plot(p, quant(p))
scatter(nods, points, s=10, color = "red")
title("Inverse demand function, with the nodes of the Quasi-Monte-Carlo Approximation")
xlabel("Prices")
ylabel("Quantities")

	end


	function question_2a(n)

function dem(p,t1,t2)
	exp(t1)/p .+ exp(t2).* p^-0.5 - 2
end

using FastGaussQuadrature

nodes, weights = gausshermite(n)

sigma = hcat([0.02, 0.01],[0.01,0.01])
omega = chol(sigma,Val{:U})
mu = [0.0;0.0]

nod = hcat(kron(ones(n), nodes),kron(nodes, ones(n)))
# create a grid with the nodes for the Gauss-Hermite approximation in two dimensions
# two columns correspondings to the the cartesian products of the set of nodes with itself

grid = omega* transpose(nod) + zeros(2, n*n)

wgth = kron(weights, weights) .*(1/pi)

y = []
y2 = []
for i=1:n*n
	function g(p)
		dem(p, grid[1, i], grid[2,i]) 
	end
	popt = fzero(p -> g(p), [0.001,1000])
	push!(y, popt)
	push!(y2, popt^2)
end

EP = transpose(y)*wgth
VAR = transpose(y2)*wgth - EP.^2



	end


	function question_2b(n)

function dem(p,t1,t2)
	exp(t1)/p .+ exp(t2).* p^-0.5 - 2
end

nodes = rand(n)
sigma = hcat([0.02, 0.01],[0.01,0.01])
omega = chol(sigma,Val{:U})

nod = hcat(kron(ones(n), nodes),kron(nodes, ones(n)))
grid = omega* transpose(nod) + zeros(2, n*n)

wgth = ones(n*n, 1).* 1/(n*n)

y = []
y2 = []
for i=1:n*n
	function g(p)
		dem(p, grid[1, i], grid[2,i]) 
	end
	popt = fzero(p -> g(p), [0.001,1000])
	push!(y, popt)
	push!(y2, popt^2)
end

EP = transpose(y)*wgth
VAR = transpose(y2)*wgth - EP.^2

	end	


	# function to run all questions
	function runall(n=10)
		println("running all questions of HW-integration:")
		println("results of question 1:")
		question_1b(n)	# make sure your function prints some kind of result!
		question_1c(n)
		question_1d(n)
		println("")
		println("results of question 2:")
		q2 = question_2a(n)
		println(q2)
		q2b = question_2b(n)
		println(q2b)
		println("end of HW-integration")
	end

end