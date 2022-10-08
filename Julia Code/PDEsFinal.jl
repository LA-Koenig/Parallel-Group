### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ 37d40cb2-c581-11ec-14d8-6be0f0c127a7
using SparseArrays, LinearAlgebra, BenchmarkTools, GLMakie, PlutoUI, PrettyTables, PlutoProfile

# ╔═╡ 822d116c-f1b6-42a5-990b-e690d4eaea99
md"""
# Final Exam
$$\frac{\delta u}{\delta t} = \frac{\delta^2 u}{\delta x^2} + \frac{\delta^2 u}{\delta y^2} + \frac{\delta^2 u}{\delta z^2}, t\gt0,0\lt x,y,z\lt 1$$
subject to ``u(x,y,z,0) = \sin (\pi y) \sin (\pi z), 0 \le x, y, z\le 1``; and given homogeneous Dirichlet boundary conditions ``(u\equiv0)`` on the boundary ``\delta \Omega`` of the cube.
The Fourier solution of this problem is easily found to be
$$u(x,y,z,t) = \sin (\pi y)\sin (\pi z)\sum_{n=1}^\inf [1 - (-1)^n]\frac{2}{n\pi}\sin (n\pi x)\exp(-\pi^2(n^2+2)t)$$
"""

# ╔═╡ fde5335d-c89b-4726-9302-8947795238df
md"""## Solvers"""

# ╔═╡ 9af21332-3679-4997-98ca-05911d648654
function tri_solve!(M,dprime,cprime,d)
	nd = length(d)
	
	cprime[1] = M.ev[1]/M.dv[1]
	dprime[1] = d[1]/M.dv[1]
	for i in 2:nd-1
		cprime[i] = M.ev[i]/(M.dv[i] - M.ev[i-1]*cprime[i-1])
		dprime[i] = (d[i] - M.ev[i-1]*dprime[i-1])/(M.dv[i] - M.ev[i-1]*cprime[i-1])
	end
	d[end] = (d[end] - M.ev[end-1]*dprime[end-1])/(M.dv[end] - M.ev[end-1]*cprime[end-1])

	for j in nd-1:-1:1
		d[j] = (dprime[j] - cprime[j]*d[j+1])
	end
end	

# ╔═╡ 6626a106-3066-4924-a1c0-2c708a38ed88
function back_diff_tri2(space_dims, t, init)
	numdims = length(space_dims)
	dim_lengths = [length(dim)-2 for dim in space_dims]
	max_dim_length = maximum(dim_lengths)
	lvU = reduce(*,dim_lengths)
	permutevec = [collect(2:numdims); 1]
	# Compute r
	h = Float64(space_dims[1].step)
	δt = Float64(t.step)
	r = δt/h^2

	# Initialize U
	U = zeros(Float64, dim_lengths...)
	for idx in CartesianIndices(U)
		coords = [space_dims[i][idx[i]+1] for i in 1:numdims]
		U[idx] = init(coords...)
	end
	Us = [copy(U) for i in 1:numdims]

	# Construct the matrices used in each dimension
	A = SymTridiagonal(fill(-2,max_dim_length),
					   fill(1,max_dim_length-1))
	M = (I - r*A)
	dprime = zeros(Float64, maximum(dim_lengths))
	cprime = zeros(Float64, maximum(dim_lengths)-1)

	# Iterate over each timestep
	for step in 2:length(t)
		for dim_idx in 1:numdims
			length_dim = dim_lengths[dim_idx]
			vU = vec(Us[dim_idx])
			for idx in 1:length_dim:lvU
				tri_solve!(M, dprime, cprime, view(vU, idx:idx+length_dim-1))
			end
			permutedims!(Us[(dim_idx%numdims)+1], Us[dim_idx], permutevec)
		end
	end
	Us[1]
end

# ╔═╡ 2d01c0ea-5889-4d83-806d-9d62b669930e
function extrapolation(space_dims, t, init)
	numdims = length(space_dims)
	dim_lengths = [length(dim)-2 for dim in space_dims]
	max_dim_length = maximum(dim_lengths)
	lvU = reduce(*,dim_lengths)
	permutevec = [collect(2:numdims); 1]
	# Compute r
	h = Float64(space_dims[1].step)
	δt = Float64(t.step)
	r = δt/h^2

	# Initialize U
	U = zeros(Float64, dim_lengths...)
	for idx in CartesianIndices(U)
		coords = [space_dims[i][idx[i]+1] for i in 1:numdims]
		U[idx] = init(coords...)
	end
	Us1 = [copy(U) for i in 1:numdims]
	Us2 = [copy(U) for i in 1:numdims]

	A = SymTridiagonal(fill(-2,max_dim_length),
					   fill(1,max_dim_length-1))
	M1 = (I - r*A)
	M2 = (I - 2r*A)
	dprime = zeros(Float64, max_dim_length)
	cprime = zeros(Float64, max_dim_length-1)

	# Iterate over each timestep
	for step in 3:2:length(t)
		# Apply solver twice to U1
		for i in 1:2, dim_idx in 1:numdims
			length_dim = dim_lengths[dim_idx]
			vU1 = vec(Us1[dim_idx])
			for idx in 1:length_dim:lvU
				tri_solve!(M1, dprime, cprime, view(vU1, idx:idx+length_dim-1))
			end
			permutedims!(Us1[(dim_idx%numdims)+1], Us1[dim_idx], permutevec)
		end
		# Apply solver once to U2
		for dim_idx in 1:numdims
			length_dim = dim_lengths[dim_idx]
			vU2 = vec(Us2[dim_idx])
			for idx in 1:length_dim:lvU
				tri_solve!(M2, dprime, cprime, view(vU2, idx:idx+length_dim-1))
				
			end
			permutedims!(Us2[(dim_idx%numdims)+1], Us2[dim_idx], permutevec)
		end
		# combine, apply to both arrays
		@. Us1[1] = Us2[1] = 2*Us1[1] - Us2[1]
	end
	Us1[1]
end

# ╔═╡ 8bbf5ae5-8ca5-497a-8d48-03ad18109fd3
function crank_nicholson(space_dims, t, init)
	numdims = length(space_dims)
	dim_lengths = [length(dim)-2 for dim in space_dims]
	max_dim_length = maximum(dim_lengths)
	lvU = reduce(*,dim_lengths)
	permutevec = [collect(2:numdims); 1]
	# Compute r
	h = Float64(space_dims[1].step)
	δt = Float64(t.step)
	r = δt/h^2

	# Initialize U
	U = zeros(Float64, dim_lengths...)
	for idx in CartesianIndices(U)
		coords = [space_dims[i][idx[i]+1] for i in 1:numdims]
		U[idx] = init(coords...)
	end
	Us = [copy(U) for i in 1:numdims]

	A = SymTridiagonal(fill(-2,max_dim_length),
					   fill(1,max_dim_length-1))
	Mi = (I - (1/2)*r*A)
	Me = (I + (1/2)*r*A)
	dprime = zeros(Float64, maximum(dim_lengths))
	cprime = zeros(Float64, maximum(dim_lengths)-1)

	# Iterate over each timestep
	for step in 2:length(t)
		for dim_idx in 1:numdims
			length_dim = dim_lengths[dim_idx]
			vU = vec(Us[dim_idx])
			for idx in 1:length_dim:lvU
				vU[idx:idx+length_dim-1] = view(Me, 1:length_dim, 1:length_dim) * 
										   view(vU,idx:idx+length_dim-1)
				tri_solve!(Mi, dprime, cprime, view(vU, idx:idx+length_dim-1))
			end
			permutedims!(Us[(dim_idx%numdims)+1], Us[dim_idx], permutevec)
		end
	end
	Us[1]
end

# ╔═╡ 2bc264fb-80d3-4e46-b09e-3710b5983de0
md"""## Problem specific functions"""

# ╔═╡ 847c5261-869d-4a63-8831-d033d844b4a9
function genranges()
	δhs = [0.1, 0.05, 0.025]
	δts = [0.1, 0.01]
	return [(0:δh:1, 0:δt:1) for δh ∈ δhs, δt ∈ δts]
end

# ╔═╡ c8bc65e8-6511-4ef8-8478-291ef60e8aff
function init_1(x,y,z)
	return sin(pi*y)*sin(pi*z)
end

# ╔═╡ c77b6f26-dfc9-4d78-92f6-c8f938b00c68
function init_2(x,y,z)
	return 1
end

# ╔═╡ d573929a-dc4e-4f8b-84e1-01b9c04d9618
function back_diff_wrap3((s,t), init)
	return back_diff_tri2([s,s,s], t, init)
end

# ╔═╡ cba93fa1-96da-480c-b325-7c645f3aba79
function extrap_wrap((s,t), init)
	return extrapolation([s,s,s], t, init)
end

# ╔═╡ d6174a89-0591-4cf2-bd64-1a59dcf94538
function cn_wrap((s,t), init)
	return crank_nicholson([s,s,s], t, init)
end

# ╔═╡ b82e1e04-3b06-4c2d-af31-fedc57089326
function analytical((s, t))
	ns = length(s)-2
	U = zeros(ns,ns,ns)

	for cart_idx in CartesianIndices(U)
		x = s[cart_idx[1]+1]
		y = s[cart_idx[2]+1]
		z = s[cart_idx[3]+1]
		total = 0
		for n in 1:50
			total += (1-(-1)^n)*(2/(n*pi))*sin(n*pi*x)*exp(-(pi^2)*((n^2)+2)*t[end])
		end
		U[cart_idx] = total * sin(pi*y)*sin(pi*z)
	end
	return U
end

# ╔═╡ d2d3e950-3369-4b81-aa6a-c45ac98447ea
function l2_err(U, Uan)
	return norm(vec(U-Uan))
end

# ╔═╡ 1776c0f1-36d0-4089-8a16-beb2b7775dbf
function inf_err(U, Uan)
	return norm(vec(U-Uan), Inf)
end

# ╔═╡ 9562ef4e-b33a-4dcd-a630-91fbc8e570c8
function errtable(err_bd, err_ex, err_cn, errtype)
    header = ([errtype; fill("δt = 0.1", 3); fill("δt = 0.01", 3)],
              ["Method"; repeat(["δx = 0.1", "δx = 0.05", "δx = 0.025"], 2)])

    data = [["Backward-Difference, t=1"; "Crank-Nicholson, t=1"; "Extrapolation, t=1"];;
            [reshape([err_bd[i,j] for i ∈ 1:size(err_bd, 1), j ∈ 1:size(err_bd, 2)], (1,:));
             reshape([err_cn[i,j] for i ∈ 1:size(err_cn, 1), j ∈ 1:size(err_cn, 2)], (1,:));
			 reshape([err_ex[i,j] for i ∈ 1:size(err_ex, 1), j ∈ 1:size(err_ex, 2)], (1,:))]]

    return pretty_table(String, data; backend = Val(:html), tf = tf_html_dark, header=header, formatters=ft_printf("%.3e", [2,3,4,5,6,7]))
end

# ╔═╡ 44e30008-4610-4c9f-9e45-3f3f9432908f
md"""## Part 1, initial state is sine function
### Simulations"""

# ╔═╡ 0f054650-f7f7-4d31-a6fb-0bd41c2cd100
const ranges = genranges()

# ╔═╡ c7b56b9f-5d5c-41eb-b5b4-3074e87c7aa8
const Us_bd = back_diff_wrap3.(ranges, init_1)

# ╔═╡ 415eb156-0feb-4107-a0e2-5b7edbec0d7c
const Us_ex = extrap_wrap.(ranges, init_1)

# ╔═╡ 043a6b6d-7980-4364-81aa-e59be848139e
const Us_cn = cn_wrap.(ranges, init_1)

# ╔═╡ 3b072a4b-0b7f-48f0-adca-26acdd1f55db
const Us_an = analytical.(ranges)

# ╔═╡ 9752d555-3168-46a5-8a64-c9c78343d211
md"""### Errors"""

# ╔═╡ abe8be99-a1c8-431f-a9f6-be5cb8b87861
bd_err_l2 = l2_err.(Us_bd, Us_an)

# ╔═╡ 805c3d65-849b-4108-b810-5a7ff9a118b7
bd_err_inf = inf_err.(Us_bd, Us_an)

# ╔═╡ 5c3de8ab-48d9-4409-a12a-484db456246e
ex_err_l2 = l2_err.(Us_ex, Us_an)

# ╔═╡ e3d6659f-6bcb-4791-9c8b-3eb4fd36605c
ex_err_inf = inf_err.(Us_ex, Us_an)

# ╔═╡ e9af1ed6-a43e-4a3b-b647-899dfbb8ffd5
cn_err_l2 = l2_err.(Us_cn, Us_an)

# ╔═╡ aa5e5d9d-3e73-4721-a63e-41eaa8c8a398
cn_err_inf = inf_err.(Us_cn,Us_an)

# ╔═╡ d86a3a46-90f0-4e22-b63a-ee2d735b6ac0
HTML(errtable(bd_err_l2, ex_err_l2, cn_err_l2, "L2 Error"))

# ╔═╡ 90c6a31d-d98e-4f70-9f0e-7a01d56630c2
HTML(errtable(bd_err_inf, ex_err_inf, cn_err_inf, "Inf Error"))

# ╔═╡ 2c4e1064-93cb-4545-be33-fe5d7dbcbe53
md"""### Surface Plots"""

# ╔═╡ 4d7a930b-d89f-40a3-9475-37b1ec8d5f72
zfrac = 1/4

# ╔═╡ 05bc2ec5-9cf4-49b5-926f-891fdf1ed2dc
let
	f1 = Figure(resolution = (1000,1000))
	for h in 1:3, k in 1:2
		zdx = ceil(Int, length(ranges[h,1][1])*zfrac)
		zval = ranges[h,1][1][zdx]
		hval = Float64(ranges[h,1][1].step)
		kval = Float64(ranges[1,k][2].step)
		Axis3(
			f1[h,k],
			title = "h = $(hval), k = $(kval), z = $(zval)",
			xlabel = "",
			ylabel = "",
			zlabel = "",
			elevation = 0.2pi,
			azimuth=0.2pi,
			aspect = (3,3,2))
		to_plot = zeros(fill(length(ranges[h,1][1]),3)...)
		to_plot[2:end-1,2:end-1,2:end-1] = Us_bd[h,k]
		surface!(ranges[h,k][1], ranges[h,k][1], to_plot[:,:,zdx]; shading=false, colormap=:lapaz)
	end
	supertitle = Label(f1[0,:], "Backward Difference", textsize = 30)
	f1
end

# ╔═╡ a3be3ec7-63d5-435b-af9e-048a822ca7ce
let
	f1 = Figure(resolution = (1000,1000))
	for h in 1:3, k in 1:2
		zdx = ceil(Int, length(ranges[h,1][1])*zfrac)
		zval = ranges[h,1][1][zdx]
		hval = Float64(ranges[h,1][1].step)
		kval = Float64(ranges[1,k][2].step)
		Axis3(
			f1[h,k],
			title = "h = $(hval), k = $(kval), z = $(zval)",
			xlabel = "",
			ylabel = "",
			zlabel = "",
			elevation = 0.2pi,
			azimuth=0.2pi,
			aspect = (3,3,2))
		to_plot = zeros(fill(length(ranges[h,1][1]),3)...)
		to_plot[2:end-1,2:end-1,2:end-1] = Us_ex[h,k]
		surface!(ranges[h,k][1], ranges[h,k][1], to_plot[:,:,zdx]; shading=false, colormap=:lapaz)
	end
	supertitle = Label(f1[0,:], "Extrapolation", textsize = 30)
	f1
end

# ╔═╡ 8924d816-8962-45fd-b15e-0b4e49fc5b50
let
	f1 = Figure(resolution = (1000,1000))
	for h in 1:3, k in 1:2
		zdx = ceil(Int, length(ranges[h,1][1])*zfrac)
		zval = ranges[h,1][1][zdx]
		hval = Float64(ranges[h,1][1].step)
		kval = Float64(ranges[1,k][2].step)
		Axis3(
			f1[h,k],
			title = "h = $(hval), k = $(kval), z = $(zval)",
			xlabel = "",
			ylabel = "",
			zlabel = "",
			elevation = 0.2pi,
			azimuth=0.2pi,
			aspect = (3,3,2))
		to_plot = zeros(fill(length(ranges[h,1][1]),3)...)
		to_plot[2:end-1,2:end-1,2:end-1] = Us_cn[h,k]
		surface!(ranges[h,k][1], ranges[h,k][1], to_plot[:,:,zdx]; shading=false, colormap=:lapaz)
	end
	supertitle = Label(f1[0,:], "Crank-Nicholson", textsize = 30)
	f1
end

# ╔═╡ 61a2e7ec-8d5a-4897-aa9a-526679f1936e
md"""## Part 2, initial state is 1
### Simulations"""

# ╔═╡ 99ff79a5-44d2-4885-b367-a403984eb40e
const Us_bd2 = back_diff_wrap3.(ranges, init_2)

# ╔═╡ f1413200-747a-4a64-88b2-89422c318304
const Us_ex2 = extrap_wrap.(ranges, init_2)

# ╔═╡ 8bbd7282-ce9c-4f14-9769-65150352ada7
const Us_cn2 = cn_wrap.(ranges, init_2)

# ╔═╡ ada0bd92-fb03-4e06-b177-f6e885a22436
md"""## Surface Plots"""

# ╔═╡ d3f672e1-ca18-46a7-8d9a-b92cf520984d
let
	f1 = Figure(resolution = (1000,1000))
	for h in 1:3, k in 1:2
		zdx = ceil(Int, length(ranges[h,1][1])*zfrac)
		zval = ranges[h,1][1][zdx]
		hval = Float64(ranges[h,1][1].step)
		kval = Float64(ranges[1,k][2].step)
		Axis3(
			f1[h,k],
			title = "h = $(hval), k = $(kval), z = $(zval)",
			xlabel = "",
			ylabel = "",
			zlabel = "",
			elevation = 0.2pi,
			azimuth=0.2pi,
			aspect = (3,3,2))
		to_plot = zeros(fill(length(ranges[h,1][1]),3)...)
		to_plot[2:end-1,2:end-1,2:end-1] = Us_bd2[h,k]
		surface!(ranges[h,k][1], ranges[h,k][1], to_plot[:,:,zdx]; shading=false, colormap=:lapaz)
	end
	supertitle = Label(f1[0,:], "Backward Difference", textsize = 30)
	f1
end

# ╔═╡ 61cac858-af0a-4dae-a117-e9c046a5fb58
let
	f1 = Figure(resolution = (1000,1000))
	for h in 1:3, k in 1:2
		zdx = ceil(Int, length(ranges[h,1][1])*zfrac)
		zval = ranges[h,1][1][zdx]
		hval = Float64(ranges[h,1][1].step)
		kval = Float64(ranges[1,k][2].step)
		Axis3(
			f1[h,k],
			title = "h = $(hval), k = $(kval), z = $(zval)",
			xlabel = "",
			ylabel = "",
			zlabel = "",
			elevation = 0.2pi,
			azimuth=0.2pi,
			aspect = (3,3,2))
		to_plot = zeros(fill(length(ranges[h,1][1]),3)...)
		to_plot[2:end-1,2:end-1,2:end-1] = Us_ex2[h,k]
		surface!(ranges[h,k][1], ranges[h,k][1], to_plot[:,:,zdx]; shading=false, colormap=:lapaz)
	end
	supertitle = Label(f1[0,:], "Extrapolation", textsize = 30)
	f1
end

# ╔═╡ 1d571340-0eba-4027-a348-fa7c2767b5c0
let
	f1 = Figure(resolution = (1000,1000))
	for h in 1:3, k in 1:2
		zdx = ceil(Int, length(ranges[h,1][1])*zfrac)
		zval = ranges[h,1][1][zdx]
		hval = Float64(ranges[h,1][1].step)
		kval = Float64(ranges[1,k][2].step)
		Axis3(
			f1[h,k],
			title = "h = $(hval), k = $(kval), z = $(zval)",
			xlabel = "",
			ylabel = "",
			zlabel = "",
			elevation = 0.2pi,
			azimuth=0.2pi,
			aspect = (3,3,2))
		to_plot = zeros(fill(length(ranges[h,1][1]),3)...)
		to_plot[2:end-1,2:end-1,2:end-1] = Us_cn2[h,k]
		surface!(ranges[h,k][1], ranges[h,k][1], to_plot[:,:,zdx]; shading=false, colormap=:lapaz)
	end
	supertitle = Label(f1[0,:], "Crank-Nicholson", textsize = 30)
	f1
end

# ╔═╡ 4b1b3872-839b-49be-aac3-0e436c86a960
md"""### Conclusions
The Matrix splitting method to create locally one-dimensional solvers sped up the calculation significantly. Using a sparse solver with no matrix splitting for the backwards difference method resulted in a very large computation time, ~220s, while the matrix splitting method resulted in a computation time of about ~320ms, an almost 1000x speed up. With further work, all the methods could work this well.
I made the choice to transpose my data so that at each step solving the split M for one dimension, I first transpose the data to make sure that dimension in use is on the primary (contiguous) axis. Many of the other people who presented did not make use of such a method, and I think comparing the difference of permuting the dimensions vs not on the computation time would be an excellent topic for future work.
Ultimately the Crank-Nicholson method and the Backward Difference Extrapolation method both worked far better than the plain Backward difference method, however they were both slower, as they are more computationally expensive. Crank-Nicholson also exhibited some sever oscilliations on the boundaries near discontinuities in the initial conditions, so as always care must be taken to select proper k and h values to have smooth plots.
"""