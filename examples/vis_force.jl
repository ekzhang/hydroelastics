### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ e708e17d-617e-4a3f-bca1-164d22c3dffc
begin
    import Pkg
    Pkg.activate("..")
end

# ╔═╡ b0dbee5a-9881-48be-9fe8-56cc4f09f4dd
begin
    using Revise
    using Hydroelastics
	using Polyhedra
    using MeshCat, Colors
    using GeometryBasics
    using CoordinateTransformations
    using LinearAlgebra, StaticArrays
    using PlutoUI
	using GLPK
	lib = DefaultLibrary{Float64}(GLPK.Optimizer)

end

# ╔═╡ 539ed9a3-71d4-4da6-bc46-743718018abd
begin
    local material1 = MeshPhongMaterial(color = RGBA(0.3, 0.1, 0.7, 1.0))
    #setobject!(vis[:object1], render_obj(object1), material1)
    local material2 = MeshPhongMaterial(color = RGBA(0.8, 0.3, 0.4, 1.0))
    #setobject!(vis[:object2], render_obj(object2), material2)
end

# ╔═╡ f63d1c23-78cc-489c-8acc-62c7d0aca7f3
begin
	import Hydroelastics: HalfPlane, Point, intersect_halfplanes, intersect_halfplanes_slow
	import Plots: plot, quiver!, plot!, scatter!
	
	halfplanes = [
	    HalfPlane(Point(0.05328478907262187, 0.032931810730249304), Point(-0.5257311121191336, 0.85065080835204)),
	    HalfPlane(Point(-0.1780353002260096, 0.2880671669629754), Point(0.85065080835204, 0.5257311121191336)),
	    HalfPlane(Point(-0.3302474442330619, 1.1330056073001064e-17), Point(-3.4307778215552914e-17, -1.0)),
	    HalfPlane(Point(-0.0, 1.0685339887498948), Point(-1.0, -0.0)),
	    HalfPlane(Point(0.007597871376374749, -0.004695742752749542), Point(-0.5257311121191337, -0.8506508083520398)),
	    HalfPlane(Point(0.1547801932600118, 0.2504396134799765), Point(0.85065080835204, -0.5257311121191336)),
	    HalfPlane(Point(0.3501096768180673, -1.2011485143393161e-17), Point(3.4307778215552914e-17, 1.0)),
	    HalfPlane(Point(1.919491578940921e-16, 1.0685339887498948), Point(-1.0, 1.7963785889362146e-16)),
	]

	halfplanes = [
	HalfPlane(Point(0,0), Point(1,0)),
	HalfPlane(Point(1,0), Point(0,1)),
	HalfPlane(Point(1,1), Point(-1,0)),
	HalfPlane(Point(0,1), Point(0,-1))
]

	
	plot()
	for hp in halfplanes
		x, y = hp.p
		x1, y1 = hp.pq
		quiver!([x - 0.5 * x1], [y - 0.5 * y1]; quiver=([x1], [y1]))
	end

	pts = intersect_halfplanes(halfplanes)
	pts2 = intersect_halfplanes_slow(halfplanes)
	println("pts", pts)
	println("pts2", pts2)
	if !isnothing(pts)
		scatter!(pts[1, :], pts[2, :])
	end

	plot!(legend=:none)
end

# ╔═╡ 23b30021-68a8-4606-9a00-c8fb16ae65bf
begin 
	A3 = [[1.0,0.0,-1.0,0.0] [0.0,1.0,0.0,-1.0]]
	b3 = [1.0,1.0,1.0,1.0]
	res3 = polyhedron(hrep(A3, b3, BitSet()))
	plot(res3, ratio=:equal)
end

# ╔═╡ 71ded04f-7bef-4496-a4d3-d57bca6c6a10


# ╔═╡ Cell order:
# ╠═e708e17d-617e-4a3f-bca1-164d22c3dffc
# ╠═b0dbee5a-9881-48be-9fe8-56cc4f09f4dd
# ╠═23b30021-68a8-4606-9a00-c8fb16ae65bf
# ╠═539ed9a3-71d4-4da6-bc46-743718018abd
# ╠═f63d1c23-78cc-489c-8acc-62c7d0aca7f3
# ╠═71ded04f-7bef-4496-a4d3-d57bca6c6a10
