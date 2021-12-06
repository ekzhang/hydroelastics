### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ e708e17d-617e-4a3f-bca1-164d22c3dffc
begin
    import Pkg
    Pkg.activate("..")
end

# ╔═╡ b0dbee5a-9881-48be-9fe8-56cc4f09f4dd
begin
    using Revise
    using Hydroelastics

    using MeshCat, Colors
    using GeometryBasics
    using CoordinateTransformations
    using LinearAlgebra, StaticArrays
    using PlutoUI
end

# ╔═╡ 8bebb189-d470-4077-88dd-a9ad97b94b76
@bind d Slider(0.0:0.01:1.0; show_value = true, default = 0.5)

# ╔═╡ 27b6e401-897c-459a-b6f6-575f3f89b0b5
begin
    vis = Visualizer()
    IJuliaCell(vis)
end

# ╔═╡ 90ff5a44-0464-4e10-8494-99a91b301ed4
sphere = make_icosphere(2)

# ╔═╡ 51f5ad18-2172-41b9-ba71-96e0e1871eff
cube = make_cube()

# ╔═╡ c9706e48-450e-4d1f-8541-737d7e0aa470
begin
    object1 = translate(sphere, @SVector [0.031, -d, 0.052])
    object2 = translate(sphere, @SVector [0, d, 0])
end

# ╔═╡ e32777f6-d54f-4ef1-be3f-0c57c1d01cd0
compute_force(object1, object2)

# ╔═╡ 539ed9a3-71d4-4da6-bc46-743718018abd
begin
    local material1 = MeshPhongMaterial(color = RGBA(0.3, 0.1, 0.7, 1.0))
    setobject!(vis[:object1], render_obj(object1), material1)
    local material2 = MeshPhongMaterial(color = RGBA(0.8, 0.3, 0.4, 1.0))
    setobject!(vis[:object2], render_obj(object2), material2)
end

# ╔═╡ f63d1c23-78cc-489c-8acc-62c7d0aca7f3
begin
	import Hydroelastics: HalfPlane, Point, intersect_halfplanes
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
	    HalfPlane(Point(-0.16469227056401645, -0.2664776914569725), Point(-0.8506508083520399, 0.5257311121191335)),
	    HalfPlane(Point(0.028287259696359788, -0.040547756753354335), Point(0.8201447309019166, 0.5721561153861965)),
	    HalfPlane(Point(-0.16620458460952442, 0.09547326678562677), Point(-0.49810107091297723, -0.8671189786617204)),
	    HalfPlane(Point(0.2970629430715192, -0.23999299562154927), Point(0.6284282897093981, 0.77786752387082)),
	    HalfPlane(Point(-0.05197402131154934, 0.06002998491806123), Point(-0.7560128937297079, -0.6545567236797993)),
	    HalfPlane(Point(0.1879473775300143, 0.30410524493997143), Point(0.85065080835204, -0.5257311121191336)),
	    HalfPlane(Point(0.1619918586413753, -0.0851641599977757), Point(0.46534112719498627, 0.8851314226378474)),
	    HalfPlane(Point(0.2970629430715192, -0.23999299562154938), Point(0.6284282897093982, 0.7778675238708199)),
	]

	plot()
	for hp in halfplanes
		x, y = hp.p
		x1, y1 = hp.pq
		quiver!([x - 0.5 * x1], [y - 0.5 * y1]; quiver=([x1], [y1]))
	end

	points = intersect_halfplanes(halfplanes)
	if !isnothing(points)
		scatter!(points[1, :], points[2, :])
	end

	plot!(legend=:none)
end

# ╔═╡ Cell order:
# ╠═e708e17d-617e-4a3f-bca1-164d22c3dffc
# ╠═b0dbee5a-9881-48be-9fe8-56cc4f09f4dd
# ╠═8bebb189-d470-4077-88dd-a9ad97b94b76
# ╠═e32777f6-d54f-4ef1-be3f-0c57c1d01cd0
# ╠═27b6e401-897c-459a-b6f6-575f3f89b0b5
# ╠═90ff5a44-0464-4e10-8494-99a91b301ed4
# ╠═51f5ad18-2172-41b9-ba71-96e0e1871eff
# ╠═c9706e48-450e-4d1f-8541-737d7e0aa470
# ╠═539ed9a3-71d4-4da6-bc46-743718018abd
# ╠═f63d1c23-78cc-489c-8acc-62c7d0aca7f3
