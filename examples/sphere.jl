### A Pluto.jl notebook ###
# v0.17.2

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

    using MeshCat, Colors
    using GeometryBasics
    using CoordinateTransformations
    using LinearAlgebra
end

# ╔═╡ 27b6e401-897c-459a-b6f6-575f3f89b0b5
begin
    vis = Visualizer()
    IJuliaCell(vis)
end

# ╔═╡ 90ff5a44-0464-4e10-8494-99a91b301ed4
sphere = make_icosphere(3)

# ╔═╡ 539ed9a3-71d4-4da6-bc46-743718018abd
begin
	local material = MeshPhongMaterial(color=RGBA(0.3, 0.1, 0.7, 1.0))
    setobject!(vis, mesh_to_obj(sphere), material)
    settransform!(vis, Translation(0, 0, 1.0))
end

# ╔═╡ Cell order:
# ╠═e708e17d-617e-4a3f-bca1-164d22c3dffc
# ╠═b0dbee5a-9881-48be-9fe8-56cc4f09f4dd
# ╠═27b6e401-897c-459a-b6f6-575f3f89b0b5
# ╠═90ff5a44-0464-4e10-8494-99a91b301ed4
# ╠═539ed9a3-71d4-4da6-bc46-743718018abd
