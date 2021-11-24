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
    using Hydroelastics

    using MeshCat
    using GeometryBasics
    using CoordinateTransformations
end

# ╔═╡ 27b6e401-897c-459a-b6f6-575f3f89b0b5
begin
    vis = Visualizer()
    IJuliaCell(vis)
end

# ╔═╡ d4f9ca2a-1913-48a6-bed1-ab7633615752
begin
    setobject!(vis, HyperRectangle(Vec(0.0, 0, 0), Vec(1.0, 1, 1)))
    settransform!(vis, Translation(-0.5, -0.5, 0))
end

# ╔═╡ Cell order:
# ╠═e708e17d-617e-4a3f-bca1-164d22c3dffc
# ╠═b0dbee5a-9881-48be-9fe8-56cc4f09f4dd
# ╠═27b6e401-897c-459a-b6f6-575f3f89b0b5
# ╠═d4f9ca2a-1913-48a6-bed1-ab7633615752
