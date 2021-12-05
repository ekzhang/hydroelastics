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
sphere = make_icosphere(1)

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
