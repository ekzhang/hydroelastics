### A Pluto.jl notebook ###
# v0.17.3

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

# ╔═╡ c98326ca-13cf-48a6-b423-cb584b92ca5a
world = Dict(
    :sphere => translate(make_icosphere(2), @SVector [0.0, 0.0, 2.0]),
    :ground => translate(make_cube(10.0), @SVector [0.0, 0.0, -4.9]),
	:ramp => translate(rotateX(make_cube(3.0), 0.3), @SVector [0.0, 0.0, -1.0]),
)

# ╔═╡ 0be761c5-a5ed-41d8-a3f5-7e172b7dedcc
num_steps = 1000

# ╔═╡ 44f26aa2-f18e-45cc-9e76-22d00ab8c2a5
poses = simulate(world; dt = 0.001, num_steps = num_steps, forces = (a, α) -> begin
    a[:sphere] += @SVector [0, 0, -9.8] # gravity
    a[:ground] = α[:ground] = @SVector zeros(3) # ground is static
    a[:ramp] = α[:ramp] = @SVector zeros(3) # ramp is static
end)

# ╔═╡ 861e70f8-8b29-40d3-8a1f-915fd778526c
@bind step Slider(1:1:num_steps+1; show_value = true, default = 1)

# ╔═╡ 27b6e401-897c-459a-b6f6-575f3f89b0b5
begin
    vis = Visualizer()
    IJuliaCell(vis)
end

# ╔═╡ e2fe8d8d-2c82-466e-b850-02a67d70dd82
begin
    setobject!(vis[:sphere], render_obj(poses[:sphere][step]),
        MeshPhongMaterial(color = RGBA(0.3, 0.1, 0.7, 1.0)))
    setobject!(vis[:ground], render_obj(poses[:ground][step]),
        MeshPhongMaterial(color = RGBA(0.4, 0.3, 0.2, 1.0)))
    setobject!(vis[:ramp], render_obj(poses[:ramp][step]),
        MeshPhongMaterial(color = RGBA(0.8, 0.3, 0.4, 1.0)))
end

# ╔═╡ Cell order:
# ╠═e708e17d-617e-4a3f-bca1-164d22c3dffc
# ╠═b0dbee5a-9881-48be-9fe8-56cc4f09f4dd
# ╠═c98326ca-13cf-48a6-b423-cb584b92ca5a
# ╠═0be761c5-a5ed-41d8-a3f5-7e172b7dedcc
# ╠═44f26aa2-f18e-45cc-9e76-22d00ab8c2a5
# ╠═861e70f8-8b29-40d3-8a1f-915fd778526c
# ╠═27b6e401-897c-459a-b6f6-575f3f89b0b5
# ╠═e2fe8d8d-2c82-466e-b850-02a67d70dd82
