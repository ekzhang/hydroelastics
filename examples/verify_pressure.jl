### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ 74025e80-8603-4254-a4d1-9d11e656b689
begin
    import Pkg
    Pkg.activate("..")
end

# ╔═╡ daf02e8c-4cad-11ec-1211-07b7add5573d
begin
	using Revise
    using Hydroelastics
    using LinearAlgebra
	using Plots
end

# ╔═╡ 36ac084d-c85c-488c-bd20-7a490fa471c9
function get_cube(com, sz)
    """
    returns a cube with center at com
    """
    cube_verts =
        com .+ sz * [
            -1.0 -1.0 -1.0 -1.0 1.0 1.0 1.0 1.0 0.0
            -1.0 -1.0 1.0 1.0 -1.0 -1.0 1.0 1.0 0.0
            -1.0 1.0 -1.0 1.0 -1.0 1.0 -1.0 1.0 0.0
        ]
    cube_tets = [
        1 4 1 1 1 1 8 4 8 7 4 8
        2 2 2 5 3 5 4 3 7 6 8 6
        3 3 6 6 7 7 7 7 6 5 2 2
        9 9 9 9 9 9 9 9 9 9 9 9
    ]
    cube_pots = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    cube = Mesh(cube_verts, cube_tets, cube_pots)
    cube
end


# ╔═╡ 22c3247a-7cf8-4cc3-b84f-8a510942d399
begin
	x = []
	y = []
	for i in LinRange(-1.0, 1.0, 300)
		cu1 = get_cube([0.0, 0.0, 0.0], 1.0)
	    cu2 = get_cube([i,1e-5,-2.4e-5],1.0)
		append!(x,i)
		append!(y,norm(mesh_force(cu1,cu2)))
	end
	plot(x,y)
end

# ╔═╡ Cell order:
# ╠═74025e80-8603-4254-a4d1-9d11e656b689
# ╠═daf02e8c-4cad-11ec-1211-07b7add5573d
# ╠═36ac084d-c85c-488c-bd20-7a490fa471c9
# ╠═22c3247a-7cf8-4cc3-b84f-8a510942d399
