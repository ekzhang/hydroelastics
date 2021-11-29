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
	using Polyhedra
	#using PlutoUI
	using MeshCat
end

# ╔═╡ 36ac084d-c85c-488c-bd20-7a490fa471c9
function get_cube(com, sz)
    """
    returns a cube with center at com
    """
    cube_verts =
        com .+
        sz * [
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
    cube = Hydroelastics.Mesh(cube_verts, cube_tets, cube_pots)
    cube
end


# ╔═╡ 22c3247a-7cf8-4cc3-b84f-8a510942d399
begin
    x = []
    y = []
	all_forces = []
	specific_tet_force = []
	#vis = Visualizer()
	#m1 = polyhedron(vrep([-1.0 -1.0 -1.0 0.0; 1.0 -1.0 1.0 0.0; 1.0 1.0 -1.0 0.0]'))
	#m2 = polyhedron(vrep([-1.0 1.0 -1.0 0.0; 1.0 1.0 -1.0 0.0; 1.0 1.0 1.0 0.0]'))
	#setobject!(vis, Polyhedra.Mesh(m1))
	#setobject!(vis, Polyhedra.Mesh(m2))

	#IJuliaCell(vis)
	#plot([1,2,3], [4,5,6])
	#with_terminal() do 
	    for i in LinRange(-3.0, 3.0, 301)
	        #for i in LinRange(-0.1, 0.1, 9)
	        cu1 = get_cube([0.0, 0.0, 0.0], 1.0)
	        cu2 = get_cube([i, 0.0, 0.0], 1.0)
			push!(x, i)
			
	        #cu2 = get_cube([i,1e-5,-2.4e-5],1.0)
	        push!(all_forces, norm(mesh_force(cu1, cu2)))
			if abs(i) < 1e-5
				for i in 1:12
					for j in 1:12
						if i == j
							@assert (norm(tet_force(cu1, cu2, i, i)) < 1e-5)
						else
							@assert (norm(tet_force(cu1, cu2, i, j) + tet_force(cu1, cu2, j, i)) < 1e-5)
						end
					end
				end
				println("intersect tets ", tet_force(cu1, cu2, 2, 8))
				println("intersect tets ", tet_force(cu1, cu2, 8, 2))
			end
			#push!(specific_tet_force, tet_force(cu1, cu2, 2, 11)[1])
	    end
		z = []
		#for i in 5:5
		#	for j in 10:10
		#		ys = [ele[12 * i + j][1] for ele in all_forces]
				#print(ys)
		#		print(size(ys))
		#		z = ys
		#		print(size(x))
		#	end
		#end
		#plot(x, specific_tet_force)
	#end
	plot(x, all_forces)
end

# ╔═╡ d15ae364-c519-48db-838f-4f58226c1197
plot([1,2,3], [4,5,6])

# ╔═╡ Cell order:
# ╠═74025e80-8603-4254-a4d1-9d11e656b689
# ╠═daf02e8c-4cad-11ec-1211-07b7add5573d
# ╠═36ac084d-c85c-488c-bd20-7a490fa471c9
# ╠═22c3247a-7cf8-4cc3-b84f-8a510942d399
# ╠═d15ae364-c519-48db-838f-4f58226c1197
