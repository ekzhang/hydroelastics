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
    #using PlutoUI
    using MeshCat
end

# ╔═╡ 36ac084d-c85c-488c-bd20-7a490fa471c9
cu1 = make_cube()

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
    for i in LinRange(-1.0, 1.0, 301)
        #for i in LinRange(-0.1, 0.1, 9)
        cu1 = make_cube()
        cu2 = make_cube([i, 0, 0])

        push!(x, i)

        push!(all_forces, norm(compute_force(cu1, cu2).F_AB))

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

# ╔═╡ 431997b4-0cb0-47c9-8ffb-38ae5615593c
begin
    i = 0.05
    cu2 = make_cube([i, 0, 0])
    transform(cu2.mesh.verts[:, 1:9], cu2.pose)
end

# ╔═╡ d15ae364-c519-48db-838f-4f58226c1197
plot([1, 2, 3], [4, 5, 6])

# ╔═╡ Cell order:
# ╠═74025e80-8603-4254-a4d1-9d11e656b689
# ╠═daf02e8c-4cad-11ec-1211-07b7add5573d
# ╠═36ac084d-c85c-488c-bd20-7a490fa471c9
# ╠═22c3247a-7cf8-4cc3-b84f-8a510942d399
# ╠═431997b4-0cb0-47c9-8ffb-38ae5615593c
# ╠═d15ae364-c519-48db-838f-4f58226c1197
