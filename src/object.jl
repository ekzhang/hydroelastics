using Polyhedra.Core;
using LinearAlgebra;

using Polyhedra
import GLPK
lib = DefaultLibrary{Float64}(GLPK.Optimizer)


struct Mesh
    n::Int32 # number of vertices
    m::Int32 # number of tets
    verts::Matrix{Float64} # shape: [3, n]
    tets::Matrix{Int64} # shape: [4, m]
    potentials::Vector{Float64} # shape: [n]

    Mesh(verts, tets, potentials) = begin
        #println("Mesh: verts = ", verts)
        @assert(size(verts, 1) == 3, "verts should be in R3")
        @assert(size(tets, 1) == 4, "tets should be 4-tuples of indices")
        @assert(size(potentials, 1) == size(verts, 2), "potentials should be the same size as verts")
        n, m = size(verts, 2), size(tets, 2)
        new(n, m, verts, tets, potentials)
    end
end

function intersect_tets(m1, m2, a_face_idx, b_face_idx)
    # Usage: Meshes m1, m2, followed by a_fact_idx, b_face_idx. 
    # TODO: Include poses
    coords_A = m1.verts[:, m1.tets[1:4, a_face_idx]] # 3 x 4 matrix
    coords_B = m2.verts[:, m2.tets[1:4, b_face_idx]] # 3 x 4 matrix

    function get_equations(coords, potentials)
        ones_arr = ones(size(coords, 2), 1)
        mat = hcat(transpose(coords), ones_arr)
        x = mat \ potentials
        return x
    end

    a_pot = get_equations(coords_A, m1.potentials)
    b_pot = get_equations(coords_B, m2.potentials)
    intersection_eq = a_pot - b_pot  # intersection \cdot x = 0
    if norm(intersection_eq[1:3]) < 1e-6
        return []
    end
    function isect_tet_plane(intersection_eq, coords)
        ones_arr = ones(size(coords, 2), 1)
        mat = hcat(transpose(coords), ones_arr)
        dot_prods = mat * intersection_eq
        intersection_points = Vector{Float64}[]
        for i = 1:4
            for j = i+1:4
                if dot_prods[i] * dot_prods[j] < 0
                    frac = dot_prods[i] / (dot_prods[i] - dot_prods[j])
                    push!(intersection_points, (1 - frac) * coords[:, i] + frac * coords[:, j])
                end
            end
        end
        return intersection_points
    end
    isect_A = hcat(isect_tet_plane(intersection_eq, coords_A)...)
    isect_B = hcat(isect_tet_plane(intersection_eq, coords_B)...)
    begin
        xproj = false
        yproj = false
        zproj = false
        begin
            if intersection_eq[1] != 0
                twoDproj = [0.0 1.0 0.0; 0.0 0.0 1.0]
                xproj = true
            elseif intersection_eq[2] != 0
                twoDproj = [1.0 0.0 0.0; 0.0 0.0 1.0]
                yproj = true
            elseif intersection_eq[3] != 0
                twoDproj = [1.0 0.0 0.0; 0.0 1.0 0.0]
                zproj = true
            end
        end
    end

    function project_to_2D(coords, proj)
        return proj * coords
    end
    Apoly = project_to_2D(isect_A, twoDproj)
    Bpoly = project_to_2D(isect_B, twoDproj)
    PA = polyhedron(vrep(transpose(Apoly)), lib)
    PB = polyhedron(vrep(transpose(Bpoly)), lib)
    res = polyhedron(vrep(intersect(PA, PB)))
    all_points = hcat(points(res.vrep)...)
    final_res = zeros(3, size(all_points, 2))
    for i = 1:size(all_points, 2)
        if xproj
            final_res[2, i] = all_points[1, i]
            final_res[3, i] = all_points[2, i]
            final_res[1, i] = (-intersection_eq[4] - intersection_eq[2] * final_res[2, i] - intersection_eq[1] * final_res[1, i]) / intersection_eq[1]
        elseif yproj
            final_res[1, i] = all_points[1, i]
            final_res[3, i] = all_points[2, i]
            final_res[2, i] = (-intersection_eq[4] - intersection_eq[3] * final_res[3, i] - intersection_eq[1] * final_res[1, i]) / intersection_eq[2]
        elseif zproj
            final_res[1, i] = all_points[1, i]
            final_res[2, i] = all_points[2, i]
            final_res[3, i] = (-intersection_eq[4] - intersection_eq[1] * final_res[1, i] - intersection_eq[2] * final_res[2, i]) / intersection_eq[3]
        end
    end
    return final_res
end

mutable struct Object
    mesh::Mesh
    #transform::Transform3D
    #frame::CartesianFrame3D
end

