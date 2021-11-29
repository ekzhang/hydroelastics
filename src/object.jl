lib = DefaultLibrary{Float64}(diff_optimizer(GLPK.Optimizer))

struct Mesh
    n::Int32 # number of vertices
    m::Int32 # number of tets
    verts::Matrix{Float64} # shape: [3, n]
    tets::Matrix{Int64} # shape: [4, m]
    potentials::Vector{Float64} # shape: [n]
    com::Vector{Float64}

    Mesh(verts::Matrix{Float64}, tets::Matrix{Int64}, potentials::Vector{Float64}) = begin
        @assert(size(verts, 1) == 3, "verts should be in R3")
        @assert(size(tets, 1) == 4, "tets should be 4-tuples of indices")
        @assert(
            size(potentials, 1) == size(verts, 2),
            "potentials should be the same size as verts"
        )
        n, m = size(verts, 2), size(tets, 2)
        com = center_of_mass(verts, tets)
        new(n, m, verts, tets, potentials, com)
    end
end

function center_of_mass(verts::Matrix{Float64}, tets::Matrix{Int64})
    """
    compute center of mass given coords of vertices and tets
    """
    # verts [3,n]; tets [4, m]
    tet_centers = Vector{Float64}[]
    vols = Vector{Float64}()
    for i = 1:size(tets, 2)
        tet_vtxs = verts[:, tets[:, i]] # 3 x 4
        push!(tet_centers, mean(eachcol(tet_vtxs)))
        tet_vtxs_4 = [tet_vtxs; ones(Int, 1, 4)]
        append!(vols, abs(det(tet_vtxs_4) / 6))
    end
    sum(tet_centers .* vols) / sum(vols)
end

function intersect_tets(m1::Mesh, m2::Mesh, a_face_idx::Int64, b_face_idx::Int64)
    # Usage: Meshes m1, m2, followed by a_fact_idx, b_face_idx.
    # TODO: Include poses
    coords_A = m1.verts[:, m1.tets[1:4, a_face_idx]] # 3 x 4 matrix
    coords_B = m2.verts[:, m2.tets[1:4, b_face_idx]] # 3 x 4 matrix
    #println("coords_A", coords_A)
    #println("coords_B", coords_B)

    function get_equations(coords::Matrix{Float64}, potentials::Vector{Float64})
        ones_arr = ones(size(coords, 2), 1)
        mat = hcat(transpose(coords), ones_arr)
        x = mat \ potentials
        return x
    end

    a_pot = get_equations(coords_A, m1.potentials[m1.tets[1:4, a_face_idx]])
    b_pot = get_equations(coords_B, m2.potentials[m2.tets[1:4, b_face_idx]])
    intersection_eq = a_pot - b_pot  # intersection \cdot x = 0

    tet_A = polyhedron(vrep(coords_A'), lib)
    tet_B = polyhedron(vrep(coords_B'), lib)


    intersection = intersect(tet_A, tet_B)
    #println("intersection ", intersection)
    #println("intersection_eq ", intersection_eq)
    if norm(intersection_eq[1:3]) < 1e-6
        return zeros(3, 0)
    else
        affine_subspace = hrep([HyperPlane(intersection_eq[1:3], -intersection_eq[4])])
        intersection = intersect(intersection, affine_subspace)
    end
    all_points = hcat(points(polyhedron(vrep(intersection)).vrep)...)

    #intersection_2 = intersect(intersection, affine_subspace)
    #res = polyhedron(vrep(intersect(PA, PB)))
    #all_points = hcat(points(res.vrep)...)
    #println(all_points)
    if size(all_points, 2) < 3
        return zeros(3, 0)
    else
        return all_points
    end
end

#for halfspace in allhalfspaces(intersection_2)
#    println(halfspace)
#end


#intersection_final = hcat(polyhedron(vrep(intersect(intersection, HyperPlane(intersection_eq[1:3], -intersection_eq[4]))...)))


#     function isect_tet_plane(intersection_eq::Vector{Float64}, coords::Matrix{Float64})
#         ones_arr = ones(size(coords, 2), 1)
#         mat = hcat(transpose(coords), ones_arr)
#         dot_prods = mat * intersection_eq
#         intersection_points = Vector{Float64}[]
#         for i = 1:4
#             for j = i+1:4
#                 if dot_prods[i] * dot_prods[j] < 0.0
#                     frac = dot_prods[i] / (dot_prods[i] - dot_prods[j])
#                     push!(
#                         intersection_points,
#                         (1.0 - frac) * coords[:, i] + frac * coords[:, j],
#                     )
#                 end
#             end
#         end
#         return hcat(intersection_points...)
#     end

#     isect_A = isect_tet_plane(intersection_eq, coords_A)
#     isect_B = isect_tet_plane(intersection_eq, coords_B)

#     if (isempty(isect_A) || (isempty(isect_B)))
#         return zeros(3, 0)
#     end
#     xproj = false
#     yproj = false
#     zproj = false
#     if abs(intersection_eq[1]) / norm(intersection_eq[1:3]) > 1e-3
#         twoDproj = [0.0 1.0 0.0; 0.0 0.0 1.0]
#         xproj = true
#     elseif abs(intersection_eq[2]) / norm(intersection_eq[1:3]) > 1e-3
#         twoDproj = [1.0 0.0 0.0; 0.0 0.0 1.0]
#         yproj = true
#     elseif abs(intersection_eq[3]) / norm(intersection_eq[1:3]) > 1e-3
#         twoDproj = [1.0 0.0 0.0; 0.0 1.0 0.0]
#         zproj = true
#     end

#     Apoly = twoDproj * isect_A
#     Bpoly = twoDproj * isect_B
#     PA = polyhedron(vrep(Apoly'), lib)
#     PB = polyhedron(vrep(Bpoly'), lib)
#     res = polyhedron(vrep(intersect(PA, PB)))
#     all_points = hcat(points(res.vrep)...)
#     if isempty(all_points)
#         return zeros(3, 0)
#     end
#     final_res = zeros(3, size(all_points, 2))
#     if size(all_points, 2) <= 2
#         return zeros(3, 0)
#     end
#     for i = 1:size(all_points, 2)
#         if xproj
#             final_res[2, i] = all_points[1, i]
#             final_res[3, i] = all_points[2, i]
#             final_res[1, i] =
#                 (
#                     -intersection_eq[4] - intersection_eq[2] * final_res[2, i] -
#                     intersection_eq[1] * final_res[1, i]
#                 ) / intersection_eq[1]
#         elseif yproj
#             final_res[1, i] = all_points[1, i]
#             final_res[3, i] = all_points[2, i]
#             final_res[2, i] =
#                 (
#                     -intersection_eq[4] - intersection_eq[3] * final_res[3, i] -
#                     intersection_eq[1] * final_res[1, i]
#                 ) / intersection_eq[2]
#         elseif zproj
#             final_res[1, i] = all_points[1, i]
#             final_res[2, i] = all_points[2, i]
#             final_res[3, i] =
#                 (
#                     -intersection_eq[4] - intersection_eq[1] * final_res[1, i] -
#                     intersection_eq[2] * final_res[2, i]
#                 ) / intersection_eq[3]
#         end
#     end
#     final_res
# end

function triangulate_polygon(vertices::Matrix{Float64})
    """
    given out-of-order coordinates of vertices of a planar and convex polygon in a 3xN matrix, output a triangulation
    of the polygon
    """
    n = size(vertices)[2]
    com = [sum(vertices[i, :]) for i = 1:3] / n
    angles = [0.0]  # list of (angle, vtx index corresponding to angle)
    displacements = vertices - reshape(repeat(com, n), 3, n)
    mags = [sqrt(dot(displacements[:, i], displacements[:, i])) for i = 1:n]
    ind = 2
    normal = cross(displacements[:, 1], displacements[:, ind])  # compute any normal to the plane
    while norm(normal) < 1e-6 * mags[1] * mags[ind]
        ind += 1
        if ind > n
            return zeros(3, 0)
        end
        normal = cross(displacements[:, 1], displacements[:, ind])
    end
    initial_disp = displacements[:, 1]
    for ind = 2:n
        disp = displacements[:, ind]
        angle = acos(
            clamp(dot(disp, initial_disp) / (mags[1] * mags[ind]), -(1 - 1e-9), (1 - 1e-9)),
        )
        if dot(cross(disp, initial_disp), normal) < 0
            angle = 2 * pi - angle
        end
        push!(angles, angle)
    end
    order = sortperm(angles)
    # 1, 2, 3; 1, 3, 4;, 1, 4, 5 ..., 1, n-1, n
    hcat([[1, order[i], order[i+1]] for i = 2:n-1]...)::Matrix{Int64}
end

function tet_force(A::Mesh, B::Mesh, i::Int64, j::Int64)
    """
    compute overall force between two tets A[i], B[j] of objects A, B.
    return the force applied to A (in the direction of A.com - B.com)
    """
    normal = zeros(3)
    intersection_polygon = intersect_tets(A, B, i, j)
    total_force = 0.0
    if (size(intersection_polygon)[1] > 0) && (size(intersection_polygon)[2] > 0)
        triangles = triangulate_polygon(intersection_polygon)
        if isempty(triangles)
            return zeros(3)
        end
        vtx_inds = A.tets[:, i]
        vtx_coords = A.verts[:, vtx_inds]
        Acom = mean(eachcol(vtx_coords))
        Bcom = mean(eachcol(B.verts[:, B.tets[:, j]]))
        vtx_coords = vcat(vtx_coords, ones(1, 4))  # 4x4 matrix of vtxs padded w 1s
        for xyz in eachcol(triangles)
            vtxs = intersection_polygon[:, xyz]
            com = push!(mean(eachcol(vtxs)), 1)
            res = vtx_coords \ com
            pressure = sum(res .* A.potentials[vtx_inds])
            area =
                0.5 * norm(cross(vtxs[:, 1] - vtxs[:, 2], vtxs[:, 1] - vtxs[:, 3]))
            total_force += pressure * area
        end
        normal = cross(
            intersection_polygon[:, 1] - intersection_polygon[:, 2],
            intersection_polygon[:, 1] - intersection_polygon[:, 3],
        )
        normal = normal / norm(normal)
        if dot(normal, Acom - Bcom) < 0
            normal = -1 * normal
        end
    end
    total_force * normal
end

function mesh_force(A::Mesh, B::Mesh)
    """
    Computes force on mesh A due to contact with mesh B
    """
    force = zeros(3)
    for i = 1:A.m
        for j = 1:B.m
            force += tet_force(A, B, i, j)
        end
    end
    force
end

mutable struct Object
    mesh::Mesh
    #transform::Transform3D
    #frame::CartesianFrame3D
end

export Mesh, Object, intersect_tets, triangulate_polygon, tet_force, mesh_force
