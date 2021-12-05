function transform_vertices(transform, vertices)
    # concatenate ones to vertices
    ones_arr = ones(1, size(vertices, 2))
    mat = vcat(vertices, ones_arr)
    vertices = transform * mat
    return vertices[1:3, :]
end

function intersect_tets(o1::Object, o2::Object, a_face_idx::Int64, b_face_idx::Int64)
    #error("unimplemented")

    # Usage: Meshes m1, m2, followed by a_fact_idx, b_face_idx.
    m1 = o1.mesh
    m2 = o2.mesh
    coords_A = transform_vertices(o1.pose, m1.verts[:, m1.tets[1:4, a_face_idx]]) # 3 x 4 matrix
    coords_B = transform_vertices(o2.pose, m2.verts[:, m2.tets[1:4, b_face_idx]]) # 3 x 4 matrix

    function get_equations(coords::Matrix{Float64}, potentials::Vector{Float64})
        ones_arr = ones(size(coords, 2), 1)
        mat = hcat(transpose(coords), ones_arr)
        x = mat \ potentials
        return x
    end

    a_pot = get_equations(coords_A, m1.potentials[m1.tets[1:4, a_face_idx]])
    b_pot = get_equations(coords_B, m2.potentials[m2.tets[1:4, b_face_idx]])
    intersection_eq = a_pot - b_pot  # intersection \cdot x = 0
    if norm(intersection_eq[1:3]) < 1e-6
        return zeros(3, 0)
    end

    function isect_tet_plane(intersection_eq::Vector{Float64}, coords::Matrix{Float64})
        ones_arr = ones(size(coords, 2), 1)
        mat = hcat(transpose(coords), ones_arr)
        dot_prods = mat * intersection_eq
        intersection_points = Vector{Float64}[]
        for i = 1:4
            for j = i+1:4
                if dot_prods[i] * dot_prods[j] < 0.0
                    frac = dot_prods[i] / (dot_prods[i] - dot_prods[j])
                    push!(
                        intersection_points,
                        (1.0 - frac) * coords[:, i] + frac * coords[:, j],
                    )
                end
            end
        end
        for i = 1:4 # there might be a better way to deal w vtx on plane cases
            if abs(dot_prods[i]) < 1e-9 * norm(coords[:, i]) * norm(intersection_eq)
                push!(intersection_points, coords[:, i])
            end
        end
        return hcat(intersection_points...)
    end

    isect_A = isect_tet_plane(intersection_eq, coords_A)
    isect_B = isect_tet_plane(intersection_eq, coords_B)

    if (isempty(isect_A) || (isempty(isect_B)))
        return zeros(3, 0)
    end
    xproj = false
    yproj = false
    zproj = false
    if abs(intersection_eq[1]) / norm(intersection_eq[1:3]) > 1e-3
        twoDproj = [0.0 1.0 0.0; 0.0 0.0 1.0]
        xproj = true
    elseif abs(intersection_eq[2]) / norm(intersection_eq[1:3]) > 1e-3
        twoDproj = [1.0 0.0 0.0; 0.0 0.0 1.0]
        yproj = true
    elseif abs(intersection_eq[3]) / norm(intersection_eq[1:3]) > 1e-3
        twoDproj = [1.0 0.0 0.0; 0.0 1.0 0.0]
        zproj = true
    end

    Apoly = twoDproj * isect_A
    Bpoly = twoDproj * isect_B
    PA = polyhedron(vrep(Apoly'), lib)
    PB = polyhedron(vrep(Bpoly'), lib)
    res = polyhedron(vrep(intersect(PA, PB)))
    all_points = hcat(points(res.vrep)...)
    if isempty(all_points)
        return zeros(3, 0)
    end
    final_res = zeros(3, size(all_points, 2))
    if size(all_points, 2) <= 2
        return zeros(3, 0)
    end
    for i = 1:size(all_points, 2)
        if xproj
            final_res[2, i] = all_points[1, i]
            final_res[3, i] = all_points[2, i]
            final_res[1, i] =
                (
                    -intersection_eq[4] - intersection_eq[2] * final_res[2, i] -
                    intersection_eq[1] * final_res[1, i]
                ) / intersection_eq[1]
        elseif yproj
            final_res[1, i] = all_points[1, i]
            final_res[3, i] = all_points[2, i]
            final_res[2, i] =
                (
                    -intersection_eq[4] - intersection_eq[3] * final_res[3, i] -
                    intersection_eq[1] * final_res[1, i]
                ) / intersection_eq[2]
        elseif zproj
            final_res[1, i] = all_points[1, i]
            final_res[2, i] = all_points[2, i]
            final_res[3, i] =
                (
                    -intersection_eq[4] - intersection_eq[1] * final_res[1, i] -
                    intersection_eq[2] * final_res[2, i]
                ) / intersection_eq[3]
        end
    end
    final_res
end

"""
Given out-of-order coordinates of vertices of a planar and convex polygon in a
3xN matrix, output a triangulation of the polygon.
"""
function triangulate_polygon(vertices::Matrix{Float64})
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

"""
Result of a force calculation, expressed as vectors in the world frame.
"""
struct ForceResult
    F_AB::SVector{3} # The force applied to object A from B.
    F_BA::SVector{3} # The force applied to object B from A.
    τ_AB::SVector{3} # The torque applied to object A from B.
    τ_BA::SVector{3} # The torque applied to object B from A.
end

"""
Compute overall force between two tets A[i], B[j] of objects A, B.

Returns the force applied to A (in the direction of A.com - B.com), along with
the net torque vetors applied to A and B.
"""
function tet_force(A::Object, B::Object, i::Int64, j::Int64)::Vector{Float64}
    #error("unimplemented")
    normal = zeros(3)
    intersection_polygon = intersect_tets(A, B, i, j)
    total_force = 0.0
    if (size(intersection_polygon)[1] > 0) && (size(intersection_polygon)[2] > 0)
        triangles = triangulate_polygon(intersection_polygon)
        if isempty(triangles)
            return zeros(3)
        end
        vtx_inds = A.mesh.tets[:, i]
        vtx_coords = transform_vertices(A.pose, A.mesh.verts[:, vtx_inds])
        vtx_coords = vcat(vtx_coords, ones(1, 4))  # 4x4 matrix of vtxs padded w 1s
        for xyz in eachcol(triangles)
            # calculate barycentric coordinates of each point in the triangle.
            vtxs = intersection_polygon[:, xyz]
            com = push!(mean(eachcol(vtxs)), 1)
            res = vtx_coords \ com
            pressure = sum(res .* A.mesh.potentials[vtx_inds])
            area = 0.5 * norm(cross(vtxs[:, 1] - vtxs[:, 2], vtxs[:, 1] - vtxs[:, 3]))
            total_force += pressure * area
        end
        normal = cross(
            intersection_polygon[:, 1] - intersection_polygon[:, 2],
            intersection_polygon[:, 1] - intersection_polygon[:, 3],
        )
        normal = normal / norm(normal)
        if dot(normal, A.mesh.com - B.mesh.com) < 0
            normal = -1 * normal
        end
    end
    total_force * normal
end

"""
Computes force on object A due to contact with object B, along with the net
torques on both objects.
"""
function compute_force(A::Object, B::Object)::ForceResult
    #error("unimplemented")
    force = zeros(3)
    τ_AB = zeros(3)
    τ_BA = zeros(3)

    for i = 1:A.mesh.m
        for j = 1:B.mesh.m
            force += tet_force(A, B, i, j)
            #τ_AB += cross(A.mesh.com - B.mesh.com, force)
            #τ_BA += cross(B.mesh.com - A.mesh.com, force)
        end
    end
    τ_AB = cross(A.mesh.com - B.mesh.com, force)
    τ_BA = -1 * τ_AB
    #println("force: ", force)
    ForceResult(force, -1 * force, τ_AB, τ_BA)
end

export compute_force, transform_vertices
