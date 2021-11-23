lib = DefaultLibrary{Float64}(diff_optimizer(GLPK.Optimizer))

struct Mesh
    n::Int32 # number of vertices
    m::Int32 # number of tets
    verts::Matrix{Float64} # shape: [3, n]
    tets::Matrix{Int64} # shape: [4, m]
    potentials::Vector{Float64} # shape: [n]
    com::Vector{Float64}

    Mesh(verts, tets, potentials) = begin
        #println("Mesh: verts = ", verts)
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

function center_of_mass(verts, tets)
    ```
    compute center of mass given coords of vertices and tets
    ```
    # verts [3,n]; tets [4, m]
    tet_centers = Vector{Float64}()
    vols = Vector{Float64}()
    for i = 1:size(tets)[2]
        tet_vtxs = verts[:, tets[:, i]] # 3 x 4
        append!(tet_centers, [sum(tet_vtxs[i, :]) for i = 1:3] / 4)
        tet_vtxs = [tet_vtxs; ones(Int, 1, 4)]
        append!(vols, abs(det(tet_vtxs) / 6))
    end
    com = [0, 0, 0]
    tet_centers = reshape(tet_centers, 3, size(tets)[2]) # 3 x m
    for i = 1:size(tet_centers)[2]
        com = com + tet_centers[:, i] * vols[i]
    end
    com / sum(vols)
end

function intersect_tets(m1, m2, a_face_idx, b_face_idx)
    # Usage: Meshes m1, m2, followed by a_fact_idx, b_face_idx.
    # TODO: Include poses
    coords_A = m1.verts[:, m1.tets[1:4, a_face_idx]] # 3 x 4 matrix
    coords_B = m2.verts[:, m2.tets[1:4, b_face_idx]] # 3 x 4 matrix

    function get_equations(coords, potentials)
        ones_arr = ones(size(coords, 2), 1)
        mat = hcat(transpose(coords), ones_arr)
        print("mat = ", mat)
        print("potentials = ", potentials)
        x = mat \ potentials
        return x
    end

    a_pot = get_equations(coords_A, m1.potentials[m1.tets[1:4, a_face_idx]])
    b_pot = get_equations(coords_B, m2.potentials[m2.tets[1:4, a_face_idx]])
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
                    push!(
                        intersection_points,
                        (1 - frac) * coords[:, i] + frac * coords[:, j],
                    )
                end
            end
        end
        return hcat(intersection_points...)
    end

    isect_A = isect_tet_plane(intersection_eq, coords_A)
    isect_B = isect_tet_plane(intersection_eq, coords_B)

    if (isect_A == []) || (isect_B == [])
        return []
    end
    xproj = false
    yproj = false
    zproj = false
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

    Apoly = twoDproj * isect_A
    Bpoly = twoDproj * isect_B
    PA = polyhedron(vrep(Apoly'), lib)
    PB = polyhedron(vrep(Bpoly'), lib)
    res = polyhedron(vrep(intersect(PA, PB)))
    all_points = hcat(points(res.vrep)...)
    final_res = zeros(3, size(all_points, 2))
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

function triangulate_polygon(vertices)
    """
    given out-of-order coordinates of vertices of a planar and convex polygon in a 3xN matrix, output a triangulation
    of the polygon
    """
    N = size(vertices)[2]
    com = [sum(vertices[i, :]) for i = 1:3] / N
    angles = [0.0]  # list of (angle, vtx index corresponding to angle)
    displacements = vertices - reshape(repeat(com, N), 3, N)
    mags = [sqrt(dot(displacements[:, i], displacements[:, i])) for i = 1:N]
    ind = 2
    normal = cross(displacements[:, 1], displacements[:, ind])  # compute any normal to the plane
    while sqrt(dot(normal, normal)) < 0.000001 * mags[1] * mags[ind]
        ind += 1
        normal = cross(displacements[:, 1], displacements[:, ind])
    end
    initial_disp = displacements[:, 1]
    for ind = 2:N
        disp = displacements[:, ind]
        angle = acos(dot(disp, initial_disp) / (mags[1] * mags[ind]))
        if dot(cross(disp, initial_disp), normal) < 0
            angle = 2 * pi - angle
        end
        append!(angles, [angle])
    end
    order = sortperm(angles)
    # 0, 1, 2; 0, 2, 3;, 0, 3, 4 ..., 0, n-2, n-1
    hcat([[1, order[i], order[i+1]] for i = 2:N-1]...)
end

function pressure(A, B, i, j)
    """
    compute overall pressure between two tets A[i], B[j] of objects A, B
    """
    total_pressure = 0
    intersection_polygon = intersect_tets(A, B, i, j)
    if size(intersection_polygon)[1] > 0
        triangles = triangulate_polygon(intersection_polygon)
        vtx_inds = A.tets[:, i]
        vtx_coords = A.verts[:, vtx_inds]
        vtx_coords = vcat(vtx_coords, ones(4))  # 4x4 matrix of vtxs padded w 1s
        for xyz in eachcol(triangles)
            vtxs = intersection_polygon[:, xyz]
            com = push!(mean(eachcol(vtxs)), 1)
            res = vtx_coords \ com
            for k = 1:4
                total_pressure += res[k] * A.potentials[vtx_inds[k]]
            end
        end
    end
    total_pressure
end

mutable struct Object
    mesh::Mesh
    #transform::Transform3D
    #frame::CartesianFrame3D
end

export Mesh, Object, intersect_tets, pressure

function get_cube(com)
    """
    returns a cube with center at com
    """
    com = com
    cube_verts = [
        com + [-1.0, -1.0, -1.0],
        com + [-1.0, -1.0, 1.0],
        com + [-1.0, 1.0, -1.0],
        com + [-1.0, 1.0, 1.0],
        com + [1.0, -1.0, -1.0],
        com + [1.0, -1.0, 1.0],
        com + [1.0, 1.0, -1.0],
        com + [1.0, 1.0, 1.0],
        com,
    ]
    cube_tets = [
        1 4 1 1 1 1 8 4 8 7 4 8
        2 2 2 5 3 5 4 3 7 6 8 6
        3 3 6 6 7 7 7 7 6 5 2 2
        9 9 9 9 9 9 9 9 9 9 9 9
    ]
    cube_pots = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    cube = Mesh(
        verts = cube_verts,
        tets = cube_tets,
        potentials = cube_pots,
    )
    cube
end

