lib = DefaultLibrary{Float64}(diff_optimizer(GLPK.Optimizer))

function sort_polygon(polygonA::Matrix{Float64})
    # sort 2xN polygon in counter-clockwise order about center of mass
    n = size(polygonA)[2]
    com = [sum(polygonA[i, :]) for i = 1:2] / n
    angles = []  # list of (angle, vtx index corresponding to angle)
    displacements = polygonA - reshape(repeat(com, n), 2, n)
    for ind = 1:n
        disp = displacements[:, ind]
        angle = -atan(disp[2], disp[1]) # negative for counterclockwise
        push!(angles, angle)
    end
    order = sortperm(angles)
    hcat([polygonA[:, order[i]] for i = 1:n]...)::Matrix{Float64}
end

function nextidx(idx::Int64, n::Int64)
    if idx == n
        return 1
    else
        return idx + 1
    end
end

struct Point
    x::Float64
    y::Float64
end

function Base.:+(x::Point, y::Point)
    Point(x.x + y.x, x.y + y.y)
end

function Base.:-(x::Point, y::Point)
    Point(x.x - y.x, x.y - y.y)
end

function Base.:*(x::Point, y::Float64)
    Point(x.x * y, x.y * y)
end

function point_cross(x::Point, y::Point)
    x.x * y.y - x.y * y.x
end


struct HalfPlane
    p::Point
    pq::Point
    angle::Float64

    HalfPlane(p::Point, q::Point) = begin
        p = p
        pq = Point(q.x - p.x, q.y - p.y)
        angle = atan(pq.y, pq.x)
        new(p, pq, angle)
    end
end

function out(h::HalfPlane, r::Point)
    return point_cross(h.pq, r - h.p) > 1e-6
end

function Base.isless(x::HalfPlane, y::HalfPlane)
    if (abs(x.angle - y.angle) < 1e-6)
        return point_cross(x.pq, y.p - x.p) > 0
    end
    return x.angle < y.angle
end

function Base.isequal(x::HalfPlane, y::HalfPlane)
    if (abs(x.angle - y.angle) < 1e-6)
        return true
    else
        return false
    end
end

function intersect_halfplanes(s::HalfPlane, t::HalfPlane)
    alpha = point_cross((t.p - s.p), t.pq) / point_cross(s.pq, t.pq)
    return s.p + s.pq * alpha
end

function cust_unique(s::Vector{HalfPlane})
    res = Vector{HalfPlane}()
    for i = 1:size(s)[1]
        if size(res)[1] > 0
            if abs(res[size(res)[1]].angle - s[i].angle) > 1e-6
                push!(res, s[i])
            end
        else
            push!(res, s[i])
        end
    end
    res
end

function intersect_polygons(polygonA::Matrix{Float64}, polygonB::Matrix{Float64})
    #polygon A is a 2xN matrix, polygon B is a 2xM matrix. 
    #returns the intersection of polygon A and polygon B
    #if there is no intersection, returns an empty matrix
    box = [
        1e7 -1e7 -1e7 1e7
        1e7 1e7 -1e7 -1e7
    ]
    n, m = size(polygonA, 2), size(polygonB, 2)
    polygonA, polygonB = sort_polygon(polygonA), sort_polygon(polygonB)
    box = sort_polygon(box)
    halfplanes::Array{HalfPlane} = []
    for i = 1:n
        pt1 = Point(polygonA[1, i], polygonA[2, i])
        pt2 = Point(polygonA[1, nextidx(i, n)], polygonA[2, nextidx(i, n)])
        s = HalfPlane(pt1, pt2)
        push!(halfplanes, s)
    end
    for i = 1:m
        t = HalfPlane(
            Point(polygonB[1, i], polygonB[2, i]),
            Point(polygonB[1, nextidx(i, m)], polygonB[2, nextidx(i, m)]),
        )
        push!(halfplanes, t)
    end
    for i = 1:4
        t = HalfPlane(
            Point(box[1, i], box[2, i]),
            Point(box[1, nextidx(i, 4)], box[2, nextidx(i, 4)]),
        )
        push!(halfplanes, t)
    end
    sort!(halfplanes)
    halfplanes = cust_unique(halfplanes)
    dq = Vector{HalfPlane}()
    len = 0
    for i = 1:length(halfplanes)
        while ((len > 1) && out(halfplanes[i], intersect_halfplanes(dq[len], dq[len-1])))
            pop!(dq)
            len = len - 1
        end

        while ((len > 1) && out(halfplanes[i], intersect_halfplanes(dq[1], dq[2])))
            popfirst!(dq)
            len = len - 1
        end
        push!(dq, halfplanes[i])
        len = len + 1
    end

    while (len > 2 && out(dq[1], intersect_halfplanes(dq[len], dq[len-1])))
        pop!(dq)
        len = len - 1
    end
    while (len > 2 && out(dq[len], intersect_halfplanes(dq[1], dq[2])))
        popfirst!(dq)
        len = len - 1
    end

    if (len < 3)
        return Vector{Point}()
    end
    result = Vector{Point}()
    for i = 1:len
        push!(result, intersect_halfplanes(dq[i], dq[nextidx(i, len)]))
    end
    return result
end

function convertPoints(points::Vector{Point})
    if size(points)[1] == 0
        return Array{Float64}(undef, 0, 2)
    end
    hcat([[points[i].x; points[i].y] for i = 1:(size(points)[1])]...)
end

"""
Finds intersection between two tetrahedra
"""
function intersect_tets(o1::Object, o2::Object, a_face_idx::Int64, b_face_idx::Int64)
    # Usage: Meshes m1, m2, followed by a_fact_idx, b_face_idx.
    m1 = o1.mesh
    m2 = o2.mesh
    coords_A = transform(m1.verts[:, m1.tets[1:4, a_face_idx]], o1.pose) # 3 x 4 matrix
    coords_B = transform(m2.verts[:, m2.tets[1:4, b_face_idx]], o2.pose) # 3 x 4 matrix

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

    """
    Intersects Equation with plane
    """
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
    all_points = convertPoints(intersect_polygons(Apoly, Bpoly))

    # old polygon intersection code
    #PA = polyhedron(vrep(Apoly'), lib)
    #PB = polyhedron(vrep(Bpoly'), lib)
    #res = polyhedron(vrep(intersect(PA, PB)))
    #all_points = hcat(points(res.vrep)...)
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
function tet_force(
    A::Object,
    B::Object,
    i::Int64,
    j::Int64,
)::Tuple{Vector{Float64},Vector{Float64}}
    normal = zeros(3)
    intersection_polygon = intersect_tets(A, B, i, j)
    total_force = 0.0
    intersection_com = zeros(3)
    total_area = 0.0
    if (size(intersection_polygon)[1] > 0) && (size(intersection_polygon)[2] > 0)
        triangles = triangulate_polygon(intersection_polygon)
        if isempty(triangles)
            return zeros(3)
        end
        vtx_inds = A.mesh.tets[:, i]
        vtx_coords = transform(A.mesh.verts[:, vtx_inds], A.pose)
        vtx_coords = vcat(vtx_coords, ones(1, 4))  # 4x4 matrix of vtxs padded w 1s
        for xyz in eachcol(triangles)
            # calculate barycentric coordinates of each point in the triangle.
            vtxs = intersection_polygon[:, xyz]
            com = push!(mean(eachcol(vtxs)), 1)
            res = vtx_coords \ com
            pressure = sum(res .* A.mesh.potentials[vtx_inds])
            area = 0.5 * norm(cross(vtxs[:, 1] - vtxs[:, 2], vtxs[:, 1] - vtxs[:, 3]))
            total_force += pressure * area
            total_area += area
            intersection_com += area * com[1:3]
        end
        normal = cross(
            intersection_polygon[:, 1] - intersection_polygon[:, 2],
            intersection_polygon[:, 1] - intersection_polygon[:, 3],
        )
        normal = normal / norm(normal)
        if dot(normal, transform(A.mesh.com, A.pose) - transform(B.mesh.com, B.pose)) < 0
            normal = -1 * normal
        end
        intersection_com /= total_area
    end
    (total_force * normal, intersection_com)
end

"""
Computes force on object A due to contact with object B, along with the net
torques on both objects.
"""
function compute_force(A::Object, B::Object)::ForceResult
    force = zeros(3)
    τ_AB = zeros(3)
    τ_BA = zeros(3)

    for i = 1:A.mesh.m
        for j = 1:B.mesh.m
            tets_result = tet_force(A, B, i, j)
            force += tets_result[1]
            τ_AB += cross(tets_result[2] - A.mesh.com, tets_result[1])
            τ_BA += cross(tets_result[2] - B.mesh.com, tets_result[1])
        end
    end
    ForceResult(force, -1 * force, τ_AB, τ_BA)
end

export compute_force, intersect_polygons, convertPoints
