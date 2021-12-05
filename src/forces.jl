const Point = SVector{2}

struct HalfPlane
    p::Point
    pq::Point
    angle::Float64

    HalfPlane(p::Point, pq::Point) = begin
        angle = atan(pq[2], pq[1])
        new(p, pq, angle)
    end
end

function out(h::HalfPlane, r::Point)
    cross(h.pq, r - h.p) > 1e-9
end

function Base.isless(x::HalfPlane, y::HalfPlane)
    if abs(x.angle - y.angle) < 1e-12
        return cross(x.pq, y.p - x.p) > 0
    end
    x.angle < y.angle
end

function isect(s::HalfPlane, t::HalfPlane)
    alpha = cross((t.p - s.p), t.pq) / cross(s.pq, t.pq)
    s.p + s.pq * alpha
end

function remove_duplicates(s::Vector{HalfPlane})
    res = HalfPlane[]
    for hp in s
        if length(res) > 0
            if abs(res[end].angle - hp.angle) > 1e-12
                push!(res, hp)
            end
        else
            push!(res, hp)
        end
    end
    res
end

"""
Computes the polygonal intersection of two-dimensional halfplanes. If there is
no intersection, returns `nothing`.

Requires that the halfplanes being intersected have finite area.
"""
function intersect_halfplanes(halfplanes::Vector{HalfPlane})
    sort!(halfplanes)
    halfplanes = remove_duplicates(halfplanes)
    dq = HalfPlane[]
    for halfplane in halfplanes
        while length(dq) >= 2 && out(halfplane, isect(dq[end], dq[end-1]))
            pop!(dq)
        end

        while length(dq) >= 2 && out(halfplane, isect(dq[1], dq[2]))
            popfirst!(dq)
        end
        push!(dq, halfplane)
    end

    while length(dq) >= 3 && out(dq[1], isect(dq[end], dq[end-1]))
        pop!(dq)
    end
    while length(dq) >= 3 && out(dq[end], isect(dq[1], dq[2]))
        popfirst!(dq)
    end

    len = length(dq)
    if len < 3
        return nothing
    end
    hcat([vec(isect(dq[i], dq[mod1(i + 1, len)])) for i = 1:len]...)
end

"""
Returns the equi-pressure intersection between two tetrahedra as a polygon, as
well as the normal vector to the plane of intersection.

Takes in objects `o1`, `o2`, followed by the indices of the tets in each object.
"""
function intersect_tets(o1::Object, o2::Object, a_face_idx::Int64, b_face_idx::Int64)
    m1 = o1.mesh
    m2 = o2.mesh
    coords_A = transform(m1.verts[:, m1.tets[:, a_face_idx]], o1.pose) # 3 x 4 matrix
    coords_B = transform(m2.verts[:, m2.tets[:, b_face_idx]], o2.pose) # 3 x 4 matrix

    # First, do an initial check for bounding boxes of tetrahedra
    if any(max.(eachcol(coords_A)...) .< min.(eachcol(coords_B)...)) ||
       any(max.(eachcol(coords_B)...) .< min.(eachcol(coords_A)...))
        return nothing
    end

    function get_equation(coords::Matrix{Float64}, potentials::Vector{Float64})
        ones_arr = ones(size(coords, 2), 1)
        mat = [coords' ones_arr]
        mat \ potentials
    end

    a_pot = get_equation(coords_A, m1.potentials[m1.tets[:, a_face_idx]])
    b_pot = get_equation(coords_B, m2.potentials[m2.tets[:, b_face_idx]])
    plane = a_pot - b_pot  # dot(plane, [x; 1.0]) == 0.0
    if norm(plane[1:3]) < 1e-8
        return nothing
    end

    # Find a 2x3 projection from the plane onto two dimensions, along with an
    # inverse 3x2 projection that has the following properties:
    #   1. proj_mat * inv_proj_mat = I
    #   2. plane[1:3]' * (inv_proj_mat * x + inv_proj_offset) + plane[4] = 0
    if abs(plane[1]) / norm(plane[1:3]) > 1e-3
        proj_mat = [0.0 1.0 0.0; 0.0 0.0 1.0]
        inv_proj_mat = [(-plane[2]/plane[1]) (-plane[3]/plane[1]); 1.0 0.0; 0.0 1.0]
        inv_proj_offset = [-plane[4] / plane[1], 0.0, 0.0]
    elseif abs(plane[2]) / norm(plane[1:3]) > 1e-3
        proj_mat = [1.0 0.0 0.0; 0.0 0.0 1.0]
        inv_proj_mat = [1.0 0.0; (-plane[1]/plane[2]) (-plane[3]/plane[2]); 0.0 1.0]
        inv_proj_offset = [0.0, -plane[4] / plane[2], 0.0]
    elseif abs(plane[3]) / norm(plane[1:3]) > 1e-3
        proj_mat = [1.0 0.0 0.0; 0.0 1.0 0.0]
        inv_proj_mat = [1.0 0.0; 0.0 1.0; (-plane[1]/plane[3]) (-plane[2]/plane[3])]
        inv_proj_offset = [0.0, 0.0, -plane[4] / plane[3]]
    else
        error("no valid projection found, aborting...")
    end

    # Convert each tetrahedron into four 2D halfplanes, on the projection of
    # the plane to two dimensions along some orthogonal axis.
    halfplanes = HalfPlane[]
    for coords in [coords_A, coords_B]
        for i = 1:4
            hs_p1 = Point(proj_mat * coords[:, mod1(i + 1, 4)])
            hs_p2 = Point(proj_mat * coords[:, mod1(i + 2, 4)])
            hs_p3 = Point(proj_mat * coords[:, mod1(i + 3, 4)])
            hs_int = Point(proj_mat * coords[:, i])

            if norm(hs_p1 - hs_p2) > 1e-9 &&
               abs(cross(hs_int - hs_p1, hs_p1 - hs_p2)) > 1e-9
                hs_a, hs_b = hs_p1, hs_p2
            elseif norm(hs_p2 - hs_p3) > 1e-9 &&
                   abs(cross(hs_int - hs_p1, hs_p1 - hs_p2)) > 1e-9
                hs_a, hs_b = hs_p2, hs_p3
            elseif norm(hs_p3 - hs_p1) > 1e-9 &&
                   abs(cross(hs_int - hs_p1, hs_p1 - hs_p2)) > 1e-9
                hs_a, hs_b = hs_p3, hs_p1
            else
                continue
            end

            if cross(hs_int - hs_a, hs_b - hs_a) < 0
                # Does this go in the correct direction?
                hs_a, hs_b = hs_b, hs_a
            end

            push!(halfplanes, HalfPlane(hs_a, normalize(hs_b - hs_a)))
        end
    end

    all_points = intersect_halfplanes(halfplanes)
    if isnothing(all_points)
        return nothing
    end

    @assert size(all_points, 2) >= 3 "sanity check length"
    hcat([inv_proj_mat * p + inv_proj_offset for p in eachcol(all_points)]...), plane[1:3]
end

function intersect_polygons(polygonA::Matrix{Float64}, polygonB::Matrix{Float64})
    n, m = size(polygonA, 2), size(polygonB, 2)
    halfplanes::Array{HalfPlane} = []
    for i = 1:n
        push!(halfplanes, HalfPlane(
            Point(polygonA[:, i]),
            Point(polygonA[:, mod1(i + 1, n)]) - Point(polygonA[:, i]),
        ))
    end
    for i = 1:m
        push!(halfplanes, HalfPlane(
            Point(polygonB[:, i]),
            Point(polygonB[:, mod1(i + 1, m)]) - Point(polygonB[:, i]),
        ))
    end
    intersect_halfplanes(halfplanes)
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
the location of force applied to A and B.
"""
function tet_force(
    A::Object,
    B::Object,
    i::Int64,
    j::Int64,
)::Union{Tuple{Vector{Float64},Vector{Float64}},Nothing}
    result = intersect_tets(A, B, i, j)
    if isnothing(result)
        return nothing
    end

    intersection_polygon, normal = result
    normalize!(normal)
    if dot(normal, transform(A.mesh.com, A.pose) - transform(B.mesh.com, B.pose)) < 0
        normal *= -1
    end

    total_force = 0.0
    intersection_com = zeros(3)
    total_area = 0.0

    vtx_inds = A.mesh.tets[:, i]
    vtx_coords = transform(A.mesh.verts[:, vtx_inds], A.pose)
    vtx_coords = vcat(vtx_coords, ones(1, 4))  # 4x4 matrix of vtxs padded w 1s
    for i = 3:size(intersection_polygon, 2)
        # calculate barycentric coordinates of each point in the triangle.
        vtxs = intersection_polygon[:, [1, i - 1, i]]
        com = push!(mean(eachcol(vtxs)), 1)
        println("vtxs: ", vtxs)
        println("com: ", com)
        res = vtx_coords \ com
        pressure = sum(res .* A.mesh.potentials[vtx_inds])
        area = 0.5 * norm(cross(vtxs[:, 1] - vtxs[:, 2], vtxs[:, 1] - vtxs[:, 3]))
        if pressure < 0 || area < 0
            println(res, " ", pressure, " ", area)
            error("bad")
        end
        total_force += pressure * area
        total_area += area
        intersection_com += area * com[1:3]
    end

    if total_area < 1e-9
        return nothing
    end

    intersection_com /= total_area
    total_force * normal, intersection_com
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
            if !isnothing(tets_result)
                force += tets_result[1]
                τ_AB += cross(tets_result[2] - A.mesh.com, tets_result[1])
                τ_BA += cross(tets_result[2] - B.mesh.com, tets_result[1])
            end
        end
    end
    ForceResult(force, -1 * force, τ_AB, τ_BA)
end

export compute_force
