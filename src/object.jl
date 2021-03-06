"""
Wrapper type to reduce the visual noise of Pluto's `Base.show` on `RTree`.
"""
struct RTreeWrapper
    t::RTree{Float64,3}
end

"""
Adapted from https://github.com/alyst/SpatialIndexing.jl/blob/105b420e/src/rtree/show.jl.
"""
function Base.show(io::IO, rtree::RTreeWrapper)
    tree = rtree.t
    print(io, typeof(tree))
    print(io, "(variant=")
    print(io, tree.variant)
    print(io, ", tight_mbrs=")
    print(io, tree.tight_mbrs)
    print(io, ", nearmin_overlap=")
    print(io, tree.nearmin_overlap)
    print(io, ", fill_factor=")
    print(io, tree.fill_factor)
    print(io, ", split_factor=")
    print(io, tree.split_factor)
    print(io, ", reinsert_factor=")
    print(io, tree.reinsert_factor)

    print(io, ", leaf_capacity=")
    print(io, SI.capacity(SI.Leaf, tree))
    print(io, ", branch_capacity=")
    print(io, SI.capacity(SI.Branch, tree))
    println(io, ")")
    print(io, SI.length(tree))
    print(io, " element(s) in ")
    print(io, SI.height(tree))
    print(io, " level(s) (")
    print(io, join(reverse(tree.nnodes_perlevel), ", "))
    println(io, " node(s) per level)")
end

"""
Geometric data structure that stores a tetrahedral mesh with potentials.
"""
struct Mesh
    n::Int32 # number of vertices
    m::Int32 # number of tets
    verts::Matrix{Float64} # shape: [3, n]
    tets::Matrix{Int64} # shape: [4, m]
    potentials::Vector{Float64} # shape: [n]
    com::Vector{Float64} # center of mass
    rtree::RTreeWrapper # bounding box index data structure

    Mesh(verts::Matrix{Float64}, tets::Matrix{Int64}, potentials::Vector{Float64}) = begin
        @assert(size(verts, 1) == 3, "verts should be in R3")
        @assert(size(tets, 1) == 4, "tets should be 4-tuples of indices")
        @assert(
            size(potentials, 1) == size(verts, 2),
            "potentials should be the same size as verts"
        )
        n, m = size(verts, 2), size(tets, 2)
        com = center_of_mass(verts, tets)
        rtree = RTree{Float64,3}(Int64)
        for i = 1:m
            points = verts[:, tets[:, i]]
            bb_min, bb_max = bounding_box(points)
            rect = SI.Rect(
                (bb_min[1], bb_min[2], bb_min[3]),
                (bb_max[1], bb_max[2], bb_max[3]),
            )
            SI.insert!(rtree, rect, i)
        end
        new(n, m, verts, tets, potentials, com, RTreeWrapper(rtree))
    end
end

"""
Compute center of mass, given coords of vertices and tets.
"""
function center_of_mass(verts::Matrix{Float64}, tets::Matrix{Int64})
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

"""Compute the volume of a mesh."""
function volume(mesh::Mesh)
    values = mesh.verts[:, mesh.tets]
    volume = 0.0
    for tet_values in eachslice(values, dims = 3)
        a, b, c, d = eachcol(tet_values)
        # compute triple product
        volume += abs(dot(cross(b - a, c - a), d - a)) / 6.0
    end
    volume
end

"""Compute the lower and upper bounding box of a collection of points."""
function bounding_box(points::Matrix{Float64})
    min.(eachcol(points)...), max.(eachcol(points)...)
end

const Pose = SMatrix{4,4,Float64}

struct Object
    mesh::Mesh
    pose::Pose
    v::SVector{3,Float64}
    ??::SVector{3,Float64}

    """Construct a new object from a mesh, with optional pose and velocities."""
    Object(
        mesh::Mesh,
        pose::Pose = Pose(I),
        v::SVector{3,Float64} = (@SVector zeros(3)),
        ??::SVector{3,Float64} = (@SVector zeros(3)),
    ) = new(mesh, pose, v, ??)
end

"""Compute the volume of an object."""
function volume(obj::Object)
    volume(obj.mesh)
end

"""Applies transform to an object."""
function transform(obj::Object, matrix::Pose)::Object
    Object(obj.mesh, matrix * obj.pose, obj.v, obj.??)
end

"""
Applies transform to vertices.
"""
function transform(vertices::Matrix{Float64}, matrix::Pose)
    # concatenate ones to vertices
    ones_arr = ones(1, size(vertices, 2))
    mat = vcat(vertices, ones_arr)
    vertices = matrix * mat
    return vertices[1:3, :]
end

"""
Applies transform to a single vertex.
"""
function transform(vertex::Vector{Float64}, matrix::Pose)
    return (matrix*[vertex; 1.0])[1:3]
end

function translate(obj::Object, disp::SVector{3})::Object
    transform(obj, @SMatrix [
        1 0 0 disp[1]
        0 1 0 disp[2]
        0 0 1 disp[3]
        0 0 0 1
    ])
end

function rotateX(obj::Object, ??::Float64)::Object
    transform(obj, @SMatrix [
        1 0 0 0
        0 cos(??) -sin(??) 0
        0 sin(??) cos(??) 0
        0 0 0 1
    ])
end

function rotateY(obj::Object, ??::Float64)::Object
    transform(obj, @SMatrix [
        cos(??) 0 sin(??) 0
        0 1 0 0
        -sin(??) 0 cos(??) 0
        0 0 0 1
    ])
end

function rotateZ(obj::Object, ??::Float64)::Object
    transform(obj, @SMatrix [
        cos(??) -sin(??) 0 0
        sin(??) cos(??) 0 0
        0 0 1 0
        0 0 0 1
    ])
end

"""
Updates the current object's pose according to an angular velocity.

Uses the exponential map on SO(3), see https://arwilliams.github.io/so3-exp.pdf.
"""
function changeRotation(obj::Object, ??::SVector{3,Float64})::Object
    ?? = norm(??)
    if ?? < 1e-12
        return obj
    end

    # Compute the skew matrix.
    ??x = @SMatrix [
        0 -??[3] ??[2]
        ??[3] 0 -??[1]
        -??[2] ??[1] 0
    ]

    # Compute the matrix exponential / exponential map.
    ??_exp = I + (sin(??) / ??) .* ??x + ((1 - cos(??)) / ??^2) .* (??x * ??x)
    transform(obj, Pose([??_exp zeros(3); zeros(3)' 1]))
end

export Mesh,
    Pose, Object, volume, transform, translate, rotateX, rotateY, rotateZ, changeRotation
