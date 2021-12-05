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

struct Object
    mesh::Mesh
    pose::SMatrix{4,4}

    """Construct a new object from a mesh, with default pose."""
    Object(mesh::Mesh) = new(mesh, @SMatrix [
        1 0 0 0
        0 1 0 0
        0 0 1 0
        0 0 0 1
    ])

    """Construct a new object from a mesh, with a specified pose."""
    Object(mesh::Mesh, pose::SMatrix{4,4}) = new(mesh, pose)
end

"""Compute the volume of an object."""
function volume(obj::Object)
    volume(obj.mesh)
end

"""Applies transform to an object."""
function transform(obj::Object, matrix::SMatrix{4,4})::Object
    Object(obj.mesh, matrix * obj.pose)
end

"""
Applies transform to vertices.
"""
function transform(vertices::Matrix{Float64}, matrix::SMatrix{4,4})
    # concatenate ones to vertices
    ones_arr = ones(1, size(vertices, 2))
    mat = vcat(vertices, ones_arr)
    vertices = matrix * mat
    return vertices[1:3, :]
end

"""
Applies transform to a single vertex.
"""
function transform(vertex::Vector{Float64}, matrix::SMatrix{4,4})
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

function rotateX(obj::Object, angle::Float64)::Object
    transform(
        obj,
        @SMatrix [
            1 0 0 0
            0 cos(angle) -sin(angle) 0
            0 sin(angle) cos(angle) 0
            0 0 0 1
        ]
    )
end

function rotateY(obj::Object, angle::Float64)::Object
    transform(
        obj,
        @SMatrix [
            cos(angle) 0 sin(angle) 0
            0 1 0 0
            -sin(angle) 0 cos(angle) 0
            0 0 0 1
        ]
    )
end

function rotateZ(obj::Object, angle::Float64)::Object
    transform(
        obj,
        @SMatrix [
            cos(angle) -sin(angle) 0 0
            sin(angle) cos(angle) 0 0
            0 0 1 0
            0 0 0 1
        ]
    )
end

export Mesh, Object, volume, transform, translate, rotateX, rotateY, rotateZ
