struct Mesh
    n::Int32 # number of vertices
    m::Int32 # number of tets
    verts::Matrix{Float64} # shape: [3, n]
    tets::Matrix{Int64} # shape: [4, m]

    Mesh(verts, tets) = begin
        @assert size(verts, 1) == 3, "verts should be in R3"
        @assert size(tets, 1) == 4, "tets should be 4-tuples of indices"
        n, m = size(verts, 2), size(tets, 2)
        new(n, m, verts, tets)
    end
end

mutable struct Object
    mesh::Mesh
    transform::Transform3D
    frame::CartesianFrame3D
end
