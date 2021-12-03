"""
Transform a custom Mesh object into a format that is accepted by MeshCat.
"""
function mesh_to_obj(mesh::Hydroelastics.Mesh)::MeshFileGeometry
    points = Point3f0.(eachcol(mesh.verts))
    lines = ["v " * join(p, " ") for p in points]
    for (a, b, c, d) in eachcol(mesh.tets)
        for (x, y, z) in [(a, b, c), (b, c, d), (c, d, a), (d, a, b)]
            # Only include faces with potential 1, to get the outer surface.
            if mesh.potentials[[x, y, z]] ≈ [1, 1, 1]
                push!(lines, "f $x $y $z")
            end
        end
    end
    MeshFileGeometry(join(lines, "\n"), "obj")
end

export mesh_to_obj
