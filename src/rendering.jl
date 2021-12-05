"""
Transform an object into a format that is accepted by MeshCat.
"""
function render_obj(obj::Object)::MeshFileGeometry
    points = [(obj.pose*[v; 1])[1:3] for v in eachcol(obj.mesh.verts)]
    lines = ["v " * join(p, " ") for p in points]
    for (a, b, c, d) in eachcol(obj.mesh.tets)
        for (x, y, z) in [(a, b, c), (b, c, d), (c, d, a), (d, a, b)]
            # Only include faces with potential 1, to get the outer surface.
            if obj.mesh.potentials[[x, y, z]] ≈ [0, 0, 0]
                push!(lines, "f $x $y $z")
            end
        end
    end
    MeshFileGeometry(join(lines, "\n"), "obj")
end

export render_obj
