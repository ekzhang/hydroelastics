# This package contains implementations of common shapes, which can be used to
# construct `Mesh` types from those shapes.

"""
Utility function for constructing an icosphere triangle mesh.

Adapted from https://observablehq.com/@mourner/fast-icosphere-mesh
"""
function make_icosphere_mesh(order::Int64)
    # set up a 20-triangle icosahedron
    f = (1 + 5^0.5) / 2
    T = 4^order

    vertices = zeros((10 * T + 2) * 3)
    vertices[1:36] = vec(
        [
            -1 f 0
            1 f 0
            -1 -f 0
            1 -f 0
            0 -1 f
            0 1 f
            0 -1 -f
            0 1 -f
            f 0 -1
            f 0 1
            -f 0 -1
            -f 0 1
        ]',
    )
    triangles = vec(
        [
            1 12 6
            1 6 2
            1 2 8
            1 8 11
            1 11 12
            12 11 3
            6 12 5
            2 6 10
            8 2 9
            11 8 7
            4 10 5
            4 5 3
            4 3 7
            4 7 9
            4 9 10
            10 9 2
            5 10 6
            3 5 12
            7 3 11
            9 7 8
        ]',
    )

    v = 13
    midCache = (order == 0) ? nothing : Dict() # midpoint vertices cache to avoid duplicating shared vertices

    function addMidPoint(a::Int64, b::Int64)
        key = floor((a + b - 2) * (a + b - 1) / 2) + min(a, b) # Cantor's pairing function
        if midCache !== nothing && haskey(midCache, key)
            i = midCache[key]
            delete!(midCache, key)
            return i
        end
        midCache[key] = v
        for k in [1, 2, 3]
            vertices[3*(v-1)+k] = (vertices[3*(a-1)+k] + vertices[3*(b-1)+k]) / 2
        end
        i = v
        v += 1
        return i
    end

    trianglesPrev = triangles
    for i = 1:order
        # subdivide each triangle into 4 triangles
        triangles = zeros(Int64, size(trianglesPrev)[1] * 4)
        for k = 0:3:(size(trianglesPrev)[1]-3)
            v1 = trianglesPrev[k+1]
            v2 = trianglesPrev[k+2]
            v3 = trianglesPrev[k+3]
            a = addMidPoint(v1, v2)
            b = addMidPoint(v2, v3)
            c = addMidPoint(v3, v1)
            t = k * 4
            triangles[t+1] = v1
            triangles[t+2] = a
            triangles[t+3] = c
            triangles[t+4] = v2
            triangles[t+5] = b
            triangles[t+6] = a
            triangles[t+7] = v3
            triangles[t+8] = c
            triangles[t+9] = b
            triangles[t+10] = a
            triangles[t+11] = b
            triangles[t+12] = c
            t += 12
        end
        trianglesPrev = triangles
    end
    # normalize vertices
    for i = 0:3:(size(vertices)[1]-3)
        m = 1 / (vertices[i+1]^2 + vertices[i+2]^2 + vertices[i+3]^2)^0.5
        vertices[i+1] *= m
        vertices[i+2] *= m
        vertices[i+3] *= m
    end
    reshape(vertices, 3, :), reshape(triangles, 3, :)
end

"""
Create an icosphere as a tetrahedral mesh, with potential 0 at the boundary.
"""
function make_icosphere(order::Int64)
    verts, tris = make_icosphere_mesh(order)
    points::Matrix{Float64} = [verts [0; 0; 0]] # add the origin
    num_points = size(points, 2)
    tets::Matrix{Int64} = [tris; repeat([num_points], 1, size(tris, 2))]
    @assert size(tets, 1) == 4 "sanity check tets are in the right format"
    potentials = zeros(Float64, num_points)
    potentials[end] = 1.0
    Object(Mesh(points, tets, potentials))
end

function make_cube(com = [0.0, 0.0, 0.0])
    """
    returns a 2x2x2 cube with center at the origin
    """
    cube_verts = [
        -1.0 -1.0 -1.0 -1.0 1.0 1.0 1.0 1.0 0.0
        -1.0 -1.0 1.0 1.0 -1.0 -1.0 1.0 1.0 0.0
        -1.0 1.0 -1.0 1.0 -1.0 1.0 -1.0 1.0 0.0
    ]
    cube_tets = [
        1 4 1 1 1 1 8 4 8 7 4 8
        2 2 2 5 3 5 4 3 7 6 8 6
        3 3 6 6 7 7 7 7 6 5 2 2
        9 9 9 9 9 9 9 9 9 9 9 9
    ]
    cube_pots = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    Object(
        Mesh(cube_verts, cube_tets, cube_pots),
        @SMatrix [
            1 0 0 com[1]
            0 1 0 com[2]
            0 0 1 com[3]
            0 0 0 1
        ]
    )
end

export make_icosphere, make_cube
