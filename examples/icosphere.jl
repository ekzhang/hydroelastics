# adapted from https://observablehq.com/@mourner/fast-icosphere-mesh
function icosphere(order::Int64)
    # set up a 20-triangle icosahedron
    f = (1 + 5^0.5) / 2
    T = 4^order

    vertices = zeros((10 * T + 2) * 3)
    vertices[1:36] = [
        -1,
        f,
        0,
        1,
        f,
        0,
        -1,
        -f,
        0,
        1,
        -f,
        0,
        0,
        -1,
        f,
        0,
        1,
        f,
        0,
        -1,
        -f,
        0,
        1,
        -f,
        f,
        0,
        -1,
        f,
        0,
        1,
        -f,
        0,
        -1,
        -f,
        0,
        1,
    ]
    triangles = [
        1,
        12,
        6,
        1,
        6,
        2,
        1,
        2,
        8,
        1,
        8,
        11,
        1,
        11,
        12,
        12,
        11,
        3,
        6,
        12,
        5,
        2,
        6,
        10,
        8,
        2,
        9,
        11,
        8,
        7,
        4,
        10,
        5,
        4,
        5,
        3,
        4,
        3,
        7,
        4,
        7,
        9,
        4,
        9,
        10,
        10,
        9,
        2,
        5,
        10,
        6,
        3,
        5,
        12,
        7,
        3,
        11,
        9,
        7,
        8,
    ]

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
    return (vertices, triangles)
end

print(icosphere(1))
