using Test
using LinearAlgebra
using Hydroelastics

function test_isect_tet()
    tet1 = Mesh(
        [
            0.0 1.0 0.0 0.0
            0.0 0.0 1.0 0.0
            0.0 0.0 0.0 1.0
        ],
        reshape([1, 2, 3, 4], :, 1),
        [0.0, 0.0, 0.0, 1.0],
    )
    tet2 = Mesh(
        [
            0.0 1.0 0.0 0.0
            0.0 0.0 1.0 0.0
            1.8 1.8 1.8 0.8
        ],
        reshape([1, 2, 3, 4], :, 1),
        [0.0, 0.0, 0.0, 1.0],
    )

    points = intersect_tets(tet1, tet2, 1, 1)
    expected_points = [[0.0, 0.0, 0.9], [0.1, 0.0, 0.9], [0.0, 0.1, 0.9]]

    for p in eachcol(points)
        # there must exist a j in res which is close
        found = false
        for q in expected_points
            if norm(p - q) < 1e-4
                found = true
                break
            end
        end
        if !found
            error("could not find returned point $p in $expected_points")
        end
    end

    true
end

@test test_isect_tet()

function test_com()
    tet = Mesh(
        [
            0.0 3 0 0 2
            0 0 0 6 4
            0 0 9 0 6
        ],
        [
            1 2
            2 3
            3 4
            4 5
        ],
        [0.0, 0.0, 0.0, 0.0, 0.0],
    )
    @test norm([1.0, 2.0, 3.0] - tet.com) < 1e-6
end

test_com()

function test_triangulation()
    # start with (0, 0), (1, 0), (1, 1), (0, 2), (-2, 2), (-2, 1)
    # map x, y -> 2*x+1, x-y-2, 2*y-x
    # then shuffle vertices; old:new mapping is {0:0, 1:4, 2:1, 3:3, 4:5, 5:2}
    polygon = [
        1.0 3 -3 1 3 -3
        -2 -2 -5 -4 -1 -6
        0 1 4 4 -1 6
    ]
    result = [
        1 1 1 1
        3 6 4 2
        6 4 2 5
    ]
    @test triangulate_polygon(polygon) == result
end

test_triangulation()


function test_pressure()
    tet1 = Mesh(
        [
            0.0 1 0 0
            0 0 1 0
            0 0 0 1
        ],
        reshape([1, 2, 3, 4], 4, 1),
        [0.0, 0.0, 0.0, 1.0],
        #mass=1,
    )
    tet2 = Mesh(
        [
            0 1 0 0
            0 0 1 0
            1.8 1.8 1.8 0.8
        ],
        reshape([1, 2, 3, 4], 4, 1),
        [0.0, 0.0, 0.0, 1.0],
        #mass=1,
    )

    final_res = pressure(tet1, tet2, 1, 1)
    # from prev test we know the intersection is a triangle (0,0,.9), (.1, 0, .9), (0, .1, .9)
    # so the weights on the vtxs of A should be .033, .033, .033, .9
    @test abs(final_res - 0.9) < 0.1^6
end

test_pressure()
