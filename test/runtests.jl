using Test
using Hydroelastics
using LinearAlgebra
using StaticArrays

@testset "tet intersection" begin
    tet1 = Object(
        Mesh(
            [
                0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0
            ],
            reshape([1, 2, 3, 4], :, 1),
            [0.0, 0.0, 0.0, 1.0],
        ),
    )
    tet2 = Object(
        Mesh(
            [
                0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                1.8 1.8 1.8 0.8
            ],
            reshape([1, 2, 3, 4], :, 1),
            [0.0, 0.0, 0.0, 1.0],
        ),
    )

    points = Hydroelastics.intersect_tets(tet1, tet2, 1, 1)
    expected_points = [[0.0, 0.0, 0.9], [0.1, 0.0, 0.9], [0.0, 0.1, 0.9]]

    for p in eachcol(points)
        # there must exist a j in res which is close
        found = false
        for q in expected_points
            if norm(p - q) < 1e-4
                found = true
            end
        end
        @test found
    end
end

@testset "center of mass" begin
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

    @test norm(tet.com - [1.0, 2.0, 3.0]) < 1e-6
end

@testset "triangulation" begin
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
    @test Hydroelastics.triangulate_polygon(polygon) == result
end


@testset "mesh forces" begin
    tet1 = Object(
        Mesh(
            [
                0.0 1 0 0
                0 0 1 0
                0 0 0 1
            ],
            reshape([1, 2, 3, 4], 4, 1),
            [0.0, 0.0, 0.0, 1.0],
        ),
    )
    tet2 = Object(
        Mesh(
            [
                0 1 0 0
                0 0 1 0
                1.8 1.8 1.8 0.8
            ],
            reshape([1, 2, 3, 4], 4, 1),
            [0.0, 0.0, 0.0, 1.0],
        ),
    )

    force = Hydroelastics.tet_force(tet1, tet2, 1, 1)
    # from prev test we know the intersection is a triangle (0,0,.9), (.1, 0, .9), (0, .1, .9)
    # so the weights on the vtxs of A should be .033, .033, .033, .9
    expected_force = [0, 0, -0.0045]
    @test norm(force - expected_force) < 1e-6
    #cu1 = make_cube([0.0, 0.0, 0.0])
    #cu2 = make_cube([0.39103, 0.0232, 0.4312])
    #force_cubes = compute_force(cu1, cu2)
    #@test norm(force_cubes - [-1.1401553, -0.6682203, -0.3593239]) < 1e-6
end

@testset "icosphere volume" begin
    sphere = make_icosphere(5)
    @test abs(volume(sphere) - 4pi / 3) < 0.0025
end

@testset "object rotations" begin
    tet = Mesh(
        [
            0.0 1.0 0.0 0.0
            0.0 0.0 1.0 0.0
            0.0 0.0 0.0 1.0
        ],
        reshape([1, 2, 3, 4], :, 1),
        [0.0, 0.0, 0.0, 1.0],
    )
    obj = Object(tet, @SMatrix [
        1 0 0 1
        0 1 0 1
        0 0 1 1
        0 0 0 1
    ])
    new_obj = rotateX(obj, pi / 4)
    @test norm(new_obj.pose[1:3, 1:3] - @SMatrix [
        1 0 0
        0 0.7071 -.7071
        0 0.7071 0.7071
    ]) < 1e-2
    # TODO: add tests for rotateY, rotateZ, position checking
end
