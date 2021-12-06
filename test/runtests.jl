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

    @test !isnothing(points)

    for p in eachcol(points[1])
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

    force = (Hydroelastics.tet_force(tet1, tet2, 1, 1))[1]
    # from prev test we know the intersection is a triangle (0,0,.9), (.1, 0, .9), (0, .1, .9)
    # so the weights on the vtxs of A should be .033, .033, .033, .9
    expected_force = [0, 0, -0.0045]
    @test norm(force - expected_force) < 1e-6

    sphere = make_icosphere(1)
    object1 = translate(sphere, @SVector [0.031, -0.5, 0.052])
    object2 = translate(sphere, @SVector [0, 0.5, 0])
    forces = compute_force(object1, object2)

    @test norm(forces.F_AB - [0.012750, -0.444069, 0.025010]) < 1e-6
    @test norm(forces.F_AB + forces.F_BA) < 1e-12
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
end
