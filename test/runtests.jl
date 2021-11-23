using Test
using LinearAlgebra
using Hydroelastics

function test_isect_tet()
    tet1 = Mesh(
        [[0.0, 0.0, 0.0] [1.0, 0.0, 0.0] [0.0, 1.0, 0.0] [0.0, 0.0, 1.0]],
        reshape([1, 2, 3, 4], :, 1),
        [0.0, 0.0, 0.0, 1.0],
    )
    tet2 = Mesh(
        [[0.0, 0.0, 2 - 0.2] [1.0, 0.0, 2 - 0.2] [0.0, 1.0, 2 - 0.2] [0.0, 0.0, 1 - 0.2]],
        reshape([1, 2, 3, 4], :, 1),
        [0.0, 0.0, 0.0, 1.0],
    )

    points = intersect_tets(tet1, tet2, 1, 1)
    expected_points = [[0.0, 0.0, 0.9], [0.1, 0.0, 0.9], [0.0, 0.1, 0.9]]

    for i = 1:size(points, 2)
        # there must exist a j in res which is close
        found = false
        for j in expected_points
            if norm(points[:, i] - j) < 1e-4
                found = true
                break
            end
        end
        if !found
            return false
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
