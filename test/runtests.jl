using Test
using Hydroelastics
using LinearAlgebra

function isect_cubes()
    function get_cube(com)
        """
        returns a cube with center at com
        """
        cube_verts = com .+ [
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
        cube = Mesh(
            cube_verts,
            cube_tets,
            cube_pots,
        )
        cube
    end
    cu1 = get_cube([0.0, 0.0, 0.0])
    cu2 = get_cube([0.5, 0.5, 0.5])

    press = 0
    for i = 1:12
        for j = 1:12
            press += pressure(cu1, cu2, i, j)
        end
    end
    press
end


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

    final_res = intersect_tets(tet1, tet2, 1, 1)
    res = [[0.0, 0.0, 0.9], [0.1, 0.0, 0.9], [0.0, 0.1, 0.9]]
    works = true
    for i = 1:size(final_res, 2)
        # there must exist a j in res which is close
        found = false
        for j in res
            if norm(final_res[:, i] - j) < 1e-4
                found = true
            end
        end
        if !found
            global works = false
        end
    end
    works
end

@test isect_cubes()
@test 1 == 1
@test test_isect_tet()

@test 1 == 1

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
        #mass=1,
    )
    com = tet.com
    res = [1.0, 2.0, 3.0]
    for i = 1:3
        @test abs(res[i] - com[i] < 1e-6)
    end
end
test_com()
