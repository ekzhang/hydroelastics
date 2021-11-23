using Test
using Hydroelastics
using LinearAlgebra

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
    println(final_res)
    res = [[0.0, 0.0, 0.9], [0.1, 0.0, 0.9], [0.0, 0.1, 0.9]]
    works = true
    for i = 1:size(final_res, 2)
        # there must exist a j in res which is close
        found = false
        for j in res
            println(j)
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


@test 1 == 1
@test test_isect_tet()
