using Test
using Hydroelastics

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
