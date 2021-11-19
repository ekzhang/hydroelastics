import numpy as np
from taichi import ti

ti.init(arch=ti.cpu)

from .object import triangulate_polygon


def test_triangulation():
    # start with (0, 0), (1, 0), (1, 1), (0, 2), (-2, 2), (-2, 1)
    # map x, y -> 2*x+1, x-y-2, 2*y-x
    # then shuffle vertices; old:new mapping is {0:0, 1:4, 2:1, 3:3, 4:5, 5:2}
    polygon = np.array(
        [(1.0, -2, 0), (3, -2, 1), (-3, -5, 4), (1, -4, 4), (3, -1, -1), (-3, -6, 6)]
    )
    result = [(0, 2, 5), (0, 5, 3), (0, 3, 1), (0, 1, 4)]
    assert triangulate_polygon(polygon) == result


test_triangulation()
