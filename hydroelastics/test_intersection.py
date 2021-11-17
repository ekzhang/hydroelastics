import os

# import the object class from object.py
import taichi as ti
import numpy as np
from typing import List

ti.init(arch=ti.cpu)

from hydroelastics.object import Object, intersect


def test_isect():
    tet1 = Object(
        verts=[
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
        ],
        potentials=[0.0, 0.0, 0.0, 1.0],
        tets=[(0, 1, 2, 3)],
        mass=1,
    )
    tet2 = Object(
        verts=[
            (0, 0, 2 - 0.2),
            (1, 0, 2 - 0.2),
            (0, 1, 2 - 0.2),
            (0, 0, 1 - 0.2),
        ],
        potentials=[0.0, 0.0, 0.0, 1.0],
        tets=[(0, 1, 2, 3)],
        mass=1,
    )

    final_res = intersect(tet1, tet2, 0, 0)
    print(final_res)
    res = [
        np.array([0.0, 0.0, 0.9]),
        np.array([0.1, 0.0, 0.9]),
        np.array([0.0, 0.1, 0.9]),
    ]
    for i in final_res:
        # there must exist a j in res which is close
        close = False
        for j in res:
            if np.linalg.norm(i - j) < 1e-4:
                close = True
                break
        assert close


test_isect()
