import os

# import the object class from object.py
import taichi as ti
import numpy as np
from typing import List

ti.init(arch=ti.cpu)

from hydroelastics.object import Object, pressure


def test_pressure():
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

    final_res = pressure(tet1, tet2, 0, 0)
    # from prev test we know the intersection is a triangle (0,0,.9), (.1, 0, .9), (0, .1, .9)
    # so the weights on the vtxs of A should be .033, .033, .033, .9
    print(final_res)
    assert (final_res - 0.9) ** 2 < 1e-12


test_pressure()
