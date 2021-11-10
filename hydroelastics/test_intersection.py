import os

# import the object class from object.py
import taichi as ti
import numpy as np
from typing import List

ti.init(arch=ti.cpu)

from object import Object, intersect


def test_isect():
    os.chdir(os.path.dirname(os.path.realpath(__file__)))
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

    pts = intersect(tet1, tet2, 0, 0)
    res = [
        np.array([0.0, 0.0, 0.9]),
        np.array([0.1, 0.0, 0.9]),
        np.array([0.0, 0.1, 0.9]),
    ]
    assert len(pts) == len(res)
    for i in range(len(pts)):
        assert np.linalg.norm(pts[i] - res[i]) < 1e-6


test_isect()