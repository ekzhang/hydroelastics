import os

# import the object class from object.py
import taichi as ti
import numpy as np
from typing import List

ti.init(arch=ti.cpu)

from hydroelastics.object import Object


def test_com():
    tet = Object(
        verts=[[0, 0, 0], [3, 0, 0], [0, 0, 9], [0, 6, 0], [2, 4, 6]],
        potentials=[0.0, 0.0, 0.0, 0.0, 0.0],
        tets=[(0, 1, 2, 3), (1, 2, 3, 4)],
        mass=1,
    )
    com = tet.com
    res = np.array([1.0, 2.0, 3.0])
    for i in range(res.shape[0]):
        assert np.abs(res[i] - com[i] < 1e-6)


test_com()
