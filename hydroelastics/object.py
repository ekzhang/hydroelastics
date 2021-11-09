import taichi as ti
import numpy as np


@ti.data_oriented
class Object:
    def __init__(self, verts, potentials, tets, mass: float):
        assert len(verts) == len(potentials)
        self.n_vert = len(verts)
        self.n_tets = len(tets)
        self.mass = mass
        self.vertices = ti.Vector.field(3, dtype=ti.f32, shape=(self.n_vert,))
        self.potentials = ti.field(dtype=ti.f32, shape=(self.n_vert,))
        self.tets = ti.Vector.field(4, dtype=ti.i32, shape=(self.n_tets,))
        self.pose = ti.Matrix.field(4, 4, dtype=ti.f32, shape=())
        self.com = center_of_mass(verts, tets)

        # TODO: Initialize the fields with parameters
        # TODO: Calculate center of mass in identity pose and store in self.com

    # @ti.kernel
    def draw(self):
        # TODO: Render this object with rasterization
        raise NotImplementedError()


def center_of_mass(verts, tets):
    tet_centers = []
    vols = []
    for tet in tets:
        tet_vtxs = np.array([verts[i] for i in tet])
        tet_centers.append(np.mean(tet_vtxs, axis=0))
        tet_vtxs = np.hstack((tet_vtxs, np.ones((4, 1))))
        vols.append(np.abs(np.linalg.det(tet_vtxs) / 6))
    return np.average(np.array(tet_centers), weights=np.array(vols), axis=0)


"""
sanity check for center of mass. feel free to move/delete this; i wasn't sure where to put it
verts = [[0, 0, 0], 
[3, 0, 0], 
[0, 0, 9], 
[0, 6, 0], 
[2, 4, 6]
]
tets = [
    (0, 1, 2, 3), 
    (1, 2, 3, 4)
]
assert(np.all(center_of_mass(verts, tets) == np.array([1., 2., 3.])))
"""
