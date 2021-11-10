import taichi as ti
import numpy as np
from typing import List
import os


def center_of_mass(verts, tets):
    tet_centers = []
    vols = []
    for tet in tets:
        tet_vtxs = np.array([verts[i] for i in tet])
        tet_centers.append(np.mean(tet_vtxs, axis=0))
        tet_vtxs = np.hstack((tet_vtxs, np.ones((4, 1))))
        vols.append(np.abs(np.linalg.det(tet_vtxs) / 6))
    return np.average(np.array(tet_centers), weights=np.array(vols), axis=0)


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
        for i in range(self.n_vert):
            self.vertices[i] = verts[i]

        for i in range(self.n_tets):
            self.tets[i] = tets[i]

        for i in range(self.n_vert):
            self.potentials[i] = potentials[i]

        # TODO: Make pose not none
        self.pose = None

        # TODO: Initialize the fields with parameters
        # TODO: Calculate center of mass in identity pose and store in self.com

    # @ti.kernel
    def draw(self):
        # TODO: Render this object with rasterization
        raise NotImplementedError()


def intersect(a, b, a_face_idx: int, b_face_idx: int) -> List[ti.Vector]:
    """
    return a list of intersections between tetrahedron a_face_idx of a and
    of tetrahedron b_face_idx of b. This finds the intersection of tetrahedon A's edges
    with the equipotential surface.
    """
    self_coords = np.array([a.vertices[x] for x in a.tets[a_face_idx]])
    if a.pose is not None:
        self_coords = a.pose(self_coords)

    def get_equations(vertices, faces, potentials, pose, face_idx):
        coords = np.array([vertices[x] for x in faces[face_idx]])
        if pose != None:
            coords = pose(coords)
        coords = np.append(coords, np.ones(4).reshape(4, -1), axis=1)
        equation = np.linalg.solve(
            coords, np.array([potentials[x] for x in faces[face_idx]])
        )
        return equation

    a_pot = get_equations(a.vertices, a.tets, a.potentials, a.pose, a_face_idx)
    b_pot = get_equations(b.vertices, b.tets, b.potentials, b.pose, b_face_idx)
    intersection_points = []
    intersection = a_pot - b_pot  # intersection \cdot x = 0
    if (
        np.linalg.norm(intersection[:3]) < 1e-6
    ):  # something degenerate -- potential functions are linear shifts.
        return []

    i_value = ti.field(dtype=ti.f32, shape=(4,))  # the inner products of each vertex
    for i in range(4):
        i_value[i] = np.dot(intersection, np.append(self_coords[i], 1.0))

    for i in range(4):
        for j in range(i, 4):
            if i_value[i] * i_value[j] < 0:  # must be different signs
                print(i_value[i], i_value[j])
                frac = i_value[i] / (i_value[i] - i_value[j])  # must be between 0 and 1
                intersection = (1 - frac) * self_coords[i] + frac * self_coords[j]
                intersection_points.append(intersection)
    return intersection_points


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
