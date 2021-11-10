import taichi as ti
import numpy as np
from typing import List, Tuple
import os


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

    # @ti.kernel
    def draw(self):
        # TODO: Render this object with rasterization
        raise NotImplementedError()


def center_of_mass(verts, tets):
    """
    compute center of mass of an object given vertex positions and tet shape
    """
    tet_centers = []
    vols = []
    for tet in tets:
        tet_vtxs = np.array([verts[i] for i in tet])
        tet_centers.append(np.mean(tet_vtxs, axis=0))
        tet_vtxs = np.hstack((tet_vtxs, np.ones((4, 1))))
        vols.append(np.abs(np.linalg.det(tet_vtxs) / 6))
    return np.average(np.array(tet_centers), weights=np.array(vols), axis=0)


def intersect(
    a, b, a_face_idx: int, b_face_idx: int
) -> Tuple[ti.Vector.field(3, dtype=ti.f32, shape=(6,)), int]:
    """
    return a list of intersections between tetrahedron a_face_idx of a and
    of tetrahedron b_face_idx of b. This finds the intersection of tetrahedon A's edges
    with the equipotential surface.
    """
    self_coords = ti.Vector.field(3, dtype=ti.f32, shape=(4,))
    for i in range(4):
        self_coords[i] = a.vertices[a.tets[a_face_idx][i]]
    print("self_coords", self_coords)

    # if a.pose is not None:
    #    self_coords = a.pose(self_coords)

    def get_equations(vertices, tets, potentials, pose, face_idx):
        self_coords = ti.Vector.field(3, dtype=ti.f32, shape=(4,))
        for i in range(4):
            self_coords[i] = vertices[tets[face_idx][i]]
        print("self_coords", self_coords)
        coord_matrix = ti.Matrix.field(4, 4, dtype=ti.f32, shape=())
        for i in range(4):
            for j in range(3):
                coord_matrix[None][i, j] = self_coords[i][j]
            coord_matrix[None][i, 3] = 1
        print("coord_matrix", coord_matrix)
        coords = coord_matrix[None].to_numpy()
        # print("coords", x[None].inverse())
        potential_vector = ti.Vector.field(4, dtype=ti.f32, shape=())
        for i in range(4):
            potential_vector[None][i] = potentials[tets[face_idx][i]]
        # equation = coord_matrix[None].inverse() @ potential_vector[None]
        equation = np.linalg.solve(coords, potential_vector.to_numpy())
        print("equation", equation)
        equation_ti = ti.field(dtype=ti.f32, shape=(4,))
        equation_ti.from_numpy(equation)
        print("equation_ti", equation_ti)

        return equation_ti

    a_pot = get_equations(a.vertices, a.tets, a.potentials, a.pose, a_face_idx)
    b_pot = get_equations(b.vertices, b.tets, b.potentials, b.pose, b_face_idx)
    intersection_points = []
    intersection_ti = ti.Vector.field(4, dtype=ti.f32, shape=())
    for i in range(4):
        intersection_ti[None][i] = a_pot[i] - b_pot[i]
    intersection = intersection_ti.to_numpy()  # intersection \cdot x = 0
    if (
        intersection_ti[None][0] ** 2
        + intersection_ti[None][1] ** 2
        + intersection_ti[None][2] ** 2
        < 1e-6
    ):  # something degenerate -- potential functions are linear shifts.
        return ti.Vector.field(3, dtype=ti.f32, shape=(6,)), 0

    i_value = ti.field(dtype=ti.f32, shape=(4,))  # the inner products of each vertex
    for i in range(4):
        i_value[i] = 0
        for j in range(3):
            i_value[i] += intersection_ti[None][j] * self_coords[i][j]
        i_value[i] += intersection_ti[None][3]
        # i_value[i] = np.dot(intersection, np.append(self_coords[i], 1.0))
    intersection_points = ti.Vector.field(
        3, dtype=ti.f32, shape=(6,)
    )  # all intersection points
    num_intersection_points = 0
    for i in range(4):
        for j in range(i, 4):
            if i_value[i] * i_value[j] < 0:  # must be different signs
                frac = i_value[i] / (i_value[i] - i_value[j])  # must be between 0 and 1
                for k in range(3):
                    intersection_points[num_intersection_points][k] = (
                        1 - frac
                    ) * self_coords[i][k] + frac * self_coords[j][k]
                num_intersection_points += 1
    return intersection_points, num_intersection_points
