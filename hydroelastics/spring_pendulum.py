import taichi as ti
from taichi_glsl import length, sqrt, cos, sin
import numpy as np

ti.init(arch=ti.gpu)

# Double Pendulum parameters
a = 1.0
b = 1.0
PI = 3.14159265358979323846
pos = ti.Vector.field(2, dtype=float, shape=(2,))
vel = ti.Vector.field(2, dtype=float, shape=(2,))
rest_lengths = ti.field(dtype=float, shape=(2,))
cur_lengths = ti.field(dtype=float, shape=(2,))
forces = ti.Vector.field(2, dtype=float, shape=(2,))
pivot = ti.Vector.field(2, dtype=float, shape=())

gravity = 9.8
k_spring = 1e3
pivot = ti.Vector([0.5, 0.9])


@ti.kernel
def initialize():
    for i in range(2):
        rest_lengths[i] = 0.1
        cur_lengths[i] = 0.1
        pos[i] = ti.Vector([0.5 + 0.1 * float(i + 1), 0.9])
        vel[i] = ti.Vector([0.0, 0.0])
    print("hello")


@ti.kernel
def apply_pendulum_force():
    """Fill in the `forces` field with spring forces between points."""
    for i in pos:
        forces[i] = ti.Vector([0.0, -gravity])
    #alpha = PI / 2 - ti.atan2(pos[0][1], pos[0][0])
    #beta = PI / 2 - ti.atan2(pos[1][1] - pos[0][1], pos[1][0] - pos[0][0])
    #print("alpha", alpha)
    #print("beta", beta)

    for i in range(2):
        if i == 0:
            # apply force from first spring
            d = pos[0] - pivot
            forces[i] -= (
                (length(d) - rest_lengths[i]) * k_spring * d / length(d)
            )
        d2 = pos[1] - pos[0]
        second_force = (
            length(d2) - rest_lengths[1]) * k_spring * d2 / length(d2)
        if i == 0:
            forces[i] += second_force
        else:
            forces[i] -= second_force


@ti.kernel
def advance(t: float):
    """Advance the simulation by `t` units of time."""
    for i in pos:
        dv = t * forces[i]
        v = vel[i] + dv
        vel[i] += dv
        pos[i] += t * v


def run():
    gui = ti.GUI("Double Pendulum", res=(600, 600))
    initialize()
    while gui.running:
        pos_np = pos.to_numpy()
        gui.circles(pivot.to_numpy().reshape((-1, 2)), radius=5)
        print(pos_np)
        gui.circles(pos_np.reshape((-1, 2)), radius=5)
        gui.lines(
            pivot.to_numpy().reshape((-1, 2)),
            pos_np[0].reshape((-1, 2)),
            radius=1,
        )
        gui.lines(
            pos_np[0].reshape((-1, 2)),
            pos_np[1].reshape((-1, 2)),
            radius=1,
        )
        gui.show()

        for _ in range(5):
            apply_pendulum_force()
            advance(1e-3)
