import taichi as ti
from taichi_glsl import length, sqrt

ti.init(arch=ti.gpu)

# 5x9 grid of springs
n, m = 5, 9
gravity = 9.8
k_spring, k_spring_diag = 1e5, 2e4
k_resistance = 0.5

pos = ti.Vector.field(2, dtype=float, shape=(n, m))
vel = ti.Vector.field(2, dtype=float, shape=(n, m))
forces = ti.Vector.field(2, dtype=float, shape=(n, m))
rest_lengths = ti.field(dtype=float, shape=(3,))


@ti.kernel
def initialize():
    len1, len2 = 0.5 / (n - 1), 0.8 / (m - 1)
    rest_lengths[0] = len1  # i -> i + 1
    rest_lengths[1] = len2  # j -> j + 1
    rest_lengths[2] = sqrt(len1 ** 2 + len2 ** 2)  # (i, j) -> (i + 1, j + 1)
    for i, j in pos:
        pos[i, j] = ti.Vector([0.1 + 0.8 * j / (m - 1), 0.9 - 0.5 * i / (n - 1)])
        vel[i, j] = ti.Vector([0.0, 0.0])


@ti.kernel
def apply_spring_force():
    """Fill in the `forces` field with spring forces between points."""
    for i, j in pos:
        forces[i, j] = ti.Vector([0.0, -gravity])
        forces[i, j] -= k_resistance * vel[i, j]
    for i, j in pos:
        if i < n - 1:
            d = pos[i, j] - pos[i + 1, j]
            force = k_spring * (length(d) - rest_lengths[0])
            forces[i, j] += -force * d
            forces[i + 1, j] += force * d
        if j < m - 1:
            d = pos[i, j] - pos[i, j + 1]
            force = k_spring * (length(d) - rest_lengths[1])
            forces[i, j] += -force * d
            forces[i, j + 1] += force * d
        if i < n - 1 and j < m - 1:
            d = pos[i, j] - pos[i + 1, j + 1]
            force = k_spring_diag * (length(d) - rest_lengths[2])
            forces[i, j] += -force * d
            forces[i + 1, j + 1] += force * d

            d = pos[i + 1, j] - pos[i, j + 1]
            force = k_spring_diag * (length(d) - rest_lengths[2])
            forces[i + 1, j] += -force * d
            forces[i, j + 1] += force * d


@ti.kernel
def advance(t: float):
    """Advance the simulation by `t` units of time using Euler's method."""
    for i, j in pos:
        if i != 0 or (j != 0 and j != m - 1):
            pos[i, j] += t * vel[i, j]
            vel[i, j] += t * forces[i, j]


def run():
    gui = ti.GUI("Cloth Simulation", res=(600, 600))
    initialize()
    while gui.running:
        pos_np = pos.to_numpy()
        gui.circles(pos_np.reshape((-1, 2)), radius=5)
        gui.lines(
            pos_np[:-1, :].reshape((-1, 2)),
            pos_np[1:, :].reshape((-1, 2)),
            radius=1,
        )
        gui.lines(
            pos_np[:, :-1].reshape((-1, 2)),
            pos_np[:, 1:].reshape((-1, 2)),
            radius=1,
        )
        gui.lines(
            pos_np[:-1, :-1].reshape((-1, 2)),
            pos_np[1:, 1:].reshape((-1, 2)),
            radius=1,
        )
        gui.lines(
            pos_np[1:, :-1].reshape((-1, 2)),
            pos_np[:-1, 1:].reshape((-1, 2)),
            radius=1,
        )
        gui.show()

        for _ in range(20):
            apply_spring_force()
            advance(1e-4)
