from manim import *
from manim_slides import Slide
import symfem

# ctrl + p to run in presentation
# this is an own shortcut made with the tasks.json configuration in VS Code
# and the shortcut is made in the ctrl+shift+p "Preferences: open Keyboard Shortcuts (JSON)" menu

# below command in terminal in correct directory converts the Manim slides to an HTML presentation
# manim-slides convert --to html -v DEBUG FiniteElementPresentation presentation.html 

n_cells = 5
p       = 3

# Create reference finite element
element = symfem.create_element(
    "interval",
    "Lagrange",
    p,
)
import numpy as np
import sympy as sp
from symfem.symbols import x


def build_reconstruction_matrix(
    element,
    cell_maps,
    points_per_cell=100,
):
    """
    Build R such that

        u_plot = R @ u

    for a 1D conforming scalar Symfem element.

    Parameters
    ----------
    element
        Symfem element on an interval.

    cell_maps
        List of functions F_e(xi), xi in [0,1],
        mapping the unit reference interval to each physical cell.

    points_per_cell
        Number of plotting points per cell.

    Returns
    -------
    x_plot : ndarray
        Physical coordinates.

    R : ndarray
        Reconstruction matrix.
        R[:, i] is global basis function i.

    dofmap : ndarray
        dofmap[e, a] = global DOF corresponding to
        local basis function a on element e.
    """

    n_cells = len(cell_maps)
    n_local = element.space_dim

    # --------------------------------------------------
    # 1. Find which LOCAL basis functions belong to
    #    left vertex, right vertex, and cell interior
    # --------------------------------------------------

    left_dofs = list(element.entity_dofs(0, 0))
    right_dofs = list(element.entity_dofs(0, 1))
    interior_dofs = list(element.entity_dofs(1, 0))

    if len(left_dofs) != len(right_dofs):
        raise ValueError(
            "This routine expects the same number of DOFs "
            "on the left and right endpoint."
        )

    n_vertex_dofs = len(left_dofs)

    # --------------------------------------------------
    # 2. Automatically construct global DOF numbering
    #
    # Every interface shares its vertex DOFs.
    # Interior DOFs remain local to a cell.
    # --------------------------------------------------

    n_vertices = n_cells + 1

    vertex_dofs = np.arange(
        n_vertices * n_vertex_dofs
    ).reshape(n_vertices, n_vertex_dofs)

    next_dof = n_vertices * n_vertex_dofs

    dofmap = np.full(
        (n_cells, n_local),
        -1,
        dtype=int
    )

    for e in range(n_cells):

        # left endpoint
        for k, a in enumerate(left_dofs):
            dofmap[e, a] = vertex_dofs[e, k]

        # right endpoint
        for k, a in enumerate(right_dofs):
            dofmap[e, a] = vertex_dofs[e + 1, k]

        # cell-internal DOFs
        for a in interior_dofs:
            dofmap[e, a] = next_dof
            next_dof += 1

    if np.any(dofmap < 0):
        raise ValueError(
            "Element contains DOFs that are neither vertex nor "
            "cell-interior DOFs."
        )

    n_global = next_dof

    # --------------------------------------------------
    # 3. Build physical plotting points and R
    # --------------------------------------------------

    xi_plot = np.linspace(0.0, 1.0, points_per_cell)

    x_blocks = []
    R_blocks = []

    for e, F in enumerate(cell_maps):

        # User's reference -> physical coordinate map
        x_cell = np.asarray(
            [F(xi) for xi in xi_plot],
            dtype=float
        )

        # Physical endpoints
        a = float(F(0.0))
        b = float(F(1.0))

        # ------------------------------------------------
        # Ask Symfem for the CORRECT physical-cell basis.
        #
        # This is the important part. We are NOT manually
        # constructing Lagrange/Hermite/etc. functions.
        # ------------------------------------------------

        physical_basis = element.map_to_cell(
            [(a,), (b,)]
        )

        R_cell = np.zeros(
            (points_per_cell, n_global)
        )

        # Evaluate each mapped local basis function
        for local_i, phi in enumerate(physical_basis):

            global_i = dofmap[e, local_i]

            expr = phi.as_sympy()

            # Turn symbolic expression into a normal NumPy function
            phi_fun = sp.lambdify(
                x[0],
                expr,
                modules="numpy"
            )

            values = np.asarray(
                phi_fun(x_cell),
                dtype=float
            )

            # Constant functions can produce a scalar
            if values.ndim == 0:
                values = np.full_like(
                    x_cell,
                    float(values)
                )

            R_cell[:, global_i] += values

        # Don't duplicate interface points in the returned plot
        if e > 0:
            x_cell = x_cell[1:]
            R_cell = R_cell[1:, :]

        x_blocks.append(x_cell)
        R_blocks.append(R_cell)

    x_plot = np.concatenate(x_blocks)
    R = np.vstack(R_blocks)

    return x_plot, R, dofmap

element = symfem.create_element(
    "interval",
    "Hermite",
    3
)


L = 5.0
n_cells = 4

cell_maps = [
    lambda xi, e=e: L * (e + xi) / n_cells
    for e in range(n_cells)
]

x_plot, R, dofmap = build_reconstruction_matrix(
    element,
    cell_maps
)


dof_positions = np.zeros(R.shape[1])
counts = np.zeros(R.shape[1])

for e, F in enumerate(cell_maps):
    for local_i, global_i in enumerate(dofmap[e]):

        # Symfem DOF point on reference cell
        dof = element.dofs[local_i]

        # For point-evaluation Lagrange DOFs
        xi = float(dof.dof_point()[0])

        x_phys = F(xi)

        dof_positions[global_i] += x_phys
        counts[global_i] += 1

dof_positions /= counts


class FiniteElementPresentation(Slide):

    def construct(self):

        # Axes
        axes = Axes(
            x_range=[0, 5, 1],
            y_range=[
                min(-0.5, 1.2),
                max(0.5, 1.2),
                1
            ],
            # x_length=10,
            # y_length=5,
            tips=False,
        )

        labels = axes.get_axis_labels(
            MathTex("x"),
            MathTex("u_h(x)")
        )

        self.play(
            Create(axes),
            Write(labels),
        )

        # Zero line
        zero_line = axes.plot(
            lambda x: 0,
            x_range=[0, 5],
        )

        self.play(Create(zero_line))

        # Build unit vectors e_i
        us = [np.zeros(n_dofs) for _ in range(n_dofs)]

        for i in range(n_dofs):
            us[i][i] = 1.0


        solutions = []
        zero_solutions = []

        for i in range(n_dofs):
            x_plot, u_plot = reconstruct_u(us[i])

            # Actual basis function
            solution = VMobject()
            solution.set_points_as_corners(
                [axes.c2p(x, y) for x, y in zip(x_plot, u_plot)]
            )
            solution.set_stroke(width=2)

            # Same x coordinates, but y = 0
            zero_solution = VMobject()
            zero_solution.set_points_as_corners(
                [axes.c2p(x, 0) for x in x_plot]
            )
            zero_solution.set_stroke(width=2)

            solutions.append(solution)
            zero_solutions.append(zero_solution)


        # Start with all basis functions collapsed onto the zero line
        self.add(*zero_solutions)

        self.next_slide()

        # Transform all of them at the same time
        self.play(
            *[
                Transform(zero_solutions[i], solutions[i])
                for i in range(n_dofs)
            ],
            run_time=3
        )

        self.wait()

        # # Start from the zero line:
        # solution_zero = VMobject()

        # zero_points = [
        #     axes.c2p(x, 0)
        #     for x in x_plot
        # ]

        # solution_zero.set_points_as_corners(zero_points)

        # # Add zero-version first, then morph into FEM solution
        # self.add(solution_zero)

        # self.next_slide()

        # self.play(
        #     Transform(solution_zero, solution),
        #     run_time=3
        # )

        # self.wait()