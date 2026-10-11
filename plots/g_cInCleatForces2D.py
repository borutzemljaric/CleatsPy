import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

"""
the code in this file serve to extract and draw axial forces in triangle cleat

"""


def axial_state(F_left, F_right, tol=1e-9):
    """Signed difference F_right - F_left, classified per time step."""
    F_L = np.asarray(F_left,  dtype=float).ravel()
    F_R = np.asarray(F_right, dtype=float).ravel()

    T = F_R - F_L
    tension     = T >  tol
    compression = T < -tol
    neutral     = ~(tension | compression)
    return T, tension, compression, neutral


def masked(values, mask):
    """Values with NaN outside the mask (so plots break there)."""
    values = np.asarray(values, dtype=float).ravel()
    out = np.full_like(values, np.nan, dtype=float)
    out[mask] = values[mask]
    return out


def compression_resultant(Tx, Ty, Tz):
    """
    Magnitude of the resultant of the compression components only.

    For each direction keep only the negative part (compression) and set
    the tension part to zero, then take the Euclidean length:

        T_comp_res = sqrt( cx^2 + cy^2 + cz^2 )

        cx = -Tx  if Tx < 0  else 0
        cy = -Ty  if Ty < 0  else 0
        cz = -Tz  if Tz < 0  else 0

    So the curve is zero (or NaN, if you prefer) when no direction is in
    compression, and positive otherwise.
    """
    Tx = np.asarray(Tx, dtype=float).ravel()
    Ty = np.asarray(Ty, dtype=float).ravel()
    Tz = np.asarray(Tz, dtype=float).ravel()

    cx = np.where(Tx < 0, -Tx, 0.0)
    cy = np.where(Ty < 0, -Ty, 0.0)
    cz = np.where(Tz < 0, -Tz, 0.0)

    return np.sqrt(cx * cx + cy * cy + cz * cz)


def tension_resultant(Tx, Ty, Tz):
    """
    Magnitude of the resultant of the tension components only.

    For each direction keep only the positive part (tension), set the
    compression part to zero, then take the Euclidean length:

        T_ten_res = sqrt( tx^2 + ty^2 + tz^2 )

        tx = Tx  if Tx > 0  else 0
        ty = Ty  if Ty > 0  else 0
        tz = Tz  if Tz > 0  else 0
    """
    Tx = np.asarray(Tx, dtype=float).ravel()
    Ty = np.asarray(Ty, dtype=float).ravel()
    Tz = np.asarray(Tz, dtype=float).ravel()

    tx = np.where(Tx > 0, Tx, 0.0)
    ty = np.where(Ty > 0, Ty, 0.0)
    tz = np.where(Tz > 0, Tz, 0.0)

    return np.sqrt(tx * tx + ty * ty + tz * tz)


def plot_axial_state(t, x1, x2, y1, y2, z1, z2, tol=1e-9,
                     title="Effective axial force per direction   "
                           "T = F_right - F_left"):
    """
    Plot T_x, T_y, T_z and the resultant compression magnitude on one graph.

    Line style encodes the state:
        solid   -> tension      (T >  tol)
        dashed  -> compression  (T < -tol)
        dotted  -> neutral      (|T| <= tol)

    Colour encodes the direction:
        red     -> X
        green   -> Y
        blue    -> Z
        black   -> |compression resultant|  (4th curve)
    """
    t = np.asarray(t, dtype=float).ravel()

    Tx, tx_ten, tx_com, tx_neu = axial_state(x1, x2, tol=tol)
    Ty, ty_ten, ty_com, ty_neu = axial_state(y1, y2, tol=tol)
    Tz, tz_ten, tz_com, tz_neu = axial_state(z1, z2, tol=tol)

    # fourth curve: resultant magnitude of the tension parts
    Tt = tension_resultant(Tx, Ty, Tz)
    # show it only where at least one direction is actually in tension
    any_tens = tx_ten | ty_ten | tz_ten
    Tt_masked = masked(Tt, any_tens)

    fig, ax = plt.subplots(figsize=(12, 5))

    dirs = [
        ("X", Tx, tx_ten, tx_com, tx_neu, "tab:red"),
        ("Y", Ty, ty_ten, ty_com, ty_neu, "tab:green"),
        ("Z", Tz, tz_ten, tz_com, tz_neu, "tab:blue"),
    ]

    legend_handles = []
    for label, T, ten, com, neu, color in dirs:
        ax.plot(t, masked(T, ten), color=color, lw=1.8, ls="-")
        ax.plot(t, masked(T, com), color=color, lw=1.8, ls="--")
        ax.plot(t, masked(T, neu), color=color, lw=1.4, ls=":")
        legend_handles.append(
            Line2D([0], [0], color=color, lw=2, label=f"T_{label}")
        )

    # ---- fourth curve: tension resultant magnitude ----------------------
    ax.plot(t, Tt_masked, color="k", lw=2.2, ls="-", zorder=5)
    legend_handles.append(
        Line2D([0], [0], color="k", lw=2.2, ls="-",
               label="|T_tension| resultant")
    )

    # style legend (shared by X/Y/Z curves)
    legend_handles += [
        Line2D([0], [0], color="k", lw=2, ls="-",  label="tension  (T>0)"),
        Line2D([0], [0], color="k", lw=2, ls="--", label="compression  (T<0)"),
        Line2D([0], [0], color="k", lw=2, ls=":",  label="neutral  (T≈0)"),
    ]

    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlim(t[0], t[-1])

    ax.set_xlabel("Time [s]")
    ax.set_ylabel("T  [N]   (+ tension, − compression)")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(handles=legend_handles, loc="best", ncol=3, fontsize=9)

    plt.tight_layout()
    return fig, (Tx, Ty, Tz, Tt), \
        (tx_ten, ty_ten, tz_ten), \
        (tx_com, ty_com, tz_com)

def solve_triangle_axial_time(F1, F2, F3, a=1.0):
    """
    Time-domain axial-force solver for an equilateral triangle
    lying in the **X–Z plane** (Y is ignored / treated as out-of-plane).

    Triangle vertices (in the X–Z plane):
        L1 = (0,   0, 0)
        L2 = (a,   0, 0)
        L3 = (a/2, 0, a*sqrt(3)/2)

    Members:
        member 1 : L1-L2     member 2 : L2-L3     member 3 : L3-L1

    Sign convention:  N > 0 = tension, N < 0 = compression.
    """
    F1 = np.atleast_2d(np.asarray(F1, dtype=float))
    F2 = np.atleast_2d(np.asarray(F2, dtype=float))
    F3 = np.atleast_2d(np.asarray(F3, dtype=float))
    assert F1.shape[1] == F2.shape[1] == F3.shape[1] == 3, \
        "Each force must be (n, 3)"
    n = F1.shape[0]

    # --- geometry: equilateral triangle in the X–Z plane ---------------
    L1 = np.array([0.0,      0.0, 0.0])
    L2 = np.array([a,        0.0, 0.0])
    L3 = np.array([a / 2.0,  0.0, a * np.sqrt(3) / 2.0])

    def unit(p, q):
        v = q - p
        return v / np.linalg.norm(v)

    u12 = unit(L1, L2)     # ( 1,     0,  0)
    u23 = unit(L2, L3)     # (-1/2,   0, +sqrt3/2)
    u31 = unit(L3, L1)     # (-1/2,   0, -sqrt3/2)

    # --- 3x3 system in the X–Z plane -----------------------------------
    #   Node L1 :  N1*u12 - N3*u31 = -F1  (take x and z rows)
    #   Node L2 : -N1*u12 + N2*u23 = -F2  (take x row only)
    A = np.array([
        [ u12[0],  0.0,     -u31[0]],   # L1 x
        [ u12[2],  0.0,     -u31[2]],   # L1 z
        [-u12[0],  u23[0],   0.0    ],  # L2 x
    ])

    # --- right-hand side, vectorised over time -------------------------
    b = np.empty((n, 3))
    b[:, 0] = -F1[:, 0]     # L1 x
    b[:, 1] = -F1[:, 2]     # L1 z   <-- was -F1[:, 1]
    b[:, 2] = -F2[:, 0]     # L2 x

    Ainv = np.linalg.inv(A)
    N    = b @ Ainv.T       # shape (n, 3)

    N1, N2, N3 = N[:, 0], N[:, 1], N[:, 2]

    # --- diagnostics: full 3-D residuals per node ----------------------
    def resid(F, terms):
        R = F.copy()
        for coeff, u in terms:
            R = R + coeff[:, None] * u[None, :]
        return R

    R1 = resid(F1, [( N1, u12), (-N3, u31)])
    R2 = resid(F2, [(-N1, u12), ( N2, u23)])
    R3 = resid(F3, [(-N2, u23), ( N3, u31)])

    info = {
        "L1": L1, "L2": L2, "L3": L3,
        "u12": u12, "u23": u23, "u31": u31,
        "N": N,
        "residual_L1": R1, "residual_L2": R2, "residual_L3": R3,
        "F_sum":   F1 + F2 + F3,
        "Fz_total": F1[:, 2] + F2[:, 2] + F3[:, 2],
    }
    return N1, N2, N3, info

def classify_triangle(N, tol=1e-9):
    """Bool masks for tension / compression / neutral of an axial force series."""
    N = np.asarray(N, dtype=float).ravel()
    tension     = N >  tol
    compression = N < -tol
    neutral     = ~(tension | compression)
    return N, tension, compression, neutral


def plot_triangle_axial(t, N1, N2, N3, tol=1e-9,
                        title="Axial forces in the triangle members"):
    """
    Plot N1, N2, N3 vs time.  Line style encodes the state:
        solid   -> tension      (N >  tol)
        dashed  -> compression  (N < -tol)
        dotted  -> neutral      (|N| <= tol)

    Colour encodes the member:
        red     -> N1  (L1-L2)
        green   -> N2  (L2-L3)
        blue    -> N3  (L3-L1)
    """
    t  = np.asarray(t,  dtype=float).ravel()
    N1 = np.asarray(N1, dtype=float).ravel()
    N2 = np.asarray(N2, dtype=float).ravel()
    N3 = np.asarray(N3, dtype=float).ravel()

    members = [
        ("N1 (L1-L2)", N1, "tab:red"),
        ("N2 (L2-L3)", N2, "tab:green"),
        ("N3 (L3-L1)", N3, "tab:blue"),
    ]

    fig, ax = plt.subplots(figsize=(12, 5))

    legend_handles = []
    for label, N, color in members:
        N, ten, com, neu = classify_triangle(N, tol=tol)
        ax.plot(t, masked(N, ten), color=color, lw=1.8, ls="-")
        ax.plot(t, masked(N, com), color=color, lw=1.8, ls="--")
        ax.plot(t, masked(N, neu), color=color, lw=1.4, ls=":")
        legend_handles.append(
            Line2D([0], [0], color=color, lw=2, label=label)
        )

    legend_handles += [
        Line2D([0], [0], color="k", lw=2, ls="-",  label="tension  (N>0)"),
        Line2D([0], [0], color="k", lw=2, ls="--", label="compression  (N<0)"),
        Line2D([0], [0], color="k", lw=2, ls=":",  label="neutral  (N≈0)"),
    ]

    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlim(t[0], t[-1])

    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Axial force  N  [N]   (+ tension, − compression)")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(handles=legend_handles, loc="best", ncol=3, fontsize=9)

    plt.tight_layout()
    return fig, ax


def axial_state_T(N, tol=1e-9):
    """Classify an axial-force time series (reuse for plotting)."""
    N = np.asarray(N, dtype=float).ravel()
    tension     = N >  tol
    compression = N < -tol
    neutral     = ~(tension | compression)
    return N, tension, compression, neutral


def g_cIncleatForces2D(point_idx, data):
    """Plot effective axial force (X, Y, Z) between cleats vs time."""
    point_idx = point_idx - 1   # to python indexing

    t = np.asarray(data["t"], dtype=float)

    keys = ["FXcleave_1", "FYcleave_1", "FZcleave_1"]
    split_parts = {k: np.array_split(data.get(k), 3, axis=1) for k in keys}
    parts_X, parts_Y, parts_Z = [split_parts[k] for k in keys]

    # forces for each cable by components
    x1 = parts_X[0][:, point_idx]
    x2 = parts_X[1][:, point_idx]
    x3 = parts_X[2][:, point_idx]

    # trick, the longitudinal part is removed, only 2D part is calculated
    y1 = parts_Y[0][:, point_idx] * 0  #manualy set to zero
    y2 = parts_Y[1][:, point_idx] * 0
    y3 = parts_Y[2][:, point_idx] * 0

    z1 = parts_Z[0][:, point_idx]
    z2 = parts_Z[1][:, point_idx]
    z3 = parts_Z[2][:, point_idx]

    # fig, Ts, tens, coms = plot_axial_state(
    #     t, x1, x2, y1, y2, z1, z2, tol=1e-6,
    #     title="Tension / compression between cleats   "
    #           "(cable 1 vs cable 2, X/Y/Z)"
    # )
    #plt.show()


    F1 = np.column_stack([x1, y1, z1])
    F2 = np.column_stack([x2, y2, z2])
    F3 = np.column_stack([x3, y3, z3])

    N1, N2, N3, info = solve_triangle_axial_time(F1, F2, F3, a=1)  # function call
    print("N1 range:", N1.min(), "..", N1.max())
    print("N2 range:", N2.min(), "..", N2.max())
    print("N3 range:", N3.min(), "..", N3.max())

    # ---- sanity check: Y residuals should be ~0 (triangle is in X–Z) --
    # print("max |residual Y|:",
    #       max(np.abs(info["residual_L1"][:, 1]).max(),
    #           np.abs(info["residual_L2"][:, 1]).max(),
    #           np.abs(info["residual_L3"][:, 1]).max()))

    # single-panel plot
    fig, ax = plot_triangle_axial(t, N1, N2, N3, tol=1e-6,
                                  title="Axial forces in the triangle")

    plt.show()


    #return fig, Ts, tens, coms