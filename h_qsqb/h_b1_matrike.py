import numpy as np

TEMPLATE_B1 = [
# row 0
[ 'A',  0,  0, '-B',  0,  0, '-A',  0,  0,  'D',  0,  0],
# row 1
[  0, 'A',  0,   0, '-B',  0,   0, '-A',  0,   0,  'D',  0],
# row 2
[  0,  0, 'A',   0,   0, '-B',  0,   0, '-A',  0,   0, 'D'],
# row 3
['-B',  0,  0,  'C',  0,  0,  'B',  0,  0, '-F',  0,  0],
# row 4
[  0, '-B', 0,   0,  'C',  0,   0,  'B',  0,   0, '-F',  0],
# row 5
[  0,  0, '-B',  0,   0,  'C',  0,   0,  'B',  0,   0, '-F'],
# row 6
['-A',  0,  0,  'B',  0,  0,  'A',  0,  0, '-D',  0,  0],
# row 7
[  0, '-A', 0,   0,  'B',  0,   0,  'A',  0,   0, '-D',  0],
# row 8
[  0,  0, '-A',  0,   0,  'B',  0,   0,  'A',  0,   0, '-D'],
# row 9
[ 'D',  0,  0, '-F',  0,  0, '-D',  0,  0,  'E',  0,  0],
# row 10
[  0, 'D',  0,   0, '-F',  0,   0, '-D',  0,   0,  'E',  0],
# row 11
[  0,  0, 'D',   0,   0, '-F',  0,   0, '-D',  0,   0, 'E'],
]


def h_b1_matrike(k, L):
    """
    B matrices for strain case. Returns b1 : (12, 12, n).
    """
    k = np.asarray(k, dtype=float)
    n = k.size
    b1 = np.zeros((12, 12, n))

    for idx in range(n):
        ki = k[idx]
        t1 = 6*ki - 6*ki**2
        t2 = 3*ki**2 - 4*ki + 1
        t3 = 2*ki - 3*ki**2

        values = {
            'A': (t1*t1) / L**2,
            'B': (t1*t2) / L,
            'C':  t2*t2,
            'D': (t3*t1) / L,
            'E':  t3*t3,
            'F':  t3*t2,
        }

        for i in range(12):
            for j in range(12):
                cell = TEMPLATE_B1[i][j]
                if cell == 0:
                    continue
                if cell.startswith('-'):
                    b1[i, j, idx] = -values[cell[1:]]
                else:
                    b1[i, j, idx] =  values[cell]

    return b1


#below is original version, ugraded vith upper one for speed
# def _orig_b1(k, L):
#     """
#     Calculate B matrices for strain case, B matrices are 3D arrays with depth of numel(k)
#     k is a column vector
#     """
#     n = len(k)  # get number of integration points
#     b1 = np.zeros((12, 12, n))
#
#     for i in range(n):  # create index
#         ki = k[i]
#
#         b1[:, :, i] = np.array([
#             [(6 * ki - 6 * ki ** 2) ** 2 / L ** 2, 0, 0, -((6 * ki - 6 * ki ** 2) * (3 * ki ** 2 - 4 * ki + 1)) / L, 0,
#              0, -(6 * ki - 6 * ki ** 2) ** 2 / L ** 2, 0, 0, ((2 * ki - 3 * ki ** 2) * (6 * ki - 6 * ki ** 2)) / L, 0,
#              0],
#             [0, (6 * ki - 6 * ki ** 2) ** 2 / L ** 2, 0, 0, -((6 * ki - 6 * ki ** 2) * (3 * ki ** 2 - 4 * ki + 1)) / L,
#              0, 0, -(6 * ki - 6 * ki ** 2) ** 2 / L ** 2, 0, 0, ((2 * ki - 3 * ki ** 2) * (6 * ki - 6 * ki ** 2)) / L,
#              0],
#             [0, 0, (6 * ki - 6 * ki ** 2) ** 2 / L ** 2, 0, 0,
#              -((6 * ki - 6 * ki ** 2) * (3 * ki ** 2 - 4 * ki + 1)) / L, 0, 0, -(6 * ki - 6 * ki ** 2) ** 2 / L ** 2, 0,
#              0, ((2 * ki - 3 * ki ** 2) * (6 * ki - 6 * ki ** 2)) / L],
#             [-((6 * ki - 6 * ki ** 2) * (3 * ki ** 2 - 4 * ki + 1)) / L, 0, 0, (3 * ki ** 2 - 4 * ki + 1) ** 2, 0, 0,
#              ((6 * ki - 6 * ki ** 2) * (3 * ki ** 2 - 4 * ki + 1)) / L, 0, 0,
#              -(2 * ki - 3 * ki ** 2) * (3 * ki ** 2 - 4 * ki + 1), 0, 0],
#             [0, -((6 * ki - 6 * ki ** 2) * (3 * ki ** 2 - 4 * ki + 1)) / L, 0, 0, (3 * ki ** 2 - 4 * ki + 1) ** 2, 0, 0,
#              ((6 * ki - 6 * ki ** 2) * (3 * ki ** 2 - 4 * ki + 1)) / L, 0, 0,
#              -(2 * ki - 3 * ki ** 2) * (3 * ki ** 2 - 4 * ki + 1), 0],
#             [0, 0, -((6 * ki - 6 * ki ** 2) * (3 * ki ** 2 - 4 * ki + 1)) / L, 0, 0, (3 * ki ** 2 - 4 * ki + 1) ** 2, 0,
#              0, ((6 * ki - 6 * ki ** 2) * (3 * ki ** 2 - 4 * ki + 1)) / L, 0, 0,
#              -(2 * ki - 3 * ki ** 2) * (3 * ki ** 2 - 4 * ki + 1)],
#             [-(6 * ki - 6 * ki ** 2) ** 2 / L ** 2, 0, 0, ((6 * ki - 6 * ki ** 2) * (3 * ki ** 2 - 4 * ki + 1)) / L, 0,
#              0, (6 * ki - 6 * ki ** 2) ** 2 / L ** 2, 0, 0, -((2 * ki - 3 * ki ** 2) * (6 * ki - 6 * ki ** 2)) / L, 0,
#              0],
#             [0, -(6 * ki - 6 * ki ** 2) ** 2 / L ** 2, 0, 0, ((6 * ki - 6 * ki ** 2) * (3 * ki ** 2 - 4 * ki + 1)) / L,
#              0, 0, (6 * ki - 6 * ki ** 2) ** 2 / L ** 2, 0, 0, -((2 * ki - 3 * ki ** 2) * (6 * ki - 6 * ki ** 2)) / L,
#              0],
#             [0, 0, -(6 * ki - 6 * ki ** 2) ** 2 / L ** 2, 0, 0,
#              ((6 * ki - 6 * ki ** 2) * (3 * ki ** 2 - 4 * ki + 1)) / L, 0, 0, (6 * ki - 6 * ki ** 2) ** 2 / L ** 2, 0,
#              0, -((2 * ki - 3 * ki ** 2) * (6 * ki - 6 * ki ** 2)) / L],
#             [((2 * ki - 3 * ki ** 2) * (6 * ki - 6 * ki ** 2)) / L, 0, 0,
#              -(2 * ki - 3 * ki ** 2) * (3 * ki ** 2 - 4 * ki + 1), 0, 0,
#              -((2 * ki - 3 * ki ** 2) * (6 * ki - 6 * ki ** 2)) / L, 0, 0, (2 * ki - 3 * ki ** 2) ** 2, 0, 0],
#             [0, ((2 * ki - 3 * ki ** 2) * (6 * ki - 6 * ki ** 2)) / L, 0, 0,
#              -(2 * ki - 3 * ki ** 2) * (3 * ki ** 2 - 4 * ki + 1), 0, 0,
#              -((2 * ki - 3 * ki ** 2) * (6 * ki - 6 * ki ** 2)) / L, 0, 0, (2 * ki - 3 * ki ** 2) ** 2, 0],
#             [0, 0, ((2 * ki - 3 * ki ** 2) * (6 * ki - 6 * ki ** 2)) / L, 0, 0,
#              -(2 * ki - 3 * ki ** 2) * (3 * ki ** 2 - 4 * ki + 1), 0, 0,
#              -((2 * ki - 3 * ki ** 2) * (6 * ki - 6 * ki ** 2)) / L, 0, 0, (2 * ki - 3 * ki ** 2) ** 2]
#         ])
#
#     return b1
#
#
#
#
# rng = np.random.default_rng(0)
# for trial in range(5):
#     L = rng.uniform(0.5, 5.0)
#     k = rng.uniform(0.05, 0.95, size=4)
#
#     old = _orig_b1(k, L)
#     new = h_b1_matrike(k, L)
#
#     assert np.allclose(old, new, rtol=1e-12, atol=1e-14), \
#         f"Mismatch at trial={trial}, L={L}, k={k}"
#
# print("B1 template verified ✔")

