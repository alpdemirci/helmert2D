# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 14:04:56 2021

@author: Alp
"""
import numpy as np
import pandas as pd

ctr = []
err = []
# Transformation parameters arrays
# ------------------------------
DYDXp = []
DXDYp = []
DYDYp = []
DXDXp = []
Sp2 = []
S2 = []
# ------------------------------
# Error calculation arrays
# ------------------------------
VV = []
non = []
# ------------------------------

# Input values
# -----------------------------
# Coordinates of the common stations
first_coor = pd.DataFrame([[103, 5610.789, 34923.955],
                           [139, 4960.700, 34937.870],
                           [84,  4137.970, 35087.280],
                           [85,  2678.804, 36319.719]],
                          columns=['P_Id', 'Y', 'X'])
secon_coor = pd.DataFrame([[103, 32155.825, 28768.351],
                           [139, 31505.671, 28779.647],
                           [84,  30682.347, 28925.755],
                           [85,  29218.210, 30152.331]],
                          columns=['P_Id', 'Y', 'X'])
# Coordinates of the remainig non-common stations
non_common = pd.DataFrame([[1405, 4858.600, 34987.200],
                           [1406, 4781.160, 35006.580],
                           [1407, 4584.770, 34965.960],
                           [1408, 4437.700, 35099.660],
                           [1409, 4311.020, 35085.460]],
                          columns=['P_Id', 'Y', 'X'])
# -------------------------------
# Coordinates of center of gravity of the common stations
Ys_first = sum(first_coor['Y']) / len(first_coor)
Xs_first = sum(first_coor['X']) / len(first_coor)
Ys_secon = sum(secon_coor['Y']) / len(secon_coor)
Xs_secon = sum(secon_coor['X']) / len(secon_coor)

# The coordinates of common points with respect to the center of gravity
for i in range(len(first_coor)):
    ctr.append([first_coor['Y'][i] - Ys_first, first_coor['X'][i] - Xs_first,
                secon_coor['Y'][i] - Ys_secon, secon_coor['X'][i] - Ys_secon])
ctr = pd.DataFrame(ctr, columns=['DYp', 'DXp', 'DY', 'DX'])

for i in range(len(ctr)):
    Sp2.append(ctr['DXp'][i] ** 2 + ctr['DYp'][i] ** 2)
    S2.append(ctr['DX'][i] ** 2 + ctr['DY'][i] ** 2)

# Generating transformation parameters
for i in range(len(first_coor)):
    DYDXp.append(ctr['DY'][i] * ctr['DXp'][i])
    DXDYp.append(ctr['DX'][i] * ctr['DYp'][i])
    DYDYp.append(ctr['DY'][i] * ctr['DYp'][i])
    DXDXp.append(ctr['DX'][i] * ctr['DXp'][i])
k11 = (sum(DXDXp) + sum(DYDYp)) / sum(Sp2)
k12 = (sum(DYDXp) - sum(DXDYp)) / sum(Sp2)
k01 = Xs_secon - Xs_first * k11 + Ys_first * k12
k02 = Ys_secon - Ys_first * k11 - Xs_first * k12
lamda = (k11 ** 2 + k12 ** 2) ** 0.5
epsilon = np.arctan2(k12, k11) * 200 / np.pi

# Error computations
for i in range(len(first_coor)):
    err.append([secon_coor['P_Id'][i], secon_coor['Y'][i], secon_coor['X'][i],
                k02 + first_coor['X'][i] * k12 + first_coor['Y'][i] * k11,
                k01 + first_coor['X'][i] * k11 - first_coor['Y'][i] * k12,
                k02 + first_coor['X'][i] * k12 + first_coor['Y'][i] * k11 - secon_coor['Y'][i],
                k01 + first_coor['X'][i] * k11 - first_coor['Y'][i] * k12 - secon_coor['X'][i]])
err = pd.DataFrame(err, columns=['P_Id', 'Y', 'X', 'Yp', 'Xp', 'Vy', 'Vx'])

for i in range(len(err)):
    VV.append(err['Vx'][i] ** 2)
    VV.append(err['Vy'][i] ** 2)

m0 = (sum(VV) / (2 * len(err) - 4)) ** 0.5
mk11 = m0 * (1 / sum(Sp2)) ** 0.5
mk01 = m0 * (1 / len(err) + (Ys_first ** 2 + Xs_first ** 2) / sum(Sp2)) ** 0.5
# Calculation of the 2nd system coordinate and accuracies of the non-common
for i in range(len(non_common)):
    non.append([non_common['P_Id'][i], non_common['Y'][i], non_common['X'][i],
                k02 + non_common['X'][i] * k12 + non_common['Y'][i] * k11,
                k01 + non_common['X'][i] * k11 - non_common['Y'][i] * k12])
non = pd.DataFrame(non, columns=['P_Id', 'Yp', 'Xp', 'Y', 'X'])
print(non)
