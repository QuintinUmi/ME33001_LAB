##########################################################################################################
# Copyright: QIN Qijun
# 
# This document is the property of Qin Qijun, please cite the source if you need to quote it.
# Workshop: https://github.com/QuintinUmi
# Contact:  qqj030212@gmail.com
##########################################################################################################


import numpy as np
from lab_report_tool_package import h_beam_package as hb

LENGTH = 100.
MAX_AREA = 500.
ACCURACY = 2
LENGTH_LIMIT_COE = 500
LL = LENGTH_LIMIT_COE * MAX_AREA**0.5
HOR_MAX = 39
DENSITY = 1.05 * 10**(-6)
DELH_COE = 0.5

POISSON_RATIO = 0.4
ELASTIC_MODULUS = 2.2 * 10**3
SHEAR_MODULUS = ELASTIC_MODULUS / 2 / (1 + POISSON_RATIO)
YIELD_STRESS = 31
FRACTURE_STRESS = 33
SHEAR_STRENGTH = 1 / (3**(1/2)) * YIELD_STRESS
F = 0.
M = 0.
K = 0.5

EXPERIMENTAL_MAX_FORCE = 1325.

b = 38.5
s = 2
h = 34
t = 2

C1 = 1.348
C2 = 0.630
k = 0.8

f_buckling = hb.buckling(ELASTIC_MODULUS * 10**(6), LENGTH, t, h, K) * 10**(-6)
i_x = hb.inertia_x(b, s, h, t)
i_z = hb.inertia_z(b, s, h, t)
i_w = hb.warping_constant(b, s, h, t)
i_t = hb.torsion_constant(b, s, h, t)
mcr = hb.mcr_cal(C1, C2, ELASTIC_MODULUS * 10**(-3), SHEAR_MODULUS * 10**(-3), 
                i_z, i_w, i_t, 
                (LENGTH - 20), (h / 2 + s), k, 1.0) * 10**3


def force_to_moment(force, deflection, position):

    moment = []
    for i in range(len(deflection)):
        if(10 <= position <= 50):
            moment.append(1/2 * force[i] * (position - 10))
        elif(50 < position <= 90):
            moment.append(1/2 * force[i] * (90 - position))
        else:
            moment.append(0)
    return moment


def force_v(force, deflection, position):

    v = []
    for i in range(len(deflection)):
        if(10 <= position <= 50):
            v.append(1/2 * force[i])
        elif(50 < position <= 90):
            v.append(-1/2 * force[i])
        else:
            v.append(0)
    return v