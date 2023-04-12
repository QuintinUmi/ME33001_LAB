##########################################################################################################
# Copyright: QIN Qijun
# 
# This document is the property of Qin Qijun, please cite the source if you need to quote it.
# Workshop: https://github.com/QuintinUmi
# Contact:  qqj030212@gmail.com
##########################################################################################################


import numpy as np
import math

PI = 3.1415926535897932

def area_inertia(b, s, h, t):
    return (b * (2 * s + h)**3 - h**3 * (b - t)) / 12

def static_moment(b, s, h, t):
    return s*b*(h+s)/2 + h*h*t/8

def buckling(ELASTIC_MODULUS, LENGTH, t, h, K):

    f_buckling = PI**2 * ELASTIC_MODULUS * (1/12 * t * h**3) / (K * LENGTH)**2
    return f_buckling

def inertia_x(b, s, h, t):
    i_x = (b * (2 * s + h)**3 - h**3 * (b - t)) / 12
    return i_x

def inertia_z(b, s, h, t):
    i_z = 2 * 1/12 * s * b**3 + 1/12 * h * t**3
    return i_z

def torsion_constant(b, s, h, t):
    i_t = (2 * b * s**3 + (h + s) * t**3) / 3
    return i_t
    
def  warping_constant(b, s, h, t):
    # i_w = (h + s)**2 * b**3 * s / 24
    i_w = inertia_z(b, s, h, t) * (h - s)**2 / 4
    return i_w

def mcr_cal(C1, C2, E, G, i_z, i_w, i_t, l_lt, z_g, k=1.0, k_w=1.0):

    # i_z = i_z * 10**(-12)
    # i_w = i_w * 10**(-18)
    # i_t = i_t * 10**(-12)
    # l_lt = l_lt * 10**(-3)
    # z_g = z_g * 10**(-3)

    mcr = C1 * PI**2 * E * i_z / (k * l_lt)**2 * (((k / k_w)**2 * (i_w / i_z) + (k * l_lt)**2 * G * i_t / (PI**2 * E * i_z ) + (C2 * z_g)**2)**0.5 - C2 * z_g)
    # mcr = PI / l_lt * (E * i_z * G * i_t)**(0.5)
    
    return mcr
    
def experiment_theta(b, s, h, t, FRACTURE_STRESS, EXPERIMENTAL_MAX_FORCE):

    i_y = 1/12 * h * t**3
    M_y = FRACTURE_STRESS * i_y / (t/2)
    F_y = M_y / (h * 0.75) * 10**(-6)
    theta = math.asin(F_y / EXPERIMENTAL_MAX_FORCE)

    return theta