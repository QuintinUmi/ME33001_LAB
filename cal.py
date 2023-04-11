import os
import math

LENGTH = 100.
MAX_AREA = 500.
ACCURACY = 2
LENGTH_LIMIT_COE = 500
LL = LENGTH_LIMIT_COE * MAX_AREA**0.5
HOR_MAX = 39
DENSITY = 1.05 * 10**(-6)
DELH_COE = 0.5

POISSON_RATIO = 0.4
ELASTIC_MODULUS = 2.2 * 10**9
SHEAR_MODULUS = ELASTIC_MODULUS / 2 / (1 + POISSON_RATIO)
YIELD_STRESS = 31 * 10**6
FRACTURE_STRESS = 33 * 10**6
SHEAR_STRENGTH = 1 / (3**(1/2)) * YIELD_STRESS
F = 0.
M = 0.
K = 0.65

POS_BASE_FIXTURES = 10

EXPERIMENTAL_MAX_FORCE = 1325.

PI = 3.1415926535897932

d = t = a = s = n = ACCURACY #I beam


def mcr_cal(C1, C2, E, G, i_z, i_w, i_t, l_lt, z_g, k=1.0, k_w=1.0):

    i_z = i_z * 10**(-12)
    i_w = i_w * 10**(-18)
    i_t = i_t * 10**(-12)
    l_lt = l_lt * 10**(-3)
    z_g = z_g * 10**(-3)

    mcr = C1 * PI**2 * E * i_z / (k * l_lt)**2 * (((k / k_w)**2 * (i_w / i_z) + (k * l_lt)**2 * G * i_t / (PI**2 * E * i_z ) + (C2 * z_g )**2)**0.5 - C2 * z_g)
    # mcr = PI / l_lt * (E * i_z * G * i_t)**(0.5)
    
    return mcr

def stability_factor(mcr, c, i_x, yield_strength):

    mcr = mcr * 10**(-9)
    c = c * 10**(-3)
    i_x = i_x * 10**(-12)

    rou_b = mcr * c / i_x / yield_strength
    return rou_b


class HollowRectangular:
    def __init__(self):
        self.area = MAX_AREA
        self.i_hr_max = 0.
        self.M_hr_max = 0.
        self.F_hr_max = 0.
        self.bf = self.df = self.hf = self.kf = 0.

        self.fm = 0.

    def cal(self):

        b = d = h = k = ACCURACY #hallow squre

        while(self.area > 4 * ACCURACY**2):

            area_temp = self.area

            d = 2 * ACCURACY
            b = 2 * ACCURACY
            h = 0
            k = 0
            
            while(area_temp <= self.area):
                    
                while(area_temp <= self.area):

                    h = b - 2 * ACCURACY
                    k = d - 2 * ACCURACY
                    area_temp = b * d - h * k
                    while(area_temp <= self.area and h >= 0.):
                            
                        while(area_temp <= self.area and k >= 0.):
                        
                            mass = 100 * area_temp * DENSITY

                            i_temp = (b * d**3 - h * k**3) / 12
                            M_temp = FRACTURE_STRESS * i_temp / (d/2)
                            F_normal_temp = M_temp / 20 *10**(-6)

                            # Q_shear = (b*d/2 - k*h/2) * ((b)*(d/2)*(d/4) - (k/2)*h*(k/4)) / (b*d/2 - k*h/2)
                            Q_shear = 2 * (d/4) * ((d/2) * (b - h)) + ((d/2) + (k/2)) / 2 * (h * (d - k))
                            F_shear_temp = 4 * SHEAR_STRENGTH * (b - h) * i_temp / Q_shear *10**(-6)

                            M_y = YIELD_STRESS / (d/2) * i_temp * 10**(-9)
                            F_y = M_y / (math.sin(H_Beam().theta) * (d/2) * 10**(-3))

                            F_temp = min(F_normal_temp, F_shear_temp, F_y)
                            # print(F_normal_temp > F_shear_temp)
                            # print(F_shear_temp, F_normal_temp)

                            # if(b == 38 and d == 38 and h == 34 and k == 30):
                            #     print(F_normal_temp, F_shear_temp, F_y, mass, self.area)

                            # if(i_temp > self.i_hr_max):
                            if(F_temp/mass > self.fm):

                                self.fm = F_temp/mass
                                        
                                self.i_hr_max = i_temp
                                self.M_hr_max = M_temp
                                self.F_hr_max = F_temp
                                        
                                self.bf = b
                                self.df = d
                                self.hf = h
                                self.kf = k

                                self.areaf = b * d - h * k

                            k -= ACCURACY
                            area_temp = b * d - h * k

                        h -= ACCURACY
                        k = d - 2 * ACCURACY
                        area_temp = b * d - h * k

                    d += ACCURACY
                    h = b - 2 * ACCURACY
                    k = d - 2 * ACCURACY
                    area_temp = b * d - h * k
                    if(b >= HOR_MAX or d >= HOR_MAX):
                        break
                
                b += ACCURACY
                d = 2 * ACCURACY
                h = b - 2 * ACCURACY
                k = d - 2 * ACCURACY
                area_temp = b * d - h * k
                if(b >= HOR_MAX or d >= HOR_MAX):
                    break

            self.area -= ACCURACY


        print("bf = {}\ndf = {}\nhf = {}\nkf = {}".format(self.bf, self.df, self.hf, self.kf))
        print("Inertia = {}".format(self.i_hr_max))
        # print("--Normal Stress Max = {}\n---Shear Stress Max = {}".format(self.normal_stress_max, self.shear_stress_max))
        print("Area =", self.areaf)
        print("Force_max = ", self.F_hr_max)
        print("fm_max = ", self.fm) # N/kg
        print("hollow Rectangular finished ----------------------------------------------")




class HollowSquare:
    def __init__(self):
        self.area = MAX_AREA
        self.i_hs_max = 0.
        self.M_hs_max = 0.
        self.F_hs_max = 0.
        self.af = self.bf = 0.

        self.fm = 0.

    def cal(self):

        a = b = ACCURACY #hallow squre

        while(self.area > 4 * ACCURACY**2):

            area_temp = self.area

            a = (self.area + 4 * ACCURACY**2) / (4 * ACCURACY)
            
            while(a >= 2):

                b = a - 2 * ACCURACY
                # print(area_temp, a, b)
                area_temp = a**2 - b**2
                while(area_temp <= self.area and b >= 0.):
                    
                    mass = 100 * area_temp * DENSITY

                    i_temp = (a**4 - b**4) / 12
                    M_temp = FRACTURE_STRESS * i_temp / (a/2)
                    F_normal_temp = M_temp / 20 *10**(-6)

                    Q_shear = (a*a/2 - b*b/2) * ((a)*(a/2)*(a/4) - (b/2)*b*(b/4)) / (a*a/2 - b*b/2)
                    F_shear_temp = SHEAR_STRENGTH * (a - b) * i_temp / Q_shear *10**(-6)

                    F_temp = min(F_normal_temp, F_shear_temp)
                    # print(F_normal_temp > F_shear_temp)

                    # print(F_normal_temp, F_shear_temp)
                    
                    # if(i_temp > self.i_hr_max):
                    if(F_temp/mass > self.fm and a < HOR_MAX):
                        
                        # print(mass)
                        self.fm = F_temp/mass
                                
                        self.i_hs_max = i_temp
                        self.M_hs_max = M_temp
                        self.F_hs_max = F_temp
                                
                        self.af = a
                        self.bf = b

                        self.areaf = a**2 - b**2

                    b -= ACCURACY
                    area_temp = a**2 - b**2

                a -= ACCURACY
                if(a >= HOR_MAX):
                    break

            self.area -= ACCURACY


        print("af = {}\nbf = {}".format(self.af, self.bf))
        print("Inertia = {}".format(self.i_hs_max))
        # print("--Normal Stress Max = {}\n---Shear Stress Max = {}".format(self.normal_stress_max, self.shear_stress_max))
        print("Area =", self.areaf)
        print("Force_max = ", self.F_hs_max)
        print("fm_max = ", self.fm) # N/kg
        print("HollowSquare finished ----------------------------------------------")


class H_Beam:
    def __init__(self):
        self.C1 = 1.348
        self.C2 = 0.630
        self.k = 0.8

        self.area = MAX_AREA
        self.i_hb_max = 0.
        self.M_hb_max = 0.
        self.F_hb_max = 0.
        self.bf = self.sf = self.hf = self.tf = 0.

        self.theta = self.experiment_theta(38.5, 2, 34, 2)

        self.areaf = MAX_AREA
        self.normal_stress_max = 0.
        self.shear_stress_max = 0.

        self.fm = 0.

    def cal(self):
        f = open("fm log.txt", "w")
        b = s = h = t = ACCURACY #H beam
        count = 0
        while(self.area > 4 * ACCURACY**2):
            
            h = ACCURACY
            t = ACCURACY
            s = ACCURACY
            b = (self.area - h * t) / (2 * s + h - h)
            mass = LENGTH * self.area * DENSITY
            while(b > t + ACCURACY):
                while(b > t + ACCURACY):
                    while(b > t + ACCURACY):

                        
                        count += 1
                                                                     
                        i_temp = (b * (2 * s + h)**3 - h**3 * (b - t)) / 12
                        M_temp = YIELD_STRESS * i_temp / (h/2 + s)
                        F_normal_temp = M_temp / (20 * 10**(-3)) *10**(-9)

                        
                        Q_shear = s*b*(h+s)/2 + h*h*t/8
                        F_shear_temp = 2 * SHEAR_STRENGTH * t * i_temp / Q_shear *10**(-6)

                        # print(F_normal_temp, F_shear_temp, F_buckling)


                        mcr_temp = mcr_cal(self.C1, self.C2, ELASTIC_MODULUS, SHEAR_MODULUS, self.inertia_z(b, s, h, t), 
                                self.warping_constant(b, s, h, t), self.inertia_t(b, s, h, t), LENGTH - 20, 
                                h / 2 + s, self.k, 1.0)
                        F_buckling_temp = mcr_temp / (20 * 10**(-3))
                        rou_b = stability_factor(mcr_temp, h / 2 + s, i_temp, YIELD_STRESS)

                        M_y = YIELD_STRESS / (h/2) * (h * t**3) * 10**(-9)
                        F_y = M_y / (math.sin(self.theta) * (h/2) * 10**(-3))
                        # M_y = F_temp * math.sin(self.theta) * (h/2) * 10**(-3)
                        # stress_y = M_y * (h/2) / (h * t**3) * 10**(9)

                        ###################################################################################
                        
                        
                        Pcr = PI**2 * ELASTIC_MODULUS * self.inertia_z(b, s, h, t) / (K * (LENGTH - 2*POS_BASE_FIXTURES))**2 *10**(-6)

                        Pt = (PI**2 * ELASTIC_MODULUS * self.inertia_t(b, s, h, t)) / ((LENGTH - 2 * POS_BASE_FIXTURES)**2) *10**(-6)
                        # Pt = (SHEAR_MODULUS * self.torsion_constant(b, s, h, t)) / (K * (LENGTH - 2*POS_BASE_FIXTURES))**2 *10**(-6)

                        # Pj = self.c(b, s, h, t) * ELASTIC_MODULUS * self.inertia_z(b, s, h, t) / ((LENGTH - 2*POS_BASE_FIXTURES)**2) *10**(-9)

                        F_temp = min(F_normal_temp, F_shear_temp, Pcr, Pt)

                        if(b == 38.5 and s == 2 and h == 34 and t == 2):
                            print(F_normal_temp, F_shear_temp, Pcr, Pt, F_y, mass, self.area)

                        
                        # print(stress_y >= YIELD_STRESS)

                        
                        # if(F_temp == F_normal_temp):
                        #     print(1)
                        # if(F_temp == F_shear_temp):
                        #     print(2)
                        # if(F_temp == F_buckling):
                        #     print(3)
                        # if(F_temp == F_buckling_temp):
                        #     print("b = ", b, "s = ", s, "h = ", h ,"t = ", t)
                        #     print(F_buckling_temp)
                        # if(F_temp == Pcr):
                        #     print("b = ", b, "s = ", s, "h = ", h ,"t = ", t)
                        #     print(Pcr)
                        
                        # print("b = ", b, "s = ", s, "h = ", h ,"t = ", t, "fm = ",F_temp/mass)
                        
                        # print("b = ", b, "s = ", s, "h = ", h ,"t = ", t, "fm = ", F_temp/mass, file=f, flush=False)
                        
                        # print(M_temp, mcr_temp)
                        # print(F_buckling_temp)
                        
                        if(F_temp/mass > self.fm and b > t + ACCURACY and b < HOR_MAX and h + 2 * s < HOR_MAX):
                            # print(stress_y)
                            # print(b, d, h, k, F_temp, mass, F_temp/mass)
                            # print(F_normal_temp, F_shear_temp, F_buckling_temp, Pcr, Pt, Pj)
                            # print(M_temp, mcr_temp)
                            # print("b = ", b, "s = ", s, "h = ", h ,"t = ", t)
                            # print(M_temp, mcr_temp)
                            # print(F_buckling_temp)

                            # if(F_temp == Pcr):
                            #     print("Pcr b = ", b, "s = ", s, "h = ", h ,"t = ", t)
                            #     print(Pcr)
                            # if(F_temp == Pt):
                            #     print("Pt b = ", b, "s = ", s, "h = ", h ,"t = ", t)
                            #     print(Pt)
                            # if(F_temp == Pj):
                            #     print("Pj b = ", b, "s = ", s, "h = ", h ,"t = ", t)
                            #     print(Pj)
                            # if(F_temp == F_normal_temp):
                            #     print("F_normal_temp Pt b = ", b, "s = ", s, "h = ", h ,"t = ", t)
                            #     print(F_normal_temp)
                            # if(F_temp == F_shear_temp):
                            #     print("F_shear_temp Pt b = ", b, "s = ", s, "h = ", h ,"t = ", t)
                            #     print(F_shear_temp)

                            


                            self.fm = F_temp/mass
                            
                            self.i_hb_max = i_temp
                            self.M_hb_max = M_temp
                            self.F_hb_max = F_temp
                            
                            self.bf = b
                            self.sf = s
                            self.hf = h
                            self.tf = t

                            self.areaf = self.area
                            self.normal_stress_max = (20 * 10**(-3)) * F_temp * (h + s)/2 / i_temp *10**(9)
                            self.shear_stress_max = F_temp * Q_shear / (i_temp * t) *10**(6)
                        
                        t += ACCURACY
                        b = (self.area - h * t) / (2 * s + h - h)
                        if(h + 2 * s > HOR_MAX or b < t + ACCURACY):
                            break
                        
                    s += ACCURACY
                    t = ACCURACY
                    b = (self.area - h * t) / (2 * s + h - h)
                    if(h + 2 * s > HOR_MAX):
                        break

                h += ACCURACY
                s = ACCURACY
                t = ACCURACY
                b = (self.area - h * t) / (2 * s + h - h)
                if(b < t + ACCURACY or (h + 2 * s > HOR_MAX)):
                    break
            self.area -= ACCURACY
            # print(area)
            
        print("i_z =", self.inertia_z(self.bf, self.sf, self.hf, self.tf))
        print("Mcr =", mcr_cal(self.C1, self.C2, ELASTIC_MODULUS, SHEAR_MODULUS, self.inertia_z(self.bf, self.sf, self.hf, self.tf), 
                                self.warping_constant(self.bf, self.sf, self.hf, self.tf), self.inertia_t(self.bf, self.sf, self.hf, self.tf), LENGTH - 20, 
                                self.hf / 2 + self.sf, self.k, 1.0))
        print("bf = {}\nsf = {}\nhf = {}\ntf = {}".format(self.bf, self.sf, self.hf, self.tf))
        print("Inertia = {}".format(self.i_hb_max))
        print("--Normal Stress Max = {}\n---Shear Stress Max = {}".format(self.normal_stress_max, self.shear_stress_max))
        print("Area =", self.areaf)
        print("Force_max = ", self.F_hb_max)
        print("fm_max = ", self.fm)
        print("H_Beam finished ----------------------------------------------")
        f.close()

    def inertia_z(self, b, s, h, t):
        i_z = 2 * 1/12 * s * b**3 + 1/12 * h * t**3
        return i_z

    def inertia_t(self, b, s, h, t):

        i_t = (2 * b * s**3 + (h + s) * t**3) / 3
        return i_t
    
    def polar_inertia(self, b, s, h, t):

        Ixx = h**3 * t / 12 + 2 * (s**3 * b / 12 + s * b * (h + s)**2 / 4)
        Iyy = t**3 * h / 12 + 2 * (b**3 * s / 12)
        j = Ixx + Iyy
        return j

    def torsion_constant(self, b, s, h, t):

        j = (1.3/3) * (h * t**3 + 2 * b * s**3)
        return j
    
    def c(self, b, s, h, t):
        
        c = (b * h**2 * t**2)/(6 * (b/2 - t/2)) 
        return c
    
    def  warping_constant(self, b, s, h, t):

        i_w = (h + s)**2 * b**3 * s / 24
        return i_w
    
    def experiment_theta(self, be, se, he, te):

        i_y = 1/12 * he * te**3
        M_y = FRACTURE_STRESS * i_y / (te/2)
        F_y = M_y / (he * 0.75) * 10**(-6)
        theta = math.asin(F_y / EXPERIMENTAL_MAX_FORCE)

        return theta
    

HollowRectangular().cal()
# HollowSquare().cal()
H_Beam().cal()