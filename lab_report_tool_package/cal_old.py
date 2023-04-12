import os

MAX_AREA = 500.
ACCURACY = 2
LENGTH_LIMIT_COE = 500
LL = LENGTH_LIMIT_COE * MAX_AREA**0.5
HOR_MAX = 39
DENSITY = 1.05 * 10**(-6)
DELH_COE = 0.5

FRACTURE_STRESS = 31 * 10**6
SHEAR_STRENGTH = 1 / (3**(1/2)) * FRACTURE_STRESS
F = 0.
M = 0.



d = t = a = s = n = ACCURACY #I beam


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
                            F_normal_temp = M_temp / 20 *10**(-9)

                            Q_shear = (b*d/2 - k*h/2) * ((b)*(d/2)*(d/4) - (k/2)*h*(k/4)) / (b*d/2 - k*h/2)
                            F_shear_temp = SHEAR_STRENGTH * (b - h) * i_temp / Q_shear *10**(-3)

                            F_temp = min(F_normal_temp, F_shear_temp)
                            # print(F_normal_temp > F_shear_temp)
                            # print(F_shear_temp, F_normal_temp)

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
        print("fm_max = ", self.fm)
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
                area_temp = a**2 - b**2
                while(area_temp <= self.area and b >= 0.):
                    
                    mass = 100 * area_temp * DENSITY

                    i_temp = (a**4 - b**4) / 12
                    M_temp = FRACTURE_STRESS * i_temp / (a/2)
                    F_normal_temp = M_temp / 20 *10**(-9)

                    Q_shear = (a*a/2 - b*b/2) * ((a)*(a/2)*(a/4) - (b/2)*b*(b/4)) / (a*a/2 - b*b/2)
                    F_shear_temp = SHEAR_STRENGTH * (a - b) * i_temp / Q_shear *10**(-6)

                    F_temp = min(F_normal_temp, F_shear_temp)
                    # print(F_normal_temp > F_shear_temp)
                    
                    # if(i_temp > self.i_hr_max):
                    if(F_temp/mass > self.fm and a < HOR_MAX):

                        self.fm = F_temp/mass
                                
                        self.i_hs_max = i_temp
                        self.M_hs_max = M_temp
                        self.F_hs_max = F_temp
                                
                        self.af = a
                        self.bf = b

                    b -= ACCURACY
                    area_temp = a**2 - b**2

                a -= ACCURACY
                if(a >= HOR_MAX):
                    break

            self.area -= ACCURACY


        print("af = {}\nbf = {}".format(self.af, self.bf))
        print("Inertia = {}".format(self.i_hs_max))
        print("fm_max = ", self.fm)
        print("Hollow Square finished ----------------------------------------------")


class H_Beam:
    def __init__(self):
        self.area = MAX_AREA
        self.i_hb_max = 0.
        self.M_hb_max = 0.
        self.F_hb_max = 0.
        self.bf = self.sf = self.hf = self.tf = 0.

        self.areaf = MAX_AREA
        self.normal_stress_max = 0.
        self.shear_stress_max = 0.

        self.fm = 0.

    def cal(self):

        b = s = h = t = ACCURACY #H beam

        while(self.area > 4 * ACCURACY**2):
            h = ACCURACY
            t = ACCURACY
            s = ACCURACY
            b = (self.area - h * t) / (2 * s + h - h)
            mass = 100 * self.area * DENSITY
            while(b > t + ACCURACY):
                while(b > t + ACCURACY):
                    while(b > t + ACCURACY):
                                                                     
                        i_temp = (b * (2 * s + h)**3 - h**3 * (b - t)) / 12
                        M_temp = FRACTURE_STRESS * i_temp / (h/2 + s)
                        F_normal_temp = M_temp / 20 *10**(-9)
                        
                        Q_shear = s*b*(h+s)/2 + h*h*t/8
                        F_shear_temp = 2 * SHEAR_STRENGTH * t * i_temp / Q_shear *10**(-9)

                        F_temp = min(F_normal_temp, F_shear_temp) * 10**3
                        
                        
                        # if(i_temp > self.i_hr_max):
                        if(F_temp/mass > self.fm and b > t + ACCURACY and b < HOR_MAX and h + 2 * s < HOR_MAX):
                            # print(b, d, h, k, F_temp, mass, F_temp/mass)
                            print(i_temp, M_temp, F_normal_temp, F_shear_temp)
                            
                            self.fm = F_temp/mass
                            
                            self.i_hb_max = i_temp
                            self.M_hb_max = M_temp
                            self.F_hb_max = F_temp
                            
                            self.bf = b
                            self.sf = s
                            self.hf = h
                            self.tf = t

                            self.areaf = self.area
                            self.normal_stress_max = 20 * F_temp * (h + s)/2 / i_temp *10**(9)
                            self.shear_stress_max = F_temp * Q_shear / (i_temp * t) *10**(6)
                        
                        s += ACCURACY
                        b = (self.area - h * t) / (2 * s + h - h)
                        if(b < HOR_MAX or h + 2 * s < HOR_MAX):
                            break
                        
                    t += ACCURACY
                    s = ACCURACY
                    b = (self.area - h * t) / (2 * s + h - h)
                    if(b < HOR_MAX or h + 2 * s < HOR_MAX):
                        break

                h += ACCURACY
                t = ACCURACY
                s = ACCURACY
                b = (self.area - h * t) / (2 * s + h - h)
                if(b < t + ACCURACY and (b < HOR_MAX or h + 2 * s < HOR_MAX)):
                    break
            self.area -= ACCURACY
            # print(area)
            

        print("bf = {}\nsf = {}\nhf = {}\ntf = {}".format(self.bf, self.sf, self.hf, self.tf))
        print("Inertia = {}".format(self.i_hb_max))
        print("--Normal Stress Max = {}\n---Shear Stress Max = {}".format(self.normal_stress_max, self.shear_stress_max))
        print("Area =", self.areaf)
        print("Force_max = ", self.F_hb_max)
        print("fm_max = ", self.fm)
        print("H_Beam finished ----------------------------------------------")

# HollowRectangular().cal()
# HollowSquare().cal()
H_Beam().cal()