import os

MAX_AREA = 500.
ACCURACY = 1.5
LENGTH_LIMIT_COE = 500
LL = LENGTH_LIMIT_COE * MAX_AREA**0.5
DENSITY = 1.05

FRACTURE_STRESS = 33 * 10**6
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

        b = d = h = k = ACCURACY #hollow rectengular

        while(self.area > 4 * ACCURACY**2):
            b = LL
            d = (self.area + h*k) / b
            h = b - 2 * ACCURACY
            k = d - 2 * ACCURACY
            mass = 100 * self.area * DENSITY
            while(b > 4):
                while(b > 4):
                    while(b > 4):
                        k -= ACCURACY
                        if(k <= 0):
                            break
                        
                        
                        b = (self.area + h*k) / d   
                        i_temp = (b * d**3 - h * k**3) / 12
                        M_temp = FRACTURE_STRESS * i_temp / (d/2)
                        F_temp = M_temp / 20
                        # print(b, d, h, k, F_temp, mass)
                        
                        if(F_temp/mass > self.fm):
                            # print(b, d, h, k, F_temp, mass, F_temp/mass)
                            self.fm = F_temp/mass
                            
                            self.i_hr_max = i_temp
                            self.M_hr_max = M_temp
                            self.F_hr_max = F_temp
                            
                            self.bf = b
                            self.df = d
                            self.hf = h
                            self.kf = k
                        
                        
                        if(b > LL or b <= 4):
                            break
                    h -= ACCURACY
                    k = d - 2 * ACCURACY
                    b = (self.area + h*k) / d  
                    if(h <= 0):
                        break
                    if(b <= 4):
                        break
                d += ACCURACY
                h = b - 2 * ACCURACY
                k = d - 2 * ACCURACY
                if(d > LL):
                    break
            self.area -= ACCURACY
            # print(area)
            

        print("bf = {}\ndf = {}\nhf = {}\nkf = {}".format(self.bf, self.df, self.hf, self.kf))
        print("Inertia = {}".format(self.i_hr_max))
        print("fm_max = ", self.fm)
        print("hollow rectengular finished ----------------------------------------------")



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
                    F_temp = M_temp / 20
                    

                    if(F_temp/mass > self.fm):

                        self.fm = F_temp/mass
                                
                        self.i_hs_max = i_temp
                        self.M_hs_max = M_temp
                        self.F_hs_max = F_temp
                                
                        self.af = a
                        self.bf = b

                    b -= ACCURACY
                    area_temp = a**2 - b**2

                a -= ACCURACY

            self.area -= ACCURACY


        print("af = {}\nbf = {}".format(self.af, self.bf))
        print("Inertia = {}".format(self.i_hs_max))
        print("fm_max = ", self.fm)
        print("hollow square finished ----------------------------------------------")


class H_Beam:
    def __init__(self):
        self.area = MAX_AREA
        self.i_hb_max = 0.
        self.M_hb_max = 0.
        self.F_hb_max = 0.
        self.bf = self.sf = self.hf = self.tf = 0.

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
                        s += ACCURACY
                        if(b < t + ACCURACY):
                            break
                                               
                        b = (self.area - h * t) / (2 * s + h - h)
                        i_temp = (b * (2 * s + h)**3 - h**3 * (b - t)) / 12
                        M_temp = FRACTURE_STRESS * i_temp / (h/2 + s)
                        F_temp = M_temp / 20
                        # print(h, b, t, i_temp, F_temp, mass)
                        # print(b, t)
                        
                        if(F_temp/mass > self.fm and b > t + ACCURACY):
                            # print(b, d, h, k, F_temp, mass, F_temp/mass)
                            self.fm = F_temp/mass
                            
                            self.i_hb_max = i_temp
                            self.M_hb_max = M_temp
                            self.F_hb_max = F_temp
                            
                            self.bf = b
                            self.sf = s
                            self.hf = h
                            self.tf = t
                        
                    t += ACCURACY
                    s = ACCURACY
                    b = (self.area - h * t) / (2 * s + h - h)
                    if(b < t + ACCURACY):
                        break

                h += ACCURACY
                t = ACCURACY
                s = ACCURACY
                b = (self.area - h * t) / (2 * s + h - h)
                if(b < t + ACCURACY):
                    break
            self.area -= ACCURACY
            # print(area)
            

        print("bf = {}\nsf = {}\nhf = {}\ntf = {}".format(self.bf, self.sf, self.hf, self.tf))
        print("Inertia = {}".format(self.i_hb_max))
        print("fm_max = ", self.fm)
        print("hollow rectengular finished ----------------------------------------------")

# HollowRectangular().cal()
HollowSquare().cal()
# H_Beam().cal()