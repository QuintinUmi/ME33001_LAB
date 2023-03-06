MAX_AREA = 500.
ACCURACY = 1.5
LENGTH_LIMIT_COE = 500
LL = LENGTH_LIMIT_COE * MAX_AREA**0.5
DENSITY = 1.05

FRACTURE_STRESS = 33 * 10**6
F = 0.
M = 0.

b = d = h = k = ACCURACY #hollow rectengular

a = b = ACCURACY #squre hallow

d = t = a = s = n = ACCURACY #I beam

b = d = h = t = ACCURACY #H beam





#hollow rectengular
area = MAX_AREA
i_hr_max = 0.
stress_max = 0.
bf = df = hf = kf = 0.

fm = 0.

b = LL
d = (area + h*k) / b
h = b - 2 * ACCURACY
k = d - 2 * ACCURACY
while(area > 16):
    mass = 100 * area * DENSITY
    while(b > 4):
        while(b > 4):
            while(b > 4):
                k -= ACCURACY
                if(k <= 0):
                    break
                
                
                b = (area + h*k) / d   
                i_temp = (b * d**3 - h * k**3) / 12
                M_temp = FRACTURE_STRESS * i_temp / (d/2)
                F_temp = M_temp / 20
                print(F_temp, mass)
                
                if(F_temp/mass > fm):
                    fm = F_temp/mass
                    
                    i_hr_max = i_temp
                    M_hr_max = M_temp
                    F_hr_max = F_temp
                    
                    bf = b
                    df = d
                    hf = h
                    kf = h
                
                
                if(b > LL or b <= 4):
                    break
            h -= ACCURACY
            k = d - 2 * ACCURACY
            if(h <= 0):
                break
            if(b <= 4):
                break
        d += ACCURACY
        h = b - 2 * ACCURACY
        k = d - 2 * ACCURACY
        if(d > LL):
            break
    area -= ACCURACY
    

print("bf = {}\ndf = {}\nhf = {}\nkf = {}".format(bf, df, hf, kf))
print("Inertia = {}".format(i_hr_max))
print("fm_max = ", fm)
print("finish")


