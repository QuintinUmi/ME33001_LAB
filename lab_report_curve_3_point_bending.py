##########################################################################################################
# Copyright: QIN Qijun
# 
# This document is the property of Qin Qijun, please cite the source if you need to quote it.
# Workshop: https://github.com/QuintinUmi
# Contact:  qqj030212@gmail.com
##########################################################################################################

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import lab_report_tool_package.read_report_file as lr
import lab_report_tool_package.curve_analyze as lc

#Please adjust your parameters according to the following instructions

x_min = -0.05   # Set the data graph display size
x_max = 5.
y_min = 0.
y_max = 75.
isSteel = False #Whether it is metal

#Set the source data path
dataLine = lr.read_file_split_data("C:\\Users\\qqj03\\Desktop\\year2sem2\\ME33001\\lab session.txt")
# print(data)

stress = np.array(lr.get_colume_data(dataLine, 4))
strain = np.array(lr.get_colume_data(dataLine, 5))

#The last parameter is the amount of data pre-fitted for the calculation of the line fit (to be adjusted manually)
young_mod, intercept, r_value, p_value, dataIndex = lc.linear_analyze(strain, stress, 0,  100)  

#The last parameter is the amount of data pre-fitted for the calculation of the line fit (to be adjusted manually)
x_elaslim, y_elaslim = lc.elastic_limit(strain, stress, young_mod * strain + intercept, 100)

x_tensile, y_tensile = lc.tensile_point(strain, stress)

if(isSteel):
    x_yield, y_yield = lc.steel_yield(strain, stress, x_elaslim, x_tensile)
else:
                        #If you get a (division by zero) error, try turning the last parameter up a bit
    x_yield, y_yield = lc.non_steel_yield(strain, young_mod * (strain - 0.2) + intercept, stress, 0.01)

x_fracture, y_fracture = lc.fracture_point(strain, stress, x_tensile)

r_mod = lc.modulus(strain, stress, x_elaslim)
t_mod = lc.modulus(strain, stress, x_fracture)

degree = 10
sIndex = lc.findIndex(strain, x_fracture)
coeff = lc.curve_fit_coeff(strain, stress, dataIndex, sIndex, degree)

#=================================================================================================================
# Above are all the calculation, following are just for drawing picture
#=================================================================================================================

x_linear = np.arange(0., x_max, 0.001)
y_linear = young_mod * x_linear + intercept
x_offset = x_linear
y_offset = young_mod * (x_linear - 0.2) + intercept
x_curve = np.arange(strain[lc.findIndex(strain, x_elaslim)], x_fracture, 0.001)
y_curve = lc.curve_fit_generate(x_curve, coeff)



# stress = stress[0: ]
# strain = strain[0: ]

print(len(stress))
print(len(strain))

plt.figure(figsize=(16, 7.5))
plt.title("G04_Acrylic")
plt.xlabel("Strain(%)")
plt.ylabel("Stress(MPa)")
plt.style.use('seaborn')

plt.grid(True)

y_min = stress[0]

plt.axis([x_min, x_max, y_min, y_max])


plt.plot(strain, stress, label="Experimental data")
plt.plot(x_linear, y_linear, c='b', linestyle='--', label="Linear Elastic")
if(not isSteel):
    plt.plot(x_offset, y_offset, c='g', linestyle='--', label="0.2% Proof Yield Strength")
plt.plot(x_curve, y_curve, c='orange', linestyle='--', label="Fit curve")

plt.fill_between(x_linear[lc.findIndex(y_linear, 0): lc.findIndex(x_linear, x_elaslim)], 0, 
                y_linear[lc.findIndex(y_linear, 0): lc.findIndex(x_linear, x_elaslim)], color="orange", alpha = 0.5)
plt.fill_between(strain[0: lc.findIndex(strain, x_fracture)], 0, 
                stress[0: lc.findIndex(strain, x_fracture)],
                    color="gray", alpha = 0.5)
plt.text(x_elaslim , y_elaslim - 20, "Resilience of Modulus\n\nUr = {} (MPa)".format(r_mod), size=10,
            bbox=dict(boxstyle='round', facecolor='#A9D7E7', alpha=0.7))
plt.text((x_min + x_max) / 2.2, (y_min + y_max) / 2, "Toughness of Modulus\n\nUt = {} (MPa)".format(t_mod), size=10,
            bbox=dict(boxstyle='round', facecolor='#A9D7E7', alpha=0.7))

plt.legend(loc='lower right', prop=None, fontsize = 12, frameon=True)

plt.scatter(x_elaslim, y_elaslim, c='black',)
plt.annotate("Elastic Limit\nStrain = " + str(round(x_elaslim, 4)) + "(%)\n" + "Stress = " + str(round(y_elaslim, 4)) + "(MPa)", 
        xy=(x_elaslim, y_elaslim), xytext=(x_elaslim + 0.2, y_elaslim - 10), arrowprops=dict(facecolor='k', headwidth=5, width=1))

plt.scatter(x_yield, y_yield, c='black',)
plt.annotate("Yield Point\nStrain = " + str(round(x_yield, 4)) + "(%)\n" + "Stress = " + str(round(y_yield, 4)) + "(MPa)", 
        xy=(x_yield, y_yield), xytext=(x_yield + 0.2, y_yield - 10), arrowprops=dict(facecolor='k', headwidth=5, width=1))

plt.scatter(x_tensile, y_tensile, c='black',)
plt.annotate("Tensile Point\nStrain = " + str(round(x_tensile, 4)) + "(%)\n" + "Stress = " + str(round(y_tensile, 4)) + "(MPa)", 
        xy=(x_tensile, y_tensile), xytext=(x_tensile, y_tensile - 10), arrowprops=dict(facecolor='k', headwidth=5, width=1))

plt.scatter(x_fracture, y_fracture, c='black',)
plt.annotate("Fracture Point\nStrain = " + str(round(x_fracture, 4)) + "(%)\n" + "Stress = " + str(round(y_fracture, 4)) + "(MPa)", 
        xy=(x_fracture, y_fracture), xytext=(x_fracture - 0.6, y_fracture - 10), arrowprops=dict(facecolor='k', headwidth=5, width=1))

plt.xticks(np.arange(0., x_max, x_max / 25))
plt.yticks(np.arange(y_min, y_max, y_max / 25))



#----------------------------------------------------------------------------------------------------

plt.figure(figsize=(12, 5.625))
plt.title("{} Report Summery".format("G04_Acryclic"))
plt.xlabel("Strain(%)")
plt.ylabel("Stress(MPa)")
plt.style.use('seaborn')
plt.grid(True)
y_min = stress[0]
plt.axis([x_min, x_max, y_min, y_max])

plt.plot(x_linear, y_linear, c='b', linestyle='--', label="Linear Elastic")
if(not isSteel):
    plt.plot(x_offset, y_offset, c='g', linestyle='--', label="0.2% Proof Strength")
plt.plot(x_curve, y_curve, c='orange', linestyle='--', label="Fit curve")
plt.legend(loc='lower right', prop=None, fontsize = 12, frameon=True)

plt.scatter(x_elaslim, y_elaslim, c='black',)
plt.annotate("Elastic Limit\nStrain = " + str(round(x_elaslim, 4)) + "(%)\n" + "Stress = " + str(round(y_elaslim, 4)) + "(MPa)", 
        xy=(x_elaslim, y_elaslim), xytext=(x_elaslim + 0.2, y_elaslim - 10), arrowprops=dict(facecolor='k', headwidth=5, width=1))

plt.scatter(x_yield, y_yield, c='black',)
plt.annotate("Yield Point\nStrain = " + str(round(x_yield, 4)) + "(%)\n" + "Stress = " + str(round(y_yield, 4)) + "(MPa)", 
        xy=(x_yield, y_yield), xytext=(x_yield + 0.2, y_yield - 10), arrowprops=dict(facecolor='k', headwidth=5, width=1))

plt.scatter(x_tensile, y_tensile, c='black',)
plt.annotate("Tensile Point\nStrain = " + str(round(x_tensile, 4)) + "(%)\n" + "Stress = " + str(round(y_tensile, 4)) + "(MPa)", 
        xy=(x_tensile, y_tensile), xytext=(x_tensile - 1, y_tensile - 10), arrowprops=dict(facecolor='k', headwidth=5, width=1))

plt.scatter(x_fracture, y_fracture, c='black',)
plt.annotate("Fracture Point\nStrain = " + str(round(x_fracture, 4)) + "(%)\n" + "Stress = " + str(round(y_fracture, 4)) + "(MPa)", 
        xy=(x_fracture, y_fracture), xytext=(x_fracture - 0.6, y_fracture - 10), arrowprops=dict(facecolor='k', headwidth=5, width=1))

opStr = \
"                        Strain(%)      Stress(MPa)\n\
Elastic Limit:    {}          {}\n\
Yield Point:      {}              {}\n\
Tensile Point    {}          {}\n\
Fracture Point: {}          {}\n\
\n\
Young's Modulus:               {} (MPa)\n\
Resilience Modulus:           {} (MPa)\n\
Toughness Modulus:          {} (MPa)".format(
            x_elaslim, y_elaslim,
            x_yield, y_yield,
            x_tensile, y_tensile,
            x_fracture, y_fracture,
            young_mod,
            r_mod, 
            t_mod
        )


cur_expr = "Curve Expression (Coefficient will save in cur_coeff.txt in your data path):\n\
            x = strain\n\
            stress = f(x) = {}x + {}, 0 <= x < {}\n\
                                = a_0*x^n + a_1*x + ... + a_(n-1)*x^1 + a_n*x^0, {} <= x <= {}\n\
            n = {}\n\
            a(0~n) = {}".format(
                young_mod, intercept, x_elaslim,
                x_elaslim, x_fracture,
                degree,
                coeff
            )

coeff_str = "\
Young's modulus: {}\n\
intercept: {}\n\
n = {}\n\
a = {}".format(young_mod, intercept, degree, coeff)
f = open("C:\\Users\\qqj03\\Desktop\\Lab Result\\G04_Acrylic_cur_coeff.txt", 'w')
f.write(opStr + '\n\n' + coeff_str)
f.close()

plt.text((x_min + x_max) / 2.2, (y_min + y_max) / 2.5, opStr, size=10,
            bbox=dict(boxstyle='round', facecolor='#A9D7E7', alpha=0.7))
plt.text((x_min + x_max) / 3, (y_min + y_max) / 7, cur_expr, size=8,
            bbox=dict(boxstyle='round', facecolor='#A9D7E7', alpha=0.7))

for i in range(0, len(strain)):
    print("strain: {}  |  stress: {}".format(strain[i], stress[i]))
plt.show()

os.system("PAUSE")