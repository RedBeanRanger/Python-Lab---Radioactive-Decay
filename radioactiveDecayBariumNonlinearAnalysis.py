# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 21:35:01 2021

@author: angel
"""

import numpy as np
from math import e
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

bkgSample, bkgCount = np.loadtxt(r"C:\Users\angel\Desktop\Docs\PHY244\RadioactiveDecay\Barium_20s_20m_background_mf2021.txt", 
                      unpack = True, skiprows = 2)
bariumSample, bariumCount = np.loadtxt(r"C:\Users\angel\Desktop\Docs\PHY244\RadioactiveDecay\Barium_20s_20m_sample_mf2021.txt", 
                      unpack = True, skiprows = 2)

noBkgCount = np.array(bariumCount - np.mean(bkgCount))
err = np.sqrt(np.add(bariumCount, bkgCount))
dt = 20
time = bariumSample * dt # convert sample number to time in seconds

plt.errorbar(time, noBkgCount, yerr = err, marker = "o", ls = "none", label = "Barium Nonlinear Decay")
plt.title("Barium Radioactive Decay Nonlinear Analysis")
plt.xlabel("Time (s)")
plt.ylabel("Radiation Count")
plt.legend()

# Model function
def calcNonlinearDecay(x, a, b):
    return b*(e**(a*x))

# Find optimal parameters
popt, pcov = curve_fit(calcNonlinearDecay, xdata = time, ydata = noBkgCount, sigma = err,
                       p0 = [-0.005, 1000])
print("The optimal parameters are: ", popt)

# Find standard deviation in halflife.
varParamA = np.sqrt(np.diag(pcov))[0] # variance in parameter a
varHalfLife = varParamA / ((popt[0])**2) # variance in half life
print("The standard deviation in half life is", np.sqrt(varHalfLife))

# Rate calculations
print("The rates are: ", noBkgCount/dt)

# Plotting the fitted curve
model_function = []

for t in time:
    model_function.append((calcNonlinearDecay(t, popt[0], popt[1])))
    
plt.errorbar(time, model_function, label = "Curve Fit")
plt.legend()

# Reduced chi squared calculations
def redChiSquared(measured_data, predicted_data, err):
    # returns chiSquared value
    # takes in arrays of measured_data, predicted_data, and err
    # assume the arrays have the same length
    sum = 0
    for x in range(0, len(measured_data)):
        #sum += ((10**(measured_data[x]) - 10**(predicted_data[x])) / 10**(err[x]))**2
        sum += (((measured_data[x]) - (predicted_data[x])) / (err[x]))**2
    return sum/(len(measured_data) - 2)

plt.show()

print("The reduced chi squared value is: ", redChiSquared(noBkgCount, model_function, err)) 