# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 23:33:02 2021

@author: angel
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

bkgSample, bkgCount = np.loadtxt(r"C:\Users\angel\Desktop\Docs\PHY244\RadioactiveDecay\Barium_20s_20m_background_mf2021.txt", 
                      unpack = True, skiprows = 2)
bariumSample, bariumCount = np.loadtxt(r"C:\Users\angel\Desktop\Docs\PHY244\RadioactiveDecay\Barium_20s_20m_sample_mf2021.txt", 
                      unpack = True, skiprows = 2)
#print(bariumCount)

noBkgCount = (bariumCount - np.mean(bkgCount))

# To prevent errors
noBkgCountWithoutInfs = [] # this array stores all the points in noBkgCount whose log will not be infinity
sampleWithoutInfs = [] # this array stores all the corresponding sample numbers of noBkgCountWithoutInfs
bariumCountWithoutInfs = []
bkgCountWithoutInfs = []

# Everything relating to valuesTakenOut is for debugging purposes
#valuesTakenOut = [] # this appends the sample number of the value that was omitted from the curve.
for x in range(0, len(noBkgCount)):
    if noBkgCount[x] != 0:
        noBkgCountWithoutInfs.append(noBkgCount[x])
        sampleWithoutInfs.append(bkgSample[x])
        bariumCountWithoutInfs.append(bariumCount[x])
        bkgCountWithoutInfs.append(bkgCount[x])
    #elif noBkgCount[x] == 0:
    #    valuesTakenOut.append(x + 1) #since x is the index number, sample number = x + 1

noBkgCountWithoutInfs = np.array(noBkgCountWithoutInfs)
sampleWithoutInfs = np.array(sampleWithoutInfs)
bariumCountWithoutInfs = np.array(bariumCountWithoutInfs)
bkgCountWithoutInfs = np.array(bkgCountWithoutInfs)

#print("The number of sample values taken out are: ", len(valuesTakenOut))
#print("The missing sample number is:", valuesTakenOut)
  
"""
For debugging purposes
print(len(noBkgCountWithoutInfs))
print(len(sampleWithoutInfs))
print(len(bariumCountWithoutInfs))
print(len(bkgCountWithoutInfs))

print(np.sqrt(np.add(bariumCountWithoutInfs, bkgCountWithoutInfs)).shape)
print(noBkgCountWithoutInfs.shape)
"""

err = np.sqrt(np.add(bariumCountWithoutInfs, bkgCountWithoutInfs))/noBkgCountWithoutInfs
logNoBkgCount = np.log(noBkgCountWithoutInfs)
dt = 20
time = sampleWithoutInfs*dt # convert sample number to time in seconds

plt.errorbar(time, logNoBkgCount, yerr = err, marker = "o", ls = "none", label = "Barium Linear Decay")
plt.title("Barium Radioactive Decay Linear Analysis")
plt.xlabel("Time (s)")
plt.ylabel("Radiation Count (log base e)")


# Model function

def calcDecayLogarithmic(x, a, b): 
    # Returns ln of I(t), radiation intensity at time t, given the following parameters:
    # b is ln of 1/2 of initial intensity of radiation (log of 1/2 * y_0)
    # a is the negative mean lifetime of the isotope (inverse of half life)
    # x is the time
    return (x*a) + b


popt, pcov = curve_fit(calcDecayLogarithmic, ydata = logNoBkgCount, xdata = time, 
                       p0 = [-0.005, 3])

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
    model_function.append((calcDecayLogarithmic(t, popt[0], popt[1])))
    #model_function.append((calcDecay(t, popt[0], popt[1])))
    
plt.errorbar(time, model_function, label = "Curve Fit")
plt.legend()

#Reduced chi squared calculations
def redChiSquared(measured_data, predicted_data, err):
    # returns chiSquared value
    # takes in arrays of measured_data, predicted_data, and err
    # assume the arrays have the same length
    sum = 0
    for x in range(0, len(measured_data)):
        #sum += ((10**(measured_data[x]) - 10**(predicted_data[x])) / 10**(err[x]))**2
        sum += (((measured_data[x]) - (predicted_data[x])) / (err[x]))**2
    return sum/(len(measured_data) - 2)

print("The reduced chi squared value is: ", redChiSquared(logNoBkgCount, model_function, err)) 
