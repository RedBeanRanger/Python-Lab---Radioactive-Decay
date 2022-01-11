# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 12:08:56 2021

@author: angel
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson
from scipy.stats import norm
from math import sqrt

# Calculate background radiation
# Assume background radiation is a constant, and that the average value of the 
# background radiation sample is a sufficient estimate. 

bkgCount = np.loadtxt(r"C:\Users\angel\Desktop\Docs\PHY244\RadioactiveDecay\Fiesta_6s_10m_background_mf2021.txt", 
                      unpack = True, skiprows = 2, usecols = 1)


avgBkgRadiation = np.mean(bkgCount) # average background radiation
#print(avgBkgRadiation) # for debugging purposes

fiestaSample, fiestaCount = np.loadtxt(r"C:\Users\angel\Desktop\Docs\PHY244\RadioactiveDecay\Fiesta_6s_10m_sample_mf2021.txt", 
                                       skiprows = 2, unpack = True)
#print(fiestaSample)
#print(fiestaCount)
fCountWithoutBkg = fiestaCount - avgBkgRadiation # fiesta radiation count without the background radiation

plt.hist(fCountWithoutBkg, bins = 12, density = True, label = "Fiesta Decay Histogram")

#uncertainty calculation
err = np.sqrt(np.add(avgBkgRadiation, bkgCount))
print("Absolute uncertainties for each data point is as follows:", err)

#make a poisson probability mass function
mu = np.mean(fCountWithoutBkg)
print("Mean radiation count is: ", mu)

poissonPlot = poisson.pmf(fiestaSample, mu)
gaussianPlot = norm.pdf(fiestaSample, mu, sqrt(mu))
#print(poissonPlot)

plt.plot(poissonPlot, label = "Poisson Plot")
plt.plot(gaussianPlot, label = "Gaussian Plot")
plt.title("Distribution of Radiation from Fiesta Plate")
plt.ylabel("Frequency of Occurrence (Normalized)")
plt.xlabel("Count of Radiation")
plt.legend()