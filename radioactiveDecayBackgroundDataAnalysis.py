# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 21:46:47 2021

@author: angel
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson
from scipy.stats import norm
from math import sqrt

bkgSample, bkgCount = np.loadtxt(r"C:\Users\angel\Desktop\Docs\PHY244\RadioactiveDecay\Fiesta_6s_10m_background_mf2021.txt", 
                      unpack = True, skiprows = 2)

avgBkgRadiation = np.mean(bkgCount) # average background radiation
print ("Mean background radiation count is: ", avgBkgRadiation)

plt.hist(bkgCount, bins = 12, density = True, label = "Background Radiation Histogram")

#uncertainty in background radiation is uncertainty in uncertainty in precision of equipment
print("Uncertainty is plus or minus", 1)

poissonPlot = poisson.pmf(bkgSample, avgBkgRadiation)
gaussianPlot = norm.pdf(bkgSample, avgBkgRadiation, sqrt(avgBkgRadiation))

plt.plot(poissonPlot, label = "Poisson Plot")
plt.plot(gaussianPlot, label = "Gaussian Plot")

ax = plt.axes()
ax.set_xlim(-2, 8)
plt.title("Distribution of Background Radiation")
plt.ylabel("Frequency of Occurrence (Normalized)")
plt.xlabel("Count of Radiation")

plt.legend()