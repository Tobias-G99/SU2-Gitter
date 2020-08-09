import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so

x = np.loadtxt('P_av_vs_beta.txt', usecols=0, skiprows=1)
y = np.loadtxt('P_av_vs_beta.txt', usecols=1, skiprows=1)
dy = np.loadtxt('P_av_vs_beta.txt', usecols=2, skiprows=1)

plt.xticks(np.arange(min(x), max(x)+0.1, 0.2))

plt.xlabel(r'$\beta$')
plt.ylabel(r'$\left\langle P \right\rangle$')
plt.ylim(min(y)-0.02,max(y)+0.02)
plt.errorbar(x,y,yerr=dy,marker='+')
plt.show()
