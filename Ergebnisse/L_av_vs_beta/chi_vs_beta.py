import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so

x = np.loadtxt('chi_vs_beta.txt', usecols=0, skiprows=1)
y = np.loadtxt('chi_vs_beta.txt', usecols=1, skiprows=1)
dy = np.loadtxt('chi_vs_beta.txt', usecols=2, skiprows=1)

plt.xticks(np.arange(min(x), max(x), 0.05))

plt.xlabel(r'$\beta$')
plt.ylabel(r'$\chi_{N_s}$')
plt.ylim(0,max(y)+3)
plt.errorbar(x,y,yerr=dy,marker='+', ls='none')
plt.legend(loc='upper right')
plt.show()
