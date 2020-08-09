import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so

data = np.loadtxt('L_av_2.050000.txt', usecols=1, skiprows=5)
data = abs(data)
plt.xticks(np.arange(0, 0.6, 0.05))

plt.xlabel(r'|$L$|')
plt.ylabel(r'P(|$L$|)')
plt.xlim(0,0.55)
plt.ylim(0,0.3)
plt.hist(x=data,bins=50,range=(0,0.5),weights=np.ones(len(data)) / len(data))
plt.legend(loc='upper right')
plt.savefig(fname='L_abs_dist_2.05.png', format='png')
plt.show()
