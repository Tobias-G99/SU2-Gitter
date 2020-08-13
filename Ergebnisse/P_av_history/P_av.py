import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so

x = np.loadtxt('P_av.txt', usecols=0, skiprows=3)
y = np.loadtxt('P_av.txt', usecols=1, skiprows=3)
c = np.zeros(x.size)
c = (c+1)*0.603945

plt.xticks(np.arange(0, max(x)+1, 5.0))

plt.xlabel(r'Nummer des Werts $P$')
plt.ylabel(r'$P$')
plt.ylim(min(y)-0.002,max(y)+0.002)
plt.scatter(x,y,marker='+')
plt.plot(x,c,color='red',label=r'$\left\langle P \right\rangle$ = 0,60395 $\pm$ 0,00018')
plt.legend(loc='upper right')
plt.show()
