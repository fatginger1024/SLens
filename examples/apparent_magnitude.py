from SLens import analyser

import numpy as np
import matplotlib.pyplot as plt



Arc = np.linspace(.1,2.,50)
Mag1 = np.zeros(len(Arc))
Mag2 = np.zeros(len(Arc))
M1 = np.zeros(len(Arc))
M2 = np.zeros(len(Arc))

for i in range(len(Arc)):
    stat = analyser(dist=Arc[i],source_mag=26)
    mag0 = stat.mag_unlensed
    mag_mu = lambda mu: mag0 - np.log10(mu)
    pop = stat.get_cross_section()
    #print(pop)
    Mag1[i] = np.abs(pop[4])
    Mag2[i] = np.abs(pop[5])
    
    M1[i] = mag_mu(np.abs(pop[4]))
    M2[i] = mag_mu(np.abs(pop[5]))

ind_mag = np.argmin(np.abs(M2-26.3))


fig,(ax1,ax2) = plt.subplots(2,1,figsize=(4,8))

ax1.plot(Arc,Mag1,label='image 1',color='lime')
ax1.plot(Arc,Mag2,label='image 2',color='orange')
ax1.set_xlabel(r'$\eta[\prime\prime]$')
ax1.set_ylabel(r'$\mu$')
ax1.legend()

ax2.plot(Arc,M1,label='image 1',color='lime',zorder=1)
ax2.plot(Arc,M2,label='image 2',color='orange',zorder=1)
ax2.scatter(Arc[ind_mag],M2[ind_mag],marker='o',s=40,facecolor='rebeccapurple',edgecolor='k',zorder=2)
ax2.fill_between(np.linspace(Arc.min(),Arc[ind_mag],100),24.8,27.4,color='plum',alpha=.8)
ax2.set_xlabel(r'$\eta[\prime\prime]$')
ax2.set_ylabel(r'$m$')
ax2.legend()
fig.tight_layout()
plt.show()
fig.savefig("./plots/mag_images.png",format='png')