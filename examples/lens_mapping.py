from SLens import analyser,gnfwSersic

import numpy as np
import matplotlib.pyplot as plt

fig,ax = plt.subplots(2,4,figsize=(8,4))
Arc = [2.1,1,.39,1e-4]
col = ['orange','lime','r']
for num,arc in enumerate(Arc):
    stat = analyser(dist=arc)
    attr = stat.get_cross_section()
    
    circle1 = plt.Circle((0,0),stat.attr[4],lw=2,facecolor='none',edgecolor='rebeccapurple')
    circle2 = plt.Circle((0,0),stat.attr[0],lw=2,facecolor='none',edgecolor='rebeccapurple')
    circle3 = plt.Circle((0,0),stat.attr[1],lw=2,facecolor='none',edgecolor='firebrick')
    xr = np.sqrt(2)/2
    sol = np.array([attr[3],attr[2],0])
    mag = np.array([attr[5],attr[4],0])
    
    def mag_str(x):
        if x < 100:
            return "%.1f" % x
        else:
            return "inf"
        
    mag_str = [mag_str(x) for x in mag]
    if num == 0:
        ax[0,num].scatter(xr*sol[1],xr*sol[1],marker='+',s=100,lw=1,color='lime')
        ax[0,num].text(xr*sol[1]-1.5,xr*sol[1]-.2,mag_str[1])
    else:   
        ax[0,num].scatter(xr*sol,xr*sol,marker='+',s=100,lw=1,c=col)
        [ax[0,num].text(xr*sol[i]+.4,xr*sol[i]-.2,mag_str[i]) for i in range(len(sol))]
        
    ax[1,num].plot(xr*arc,xr*arc,marker='+',ms=12,color='b')
    ax[0,num].add_patch(circle1)
    ax[0,num].add_patch(circle2)
    ax[1,num].add_patch(circle3)
    
    ax[0,num].set_xlim([-3.3,3.3])
    ax[0,num].set_ylim([-3.3,3.3])
    ax[1,num].set_xlim([-3.3,3.3])
    ax[1,num].set_ylim([-3.3,3.3])
    ax[0,num].set_xlabel(r'$\xi_1[\prime\prime]$')
    ax[0,num].set_ylabel(r'$\xi_2[\prime\prime]$')
    ax[1,num].set_xlabel(r'$\eta_1[\prime\prime]$')
    ax[1,num].set_ylabel(r'$\eta_2[\prime\prime]$')
    ax[0,0].text(-2.5,2.2,'image')
    ax[1,0].text(-2.5,2.22,'source')

fig.tight_layout()
plt.show()
fig.savefig("./plots/source_loc.png",format='eps')