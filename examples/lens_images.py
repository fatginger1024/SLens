from SLens import analyser

import numpy as np
import matplotlib.pyplot as plt

col = ['darkred','olive','darkcyan']
arc = .39
stat = analyser(dist=arc)
pop = stat.get_cross_section()

sol = np.array([pop[2],pop[3]])
alpha_sol = [np.sign(x)*stat.lens_alpha(np.abs(x/(stat.thetas * 206265)))[0] * stat.thetas * 206265 for x in sol]

x = np.concatenate((-stat._x[::-1],stat._x),axis=None) * stat.thetas * 206265
Y = np.concatenate((-stat.lens_alpha(stat._x)[::-1],stat.lens_alpha(stat._x)),axis=None) * stat.thetas * 206265

x1 = np.repeat(sol[0],100)
x2 = np.repeat(sol[1],100)
y1 = np.linspace(0,alpha_sol[0],100)
y2 = np.linspace(0,alpha_sol[1],100)

fig,ax = plt.subplots(1,1,figsize=(5,5))
ax.plot(x,Y,label=r'$Y=\alpha(\xi)$',color=col[0])
ax.plot(x1,y1,ls='--',color='k')
ax.plot(x2,y2,ls='--',color='k')

ax.text(sol[0],-.2,r'$\xi_1$')
ax.text(sol[1]-.7,.2,r'$-\xi_2$')
ax.text(arc-1,-.2,r'$\eta$')

        
ax.plot(x,x-arc,label=r'$Y=\xi-\eta$',zorder=1,color='crimson')
ax.plot(x,x,color='peachpuff',label=r'$Y=\xi$')
ax.plot(arc,0,ls='none',marker='o',color='indigo',alpha=.7,label='source'+r'$\, \eta$')
ax.plot(sol[0],alpha_sol[0],marker='*',ms=12,ls='none',color='lime',zorder=2,label='image 1')
ax.plot(sol[1],alpha_sol[1],marker='*',ms=12,ls='none',color='orange',zorder=2,label='image 2')


ax.axvline(x=0,ls='--',color='k')
ax.axhline(y=0,ls='--',color='k')
ax.set_xlim([-10,10])
ax.set_ylim([-2.5,2.5])
ax.set_xlabel(r'$\xi[\prime\prime]$')
ax.set_ylabel(r'$Y$')
ax.legend()
plt.show()
fig.savefig("./plots/illus.png",format='png')