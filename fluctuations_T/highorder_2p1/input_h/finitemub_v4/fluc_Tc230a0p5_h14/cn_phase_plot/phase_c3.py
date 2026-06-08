#!/usr/bin/env python
# -*- coding: utf-8 -*-
# sphinx_gallery_thumbnail_number = 3

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter  # useful for `logit` scale
import matplotlib.ticker as ticker
import matplotlib as mpl
#from scipy.interpolate import spline
from matplotlib import cm
from matplotlib import axes
from matplotlib.font_manager import FontProperties
import pylab as pl
from matplotlib.colors import LinearSegmentedColormap

mpl.style.use('classic')
# Data for plotting
c2data=np.loadtxt('./c2data.dat')
c3data=np.loadtxt('./c3data.dat')
c4data=np.loadtxt('./c4data.dat')
pb=np.loadtxt('./phaseboundary.dat')
mubdata=np.loadtxt('./freezeout/mub.dat')
FOAndro=np.loadtxt('./freezeout/Tandro.dat')
FOI=np.loadtxt('./freezeout/STARFitI.dat')
FOII=np.loadtxt('./freezeout/STARFitII.dat')
# Create figure
fig=plt.figure(figsize=(5, 3.5))
colors = ['#053061', '#134b87', '#327db7', '#6fafd2', '#c7e0ed', '#fbd2bc', '#feab88', '#b71c2c', '#8b0824']  # 深蓝-白-深红
cmap_custom = LinearSegmentedColormap.from_list('mymap', colors)
####################################################################################################
ax3=fig.add_subplot(111)
vd=0.0001
vu=0.01
im=ax3.imshow(-c3data.T, cmap=cmap_custom,interpolation='nearest',vmin=vd,vmax=vu,zorder=3,extent=[0, 555, 21, 300], origin='lower', aspect=1.7)
vnorm = mpl.colors.Normalize(vmin=vd, vmax=vu)
plt.rcParams['font.size'] = 7
cbar=plt.colorbar(im,fraction=0.029, pad=0.02,norm=vnorm)
cbar.set_label(r'$-c_3$', rotation=0,fontsize=12,labelpad=10)
#ax2.set_aspect(0.8) 
plt.scatter(620,99,color='red',marker='o',s=30,label=r'CEP',zorder=3)
pb1,=ax3.plot(pb[:,0],pb[:,1],dashes=[1,0],color='k',linewidth=1.2,alpha=0.6,label=r'Crossover',zorder=3)
star4,=ax3.plot(mubdata,FOI,dashes=[4,1,2,1],color='b',linewidth=1.2,alpha=0.6,label=r'freezeout: STAR Fit I',zorder=4)
star4,=ax3.plot(mubdata,FOII,dashes=[5,2],color='g',linewidth=1.2,alpha=0.6,label=r'freezeout: STAR Fit II',zorder=4)
And,=ax3.plot(mubdata,FOAndro,color='r',linewidth=1.2,alpha=0.6,label=r'freezeout: Andronic et al.',zorder=3)
plt.axis([0.,630.,80.,250.])
ax3.set_xlabel('$\mu_B\,[\mathrm{MeV}]$', fontsize=12, color='black')
ax3.set_ylabel(r'$T\,[\mathrm{MeV}]$', fontsize=12, color='black')
ax3.set_aspect(2.2)
for spine in ax3.spines.values():
    spine.set_zorder(10)
ax3.tick_params(which='both', zorder=10)

for label in ax3.xaxis.get_ticklabels():
    label.set_fontsize(8)
for label in ax3.yaxis.get_ticklabels():
    label.set_fontsize(8)
ax3.legend(loc=0,fontsize=7,frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
fig.subplots_adjust(top=1., bottom=0.0, left=0.112, right=0.85, hspace=0.1,
                    wspace=0.1)

fig.savefig("phasec3.pdf",dpi=300,transparent=True)