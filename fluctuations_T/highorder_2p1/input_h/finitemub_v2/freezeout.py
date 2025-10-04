#!/usr/bin/env python
# -*- coding: utf-8 -*-
# sphinx_gallery_thumbnail_number = 3

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter  # useful for `logit` scale
import matplotlib.ticker as ticker
import matplotlib as mpl
from scipy.interpolate import make_interp_spline
from matplotlib.ticker import FuncFormatter
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.cm as cm 
from matplotlib.font_manager import FontProperties
import pylab as pl
from matplotlib.ticker import FixedLocator
from decimal import Decimal, getcontext
getcontext().prec = 20

mpl.style.use('classic')
# Data for plotting
# mub=0 data
Andro=np.loadtxt('./Andro.dat')
SFI=np.loadtxt('./STARFitI.dat')
SFII=np.loadtxt('./STARFitII.dat')
################################################################################################################
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(Andro[:,0],Andro[:,1],linewidth=2,alpha=0.8,label=r'$\mathrm{Freeze-out:\,Andronic\,et\,al.}$') 
ax1.plot(SFI[:,0],SFI[:,1],linewidth=2,alpha=0.8,label=r'$\mathrm{Freeze-out:\,STAR\,Fit\,I}$') 
ax1.plot(SFII[:,0],SFII[:,1],linewidth=2,alpha=0.8,label=r'$\mathrm{Freeze-out:\,STAR\,Fit\,II}$') 
ax1.set_ylabel(r'$c_2^{1/2}(\sqrt{s_{NN}})$', fontsize=13, color='black')
ax1.set_xlabel(r'$\sqrt{s_{NN}}\,[\mathrm{GeV}]$', fontsize=13, color='black')
ax1.set_title(r'$\mathrm{QCD-assisted\,\,LEFT}$', fontsize=13, color='black')
ax1.legend(loc=0,fontsize='10',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
ax1.axis([4,100,0.15,0.3])
ax1.set_xscale('log')
for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.17, right=0.95, hspace=0.33,wspace=0.2)

fig.savefig("./c2_freezeout_model.pdf")