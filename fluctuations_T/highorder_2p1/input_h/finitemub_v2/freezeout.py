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
Andro_c2=np.loadtxt('./Andro_c2.dat')
SFI_c2=np.loadtxt('./STARFitI_c2.dat')
SFII_c2=np.loadtxt('./STARFitII_c2.dat')
ymin_c2 = np.minimum(np.minimum(Andro_c2[:,1], SFI_c2[:,1]), SFII_c2[:,1])
ymax_c2 = np.maximum(np.maximum(Andro_c2[:,1], SFI_c2[:,1]), SFII_c2[:,1])
Andro_c3=np.loadtxt('./Andro_c3.dat')
SFI_c3=np.loadtxt('./STARFitI_c3.dat')
SFII_c3=np.loadtxt('./STARFitII_c3.dat')
ymin_c3 = np.minimum(np.minimum(Andro_c3[:,1], SFI_c3[:,1]), SFII_c3[:,1])
ymax_c3 = np.maximum(np.maximum(Andro_c3[:,1], SFI_c3[:,1]), SFII_c3[:,1])
Andro_c4=np.loadtxt('./Andro_c4.dat')
SFI_c4=np.loadtxt('./STARFitI_c4.dat')
SFII_c4=np.loadtxt('./STARFitII_c4.dat')
ymin_c4 = np.minimum(np.minimum(Andro_c4[:,1], SFI_c4[:,1]), SFII_c4[:,1])
ymax_c4 = np.maximum(np.maximum(Andro_c4[:,1], SFI_c4[:,1]), SFII_c4[:,1])
################################################################################################################
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(Andro_c2[:,0],Andro_c2[:,1],linewidth=2,alpha=1,c='#A0D3FF',label=r'$\mathrm{Freeze-out:\,Andronic\,et\,al.}$') 
ax1.plot(SFI_c2[:,0],SFI_c2[:,1],linewidth=2,alpha=1,c='#FF8AD0',label=r'$\mathrm{Freeze-out:\,STAR\,Fit\,I}$') 
ax1.plot(SFII_c2[:,0],SFII_c2[:,1],linewidth=2,alpha=1,c='#A886F6',label=r'$\mathrm{Freeze-out:\,STAR\,Fit\,II}$') 
ax1.fill_between(Andro_c2[:,0],ymin_c2,ymax_c2,color='red',alpha=0.1)
ax1.set_ylabel(r'$c_2$', fontsize=13, color='black')
ax1.set_xlabel(r'$\sqrt{s_{NN}}\,[\mathrm{GeV}]$', fontsize=13, color='black')
ax1.set_title(r'$\mathrm{QCD-assisted\,\,LEFT}$', fontsize=13, color='black')
ax1.legend(loc=0,fontsize='10',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
ax1.axis([4,100,0.02,0.09])
ax1.set_xscale('log')
for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.17, right=0.95, hspace=0.33,wspace=0.2)

fig.savefig("./c2_freezeout_model.pdf")
################################################################################################################
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(Andro_c3[:,0],Andro_c3[:,1],linewidth=2,alpha=1,c='#A0D3FF',label=r'$\mathrm{Freeze-out:\,Andronic\,et\,al.}$') 
ax1.plot(SFI_c3[:,0],SFI_c3[:,1],linewidth=2,alpha=1,c='#FF8AD0',label=r'$\mathrm{Freeze-out:\,STAR\,Fit\,I}$') 
ax1.plot(SFII_c3[:,0],SFII_c3[:,1],linewidth=2,alpha=1,c='#A886F6',label=r'$\mathrm{Freeze-out:\,STAR\,Fit\,II}$') 
ax1.fill_between(Andro_c3[:,0],ymin_c3,ymax_c3,color='red',alpha=0.1)
ax1.set_ylabel(r'$c_3$', fontsize=13, color='black')
ax1.set_xlabel(r'$\sqrt{s_{NN}}\,[\mathrm{GeV}]$', fontsize=13, color='black')
ax1.set_title(r'$\mathrm{QCD-assisted\,\,LEFT}$', fontsize=13, color='black')
ax1.legend(loc=0,fontsize='10',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
ax1.axis([4,100,-0.05,-0.001])
ax1.set_xscale('log')
for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.17, right=0.95, hspace=0.33,wspace=0.2)

fig.savefig("./c3_freezeout_model.pdf")
################################################################################################################
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(Andro_c4[:,0],Andro_c4[:,1],linewidth=2,alpha=1,c='#A0D3FF',label=r'$\mathrm{Freeze-out:\,Andronic\,et\,al.}$') 
ax1.plot(SFI_c4[:,0],SFI_c4[:,1],linewidth=2,alpha=1,c='#FF8AD0',label=r'$\mathrm{Freeze-out:\,STAR\,Fit\,I}$') 
ax1.plot(SFII_c4[:,0],SFII_c4[:,1],linewidth=2,alpha=1,c='#A886F6',label=r'$\mathrm{Freeze-out:\,STAR\,Fit\,II}$') 
ax1.fill_between(Andro_c4[:,0],ymin_c4,ymax_c4,color='red',alpha=0.1)
ax1.set_ylabel(r'$c_4$', fontsize=13, color='black')
ax1.set_xlabel(r'$\sqrt{s_{NN}}\,[\mathrm{GeV}]$', fontsize=13, color='black')
ax1.set_title(r'$\mathrm{QCD-assisted\,\,LEFT}$', fontsize=13, color='black')
ax1.legend(loc=0,fontsize='10',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
ax1.axis([4,100,-0.015,0.09])
ax1.set_xscale('log')
for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.17, right=0.95, hspace=0.33,wspace=0.2)

fig.savefig("./c4_freezeout_model.pdf")
################################################################################################################
