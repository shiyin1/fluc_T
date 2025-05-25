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
mfmub0=np.loadtxt('./data/mfmub0.dat')
mfmub100=np.loadtxt('./data/mfmub100.dat')
mfmub200=np.loadtxt('./data/mfmub200.dat')
mfmub300=np.loadtxt('./data/mfmub300.dat')
mfmub400=np.loadtxt('./data/mfmub400.dat')
mfmub500=np.loadtxt('./data/mfmub500.dat')
mfmub550=np.loadtxt('./data/mfmub550.dat')
dmfdTmub0=np.loadtxt('./data/dmfdTmub0.dat')
dmfdTmub100=np.loadtxt('./data/dmfdTmub100.dat')
dmfdTmub200=np.loadtxt('./data/dmfdTmub200.dat')
dmfdTmub300=np.loadtxt('./data/dmfdTmub300.dat')
dmfdTmub400=np.loadtxt('./data/dmfdTmub400.dat')
dmfdTmub500=np.loadtxt('./data/dmfdTmub500.dat')
dmfdTmub550=np.loadtxt('./data/dmfdTmub550.dat')
l=np.loadtxt('./data/l.dat')
dldT=np.loadtxt('./data/dldTmub0.dat')
T=np.arange(1, 251, 1)
################################################################################################################
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(T,mfmub0,'g',linewidth=1.5,alpha=0.7,zorder=3,label=r'$\mu_B=0$')
ax1.plot(T,mfmub100,'r',linewidth=1.5,dashes=[5,2],alpha=0.7,zorder=4,label=r'$\mu_B=100\,\mathrm{MeV}$')
ax1.plot(T,mfmub200,'b',linewidth=1.5,dashes=[4,1,2,1],alpha=0.7,zorder=5,label=r'$\mu_B=200\,\mathrm{MeV}$')
ax1.plot(T,mfmub300,'m',linewidth=1.5,dashes=[2,2],alpha=0.7,zorder=6,label=r'$\mu_B=300\,\mathrm{MeV}$')
ax1.plot(T,mfmub400,'y',linewidth=1.5,dashes=[3,1,2,1],alpha=0.7,zorder=7,label=r'$\mu_B=400\,\mathrm{MeV}$')
ax1.plot(T,mfmub500,'c',linewidth=1.5,dashes=[2,1],alpha=0.7,zorder=8,label=r'$\mu_B=500\,\mathrm{MeV}$')
ax1.plot(T,mfmub550,'k',linewidth=1.5,dashes=[1,1],alpha=0.5,zorder=9,label=r'$\mu_B=550\,\mathrm{MeV}$')
ax1.set_ylabel(r'$m_l\,[\mathrm{MeV}]$', fontsize=13, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=13, color='black')
ax1.legend(loc=0,fontsize='9.5',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
ax1.axis([1,250,1,360])

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.135, right=0.95, hspace=0.33,wspace=0.2)

fig.savefig("./ml.pdf")
################################################################################################################
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(T,dmfdTmub0,'g',linewidth=1.5,alpha=0.7,zorder=3,label=r'$\mu_B=0$')
ax1.plot(T,dmfdTmub100,'r',linewidth=1.5,dashes=[5,2],alpha=0.7,zorder=4,label=r'$\mu_B=100\,\mathrm{MeV}$')
ax1.plot(T,dmfdTmub200,'b',linewidth=1.5,dashes=[4,1,2,1],alpha=0.7,zorder=5,label=r'$\mu_B=200\,\mathrm{MeV}$')
ax1.plot(T,dmfdTmub300,'m',linewidth=1.5,dashes=[2,2],alpha=0.7,zorder=6,label=r'$\mu_B=300\,\mathrm{MeV}$')
ax1.plot(T,dmfdTmub400,'y',linewidth=1.5,dashes=[3,1,2,1],alpha=0.7,zorder=7,label=r'$\mu_B=400\,\mathrm{MeV}$')
ax1.plot(T,dmfdTmub500,'c',linewidth=1.5,dashes=[2,1],alpha=0.7,zorder=8,label=r'$\mu_B=500\,\mathrm{MeV}$')
ax1.plot(T,dmfdTmub550,'k',linewidth=1.5,dashes=[1,1],alpha=0.5,zorder=9,label=r'$\mu_B=550\,\mathrm{MeV}$')
ax1.set_ylabel(r'$\partial\,m_l/\partial\,T$', fontsize=13, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=13, color='black')
ax1.legend(loc=0,fontsize='9',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
ax1.axis([1,250,0.,10.])

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.13, right=0.95, hspace=0.33,wspace=0.2)

fig.savefig("./dmfdT.pdf")
################################################################################################################
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(T,l,'g',linewidth=1.5,alpha=0.7,zorder=3,label=r'$\mu_B=0$')
ax1.set_ylabel(r'$l$', fontsize=13, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=13, color='black')
ax1.legend(loc=0,fontsize='9',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
ax1.axis([1,250,0.,1.])

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.13, right=0.95, hspace=0.33,wspace=0.2)

fig.savefig("./l.pdf")
################################################################################################################
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(T,dldT*100,'g',linewidth=1.5,alpha=0.7,zorder=3,label=r'$\mu_B=0$')
ax1.set_ylabel(r'$\partial\,l/\partial\,T\,[\times 10^{-2}]$', fontsize=13, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=13, color='black')
ax1.legend(loc=0,fontsize='9',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
ax1.axis([1,250,0.,0.012*100])

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.13, right=0.95, hspace=0.33,wspace=0.2)

fig.savefig("./dldT.pdf")