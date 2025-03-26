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
c0mub0=np.loadtxt('./mub0/data/chi0.dat')
c1mub0=np.loadtxt('./mub0/data/chi1.dat')
c2mub0=np.loadtxt('./mub0/data/chi2.dat')
c3mub0=np.loadtxt('./mub0/data/chi3.dat')
c4mub0=np.loadtxt('./mub0/data/chi4.dat')
c5mub0=np.loadtxt('./mub0/data/chi5.dat')
c6mub0=np.loadtxt('./mub0/data/chi6.dat')
# mub=100 data
c0mub100=np.loadtxt('./mub100/data/chi0.dat')
c1mub100=np.loadtxt('./mub100/data/chi1.dat')
c2mub100=np.loadtxt('./mub100/data/chi2.dat')
c3mub100=np.loadtxt('./mub100/data/chi3.dat')
c4mub100=np.loadtxt('./mub100/data/chi4.dat')
c5mub100=np.loadtxt('./mub100/data/chi5.dat')
c6mub100=np.loadtxt('./mub100/data/chi6.dat')
# mub=200 data
c0mub200=np.loadtxt('./mub200/data/chi0.dat')
c1mub200=np.loadtxt('./mub200/data/chi1.dat')
c2mub200=np.loadtxt('./mub200/data/chi2.dat')
c3mub200=np.loadtxt('./mub200/data/chi3.dat')
c4mub200=np.loadtxt('./mub200/data/chi4.dat')
c5mub200=np.loadtxt('./mub200/data/chi5.dat')
c6mub200=np.loadtxt('./mub200/data/chi6.dat')
# mub=300 data
c0mub300=np.loadtxt('./mub300/data/chi0.dat')
c1mub300=np.loadtxt('./mub300/data/chi1.dat')
c2mub300=np.loadtxt('./mub300/data/chi2.dat')
c3mub300=np.loadtxt('./mub300/data/chi3.dat')
c4mub300=np.loadtxt('./mub300/data/chi4.dat')
c5mub300=np.loadtxt('./mub300/data/chi5.dat')
c6mub300=np.loadtxt('./mub300/data/chi6.dat')
# mub=400 data
c0mub400=np.loadtxt('./mub400/data/chi0.dat')
c1mub400=np.loadtxt('./mub400/data/chi1.dat')
c2mub400=np.loadtxt('./mub400/data/chi2.dat')
c3mub400=np.loadtxt('./mub400/data/chi3.dat')
c4mub400=np.loadtxt('./mub400/data/chi4.dat')
c5mub400=np.loadtxt('./mub400/data/chi5.dat')
c6mub400=np.loadtxt('./mub400/data/chi6.dat')
# mub=500 data
c0mub500=np.loadtxt('./mub500/data/chi0.dat')
c1mub500=np.loadtxt('./mub500/data/chi1.dat')
c2mub500=np.loadtxt('./mub500/data/chi2.dat')
c3mub500=np.loadtxt('./mub500/data/chi3.dat')
c4mub500=np.loadtxt('./mub500/data/chi4.dat')
c5mub500=np.loadtxt('./mub500/data/chi5.dat')
c6mub500=np.loadtxt('./mub500/data/chi6.dat')
# mub=550 data
c0mub550=np.loadtxt('./mub550/data/chi0.dat')
c1mub550=np.loadtxt('./mub550/data/chi1.dat')
c2mub550=np.loadtxt('./mub550/data/chi2.dat')
c3mub550=np.loadtxt('./mub550/data/chi3.dat')
c4mub550=np.loadtxt('./mub550/data/chi4.dat')
c5mub550=np.loadtxt('./mub550/data/chi5.dat')
c6mub550=np.loadtxt('./mub550/data/chi6.dat')
# Temperature
T=np.arange(21, 301, 1)
################################################################################################################
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(T,c1mub0/T**3,linewidth=2,alpha=0.8,label=r'$\mu_B=0$') 
ax1.plot(T,c1mub100/T**3,linewidth=2,alpha=0.8,label=r'$\mu_B=100\,\mathrm{MeV}$') 
ax1.plot(T,c1mub200/T**3,linewidth=2,alpha=0.8,label=r'$\mu_B=200\,\mathrm{MeV}$') 
ax1.plot(T,c1mub300/T**3,linewidth=2,alpha=0.8,label=r'$\mu_B=300\,\mathrm{MeV}$') 
ax1.plot(T,c1mub400/T**3,linewidth=2,alpha=0.8,label=r'$\mu_B=400\,\mathrm{MeV}$') 
ax1.plot(T,c1mub500/T**3,linewidth=2,alpha=0.8,label=r'$\mu_B=500\,\mathrm{MeV}$') 
ax1.plot(T,c1mub550/T**3,linewidth=2,alpha=0.8,label=r'$\mu_B=550\,\mathrm{MeV}$') 
ax1.set_ylabel(r'$T^{-3}\partial p(T,\mu_B)/\partial T$', fontsize=13, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=13, color='black')
ax1.legend(loc=0,fontsize='8',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
#ax1.axis([0,300,-10.,-10**-2])
ax1.axis([50,300,0,22])
#ax1.set_yscale('symlog',linthresh=10**-7)
for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.17, right=0.95, hspace=0.33,wspace=0.2)

fig.savefig("./chi1.pdf")
################################################################################################################
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(T,c2mub0/T**2,linewidth=2,alpha=0.8,label=r'$\mu_B=0$') 
ax1.plot(T,c2mub100/T**2,linewidth=2,alpha=0.8,label=r'$\mu_B=100\,\mathrm{MeV}$') 
ax1.plot(T,c2mub200/T**2,linewidth=2,alpha=0.8,label=r'$\mu_B=200\,\mathrm{MeV}$') 
ax1.plot(T,c2mub300/T**2,linewidth=2,alpha=0.8,label=r'$\mu_B=300\,\mathrm{MeV}$') 
ax1.plot(T,c2mub400/T**2,linewidth=2,alpha=0.8,label=r'$\mu_B=400\,\mathrm{MeV}$') 
ax1.plot(T,c2mub500/T**2,linewidth=2,alpha=0.8,label=r'$\mu_B=500\,\mathrm{MeV}$') 
ax1.plot(T,c2mub550/T**2,linewidth=2,alpha=0.8,label=r'$\mu_B=550\,\mathrm{MeV}$') 
ax1.set_ylabel(r'$T^{-2}\partial^2p(T,\mu_B)/\partial T^2$', fontsize=13, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=13, color='black')
#ax1.legend(loc=0,fontsize='10',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
#ax1.axis([0,300,-10.,-10**-2])
ax1.axis([50,300,0,75])
#ax1.set_yscale('symlog',linthresh=10**-7)
for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.17, right=0.95, hspace=0.33,wspace=0.2)

fig.savefig("./chi2.pdf")
################################################################################################################
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(T,c3mub0/T,linewidth=2,alpha=0.8) 
ax1.plot(T,c3mub100/T,linewidth=2,alpha=0.8) 
ax1.plot(T,c3mub200/T,linewidth=2,alpha=0.8) 
ax1.plot(T,c3mub300/T,linewidth=2,alpha=0.8) 
ax1.plot(T,c3mub400/T,linewidth=2,alpha=0.8) 
ax1.plot(T,c3mub500/T,linewidth=2,alpha=0.8) 
ax1.plot(T,c3mub550/T,linewidth=2,alpha=0.8) 
ax1.set_ylabel(r'$T^{-1}\partial^3 p(T,\mu_B)/\partial T^3$', fontsize=13, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=13, color='black')
#ax1.legend(loc=0,fontsize='10',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
ax1.axis([50,300,0.,260])
#ax1.set_yscale('symlog',linthresh=10**-7)
for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.17, right=0.95, hspace=0.33,wspace=0.2)

fig.savefig("./chi3.pdf")
################################################################################################################
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(T,c4mub0,linewidth=2,alpha=0.8) 
ax1.plot(T,c4mub100,linewidth=2,alpha=0.8) 
ax1.plot(T,c4mub200,linewidth=2,alpha=0.8) 
ax1.plot(T,c4mub300,linewidth=2,alpha=0.8) 
ax1.plot(T,c4mub400,linewidth=2,alpha=0.8) 
ax1.plot(T,c4mub500,linewidth=2,alpha=0.8) 
ax1.plot(T,c4mub550,linewidth=2,alpha=0.8) 
ax1.set_ylabel(r'$\partial^4p(T,\mu_B)/\partial T^4$', fontsize=13, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=13, color='black')
#ax1.legend(loc=0,fontsize='10',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
ax1.axis([50,300,-3000,2000])
#ax1.set_yscale('symlog',linthresh=10**-10)
ax1.set_xticks([50,100,150,200,250,300])  # 先设定刻度位置
ax1.set_xticklabels(["50", "100", "150", "200", "250", "300"], fontsize=8) 
ax1.set_yticks([-3000,-2000,-1000,0,1000,2000])  # 先设定刻度位置
ax1.set_yticklabels([r'$-3\times 10^{3}$',r'$-2\times 10^{3}$',r'$-1\times 10^{3}$',r'$0$',r'$1\times 10^{3}$',r'$2\times 10^{3}$'], fontsize=8) 

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.2, right=0.95, hspace=0.33,wspace=0.2)

fig.savefig("./chi4.pdf")
################################################################################################################
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(T,c5mub0*T,linewidth=2,alpha=0.8) 
ax1.plot(T,c5mub100*T,linewidth=2,alpha=0.8) 
ax1.plot(T,c5mub200*T,linewidth=2,alpha=0.8) 
ax1.plot(T,c5mub300*T,linewidth=2,alpha=0.8) 
ax1.plot(T,c5mub400*T,linewidth=2,alpha=0.8) 
ax1.plot(T,c5mub500*T,linewidth=2,alpha=0.8) 
ax1.plot(T,c5mub550*T,linewidth=2,alpha=0.8) 
ax1.set_ylabel(r'$T\,\partial^5p(T,\mu_B)/\partial T^5$', fontsize=13, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=13, color='black')
ax1.axis([50,300,-1.1*10**5,10**5])
#ax1.set_yscale('symlog',linthresh=10**-12)
ax1.set_xticks([50,100,150,200,250,300])  # 先设定刻度位置
ax1.set_xticklabels(["50", "100", "150", "200", "250", "300"], fontsize=8) 
ax1.set_yticks([-100000,-50000,0,50000,100000])  # 先设定刻度位置
ax1.set_yticklabels([r'$-1\times 10^{5}$',r'$-5\times 10^{4}$',r'$0$',r'$5\times 10^{4}$',r'$1\times 10^{5}$'], fontsize=8) 

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.23, right=0.95, hspace=0.33,wspace=0.2)

fig.savefig("./chi5.pdf")
################################################################################################################
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(T,c6mub0*T**2,linewidth=2,alpha=0.8) 
ax1.plot(T,c6mub100*T**2,linewidth=2,alpha=0.8) 
ax1.plot(T,c6mub200*T**2,linewidth=2,alpha=0.8) 
ax1.plot(T,c6mub300*T**2,linewidth=2,alpha=0.8) 
ax1.plot(T,c6mub400*T**2,linewidth=2,alpha=0.8) 
ax1.plot(T,c6mub500*T**2,linewidth=2,alpha=0.8) 
ax1.plot(T,c6mub550*T**2,linewidth=2,alpha=0.8) 
ax1.set_ylabel(r'$T^2\,\partial^6p(T,\mu_B)/\partial T^6$', fontsize=13, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=13, color='black')
ax1.axis([50,300,-2*10**6,3*10**6])
#ax1.set_yscale('symlog',linthresh=10**-12)
ax1.set_xticks([50,100,150,200,250,300])  # 先设定刻度位置
ax1.set_xticklabels(["50", "100", "150", "200", "250", "300"], fontsize=8) 
ax1.set_yticks([-2000000,-1000000,0,1000000,2000000,3000000])  # 先设定刻度位置
ax1.set_yticklabels([r'$-2\times 10^{6}$',r'$-1\times 10^{6}$',r'$0$',r'$1\times 10^{6}$',r'$2\times 10^{6}$',r'$3\times 10^{6}$'], fontsize=8) 

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.21, right=0.95, hspace=0.33,wspace=0.2)

fig.savefig("./chi6.pdf")