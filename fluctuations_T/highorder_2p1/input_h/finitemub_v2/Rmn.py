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
chi0mub0=np.loadtxt('./mub0/data/chi0.dat')
c2mub0=np.loadtxt('./mub0/data/c2.dat')
c3mub0=np.loadtxt('./mub0/data/c3.dat')
c4mub0=np.loadtxt('./mub0/data/c4.dat')
c5mub0=np.loadtxt('./mub0/data/c5.dat')
c6mub0=np.loadtxt('./mub0/data/c6.dat')
# mub=100 data
chi0mub100=np.loadtxt('./mub100/data/chi0.dat')
c2mub100=np.loadtxt('./mub100/data/c2.dat')
c3mub100=np.loadtxt('./mub100/data/c3.dat')
c4mub100=np.loadtxt('./mub100/data/c4.dat')
c5mub100=np.loadtxt('./mub100/data/c5.dat')
c6mub100=np.loadtxt('./mub100/data/c6.dat')
# mub=200 data
chi0mub200=np.loadtxt('./mub200/data/chi0.dat')
c2mub200=np.loadtxt('./mub200/data/c2.dat')
c3mub200=np.loadtxt('./mub200/data/c3.dat')
c4mub200=np.loadtxt('./mub200/data/c4.dat')
c5mub200=np.loadtxt('./mub200/data/c5.dat')
c6mub200=np.loadtxt('./mub200/data/c6.dat')
# mub=300 data
chi0mub300=np.loadtxt('./mub300/data/chi0.dat')
c2mub300=np.loadtxt('./mub300/data/c2.dat')
c3mub300=np.loadtxt('./mub300/data/c3.dat')
c4mub300=np.loadtxt('./mub300/data/c4.dat')
c5mub300=np.loadtxt('./mub300/data/c5.dat')
c6mub300=np.loadtxt('./mub300/data/c6.dat')
# mub=400 data
chi0mub400=np.loadtxt('./mub400/data/chi0.dat')
c2mub400=np.loadtxt('./mub400/data/c2.dat')
c3mub400=np.loadtxt('./mub400/data/c3.dat')
c4mub400=np.loadtxt('./mub400/data/c4.dat')
c5mub400=np.loadtxt('./mub400/data/c5.dat')
c6mub400=np.loadtxt('./mub400/data/c6.dat')
# mub=500 data
chi0mub500=np.loadtxt('./mub500/data/chi0.dat')
c2mub500=np.loadtxt('./mub500/data/c2.dat')
c3mub500=np.loadtxt('./mub500/data/c3.dat')
c4mub500=np.loadtxt('./mub500/data/c4.dat')
c5mub500=np.loadtxt('./mub500/data/c5.dat')
c6mub500=np.loadtxt('./mub500/data/c6.dat')
# mub=550 data
chi0mub550=np.loadtxt('./mub550/data/chi0.dat')
c2mub550=np.loadtxt('./mub550/data/c2.dat')
c3mub550=np.loadtxt('./mub550/data/c3.dat')
c4mub550=np.loadtxt('./mub550/data/c4.dat')
c5mub550=np.loadtxt('./mub550/data/c5.dat')
c6mub550=np.loadtxt('./mub550/data/c6.dat')
# Temperature
T=np.arange(21, 301, 1)
# pressure
poT4mub0=(chi0mub0-chi0mub0[0])/T**4
poT4mub100=(chi0mub100-chi0mub100[0])/T**4
poT4mub200=(chi0mub200-chi0mub200[0])/T**4
poT4mub300=(chi0mub300-chi0mub300[0])/T**4
poT4mub400=(chi0mub400-chi0mub400[0])/T**4
poT4mub500=(chi0mub500-chi0mub500[0])/T**4
poT4mub550=(chi0mub550-chi0mub550[0])/T**4

blackline=np.linspace(0,300,100)
zeroline=np.zeros(100)
for i in range(1,100):
    zeroline[i]=0
################################################################################################################
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(T,c3mub0/c2mub0**2,'g',linewidth=1.5,alpha=0.7,zorder=3,label=r'$\mu_B=0$')
ax1.plot(T,c3mub100/c2mub100**2,'r',linewidth=1.5,dashes=[5,2],alpha=0.7,zorder=4,label=r'$\mu_B=100\,\mathrm{MeV}$')
ax1.plot(T,c3mub200/c2mub200**2,'b',linewidth=1.5,dashes=[4,1,2,1],alpha=0.7,zorder=5,label=r'$\mu_B=200\,\mathrm{MeV}$')
ax1.plot(T,c3mub300/c2mub300**2,'m',linewidth=1.5,dashes=[2,2],alpha=0.7,zorder=6,label=r'$\mu_B=300\,\mathrm{MeV}$')
ax1.plot(T,c3mub400/c2mub400**2,'y',linewidth=1.5,dashes=[3,1,2,1],alpha=0.7,zorder=7,label=r'$\mu_B=400\,\mathrm{MeV}$')
ax1.plot(T,c3mub500/c2mub500**2,'c',linewidth=1.5,dashes=[2,1],alpha=0.7,zorder=8,label=r'$\mu_B=500\,\mathrm{MeV}$')
ax1.plot(T,c3mub550/c2mub550**2,'k',linewidth=1.5,dashes=[1,1],alpha=0.5,zorder=9,label=r'$\mu_B=550\,\mathrm{MeV}$')
ax1.set_ylabel(r'$R_{32}$', fontsize=17, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=15, color='black')
ax1.legend(loc=0,fontsize='10',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
#ax1.axis([0,300,-10.,-10**-2])
ax1.axis([50,300,-11,0])
#ax1.set_yscale('symlog',linthresh=10**-7)
for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.17, right=0.95, hspace=0.33,wspace=0.2)

fig.savefig("./R32.pdf")
################################################################################################################
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(T,c4mub0/c2mub0**3,'g',linewidth=1.5,alpha=0.7,zorder=3,label=r'$\mu_B=0$')
ax1.plot(T,c4mub100/c2mub100**3,'r',linewidth=1.5,dashes=[5,2],alpha=0.7,zorder=4,label=r'$\mu_B=100\,\mathrm{MeV}$')
ax1.plot(T,c4mub200/c2mub200**3,'b',linewidth=1.5,dashes=[4,1,2,1],alpha=0.7,zorder=5,label=r'$\mu_B=200\,\mathrm{MeV}$')
ax1.plot(T,c4mub300/c2mub300**3,'m',linewidth=1.5,dashes=[2,2],alpha=0.7,zorder=6,label=r'$\mu_B=300\,\mathrm{MeV}$')
ax1.plot(T,c4mub400/c2mub400**3,'y',linewidth=1.5,dashes=[3,1,2,1],alpha=0.7,zorder=7,label=r'$\mu_B=400\,\mathrm{MeV}$')
ax1.plot(T,c4mub500/c2mub500**3,'c',linewidth=1.5,dashes=[2,1],alpha=0.7,zorder=8,label=r'$\mu_B=500\,\mathrm{MeV}$')
ax1.plot(T,c4mub550/c2mub550**3,'k',linewidth=1.5,dashes=[1,1],alpha=0.5,zorder=9,label=r'$\mu_B=550\,\mathrm{MeV}$')
ax1.set_ylabel(r'$R_{42}$', fontsize=17, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=15, color='black')
ax1.legend(loc=0,fontsize='10',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
ax1.axis([50,300,-20.,250])
#ax1.set_yscale('symlog',linthresh=10**-7)
for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.17, right=0.95, hspace=0.33,wspace=0.2)

fig.savefig("./R42.pdf")
################################################################################################################
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(T,c5mub0/c2mub0**4/1000,'g',linewidth=1.5,alpha=0.7,zorder=3,label=r'$\mu_B=0$')
ax1.plot(T,c5mub100/c2mub100**4/1000,'r',linewidth=1.5,dashes=[5,2],alpha=0.7,zorder=4,label=r'$\mu_B=100\,\mathrm{MeV}$')
ax1.plot(T,c5mub200/c2mub200**4/1000,'b',linewidth=1.5,dashes=[4,1,2,1],alpha=0.7,zorder=5,label=r'$\mu_B=200\,\mathrm{MeV}$')
ax1.plot(T,c5mub300/c2mub300**4/1000,'m',linewidth=1.5,dashes=[2,2],alpha=0.7,zorder=6,label=r'$\mu_B=300\,\mathrm{MeV}$')
ax1.plot(T,c5mub400/c2mub400**4/1000,'y',linewidth=1.5,dashes=[3,1,2,1],alpha=0.7,zorder=7,label=r'$\mu_B=400\,\mathrm{MeV}$')
ax1.plot(T,c5mub500/c2mub500**4/1000,'c',linewidth=1.5,dashes=[2,1],alpha=0.7,zorder=8,label=r'$\mu_B=500\,\mathrm{MeV}$')
ax1.plot(T,c5mub550/c2mub550**4/1000,'k',linewidth=1.5,dashes=[1,1],alpha=0.5,zorder=9,label=r'$\mu_B=550\,\mathrm{MeV}$')
ax1.set_ylabel(r'$R_{52}\,[\times 10^3]$', fontsize=17, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=15, color='black')
#ax1.legend(loc=0,fontsize='10',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
ax1.axis([50,300,-8000/1000,1000/1000])
#ax1.set_yscale('symlog',linthresh=10**-10)

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.15, right=0.95, hspace=0.33,wspace=0.2)

fig.savefig("./R52.pdf")
################################################################################################################
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(T,c5mub0/c2mub0**5/10000,'g',linewidth=1.5,alpha=0.7,zorder=3,label=r'$\mu_B=0$')
ax1.plot(T,c5mub100/c2mub0**5/10000,'r',linewidth=1.5,dashes=[5,2],alpha=0.7,zorder=4,label=r'$\mu_B=100\,\mathrm{MeV}$')
ax1.plot(T,c5mub200/c2mub0**5/10000,'b',linewidth=1.5,dashes=[4,1,2,1],alpha=0.7,zorder=5,label=r'$\mu_B=200\,\mathrm{MeV}$')
ax1.plot(T,c5mub300/c2mub0**5/10000,'m',linewidth=1.5,dashes=[2,2],alpha=0.7,zorder=6,label=r'$\mu_B=300\,\mathrm{MeV}$')
ax1.plot(T,c5mub400/c2mub0**5/10000,'y',linewidth=1.5,dashes=[3,1,2,1],alpha=0.7,zorder=7,label=r'$\mu_B=400\,\mathrm{MeV}$')
ax1.plot(T,c5mub500/c2mub0**5/10000,'c',linewidth=1.5,dashes=[2,1],alpha=0.7,zorder=8,label=r'$\mu_B=500\,\mathrm{MeV}$')
ax1.plot(T,c5mub550/c2mub0**5/10000,'k',linewidth=1.5,dashes=[1,1],alpha=0.5,zorder=9,label=r'$\mu_B=550\,\mathrm{MeV}$')
ax1.set_ylabel(r'$R_{62}\,[\times 10^4]$', fontsize=17, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=15, color='black')
ax1.axis([50,300,-2.1*10**5/10000,10**4/10000])
#ax1.set_yscale('symlog',linthresh=10**-12)

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.17, right=0.95, hspace=0.33,wspace=0.2)

fig.savefig("./R62.pdf")