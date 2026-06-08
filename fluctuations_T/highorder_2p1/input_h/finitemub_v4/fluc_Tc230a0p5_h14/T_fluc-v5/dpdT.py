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

blackline=np.linspace(0,300,100)
zeroline=np.zeros(100)
for i in range(1,100):
    zeroline[i]=0
################################################################################################################
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(T,c1mub0/T**3,'g',linewidth=1.5,alpha=0.7,zorder=3,label=r'$\mu_B=0$')
ax1.plot(T,c1mub100/T**3,'r',linewidth=1.5,dashes=[5,2],alpha=0.7,zorder=4,label=r'$\mu_B=100\,\mathrm{MeV}$')
ax1.plot(T,c1mub200/T**3,'b',linewidth=1.5,dashes=[4,1,2,1],alpha=0.7,zorder=5,label=r'$\mu_B=200\,\mathrm{MeV}$')
ax1.plot(T,c1mub300/T**3,'m',linewidth=1.5,dashes=[2,2],alpha=0.7,zorder=6,label=r'$\mu_B=300\,\mathrm{MeV}$')
ax1.plot(T,c1mub400/T**3,'y',linewidth=1.5,dashes=[3,1,2,1],alpha=0.7,zorder=7,label=r'$\mu_B=400\,\mathrm{MeV}$')
ax1.plot(T,c1mub500/T**3,'c',linewidth=1.5,dashes=[2,1],alpha=0.7,zorder=8,label=r'$\mu_B=500\,\mathrm{MeV}$')
#ax1.plot(T,c1mub550/T**3,'k',linewidth=1.5,dashes=[1,1],alpha=0.5,zorder=9,label=r'$\mu_B=550\,\mathrm{MeV}$')
ax1.set_ylabel(r'$\chi_1$', fontsize=17, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=15, color='black')
ax1.legend(loc=0,fontsize='8',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
#ax1.axis([0,300,-10.,-10**-2])
ax1.axis([50,300,0,21])
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
ax1.plot(T,c2mub0/T**2,'g',linewidth=1.5,alpha=0.7,zorder=3,label=r'$\mu_B=0$')
ax1.plot(T,c2mub100/T**2,'r',linewidth=1.5,dashes=[5,2],alpha=0.7,zorder=4,label=r'$\mu_B=100\,\mathrm{MeV}$')
ax1.plot(T,c2mub200/T**2,'b',linewidth=1.5,dashes=[4,1,2,1],alpha=0.7,zorder=5,label=r'$\mu_B=200\,\mathrm{MeV}$')
ax1.plot(T,c2mub300/T**2,'m',linewidth=1.5,dashes=[2,2],alpha=0.7,zorder=6,label=r'$\mu_B=300\,\mathrm{MeV}$')
ax1.plot(T,c2mub400/T**2,'y',linewidth=1.5,dashes=[3,1,2,1],alpha=0.7,zorder=7,label=r'$\mu_B=400\,\mathrm{MeV}$')
ax1.plot(T,c2mub500/T**2,'c',linewidth=1.5,dashes=[2,1],alpha=0.7,zorder=8,label=r'$\mu_B=500\,\mathrm{MeV}$')
#ax1.plot(T,c2mub550/T**2,'k',linewidth=1.5,dashes=[1,1],alpha=0.5,zorder=9,label=r'$\mu_B=550\,\mathrm{MeV}$')
ax1.set_ylabel(r'$\chi_2$', fontsize=17, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=15, color='black')
ax1.legend(loc=0,fontsize='10',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
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
ax1.plot(T,c3mub0/T,'g',linewidth=1.5,alpha=0.7,zorder=3,label=r'$\mu_B=0$')
ax1.plot(T,c3mub100/T,'r',linewidth=1.5,dashes=[5,2],alpha=0.7,zorder=4,label=r'$\mu_B=100\,\mathrm{MeV}$')
ax1.plot(T,c3mub200/T,'b',linewidth=1.5,dashes=[4,1,2,1],alpha=0.7,zorder=5,label=r'$\mu_B=200\,\mathrm{MeV}$')
ax1.plot(T,c3mub300/T,'m',linewidth=1.5,dashes=[2,2],alpha=0.7,zorder=6,label=r'$\mu_B=300\,\mathrm{MeV}$')
ax1.plot(T,c3mub400/T,'y',linewidth=1.5,dashes=[3,1,2,1],alpha=0.7,zorder=7,label=r'$\mu_B=400\,\mathrm{MeV}$')
ax1.plot(T,c3mub500/T,'c',linewidth=1.5,dashes=[2,1],alpha=0.7,zorder=8,label=r'$\mu_B=500\,\mathrm{MeV}$')
#ax1.plot(T,c3mub550/T,'k',linewidth=1.5,dashes=[1,1],alpha=0.5,zorder=9,label=r'$\mu_B=550\,\mathrm{MeV}$')
ax1.set_ylabel(r'$\chi_3$', fontsize=17, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=15, color='black')
ax1.legend(loc=0,fontsize='8',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
ax1.axis([50,300,0.,250])
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
ax1.plot(T,c4mub0,'g',linewidth=1.5,alpha=0.7,zorder=3,label=r'$\mu_B=0$')
ax1.plot(T,c4mub100,'r',linewidth=1.5,dashes=[5,2],alpha=0.7,zorder=4,label=r'$\mu_B=100\,\mathrm{MeV}$')
ax1.plot(T,c4mub200,'b',linewidth=1.5,dashes=[4,1,2,1],alpha=0.7,zorder=5,label=r'$\mu_B=200\,\mathrm{MeV}$')
ax1.plot(T,c4mub300,'m',linewidth=1.5,dashes=[2,2],alpha=0.7,zorder=6,label=r'$\mu_B=300\,\mathrm{MeV}$')
ax1.plot(T,c4mub400,'y',linewidth=1.5,dashes=[3,1,2,1],alpha=0.7,zorder=7,label=r'$\mu_B=400\,\mathrm{MeV}$')
ax1.plot(T,c4mub500,'c',linewidth=1.5,dashes=[2,1],alpha=0.7,zorder=8,label=r'$\mu_B=500\,\mathrm{MeV}$')
#ax1.plot(T,c4mub550,'k',linewidth=1.5,dashes=[1,1],alpha=0.5,zorder=9,label=r'$\mu_B=550\,\mathrm{MeV}$')
ax1.plot(blackline,zeroline,dashes=[4,2],color='k',linewidth='1',alpha=0.7)
ax1.set_ylabel(r'$\chi_4$', fontsize=17, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=15, color='black')
ax1.legend(loc=0,fontsize='8',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
ax1.axis([50,300,-2000,2000])
#ax1.set_yscale('symlog',linthresh=10**-10)
ax1.set_xticks([50,100,150,200,250,300])  # 先设定刻度位置
ax1.set_xticklabels(["50", "100", "150", "200", "250", "300"], fontsize=8) 
ax1.set_yticks([-2000,-1000,0,1000,2000])  # 先设定刻度位置
ax1.set_yticklabels([r'$-2\times 10^{3}$',r'$-1\times 10^{3}$',r'$0$',r'$1\times 10^{3}$',r'$2\times 10^{3}$'], fontsize=8)

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
ax1.plot(T,c5mub0*T,'g',linewidth=1.5,alpha=0.7,zorder=3,label=r'$\mu_B=0$')
ax1.plot(T,c5mub100*T,'r',linewidth=1.5,dashes=[5,2],alpha=0.7,zorder=4,label=r'$\mu_B=100\,\mathrm{MeV}$')
ax1.plot(T,c5mub200*T,'b',linewidth=1.5,dashes=[4,1,2,1],alpha=0.7,zorder=5,label=r'$\mu_B=200\,\mathrm{MeV}$')
ax1.plot(T,c5mub300*T,'m',linewidth=1.5,dashes=[2,2],alpha=0.7,zorder=6,label=r'$\mu_B=300\,\mathrm{MeV}$')
ax1.plot(T,c5mub400*T,'y',linewidth=1.5,dashes=[3,1,2,1],alpha=0.7,zorder=7,label=r'$\mu_B=400\,\mathrm{MeV}$')
ax1.plot(T,c5mub500*T,'c',linewidth=1.5,dashes=[2,1],alpha=0.7,zorder=8,label=r'$\mu_B=500\,\mathrm{MeV}$')
#ax1.plot(T,c5mub550*T,'k',linewidth=1.5,dashes=[1,1],alpha=0.5,zorder=9,label=r'$\mu_B=550\,\mathrm{MeV}$')
ax1.plot(blackline,zeroline,dashes=[4,2],color='k',linewidth='1',alpha=0.7)
ax1.set_ylabel(r'$\chi_5$', fontsize=17, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=15, color='black')
ax1.legend(loc=4,fontsize='8',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
ax1.axis([50,300,-7*10**4,6*10**4])
#ax1.set_yscale('symlog',linthresh=10**-12)
ax1.set_xticks([50,100,150,200,250,300])  # 先设定刻度位置
ax1.set_xticklabels(["50", "100", "150", "200", "250", "300"], fontsize=8) 
ax1.set_yticks([-50000,-25000,0,25000,50000])  # 先设定刻度位置
ax1.set_yticklabels([r'$-5\times 10^{4}$',r'$-2.5\times 10^{4}$',r'$0$',r'$2.5\times 10^{4}$',r'$5\times 10^{4}$'], fontsize=8)

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
ax1.plot(T,c6mub0*T**2,'g',linewidth=1.5,alpha=0.7,zorder=3,label=r'$\mu_B=0$')
ax1.plot(T,c6mub100*T**2,'r',linewidth=1.5,dashes=[5,2],alpha=0.7,zorder=4,label=r'$\mu_B=100\,\mathrm{MeV}$')
ax1.plot(T,c6mub200*T**2,'b',linewidth=1.5,dashes=[4,1,2,1],alpha=0.7,zorder=5,label=r'$\mu_B=200\,\mathrm{MeV}$')
ax1.plot(T,c6mub300*T**2,'m',linewidth=1.5,dashes=[2,2],alpha=0.7,zorder=6,label=r'$\mu_B=300\,\mathrm{MeV}$')
ax1.plot(T,c6mub400*T**2,'y',linewidth=1.5,dashes=[3,1,2,1],alpha=0.7,zorder=7,label=r'$\mu_B=400\,\mathrm{MeV}$')
ax1.plot(T,c6mub500*T**2,'c',linewidth=1.5,dashes=[2,1],alpha=0.7,zorder=8,label=r'$\mu_B=500\,\mathrm{MeV}$')
#ax1.plot(T,c6mub550*T**2,'k',linewidth=1.5,dashes=[1,1],alpha=0.5,zorder=9,label=r'$\mu_B=550\,\mathrm{MeV}$')
ax1.plot(blackline,zeroline,dashes=[4,2],color='k',linewidth='1',alpha=0.7)
ax1.set_ylabel(r'$\chi_6$', fontsize=17, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=15, color='black')
ax1.legend(loc=0,fontsize='8',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
ax1.axis([50,300,-1.2*10**6,2*10**6])
#ax1.set_yscale('symlog',linthresh=10**-12)
ax1.set_xticks([50,100,150,200,250,300])  # 先设定刻度位置
ax1.set_xticklabels(["50", "100", "150", "200", "250", "300"], fontsize=8) 
ax1.set_yticks([-1000000,0,1000000,2000000])  # 先设定刻度位置
ax1.set_yticklabels([r'$-1\times 10^{6}$',r'$0$',r'$1\times 10^{6}$',r'$2\times 10^{6}$'], fontsize=8)

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.21, right=0.95, hspace=0.33,wspace=0.2)

fig.savefig("./chi6.pdf")
