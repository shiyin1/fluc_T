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
c2mub0=np.loadtxt('./LEFT/mub0/data/c2.dat')
c3mub0=np.loadtxt('./LEFT/mub0/data/c3.dat')
c4mub0=np.loadtxt('./LEFT/mub0/data/c4.dat')
c5mub0=np.loadtxt('./LEFT/mub0/data/c5.dat')
c6mub0=np.loadtxt('./LEFT/mub0/data/c6.dat')
c2QCDmub0=np.loadtxt('./QCD/MUB0/data/c2.dat')
c3QCDmub0=np.loadtxt('./QCD/MUB0/data/c3.dat')
c4QCDmub0=np.loadtxt('./QCD/MUB0/data/c4.dat')
c5QCDmub0=np.loadtxt('./QCD/MUB0/data/c5.dat')
c6QCDmub0=np.loadtxt('./QCD/MUB0/data/c6.dat')
# mub=100 data
c2mub100=np.loadtxt('./LEFT/mub100/data/c2.dat')
c3mub100=np.loadtxt('./LEFT/mub100/data/c3.dat')
c4mub100=np.loadtxt('./LEFT/mub100/data/c4.dat')
c5mub100=np.loadtxt('./LEFT/mub100/data/c5.dat')
c6mub100=np.loadtxt('./LEFT/mub100/data/c6.dat')
c2QCDmub100=np.loadtxt('./QCD/MUB100/data/c2.dat')
c3QCDmub100=np.loadtxt('./QCD/MUB100/data/c3.dat')
c4QCDmub100=np.loadtxt('./QCD/MUB100/data/c4.dat')
c5QCDmub100=np.loadtxt('./QCD/MUB100/data/c5.dat')
c6QCDmub100=np.loadtxt('./QCD/MUB100/data/c6.dat')
# mub=200 data
c2mub200=np.loadtxt('./LEFT/mub200/data/c2.dat')
c3mub200=np.loadtxt('./LEFT/mub200/data/c3.dat')
c4mub200=np.loadtxt('./LEFT/mub200/data/c4.dat')
c5mub200=np.loadtxt('./LEFT/mub200/data/c5.dat')
c6mub200=np.loadtxt('./LEFT/mub200/data/c6.dat')
c2QCDmub200=np.loadtxt('./QCD/MUB200/data/c2.dat')
c3QCDmub200=np.loadtxt('./QCD/MUB200/data/c3.dat')
c4QCDmub200=np.loadtxt('./QCD/MUB200/data/c4.dat')
c5QCDmub200=np.loadtxt('./QCD/MUB200/data/c5.dat')
c6QCDmub200=np.loadtxt('./QCD/MUB200/data/c6.dat')
# mub=300 data
c2mub300=np.loadtxt('./LEFT/mub300/data/c2.dat')
c3mub300=np.loadtxt('./LEFT/mub300/data/c3.dat')
c4mub300=np.loadtxt('./LEFT/mub300/data/c4.dat')
c5mub300=np.loadtxt('./LEFT/mub300/data/c5.dat')
c6mub300=np.loadtxt('./LEFT/mub300/data/c6.dat')
c2QCDmub300=np.loadtxt('./QCD/MUB300/data/c2.dat')
c3QCDmub300=np.loadtxt('./QCD/MUB300/data/c3.dat')
c4QCDmub300=np.loadtxt('./QCD/MUB300/data/c4.dat')
c5QCDmub300=np.loadtxt('./QCD/MUB300/data/c5.dat')
c6QCDmub300=np.loadtxt('./QCD/MUB300/data/c6.dat')
# mub=400 data
c2mub400=np.loadtxt('./LEFT/mub400/data/c2.dat')
c3mub400=np.loadtxt('./LEFT/mub400/data/c3.dat')
c4mub400=np.loadtxt('./LEFT/mub400/data/c4.dat')
c2QCDmub400=np.loadtxt('./QCD/MUB400/data/c2.dat')
c3QCDmub400=np.loadtxt('./QCD/MUB400/data/c3.dat')
c4QCDmub400=np.loadtxt('./QCD/MUB400/data/c4.dat')
# mub=500 data
c2mub500=np.loadtxt('./LEFT/mub500/data/c2.dat')
c3mub500=np.loadtxt('./LEFT/mub500/data/c3.dat')
c4mub500=np.loadtxt('./LEFT/mub500/data/c4.dat')
c2QCDmub500=np.loadtxt('./QCD/MUB500/data/c2.dat')
c3QCDmub500=np.loadtxt('./QCD/MUB500/data/c3.dat')
c4QCDmub500=np.loadtxt('./QCD/MUB500/data/c4.dat')
# Temperature
T=np.arange(21, 301, 1)
TQCD=np.arange(10, 260, 2)
# pressure
pmub0=np.loadtxt('./LEFT/mub0/data/chi0.dat')
pmub100=np.loadtxt('./LEFT/mub100/data/chi0.dat')
pmub200=np.loadtxt('./LEFT/mub200/data/chi0.dat')
pmub300=np.loadtxt('./LEFT/mub300/data/chi0.dat')
pQCDmub0=np.loadtxt('./QCD/MUB0/data/chi0.dat')
pQCDmub100=np.loadtxt('./QCD/MUB100/data/chi0.dat')
pQCDmub200=np.loadtxt('./QCD/MUB200/data/chi0.dat')
pQCDmub300=np.loadtxt('./QCD/MUB300/data/chi0.dat')

poT4mub0=pmub0/T**4
poT4mub100=pmub100/T**4
poT4mub200=pmub200/T**4
poT4mub300=pmub300/T**4
poT4QCDmub0=pQCDmub0/TQCD**4
poT4QCDmub100=pQCDmub100/TQCD**4
poT4QCDmub200=pQCDmub200/TQCD**4
poT4QCDmub300=pQCDmub300/TQCD**4
################################################################################################################
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(T,c2mub0,linewidth=2,alpha=1,c='#7fb2d3') 
ax1.plot(T,c2mub100,linewidth=1.5,alpha=1,c='#8dd3c9') 
ax1.plot(T,c2mub200,linewidth=2,alpha=1,c='#ffb55f') 
ax1.plot(T,c2mub300,linewidth=2,alpha=1,c='#fc7f71') 
#ax1.plot(T,c2mub400,linewidth=2,alpha=0.8,label=r'$\mu_B=400\,\mathrm{MeV}$') 
#ax1.plot(T,c2mub500,linewidth=2,alpha=0.8,label=r'$\mu_B=500\,\mathrm{MeV}$') 

ax1.plot(TQCD,c2QCDmub0,linewidth=2,alpha=1,c='#7fb2d3',dashes=[2,1]) 
ax1.plot(TQCD,c2QCDmub100,linewidth=2,alpha=1,c='#8dd3c9',dashes=[2,1]) 
ax1.plot(TQCD,c2QCDmub200,linewidth=2,alpha=1,c='#ffb55f',dashes=[2,1]) 
ax1.plot(TQCD,c2QCDmub300,linewidth=2,alpha=1,c='#fc7f71',dashes=[2,1]) 

ax1.plot(TQCD,c2QCDmub400-10,linewidth=6,alpha=0.8,c='#7fb2d3',label=r'$\mu_B=0$') 
ax1.plot(TQCD,c2QCDmub500-10,linewidth=6,alpha=0.8,c='#8dd3c9',label=r'$\mu_B=100\,\mathrm{MeV}$') 
ax1.plot(TQCD,c2QCDmub500-10,linewidth=6,alpha=0.8,c='#ffb55f',label=r'$\mu_B=200\,\mathrm{MeV}$') 
ax1.plot(TQCD,c2QCDmub500-10,linewidth=6,alpha=0.8,c='#fc7f71',label=r'$\mu_B=300\,\mathrm{MeV}$') 

ax1.plot(TQCD,c2QCDmub500-10,linewidth=2,alpha=1,c='gray',label=r'$\mathrm{QCD-assisted\,\,LEFT}$') 
ax1.plot(TQCD,c2QCDmub500-10,linewidth=2,alpha=1,c='gray',dashes=[2,1],label=r'$\mathrm{QCD}$') 

ax1.set_ylabel(r'$c_2(T,\mu_B)$', fontsize=13, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=13, color='black')
ax1.legend(loc=0,fontsize='8',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
ax1.axis([50,250,10**-2,2])
ax1.set_yscale('symlog',linthresh=10**-7)
for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.17, right=0.95, hspace=0.33,wspace=0.2)

fig.savefig("./c2.pdf")
################################################################################################################
# Create figure
fig=plt.figure(figsize=(4.5*2, 3.5*2))
#fig=plt.figure()
ax1=fig.add_subplot(221)
ax1.plot(T,c3mub0,linewidth=2,alpha=1,c='#7fb2d3') 
ax1.plot(T,c3mub100,linewidth=1.5,alpha=1,c='#8dd3c9') 
ax1.plot(T,c3mub200,linewidth=2,alpha=1,c='#ffb55f') 
ax1.plot(T,c3mub300,linewidth=2,alpha=1,c='#fc7f71') 

ax1.plot(TQCD,c3QCDmub0,linewidth=2,alpha=1,c='#7fb2d3',dashes=[2,1]) 
ax1.plot(TQCD,c3QCDmub100,linewidth=2,alpha=1,c='#8dd3c9',dashes=[2,1]) 
ax1.plot(TQCD,c3QCDmub200,linewidth=2,alpha=1,c='#ffb55f',dashes=[2,1]) 
ax1.plot(TQCD,c3QCDmub300,linewidth=2,alpha=1,c='#fc7f71',dashes=[2,1]) 

ax1.plot(TQCD,c2QCDmub400-100,linewidth=6,alpha=0.8,c='#7fb2d3',label=r'$\mu_B=0$') 
ax1.plot(TQCD,c2QCDmub500-100,linewidth=6,alpha=0.8,c='#8dd3c9',label=r'$\mu_B=100\,\mathrm{MeV}$') 
ax1.plot(TQCD,c2QCDmub500-100,linewidth=6,alpha=0.8,c='#ffb55f',label=r'$\mu_B=200\,\mathrm{MeV}$') 
ax1.plot(TQCD,c2QCDmub500-100,linewidth=6,alpha=0.8,c='#fc7f71',label=r'$\mu_B=300\,\mathrm{MeV}$') 

ax1.plot(TQCD,c2QCDmub500-100,linewidth=2,alpha=1,c='gray',label=r'$\mathrm{QCD-assisted\,\,LEFT}$') 
ax1.plot(TQCD,c2QCDmub500-100,linewidth=2,alpha=1,c='gray',dashes=[2,1],label=r'$\mathrm{QCD}$') 

ax1.set_ylabel(r'$c_3(T,\mu_B)$', fontsize=13, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=13, color='black')
ax1.legend(loc=0,fontsize='10',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
ax1.axis([50,250,-10.,-5*10**-4])
ax1.set_xticklabels([])
ax1.set_yscale('symlog',linthresh=10**-7)
for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

ax1=fig.add_subplot(222)
ax1.plot(T,c4mub0,linewidth=2,alpha=1,c='#7fb2d3') 
ax1.plot(T,c4mub100,linewidth=1.5,alpha=1,c='#8dd3c9') 
ax1.plot(T,c4mub200,linewidth=2,alpha=1,c='#ffb55f') 
ax1.plot(T,c4mub300,linewidth=2,alpha=1,c='#fc7f71') 

ax1.plot(TQCD,c4QCDmub0,linewidth=2,alpha=1,c='#7fb2d3',dashes=[2,1]) 
ax1.plot(TQCD,c4QCDmub100,linewidth=2,alpha=1,c='#8dd3c9',dashes=[2,1]) 
ax1.plot(TQCD,c4QCDmub200,linewidth=2,alpha=1,c='#ffb55f',dashes=[2,1]) 
ax1.plot(TQCD,c4QCDmub300,linewidth=2,alpha=1,c='#fc7f71',dashes=[2,1]) 

ax1.plot(TQCD,c4QCDmub400-100,linewidth=6,alpha=0.8,c='#7fb2d3',label=r'$\mu_B=0$') 
ax1.plot(TQCD,c4QCDmub500-100,linewidth=6,alpha=0.8,c='#8dd3c9',label=r'$\mu_B=100\,\mathrm{MeV}$') 
ax1.plot(TQCD,c4QCDmub500-100,linewidth=6,alpha=0.8,c='#ffb55f',label=r'$\mu_B=200\,\mathrm{MeV}$') 
ax1.plot(TQCD,c4QCDmub500-100,linewidth=6,alpha=0.8,c='#fc7f71',label=r'$\mu_B=300\,\mathrm{MeV}$') 

ax1.plot(TQCD,c4QCDmub500-100,linewidth=2,alpha=1,c='gray',label=r'$\mathrm{QCD-assisted\,\,LEFT}$') 
ax1.plot(TQCD,c4QCDmub500-100,linewidth=2,alpha=1,c='gray',dashes=[2,1],label=r'$\mathrm{QCD}$') 

ax1.yaxis.set_label_position("right")
ax1.yaxis.tick_right()
ax1.set_ylabel(r'$c_4(T,\mu_B)$', fontsize=13, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=13, color='black')
ax1.set_xticklabels([])
ax1.axis([50,250,10**-5,1000])
ax1.set_yscale('symlog',linthresh=10**-10)

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

ax1=fig.add_subplot(223)
ax1.plot(T,c5mub0,linewidth=2,alpha=1,c='#7fb2d3') 
ax1.plot(T,c5mub100,linewidth=1.5,alpha=1,c='#8dd3c9') 
ax1.plot(T,c5mub200,linewidth=2,alpha=1,c='#ffb55f') 
ax1.plot(T,c5mub300,linewidth=2,alpha=1,c='#fc7f71') 

ax1.plot(TQCD,c5QCDmub0,linewidth=2,alpha=1,c='#7fb2d3',dashes=[2,1]) 
ax1.plot(TQCD,c5QCDmub100,linewidth=2,alpha=1,c='#8dd3c9',dashes=[2,1]) 
ax1.plot(TQCD,c5QCDmub200,linewidth=2,alpha=1,c='#ffb55f',dashes=[2,1]) 
ax1.plot(TQCD,c5QCDmub300,linewidth=2,alpha=1,c='#fc7f71',dashes=[2,1]) 

ax1.plot(TQCD,c4QCDmub400-10**10,linewidth=6,alpha=0.8,c='#7fb2d3',label=r'$\mu_B=0$') 
ax1.plot(TQCD,c4QCDmub500-10**10,linewidth=6,alpha=0.8,c='#8dd3c9',label=r'$\mu_B=100\,\mathrm{MeV}$') 
ax1.plot(TQCD,c4QCDmub500-10**10,linewidth=6,alpha=0.8,c='#ffb55f',label=r'$\mu_B=200\,\mathrm{MeV}$') 
ax1.plot(TQCD,c4QCDmub500-10**10,linewidth=6,alpha=0.8,c='#fc7f71',label=r'$\mu_B=300\,\mathrm{MeV}$') 

ax1.plot(TQCD,c4QCDmub500-10**10,linewidth=2,alpha=1,c='gray',label=r'$\mathrm{QCD-assisted\,\,LEFT}$') 
ax1.plot(TQCD,c4QCDmub500-10**10,linewidth=2,alpha=1,c='gray',dashes=[2,1],label=r'$\mathrm{QCD}$') 

ax1.set_ylabel(r'$c_5(T,\mu_B)$', fontsize=13, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=13, color='black')
#ax1.legend(loc=0,fontsize='10',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
ax1.axis([50,250,-10**4,-2*10**-6])
ax1.set_yscale('symlog',linthresh=10**-10)

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

ax1=fig.add_subplot(224)
ax1.plot(T,c6mub0,linewidth=2,alpha=1,c='#7fb2d3') 
ax1.plot(T,c6mub100,linewidth=1.5,alpha=1,c='#8dd3c9') 
ax1.plot(T,c6mub200,linewidth=2,alpha=1,c='#ffb55f') 
ax1.plot(T,c6mub300,linewidth=2,alpha=1,c='#fc7f71') 

ax1.plot(TQCD,c6QCDmub0,linewidth=2,alpha=1,c='#7fb2d3',dashes=[2,1]) 
ax1.plot(TQCD,c6QCDmub100,linewidth=2,alpha=1,c='#8dd3c9',dashes=[2,1]) 
ax1.plot(TQCD,c6QCDmub200,linewidth=2,alpha=1,c='#ffb55f',dashes=[2,1]) 
ax1.plot(TQCD,c6QCDmub300,linewidth=2,alpha=1,c='#fc7f71',dashes=[2,1]) 

ax1.plot(TQCD,c4QCDmub400-100,linewidth=6,alpha=0.8,c='#7fb2d3',label=r'$\mu_B=0$') 
ax1.plot(TQCD,c4QCDmub500-100,linewidth=6,alpha=0.8,c='#8dd3c9',label=r'$\mu_B=100\,\mathrm{MeV}$') 
ax1.plot(TQCD,c4QCDmub500-100,linewidth=6,alpha=0.8,c='#ffb55f',label=r'$\mu_B=200\,\mathrm{MeV}$') 
ax1.plot(TQCD,c4QCDmub500-100,linewidth=6,alpha=0.8,c='#fc7f71',label=r'$\mu_B=300\,\mathrm{MeV}$') 

ax1.plot(TQCD,c4QCDmub500-100,linewidth=2,alpha=1,c='gray',label=r'$\mathrm{QCD-assisted\,\,LEFT}$') 
ax1.plot(TQCD,c4QCDmub500-100,linewidth=2,alpha=1,c='gray',dashes=[2,1],label=r'$\mathrm{QCD}$') 

ax1.yaxis.set_label_position("right")
ax1.yaxis.tick_right()
ax1.set_ylabel(r'$c_6(T,\mu_B)$', fontsize=13, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=13, color='black')
ax1.axis([50,250,10**-7,9*10**5])
ax1.set_yscale('symlog',linthresh=10**-10)

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.95, bottom=0.1, left=0.1, right=0.92, hspace=0.,wspace=0.)

fig.savefig("./c3toc6.pdf")
################################################################################################################
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(T,poT4mub0,linewidth=2,alpha=1,c='#7fb2d3') 
ax1.plot(T,poT4mub100,linewidth=1.5,alpha=1,c='#8dd3c9') 
ax1.plot(T,poT4mub200,linewidth=2,alpha=1,c='#ffb55f') 
ax1.plot(T,poT4mub300,linewidth=2,alpha=1,c='#fc7f71') 

ax1.plot(TQCD,poT4QCDmub0,linewidth=2,alpha=1,c='#7fb2d3',dashes=[2,1]) 
ax1.plot(TQCD,poT4QCDmub100,linewidth=2,alpha=1,c='#8dd3c9',dashes=[2,1]) 
ax1.plot(TQCD,poT4QCDmub200,linewidth=2,alpha=1,c='#ffb55f',dashes=[2,1]) 
ax1.plot(TQCD,poT4QCDmub300,linewidth=2,alpha=1,c='#fc7f71',dashes=[2,1]) 

ax1.plot(TQCD,poT4QCDmub0-1000,linewidth=6,alpha=0.8,c='#7fb2d3',label=r'$\mu_B=0$') 
ax1.plot(TQCD,poT4QCDmub0-1000,linewidth=6,alpha=0.8,c='#8dd3c9',label=r'$\mu_B=100\,\mathrm{MeV}$') 
ax1.plot(TQCD,poT4QCDmub0-1000,linewidth=6,alpha=0.8,c='#ffb55f',label=r'$\mu_B=200\,\mathrm{MeV}$') 
ax1.plot(TQCD,poT4QCDmub0-1000,linewidth=6,alpha=0.8,c='#fc7f71',label=r'$\mu_B=300\,\mathrm{MeV}$') 

ax1.plot(TQCD,poT4QCDmub0-1000,linewidth=2,alpha=1,c='gray',label=r'$\mathrm{QCD-assisted\,\,LEFT}$') 
ax1.plot(TQCD,poT4QCDmub0-1000,linewidth=2,alpha=1,c='gray',dashes=[2,1],label=r'$\mathrm{QCD}$') 

ax1.set_ylabel(r'$p(T,\mu_B)/T^4$', fontsize=13, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=13, color='black')
ax1.legend(loc=0,fontsize='8',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
ax1.axis([0,250,0,3])
#ax1.set_yscale('symlog',linthresh=10**-12)
for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.17, right=0.95, hspace=0.33,wspace=0.2)

fig.savefig("./poT4.pdf")