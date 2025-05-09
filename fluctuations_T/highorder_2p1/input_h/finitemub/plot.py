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
################################################################################################################
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(T,c2mub0,linewidth=2,alpha=0.8,label=r'$\mu_B=0$') 
ax1.plot(T,c2mub100,linewidth=2,alpha=0.8,label=r'$\mu_B=100\,\mathrm{MeV}$') 
ax1.plot(T,c2mub200,linewidth=2,alpha=0.8,label=r'$\mu_B=200\,\mathrm{MeV}$') 
ax1.plot(T,c2mub300,linewidth=2,alpha=0.8,label=r'$\mu_B=300\,\mathrm{MeV}$') 
ax1.plot(T,c2mub400,linewidth=2,alpha=0.8,label=r'$\mu_B=400\,\mathrm{MeV}$') 
ax1.plot(T,c2mub500,linewidth=2,alpha=0.8,label=r'$\mu_B=500\,\mathrm{MeV}$') 
ax1.plot(T,c2mub550,linewidth=2,alpha=0.8,label=r'$\mu_B=550\,\mathrm{MeV}$') 
ax1.set_ylabel(r'$c_2(T,\mu_B)$', fontsize=13, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=13, color='black')
ax1.legend(loc=0,fontsize='10',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
#ax1.axis([0,300,-10.,-10**-2])
ax1.axis([50,300,10**-2,2])
ax1.set_yscale('symlog',linthresh=10**-7)
for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.17, right=0.95, hspace=0.33,wspace=0.2)

fig.savefig("./c2.pdf")
################################################################################################################
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(T,c3mub0,linewidth=2,alpha=0.8) 
ax1.plot(T,c3mub100,linewidth=2,alpha=0.8) 
ax1.plot(T,c3mub200,linewidth=2,alpha=0.8) 
ax1.plot(T,c3mub300,linewidth=2,alpha=0.8) 
ax1.plot(T,c3mub400,linewidth=2,alpha=0.8) 
ax1.plot(T,c3mub500,linewidth=2,alpha=0.8) 
ax1.plot(T,c3mub550,linewidth=2,alpha=0.8) 
ax1.set_ylabel(r'$c_3(T,\mu_B)$', fontsize=13, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=13, color='black')
#ax1.legend(loc=0,fontsize='10',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
ax1.axis([50,300,-10.,-10**-4])
ax1.set_yscale('symlog',linthresh=10**-7)
for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.17, right=0.95, hspace=0.33,wspace=0.2)

fig.savefig("./c3.pdf")
################################################################################################################
c4_mub0 = np.ma.masked_where(c4mub0 <= 0, c4mub0)
c4_mub100 = np.ma.masked_where(c4mub100 <= 0, c4mub100)
c4_mub200 = np.ma.masked_where(c4mub200 <= 0, c4mub200)
c4_mub300 = np.ma.masked_where(c4mub300 <= 0, c4mub300)
c4_mub400 = np.ma.masked_where(c4mub400 <= 0, c4mub400)
c4_mub500 = np.ma.masked_where(c4mub500 <= 0, c4mub500)
c4_mub550 = np.ma.masked_where(c4mub550 <= 0, c4mub550)
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(T,c4_mub0,linewidth=2,alpha=0.8) 
ax1.plot(T,c4_mub100,linewidth=2,alpha=0.8) 
ax1.plot(T,c4_mub200,linewidth=2,alpha=0.8) 
ax1.plot(T,c4_mub300,linewidth=2,alpha=0.8) 
ax1.plot(T,c4_mub400,linewidth=2,alpha=0.8) 
ax1.plot(T,c4_mub500,linewidth=2,alpha=0.8) 
ax1.plot(T,c4_mub550,linewidth=2,alpha=0.8) 
ax1.set_ylabel(r'$c_4(T,\mu_B)$', fontsize=13, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=13, color='black')
#ax1.legend(loc=0,fontsize='10',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
ax1.axis([50,300,10**-5,1000])
ax1.set_yscale('symlog',linthresh=10**-10)

ax_inset = fig.add_axes([0.55, 0.5, 0.35, 0.35])  # 相对于 figure 的坐标
ax_inset.plot(T,c4mub0,linewidth=2,alpha=0.8) 
ax_inset.plot(T,c4mub100,linewidth=2,alpha=0.8) 
ax_inset.plot(T,c4mub200,linewidth=2,alpha=0.8) 
ax_inset.plot(T,c4mub300,linewidth=2,alpha=0.8) 
ax_inset.plot(T,c4mub400,linewidth=2,alpha=0.8) 
ax_inset.plot(T,c4mub500,linewidth=2,alpha=0.8) 
ax_inset.plot(T,c4mub550,linewidth=2,alpha=0.8) 
ax_inset.axis([50,300,-5*10**-4,10**-3])

ax_inset.set_xticks([50,100,150,200,250,300])  # 先设定刻度位置
ax_inset.set_xticklabels(["50", "100", "150", "200", "250", "300"], fontsize=8) 
ax_inset.set_yticks([-4*10**-4,-2*10**-4,0,2*10**-4,4*10**-4,6*10**-4,8*10**-4,10**-3])  # 先设定刻度位置
ax_inset.set_yticklabels([r'$-4\times 10^{-4}$', r'$-2\times 10^{-4}$', r'0', r'$2\times 10^{-4}$', r'$4\times 10^{-4}$', r'$6\times 10^{-4}$',r'$8\times 10^{-4}$',r'$1\times 10^{-3}$'], fontsize=8) 

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.17, right=0.95, hspace=0.33,wspace=0.2)

fig.savefig("./c4.pdf")
################################################################################################################
c5_mub0 = np.ma.masked_where(c5mub0 >= 0, c5mub0)
c5_mub100 = np.ma.masked_where(c5mub100 >= 0, c5mub100)
c5_mub200 = np.ma.masked_where(c5mub200 >= 0, c5mub200)
c5_mub300 = np.ma.masked_where(c5mub300 >= 0, c5mub300)
c5_mub400 = np.ma.masked_where(c5mub400 >= 0, c5mub400)
c5_mub500 = np.ma.masked_where(c5mub500 >= 0, c5mub500)
c5_mub550 = np.ma.masked_where(c5mub550 >= 0, c5mub550)
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(T,c5_mub0,linewidth=2,alpha=0.8) 
ax1.plot(T,c5_mub100,linewidth=2,alpha=0.8) 
ax1.plot(T,c5_mub200,linewidth=2,alpha=0.8) 
ax1.plot(T,c5_mub300,linewidth=2,alpha=0.8) 
ax1.plot(T,c5_mub400,linewidth=2,alpha=0.8) 
ax1.plot(T,c5_mub500,linewidth=2,alpha=0.8) 
ax1.plot(T,c5_mub550,linewidth=2,alpha=0.8) 
ax1.set_ylabel(r'$c_5(T,\mu_B)$', fontsize=13, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=13, color='black')
ax1.axis([50,300,-10**4,-10**-6])
ax1.set_yscale('symlog',linthresh=10**-12)

ax_inset = fig.add_axes([0.56, 0.2, 0.35, 0.35])  # 相对于 figure 的坐标
ax_inset.plot(T,c5mub0,linewidth=2,alpha=0.8) 
ax_inset.plot(T,c5mub100,linewidth=2,alpha=0.8) 
ax_inset.plot(T,c5mub200,linewidth=2,alpha=0.8) 
ax_inset.plot(T,c5mub300,linewidth=2,alpha=0.8) 
ax_inset.plot(T,c5mub400,linewidth=2,alpha=0.8) 
ax_inset.plot(T,c5mub500,linewidth=2,alpha=0.8) 
ax_inset.plot(T,c5mub550,linewidth=2,alpha=0.8) 
ax_inset.axis([50,300,-2*10**-3,5*10**-4])

ax_inset.set_xticks([50,100,150,200,250,300])  # 先设定刻度位置
ax_inset.set_xticklabels(["50", "100", "150", "200", "250", "300"], fontsize=8) 
ax_inset.set_yticks([-2*10**-3,-1.5*10**-3,-1*10**-3,-5*10**-4,0,5*10**-4])  # 先设定刻度位置
ax_inset.set_yticklabels([r'$-2\times 10^{-3}$', r'$-1.5\times 10^{-3}$', r'$-1\times 10^{-3}$', r'$-0.5\times 10^{-3}$', r'$0$',r'$0.5\times 10^{-3}$'], fontsize=8) 

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.17, right=0.95, hspace=0.33,wspace=0.2)

fig.savefig("./c5.pdf")
################################################################################################################
c6_mub0 = np.ma.masked_where(c6mub0 <= 0, c6mub0)
c6_mub100 = np.ma.masked_where(c6mub100 <= 0, c6mub100)
c6_mub200 = np.ma.masked_where(c6mub200 <= 0, c6mub200)
c6_mub300 = np.ma.masked_where(c6mub300 <= 0, c6mub300)
c6_mub400 = np.ma.masked_where(c6mub400 <= 0, c6mub400)
c6_mub500 = np.ma.masked_where(c6mub500 <= 0, c6mub500)
c6_mub550 = np.ma.masked_where(c6mub550 <= 0, c6mub550)
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(T,c6_mub0,linewidth=2,alpha=0.8) 
ax1.plot(T,c6_mub100,linewidth=2,alpha=0.8) 
ax1.plot(T,c6_mub200,linewidth=2,alpha=0.8) 
ax1.plot(T,c6_mub300,linewidth=2,alpha=0.8) 
ax1.plot(T,c6_mub400,linewidth=2,alpha=0.8) 
ax1.plot(T,c6_mub500,linewidth=2,alpha=0.8) 
ax1.plot(T,c6_mub550,linewidth=2,alpha=0.8) 
ax1.set_ylabel(r'$c_6(T,\mu_B)$', fontsize=13, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=13, color='black')
ax1.axis([50,300,10**-7,10**6])
ax1.set_yscale('symlog',linthresh=10**-12)

ax_inset = fig.add_axes([0.56, 0.51, 0.35, 0.35])  # 相对于 figure 的坐标
ax_inset.plot(T,c6mub0,linewidth=2,alpha=0.8) 
ax_inset.plot(T,c6mub100,linewidth=2,alpha=0.8) 
ax_inset.plot(T,c6mub200,linewidth=2,alpha=0.8) 
ax_inset.plot(T,c6mub300,linewidth=2,alpha=0.8) 
ax_inset.plot(T,c6mub400,linewidth=2,alpha=0.8) 
ax_inset.plot(T,c6mub500,linewidth=2,alpha=0.8) 
ax_inset.plot(T,c6mub550,linewidth=2,alpha=0.8) 
ax_inset.axis([50,300,-1.5*10**-3,1*10**-3])

ax_inset.set_xticks([50,100,150,200,250,300])  # 先设定刻度位置
ax_inset.set_xticklabels(["50", "100", "150", "200", "250", "300"], fontsize=8) 
ax_inset.set_yticks([-1.5*10**-3,-1.0*10**-3,-0.5*10**-3,0,0.5*10**-3,1.0*10**-3])  # 先设定刻度位置
ax_inset.set_yticklabels([r'$-1.5\times 10^{-3}$', r'$-1.0\times 10^{-3}$', r'$-0.5\times 10^{-3}$', r'$0$',r'$0.5\times 10^{-3}$',r'$1.0\times 10^{-3}$'], fontsize=8) 

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.17, right=0.95, hspace=0.33,wspace=0.2)

fig.savefig("./c6.pdf")
################################################################################################################
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(T,poT4mub0,linewidth=2,alpha=0.8,label=r'$\mu_B=0$') 
ax1.plot(T,poT4mub100,linewidth=2,alpha=0.8,label=r'$\mu_B=100\,\mathrm{MeV}$') 
ax1.plot(T,poT4mub200,linewidth=2,alpha=0.8,label=r'$\mu_B=200\,\mathrm{MeV}$') 
ax1.plot(T,poT4mub300,linewidth=2,alpha=0.8,label=r'$\mu_B=300\,\mathrm{MeV}$') 
ax1.plot(T,poT4mub400,linewidth=2,alpha=0.8,label=r'$\mu_B=400\,\mathrm{MeV}$') 
ax1.plot(T,poT4mub500,linewidth=2,alpha=0.8,label=r'$\mu_B=500\,\mathrm{MeV}$') 
ax1.plot(T,poT4mub550,linewidth=2,alpha=0.8,label=r'$\mu_B=550\,\mathrm{MeV}$') 
ax1.set_ylabel(r'$p(T,\mu_B)/T^4$', fontsize=13, color='black')
ax1.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=13, color='black')
ax1.legend(loc=0,fontsize='10',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
ax1.axis([0,300,0,4.5])
#ax1.set_yscale('symlog',linthresh=10**-12)
for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.17, right=0.95, hspace=0.33,wspace=0.2)

fig.savefig("./poT4.pdf")