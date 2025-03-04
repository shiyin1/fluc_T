#!/usr/bin/env python
# -*- coding: utf-8 -*-
# sphinx_gallery_thumbnail_number = 3

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter  # useful for `logit` scale
import matplotlib.ticker as ticker
import matplotlib as mpl

mpl.style.use('classic')


# Data for plotting

T=np.loadtxt('./T.dat')
B2=np.loadtxt('../dataB/chi2.dat')
S2=np.loadtxt('../dataS/chi2.dat')
BS=np.loadtxt('../dataBS/chi2.dat')
CBS=(BS-S2-B2)/2.

BS1910=np.loadtxt('./BS_1910.dat')
BS2107=np.loadtxt('./BS_2107.dat')
# Lattice QCD Data

# Create figure
fig1=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig1.add_subplot(111)
#ax1=plt.subplot(111)

ax1.errorbar(BS1910[:,0],-BS1910[:,1],yerr=BS1910[:,2],color='red',ecolor='red',fmt='s', linewidth=2, markersize=3,fillstyle='none',label='Latt. QCD 1910.14592')
ax1.errorbar(BS2107[:,0],BS2107[:,7],yerr=BS2107[:,8],color='green',ecolor='green',fmt='s', linewidth=2, markersize=3,fillstyle='none',label='Latt. QCD 2107.10011')
ax1.plot(T,-CBS/S2,'-',color='blue',label='FRG,PQM')

ax1.axis([100,200,0,0.35])
#ax1.set_xscale('log')

ax1.set_xlabel(r'$T[\mathrm{MeV}^2]$', fontsize=14, color='black')
ax1.set_ylabel(r'$\chi_{BS}^{11}/\chi_S^2$', fontsize=14, color='black')

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

ax1.legend(loc=4,ncol=1,fontsize=8,frameon=True,shadow=True,handlelength=3.,borderpad=0.3,borderaxespad=0.5,numpoints=1,mode=None,framealpha=1.,labelspacing=0.37)

fig1.subplots_adjust(top=0.9, bottom=0.15, left=0.17, right=0.95, hspace=0.35,
                    wspace=0.35)

fig1.savefig("./BS_S2.pdf")


# Create figure
fig2=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax2=fig2.add_subplot(111)

ax2.errorbar(BS2107[:,0],BS2107[:,1],yerr=BS2107[:,2],color='green',ecolor='green',fmt='s', linewidth=2, markersize=3,fillstyle='none',label='Latt. QCD 2107.10011')
ax2.plot(T,B2,'-',color='blue',label='FRG,PQM')

ax2.axis([100,200,0,0.35])
#ax1.set_xscale('log')

ax2.set_xlabel(r'$T[\mathrm{MeV}^2]$', fontsize=14, color='black')
ax2.set_ylabel(r'$\chi_B^2$', fontsize=14, color='black')

for label in ax2.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax2.yaxis.get_ticklabels():
    label.set_fontsize(10)

ax2.legend(loc=2,ncol=1,fontsize=8,frameon=True,shadow=True,handlelength=3.,borderpad=0.3,borderaxespad=0.5,numpoints=1,mode=None,framealpha=1.,labelspacing=0.37)

fig2.subplots_adjust(top=0.9, bottom=0.15, left=0.17, right=0.95, hspace=0.35,
                    wspace=0.35)

fig2.savefig("./B2.pdf")

# Create figure
fig3=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax3=fig3.add_subplot(111)

ax3.errorbar(BS2107[:,0],BS2107[:,3],yerr=BS2107[:,4],color='green',ecolor='green',fmt='s', linewidth=2, markersize=3,fillstyle='none',label='Latt. QCD 2107.10011')
ax3.plot(T,S2,'-',color='blue',label='FRG,PQM')

ax3.axis([100,200,0,0.7])
#ax1.set_xscale('log')

ax3.set_xlabel(r'$T[\mathrm{MeV}^2]$', fontsize=14, color='black')
ax3.set_ylabel(r'$\chi_S^2$', fontsize=14, color='black')

for label in ax3.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax3.yaxis.get_ticklabels():
    label.set_fontsize(10)

ax3.legend(loc=2,ncol=1,fontsize=8,frameon=True,shadow=True,handlelength=3.,borderpad=0.3,borderaxespad=0.5,numpoints=1,mode=None,framealpha=1.,labelspacing=0.37)

fig3.subplots_adjust(top=0.9, bottom=0.15, left=0.17, right=0.95, hspace=0.35,
                    wspace=0.35)

fig3.savefig("./S2.pdf")

# Create figure
fig4=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax4=fig4.add_subplot(111)

ax4.errorbar(BS2107[:,0],-BS2107[:,5],yerr=BS2107[:,6],color='green',ecolor='green',fmt='s', linewidth=2, markersize=3,fillstyle='none',label='Latt. QCD 2107.10011')
ax4.plot(T,-CBS,'-',color='blue',label='FRG,PQM')

ax4.axis([100,200,0,0.3])
#ax1.set_xscale('log')

ax4.set_xlabel(r'$T[\mathrm{MeV}^2]$', fontsize=14, color='black')
ax4.set_ylabel(r'$\chi_{BS}^{11}$', fontsize=14, color='black')

for label in ax4.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax4.yaxis.get_ticklabels():
    label.set_fontsize(10)

ax4.legend(loc=2,ncol=1,fontsize=8,frameon=True,shadow=True,handlelength=3.,borderpad=0.3,borderaxespad=0.5,numpoints=1,mode=None,framealpha=1.,labelspacing=0.37)

fig4.subplots_adjust(top=0.9, bottom=0.15, left=0.17, right=0.95, hspace=0.35,
                    wspace=0.35)

fig4.savefig("./BS.pdf")
