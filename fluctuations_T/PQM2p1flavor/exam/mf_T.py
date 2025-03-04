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

T=np.loadtxt('BUFFER/T.DAT')
MF=np.loadtxt('BUFFER/MF.DAT')
MBO=np.loadtxt('BUFFER/MBO.DAT')

Tmid=(T[1:250]+T[0:249])/2.
dmF=MF[1:250,1]-MF[0:249,1]
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
#ax1=plt.subplot(111)

ax1.plot(T,MF[:,0],'-',color='green',markeredgecolor='none',linewidth=2.,markersize=4)
ax1.plot(T,MF[:,1],'-',color='red',markeredgecolor='none',linewidth=2.,markersize=4)

ax1.axis([100,200,0,530.])
#ax1.set_xscale('log')

ax1.set_xlabel(r'$T\,\mathrm{MeV}]$', fontsize=14, color='black')
ax1.set_ylabel(r'$m\,[\mathrm{MeV}]$', fontsize=14, color='black')

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

#ax1.legend(loc=4,ncol=1,fontsize=8,frameon=True,shadow=True,handlelength=3.,borderpad=0.3,borderaxespad=0.5,numpoints=1,mode=None,framealpha=1.,labelspacing=0.37)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.17, right=0.95, hspace=0.35,
                    wspace=0.35)


fig.savefig("BUFFER/MF.pdf")

fig2=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax2=fig2.add_subplot(111)
#ax1=plt.subplot(111)

ax2.plot(Tmid,dmF,'-',color='green',markeredgecolor='none',linewidth=2.,markersize=4)

ax2.axis([100,200,-7.,1.])
#ax1.set_xscale('log')

ax2.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=14, color='black')
ax2.set_ylabel(r'$dm/dT\,[\mathrm{MeV}]$', fontsize=14, color='black')

for label in ax2.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax2.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig2.subplots_adjust(top=0.9, bottom=0.15, left=0.17, right=0.95, hspace=0.35,
                    wspace=0.35)


fig2.savefig("BUFFER/dmldT.pdf")

fig3=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax3=fig3.add_subplot(111)

ax3.plot(T,MBO[:,5],'-',color='green',markeredgecolor='none',linewidth=2.,markersize=4)
ax3.plot(T,MBO[:,6],'-',color='red',markeredgecolor='none',linewidth=2.,markersize=4)

ax3.axis([100,200,0,900.])

ax3.set_xlabel(r'$T\,[\mathrm{MeV}]$', fontsize=14, color='black')
ax3.set_ylabel(r'$m\,[\mathrm{MeV}]$', fontsize=14, color='black')

for label in ax3.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax3.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig3.subplots_adjust(top=0.9, bottom=0.15, left=0.17, right=0.95, hspace=0.35,
                    wspace=0.35)


fig3.savefig("BUFFER/mbo.pdf")
