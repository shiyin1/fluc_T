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

S="./T.dat"
T=np.loadtxt(S)
print(T)
print(T[0:299]+T[1:300])

T1=T[0:299]
T2=T[1:300]

np.add(T1,T2)
print(np.add(T[0:299],T[1:300]))

for i in range(10):
    S="Top"+str(i+1)
    print(S)

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
#ax1.plot(T,-CBS/S2,'-',color='blue',label='FRG,PQM')

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

