import numpy as np
import pylab as pl
import matplotlib.ticker as ticker
import os
import pickle


#==============================================================================#
# Script to generate bifurcation diagram plots from Auto07p data, used in the  #
# paper "Notch-Mediated Determination of Hair-Bundle Polarity in Mechanosensory#
# Hair Cells of the Zebrafish's Lateral Line" published in Current Biology.    #
# The code reads the bifurcation diagrams for NICD (I), calculates the         #
# corresponding to Emx2 (E) and the Polarity Effector (P) and plots them.      #
# Reffer to the methods section of the paper for more details and the          #
# corresponding equations.                                                     #
# Author: Adrian Jacobo (2018)                                                 #
#==============================================================================#

#==============================================================================#
pl.rc('axes', linewidth=1)
pl.rc('lines', markeredgewidth=1,antialiased=True)
pl.rcParams['xtick.major.pad']='1'
pl.rcParams['ytick.major.pad']='1'
pl.rcParams['axes.labelpad']='1'
pl.rcParams['xtick.labelsize']='7'
pl.rcParams['ytick.labelsize']='7'
pl.rcParams['font.size']='7'
pl.rcParams['font.sans-serif'] = "Arial"
pl.rcParams['font.family']='sans-serif'
pl.rcParams['pdf.fonttype']=42
pl.rcParams['ps.fonttype'] = 42
#==============================================================================#

def plot_var(branch_list,var,ax,color):
    x_scale=1e5
    for branch in branch_list:
        idx = branch['varnames'].index(var)
        if idx == -1:
            print 'Variable '+ var+ ' not found'
        mask_stable = branch['stability'] < 0
        stable_plot = ax.plot(branch['data'][0,mask_stable]*x_scale,branch['data'][idx,mask_stable],color=color)
        mask_unstable = branch['stability'] > 0
        ax.plot(branch['data'][0,mask_unstable]*x_scale,branch['data'][idx,mask_unstable],linestyle='--'
                ,color=stable_plot[0].get_color())

branch_list= pickle.load( open( "ktrun.p", "rb" ) )

tscale=100.*1e-6
rostrad_color='#29ABE2'
caudad_color='#fa3200'
fig=pl.figure(facecolor='white',figsize=(6.4,1.64),dpi=200)
ax = fig.add_subplot(131)
#fig.subplots_adjust(bottom=0.17,top=0.95,left=0.15,right=0.95)
#plot_var(branch_list,'I1',ax,'#1f77b4')
#plot_var(branch_list,'I2',ax,'#ff7f0e')
for branch in branch_list:
    mask_stable = branch['stability'] < 0
    idx=branch['varnames'].index('I1')
    I1 = branch['data'][idx,mask_stable]
    idx=branch['varnames'].index('I2')
    I2 = branch['data'][idx,mask_stable]
    ax.plot(branch['data'][0,mask_stable]/tscale,I1,color=caudad_color)
    ax.plot(branch['data'][0,mask_stable]/tscale,I2,color=rostrad_color)

    mask_unstable = branch['stability'] > 0
    idx=branch['varnames'].index('I1')
    I1 = branch['data'][idx,mask_unstable]
    idx=branch['varnames'].index('I2')
    I2 = branch['data'][idx,mask_unstable]
    ax.plot(branch['data'][0,mask_unstable]/tscale,I1,color=caudad_color,linestyle='--')
    ax.plot(branch['data'][0,mask_unstable]/tscale,I2,color=rostrad_color,linestyle='--')

ax.set_ylabel("NICD")
ax.set_xlabel(r"k$_t$ (s$^{-1}$)")

ax.set_ylim([-1100*0.05,1100*1.05])
ax.set_xlim([0,2.])
ax.xaxis.set_major_locator(pl.MaxNLocator(5))
ax.yaxis.set_major_locator(pl.MaxNLocator(5))

ax1 = fig.add_subplot(132)
se=20.#100.
E0=200.
for branch in branch_list:
    mask_stable = branch['stability'] < 0
    idx=branch['varnames'].index('I1')
    I1 = branch['data'][idx,mask_stable]
    E1 = E0*(1./(1.+(I1/se)**4.))
    idx=branch['varnames'].index('I2')
    I2 = branch['data'][idx,mask_stable]
    E2 = E0*(1./(1.+(I2/se)**4.))
    ax1.plot(branch['data'][0,mask_stable]/tscale,E1,color=caudad_color)
    ax1.plot(branch['data'][0,mask_stable]/tscale,E2,color=rostrad_color)
ax1.set_ylabel("Emx2")
ax1.set_xlabel(r"k$_t$ (s$^{-1}$)")

ax1.set_ylim([-200*0.05,200*1.05])
ax1.set_xlim([0,2.])
ax1.xaxis.set_major_locator(pl.MaxNLocator(5))
ax1.yaxis.set_major_locator(pl.MaxNLocator(5))


ax3 = fig.add_subplot(133)
#fig.subplots_adjust(bottom=0.17,top=0.95,left=0.15,right=0.95)
kn= 40.
ke= 10.
P0= 100.

for branch in branch_list:
    for Ex in [0.]:
        mask_stable = branch['stability'] < 0
        idx=branch['varnames'].index('I1')
        I1 = branch['data'][idx,mask_stable]
        E1 = E0*(1./(1.+(I1/se)**4.))+Ex
        idx=branch['varnames'].index('I2')
        I2 = branch['data'][idx,mask_stable]
        E2 = E0*(1./(1.+(I2/se)**4.))+Ex
        P1 = P0*kn*E1/(ke*I1+kn*E1+kn*ke)
        P2 = P0*kn*E2/(ke*I2+kn*E2+kn*ke)
        ax3.plot(branch['data'][0,mask_stable]/tscale,P1,color=caudad_color)
        ax3.plot(branch['data'][0,mask_stable]/tscale,P2,color=rostrad_color)
ax3.set_ylabel(r"PE")
ax3.set_xlabel(r"k$_t$ (s$^{-1}$)")

ax3.set_ylim([-100*0.05,100*1.05])
ax3.set_xlim([0,2.])
ax3.xaxis.set_major_locator(pl.MaxNLocator(5))
ax3.yaxis.set_major_locator(pl.MaxNLocator(5))

pl.tight_layout(pad=0.,w_pad=0.,h_pad=0.)
pl.show()
name='kt_Emx2_PE_bifdiag'
fig.savefig(name+'.pdf',dpi=600)
