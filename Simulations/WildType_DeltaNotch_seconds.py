import PyDSTool as dst
import numpy as np
from matplotlib import pyplot as pl

#==============================================================================#
# Python  script to simulate the time traces plotted in the paper              #
# "Notch-Mediated Determination of Hair-Bundle Polarity in Mechanosensory Hair #
# Cells of the Zebrafish's Lateral Line" published in Current Biology.         #
# Refer to the methods section of the paper for more details on the equations. #
# This code uses PyDStool to integrate the differential equations of the       #
# system.                                                                       #
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

# we must give a name
DSargs = dst.args(name='Two Cells Delta-Notch Model')

sx=2.
sy=1.5
tscale=100.
DSargs.pars = { 'g': 1e-1,
               'gi': 5e-1,
               'kc': 5e-4,
               'kt': 2e-4,
               's0': 2e2,
               'n0': 5e2,
               'd0': 1e3,
               'i0': 0.,
               'ln': 2.,
               'ld': 0.,
               'de': 0.}
for key in DSargs.pars.keys():
    print key,DSargs.pars[key]/tscale
# auxiliary helper function(s) -- function name: ([func signature], definition)
DSargs.fnspecs  = {'hs1n': (['i1'], '(1./(1.+(i1/s0)**2))+ln*((i1/s0)**2)/(1.+(i1/s0)**2)'),
                   'hs1d': (['i1'], '(1./(1.+(i1/s0)**2))+ld*((i1/s0)**2)/(1.+(i1/s0)**2)'),
                   'hs2n': (['i2'], '(1./(1.+(i2/s0)**2))+ln*((i2/s0)**2)/(1.+(i2/s0)**2)'),
                   'hs2d': (['i2'], '(1./(1.+(i2/s0)**2))+ld*((i2/s0)**2)/(1.+(i2/s0)**2)')}

# rhs of the differential equations
DSargs.varspecs = {'n1': 'n0*hs1n(i1)-kc*n1*d1-kt*n1*d2-g*n1',
                   'd1': 'd0*hs1d(i1)-kc*d1*n1-kt*d1*n2-g*d1+de',
                   'i1': 'kt*n1*d2-gi*i1+i0',
                   'n2': 'n0*hs2n(i2)-kc*n2*d2-kt*n2*d1-g*n2',
                   'd2': 'd0*hs2d(i2)-kc*d2*n2-kt*d2*n1-g*d2+de',
                   'i2': 'kt*n2*d1-gi*i2+i0'}
# initial conditions
DSargs.ics      = {'n1': 1., 'd1': 1., 'i1':1.,
                   'n2': 10., 'd2': 10., 'i2':10. }

DSargs.tdomain = [0,100]                         # set the range of integration.
ode  = dst.Generator.Vode_ODEsystem(DSargs)     # an instance of the 'Generator' class.
traj = ode.compute('Lateral Inhibition')        # integrate ODE
pts  = traj.sample(dt=1.)                      # Data for plotting

# PyPlot commands
fig=pl.figure(facecolor='none',figsize=(sx,sy))
#fig.subplots_adjust(bottom=0.15,left=0.08,top=0.95,right=0.98,wspace = 0.35)
ax2 = fig.add_subplot(111)
d1 = pts['d1']
d2 = pts['d2']
ax2.plot(pts['t']*tscale, d1,color='#fa3200')
ax2.plot(pts['t']*tscale, d2,color='#29ABE2')
ax2.set_xlabel('Time (s)')                              # Axes labels
ax2.set_ylabel('Delta')                           # ...
ax2.set_xlim([DSargs.tdomain[0],DSargs.tdomain[1]*tscale])                           # ...
#ax4.set_ylim([-800*0.015,800*1.05])                           # ...
pl.tight_layout()
pl.savefig('WildType_Delta.pdf',transparent=True)



fig2=pl.figure(facecolor='none',figsize=(sx,sy))
#fig.subplots_adjust(bottom=0.15,left=0.08,top=0.95,right=0.98,wspace = 0.35)
ax4 = fig2.add_subplot(111)
I1 = pts['i1']
I2 = pts['i2']
ax4.plot(pts['t']*tscale, I1,color='#fa3200')
ax4.plot(pts['t']*tscale, I2,color='#29ABE2')
ax4.set_xlabel('Time (s)')                              # Axes labelsa
ax4.set_ylabel('NICD')
ax4.set_xlim([DSargs.tdomain[0],DSargs.tdomain[1]*tscale])                           # ...
#ax4.set_ylim([-800*0.015,800*1.05])                           # ...
pl.tight_layout()
pl.savefig('WildType_NICD.pdf',transparent=True)

fig3=pl.figure(facecolor='none',figsize=(sx,sy))

se=20.#100.
E0=200.

E1 = E0*(1./(1.+(I1/se)**4.))
E2 = E0*(1./(1.+(I2/se)**4.))

ax5 = fig3.add_subplot(111)
ax5.plot(pts['t']*tscale, E1,color='#fa3200')
ax5.plot(pts['t']*tscale, E2,color='#29ABE2')
ax5.set_xlabel('Time (s)')                              # Axes labelsa
ax5.set_ylabel('Emx2')
ax5.set_xlim([DSargs.tdomain[0],DSargs.tdomain[1]*tscale])                           # ...
ax5.set_ylim([-200*0.05,200*1.05])                           # ...
pl.tight_layout()
pl.savefig('WildType_Emx2.pdf',dpi=600,transparent=True)

fig4=pl.figure(facecolor='none',figsize=(sx,sy))

kn= 40.
ke= 10.
P0= 100.

P1 = P0*kn*E1/(ke*I1+kn*E1+kn*ke)
P2 = P0*kn*E2/(ke*I2+kn*E2+kn*ke)
ax6 = fig4.add_subplot(111)
ax6.plot(pts['t']*tscale, P1,color='#fa3200')
ax6.plot(pts['t']*tscale, P2,color='#29ABE2')
ax6.plot(DSargs.tdomain,[0.5,0.5],'k--')
ax6.set_xlabel('Time (s)')                              # Axes labelsa
ax6.set_ylabel('P')
ax6.set_xlim([DSargs.tdomain[0],DSargs.tdomain[1]*tscale])                           # ...
ax6.set_ylim([-100*0.05,100*1.05])                           # ...
pl.tight_layout()
pl.savefig('WildType_P.pdf',transparent=True)

print I1[-1],I2[-1]
print E1[-1],E2[-1]
print P1[-1],P2[-1]
