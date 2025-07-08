from matplotlib import gridspec
from functions import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.patches import ConnectionPatch
from scipy import *
from scipy.integrate import odeint
import sys

#######################################################
# set up figures
#######################################################

fig1, axs1 = plt.subplots(1, 2, figsize=(12, 4), dpi=300)
fig1.subplots_adjust(wspace=0.3)


# Example: Blank panels (no data, no labels)
for ax in axs1:
    ax.set_frame_on(True)

# --- NEW FIGURE 2: One large on top, three taller, multi-row small panels ---
fig2 = plt.figure(figsize=(14, 14), dpi=300)
fig3, axc = plt.subplots(2,2,figsize=[12,10])
axc = axc.flatten()

# Custom gridspec: 2 rows, but top row (single panel) has 1/3 height, bottom row (three panels) has 2/3 height
# We'll manually create inner grids for the bottom panels
outer_gs = gridspec.GridSpec(2, 1, height_ratios=[1, 2], hspace=0.2)

# Top: Large square-like panel
top_gs = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer_gs[0], wspace=0.1)
ax_top = fig2.add_subplot(top_gs[1])  # Use center 1/3 only
ax_top.set_box_aspect(1)              # Make it square in terms of figure space, but do NOT distort axes
ax_top.set_frame_on(True)

# Bottom: 3 taller panels each with 3 subrows
bottom_gs = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer_gs[1], wspace=0.3)

# Create 3 columns, each will have its own 3 stacked rows (no space between)
inner_axes = []
for i in range(3):
    inner_gs = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=bottom_gs[i], hspace=0)
    col_axes = []
    for j in range(3):
        ax = fig2.add_subplot(inner_gs[j])
        ax.set_frame_on(True)
        col_axes.append(ax)
    inner_axes.append(col_axes)

#fig2.tight_layout()

################################################################################
# parameter and variable Set UP 
##################################################################################

step = 0.001
ndays = 200
mtimes = np.linspace(0,ndays,int(ndays/step))

P = np.array([])
S  = np.array([])
N = np.array([])
H = np.array([])

y = [P,S,N,H]

#initial values to be used for odeint start 

P0 = 1e6
S0 = 1e6
N0 = 0.1 #nM 
H0 = 1    #nM

inits = (P0,S0,N0,H0)

scale = 1e+9

# for contours
SNs = np.linspace(1.0e5, 2.0e6, num = 10)/scale
Shs = np.linspace(0, 500, num = 10)
Z = np.zeros((int(SNs.shape[0]),int(Shs.shape[0])),float)
Nc,Pc,Sc,Hc = Z.copy(),Z.copy(),Z.copy(),Z.copy()

# axis limits
cmin,cmax = 1e+5,1e+8
nmin,nmax = 4e-4,1e-0
hmin,hmax = 1e-2,3e+2

SNlab = r'Nitrogen supply rate ($\mu$M day$^{-1}$)'
SHlab = 'H$_2$O$_2$ supply rate (nM day$^{-1}$)'

###########################################################################
# contours
###########################################################################

# colorbar labels
nlab = r'Nitrogen concentration ($\mu$M)'
plab = r'$Prochlorococcus$ (cells mL$^{-1}$)'
slab = r'$Synechococcus$ (cells mL$^{-1}$)'
clab = r'Cells (mL$^{-1}$)'
hlab = r'H$_2$O$_2$ concentration (nM)'
cmlabs = [nlab,plab,slab,hlab]

cbarlab = 'Fraction of cyanobacteria cells \n that are $Prochlorococcus$'

#dynamics, leaky syn
for (i,SN) in zip(range(SNs.shape[0]),SNs):
    for (j,Sh) in zip(range(Shs.shape[0]),Shs):
        params = [ksp,kss,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh,Qnp,Qns,0]
        leaky  = odeint(leak, inits, mtimes, args = (params,))
        Psc = leaky[:,0]
        Ssc = leaky[:,1]
        Nsc = leaky[:,2]
        Hsc = leaky[:,3]
        Ssc_av = np.mean(Ssc[-50:])
        Psc_av = np.mean(Psc[-50:])
        ratio = (Ssc_av/( Ssc_av+ Psc_av))
        Z[i,j],Nc[i,j],Pc[i,j],Sc[i,j],Hc[i,j] = ratio,np.mean(Nsc[-50:]),np.mean(Psc[-50:]),np.mean(Ssc[-50:]),np.mean(Hsc[-50:])
for (f,ax) in zip([fig1,fig2],[axs1[1],ax_top]):
    grid = ax.pcolormesh( Shs, SNs, Z, vmin=0, vmax=np.max(Z), cmap = 'summer', shading= 'auto'  )  #'gouraud'
    ax.set(xlabel=SHlab)
    ax.set(ylabel=SNlab)
    f.colorbar(grid,ax=ax,label=cbarlab)

for (ax,C,cm,cmlab) in zip(axc,[Nc,Pc,Sc,Hc],['Purples','Blues','Reds','Greens'],cmlabs):
    grid = ax.pcolormesh( Shs, SNs, C,norm=colors.LogNorm(), cmap = cm, shading= 'auto'  )  #'gouraud'
    ax.set(xlabel=SHlab)
    ax.set(ylabel=SNlab)
    fig3.colorbar(grid,ax=ax,label=cmlab)

for (ax,C,levs) in zip(axc[1:],[Pc,Sc,Hc],[[1e+4,1e+5],[1e+4,1e+5],[0,200]]):
    contour = ax.contour(Shs, SNs, C, levels=levs, colors='black')
    plt.clabel(contour, inline=True, fontsize=12)

phi = 0.0
#dynamics in non leaky  - only phio froom bacteria - none from sYn
for (i,SN) in zip(range(SNs.shape[0]),SNs):
    for (j,Sh) in zip(range(Shs.shape[0]),Shs):
        params = [ksp,kss,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh,Qnp,Qns,0]
        leaky  = odeint(leak, inits, mtimes, args = (params,))
        Psc = leaky[:,0]
        Ssc = leaky[:,1]
        Nsc = leaky[:,2]
        Hsc = leaky[:,3]
        Ssc_av = np.mean(Ssc[-50:])
        Psc_av = np.mean(Psc[-50:])
        ratio = (Ssc_av/( Ssc_av+ Psc_av))
        Z[i,j] = ratio
grid = axs1[0].pcolormesh( Shs, SNs, Z, vmin=0, vmax=np.max(Z), cmap = 'summer', shading= 'auto'  )  #'gouraud'
axs1[0].set(xlabel=SHlab)
axs1[0].set(ylabel=SNlab)
fig1.colorbar(grid, cmap= 'summer',label = cbarlab)
phi = 1.7e-6

ax_top.text(-0.3,1.0,'a',ha='center',va='center',color='k',transform=ax_top.transAxes,fontsize=14)

###########################################################################
# time-dependent dynamics
###########################################################################

#params for P to Win 
Sh = Shs[0]
SN  = SNs[2]
params = [ksp,kss,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh,Qnp,Qns,0]

#get equilibria with function 
Nstar, Pstar, Sstar, Hstar = Pwins(params)

#run model 
competition  = odeint(leak, inits, mtimes, args = (params,))

#grab values to graph 
Ps = competition[:,0]
Ss = competition[:,1]
Ns = competition[:,2]
Hs = competition[:,3]

#graaph dynamics where P dominates 
ax1, ax2,ax3 = inner_axes[0]
plt.subplots_adjust(right=0.95, wspace = 0.45, left = 0.10, hspace = 0.15, bottom = 0.10)
ax1.set(xlabel='Time (days)', ylabel=clab)
ax2.set(xlabel='Time (days)', ylabel=nlab)
ax3.set(xlabel='Time (days)', ylabel=hlab)

ax1.plot(mtimes, Ps , linewidth = 3, color = 'g', label = 'Pro abundance') #'k1 =' + str(k1p))
ax1.plot(mtimes, Ss , linewidth = 3, color = 'orange', label = 'Syn abundance') #'k1 =' + str(k1s))
ax2.plot(mtimes, Ns, linewidth = 3, color = 'purple', label = "Nutrient Concentration ")
ax3.plot(mtimes, Hs,  linewidth = 3, color = 'red', label = "H$_2$O$_2$ concentration ")

##### graphing stars equations ############### 
ax1.axhline(Sstar,color = 'orangered', linestyle = "-.",label = 'S*')
ax1.axhline(Pstar,color = 'darkgreen', linestyle = ":",label = 'P*')
ax2.axhline(Nstar,color = 'purple', linestyle = "-.",label = 'N*')
ax3.axhline(Hstar,color = 'red', linestyle = "-.",label = 'H*')

ax1.set_ylim([cmin,cmax])
ax2.set_ylim([nmin,nmax])
ax3.set_ylim([hmin,hmax])

ax1.semilogy()
ax2.semilogy()
ax3.semilogy()

# connecting arrow
xyA = (Sh, SN)
xyB = (0.75, 1.2)
ax_top.plot(*xyA, 'o', color='blue')
con = ConnectionPatch(xyA, xyB, axesA = ax_top, axesB = ax1, coordsA=ax_top.transData, coordsB='axes fraction',
                      color='black', arrowstyle='->', linewidth=2, connectionstyle='arc3,rad=0.3')
fig2.add_artist(con)

# letter labels
ax1.text(-0.2,1.2,'b',ha='center',va='center',color='k',transform=ax1.transAxes,fontsize=14)

#####################################################
#params for S to Win 
Sh = Shs[-3]
SN  = SNs[4]
params = [ksp,kss,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh,Qnp,Qns,0]

#get equilibria with function 
Nstar, Pstar, Sstar, Hstar = Swins(params)

#run model 
competition  = odeint(leak, inits, mtimes, args = (params,))

#grab values to graph 
Ps = competition[:,0]
Ss = competition[:,1]
Ns = competition[:,2]
Hs = competition[:,3]

ax1, ax2,ax3 = inner_axes[1]
ax1.set(xlabel='Time (days)', ylabel=clab)
ax2.set(xlabel='Time (days)', ylabel=nlab)
ax3.set(xlabel='Time (days)', ylabel=hlab)

ax1.plot(mtimes, Ps , linewidth = 3, color = 'g', label = 'Pro abundance') #'k1 =' + str(k1p))
ax1.plot(mtimes, Ss , linewidth = 3, color = 'orange', label = 'Syn abundance ') #'k1 =' + str(k1s))
ax2.plot(mtimes, Ns, linewidth = 3, color = 'purple', label = "Nutrient Concentration ")
ax3.plot(mtimes, Hs,  linewidth = 3, color = 'red', label = "H$_2$O$_2$ concentration ")

ax1.axhline(Sstar,color = 'orangered', linestyle = "-.",label = 'S*')
ax1.axhline(Pstar,color = 'darkgreen', linestyle = ":",label = 'P*')
ax2.axhline(Nstar,color = 'purple', linestyle = "-.",label = 'N*')
ax3.axhline(Hstar,color = 'red', linestyle = "-.",label = 'H*')

ax1.set_ylim([cmin,cmax])
ax2.set_ylim([nmin,nmax])
ax3.set_ylim([hmin,hmax])

ax1.semilogy()
ax2.semilogy()
ax3.semilogy()

# connecting arrow
xyA = (Sh, SN)
xyB = (0.75, 1.2)
ax_top.plot(*xyA, 'o', color='blue')
con = ConnectionPatch(xyA, xyB, axesA = ax_top, axesB = ax1, coordsA=ax_top.transData, coordsB='axes fraction',
                      color='black', arrowstyle='->', linewidth=2, connectionstyle='arc3,rad=-0.3')
fig2.add_artist(con)

# letter labels
ax1.text(-0.2,1.2,'c',ha='center',va='center',color='k',transform=ax1.transAxes,fontsize=14)

####################################################################

#params for coexists
Sh = Shs[2]
SN = SNs[5]
params = [ksp,kss,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh,Qnp,Qns,0]

#run model 
competition  = odeint(leak, inits, mtimes, args = (params,))

#get equilibria with function 
Nstar, Pstar, Sstar, Hstar = Coexist(params)

#grab values to graph 
Ps = competition[:,0]
Ss = competition[:,1]
Ns = competition[:,2]
Hs = competition[:,3]

#graph coexistence dynamics 
ax1,ax2,ax3 = inner_axes[2]
ax1.set(xlabel='Time (days)', ylabel=clab)
ax2.set(xlabel='Time (days)', ylabel=nlab)
ax3.set(xlabel='Time (days)', ylabel=hlab)

ax1.plot(mtimes, Ps , linewidth = 3, color = 'g', label = 'Pro abundance') #'k1 =' + str(k1p))
ax1.plot(mtimes, Ss , linewidth = 3, color = 'orange', label = 'Syn abundance') #'k1 =' + str(k1s))
ax2.plot(mtimes, Ns, linewidth = 3, color = 'purple', label = "Nutrient Concentration ")
ax3.plot(mtimes, Hs,  linewidth = 3, color = 'red', label = "H$_2$O$_2$ concentration ")

##### graphing stars equations ############### 
ax1.axhline(Sstar,color = 'orangered', linestyle = "-.",label = 'S*')
ax1.axhline(Pstar,color = 'darkgreen', linestyle = ":",label = 'P*')
ax2.axhline(Nstar,color = 'purple', linestyle = "-.",label = 'N*')
ax3.axhline(Hstar,color = 'red', linestyle = "-.",label = 'H*')

ax1.set_ylim([cmin,cmax])
ax2.set_ylim([nmin,nmax])
ax3.set_ylim([hmin,hmax])

ax1.semilogy()
ax2.semilogy()
ax3.semilogy()

# connecting arrow
xyA = (Sh, SN)
xyB = (0.5, 1.2)
ax_top.plot(*xyA, 'o', color='blue')
con = ConnectionPatch(xyA, xyB, axesA = ax_top, axesB = ax1, coordsA=ax_top.transData, coordsB='axes fraction',
                      color='black', arrowstyle='->', linewidth=2, connectionstyle='arc3,rad=-0.3')
fig2.add_artist(con)

# letter labels
ax1.text(-0.2,1.2,'d',ha='center',va='center',color='k',transform=ax1.transAxes,fontsize=14)
for (ax,l) in zip(axs1,'ab'):
    ax.text(0.05,0.9,l,ha='center',va='center',c='w',transform=ax.transAxes,fontsize=14)

for (ax,l,c) in zip(axc,'abcd',['w','w','w','k']):
    ax.text(0.12,0.9,l,ha='center',va='center',color=c,transform=ax.transAxes)

# savefigs
fig1.savefig("../figures/dynamic_vs_equilibrium", bbox_inches='tight')
fig2.savefig("../figures/figure1.png", bbox_inches='tight')
fig3.savefig("../figures/figure2.png", bbox_inches='tight')

###########################

##############################
#Contour and model calculations
##############################

# Create the figure
fig4 = plt.figure(figsize=(12, 12))

# Define the GridSpec for precise layout control
gs = gridspec.GridSpec(
    4, 6,
    height_ratios=[2, 0.1, 1, 1],  # Top row is taller; row 1 controls vertical space
    width_ratios=[1, 1, 1, 1, 1, 1],  # Equal width columns for symmetry
    hspace=0.3,  # Vertical space between top and middle-bottom
    wspace=1.0   # Horizontal space between top row plots
)

# Top row: Two plots taking half the width each
ax1 = fig4.add_subplot(gs[0, 0:3])
ax2 = fig4.add_subplot(gs[0, 3:6])

# Middle and bottom rows: Centered, stacked vertically, about 2/3 of the width
ax3 = fig4.add_subplot(gs[2, 1:5])
ax4 = fig4.add_subplot(gs[3, 1:5])

EZ55s = np.linspace(1e+5, 1e+6, num = 10)
Shs = np.linspace(0, 500, num = 10)
Z = np.zeros((int(EZ55s.shape[0]),int(Shs.shape[0])),float)
HOOH = np.zeros((int(EZ55s.shape[0]),int(Shs.shape[0])),float)
SN = SNs[5]
inits = (1e+4,1e+4,0.1,10)

for (i,EZ55) in zip(range(EZ55s.shape[0]),EZ55s):
    for (j,Sh) in zip(range(Shs.shape[0]),Shs):
        params = [ksp,kss,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh,Qnp,Qns,EZ55]
        leaky  = odeint(leak, inits, mtimes, args = (params,))
        Nstar, Pstar, Sstar, Hstar = Coexist(params)
        Psc = leaky[:,0]
        Ssc = leaky[:,1]
        Nsc = leaky[:,2]
        Hsc = leaky[:,3]
        Ssc_av = np.mean(Ssc[-10:])
        Psc_av = np.mean(Psc[-10:])
        ratio = (Ssc_av/( Ssc_av+ Psc_av))
        Z[i,j] = ratio
        HOOH[i,j]=Hsc[-1]

#######################################
# Graphing Cotour plots
######################################

grid = ax1.pcolormesh( Shs, EZ55s/1e+6, np.where(Z == -1, np.nan, Z), vmin=0, vmax=np.max(Z), cmap = 'summer', shading= 'auto'  )  #'gouraud'
ax1.set(xlabel='Supply H$_2$O$_2$ (nM day$^{-1}$)')
ax1.set(ylabel='EZ55 cell density (x10$^6$ mL$^{-1}$)')
fig4.colorbar(grid, cmap= 'summer',label = cbarlab,ax=ax1)

grid = ax2.pcolormesh( Shs, EZ55s/1e+6, HOOH, vmin=0, cmap = 'Reds', shading= 'auto'  )  #'gouraud'
ax2.set(xlabel='Supply H$_2$O$_2$ (nM day$^{-1}$)')
ax2.set(ylabel='EZ55 cell density (x10$^6$ mL$^{-1}$)')
fig4.colorbar(grid, cmap= 'summer',label = 'H$_2$O$_2$ concentration (nM)',ax=ax2)

# Create second y-axis sharing the same x-axis
ax1t = ax1.twinx()

# Move ax2's y-axis to the left side
ax1t.spines["left"].set_position(("outward", 60))  # 60 points left from ax1
ax1t.spines["left"].set_visible(True)
ax1t.yaxis.set_label_position("left")
ax1t.yaxis.set_ticks_position("left")
ax1t.spines["right"].set_visible(False)

# Plot on secondary y-axis
#y2 = [100, 200, 300, 400]
grid = ax1t.pcolormesh( Shs, EZ55s*13e-6, np.where(Z == -1, np.nan, Z), vmin=0, vmax=np.max(Z), cmap = 'summer', shading= 'auto'  )  #'gouraud'
ax1t.set_ylabel(r'abiotic decay rate (day$^{-1}$)')

####################
# P WINS
####################
# connecting arrows
Sh,EZ55 = Shs[3], EZ55s[8]
xyA = (Sh,EZ55/1e+6)
xyB = (-0.2, 0.5)
ax1.plot(*xyA, 'o', color='blue')
con = ConnectionPatch(xyA, xyB, axesA = ax1, axesB = ax4, coordsA=ax1.transData, coordsB='axes fraction',
                      color='black', arrowstyle='->', linewidth=2, connectionstyle='arc3,rad=0.2')
fig4.add_artist(con)

# dynamics
params = [ksp,kss,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh,Qnp,Qns,EZ55]
competition  = odeint(leak, inits, mtimes, args = (params,))
Nstar, Pstar, Sstar, Hstar = Pwins(params)
Ps = competition[:,0]
Ss = competition[:,1]
ax3.plot(mtimes, Ps , linewidth = 3, color = 'g',label='$Prochlorococcus$')
ax3.plot(mtimes, Ss , linewidth = 3, color = 'orange',label='$Synechococcus$')
ax3.axhline(Sstar,color = 'orangered', linestyle = "-.")
ax3.axhline(Pstar,color = 'darkgreen', linestyle = ":")

# legend
l3 = ax3.legend()
l3.draw_frame(False)

####################
# S WINS
####################
# connecting arrows
Sh,EZ55 = Shs[7], EZ55s[7]
xyA = (Sh,EZ55/1e+6)
xyB = (0.25, 1.2)
ax1.plot(*xyA, 'o', color='blue')
con = ConnectionPatch(xyA, xyB, axesA = ax1, axesB = ax3, coordsA=ax1.transData, coordsB='axes fraction',
                      color='black', arrowstyle='->', linewidth=2, connectionstyle='arc3,rad=-0.3')
fig4.add_artist(con)

# dynamics
params = [ksp,kss,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh,Qnp,Qns,EZ55]
competition  = odeint(leak, inits, mtimes, args = (params,))
Nstar, Pstar, Sstar, Hstar = Coexist(params)
Ps = competition[:,0]
Ss = competition[:,1]
ax4.plot(mtimes, Ps , linewidth = 3, color = 'g')
ax4.plot(mtimes, Ss , linewidth = 3, color = 'orange')
ax4.axhline(Sstar,color = 'orangered', linestyle = "-.",label = 'S*')
ax4.axhline(Pstar,color = 'darkgreen', linestyle = ":",label = 'P*')

ax3.set(xlabel='Time (days)', ylabel=clab)
ax4.set(xlabel='Time (days)', ylabel=clab)

ax3.semilogy()
ax4.semilogy()

ax4.set_ylim([1e-2,1e+8])

# labels
ax1.text(0.1,0.9,'a',ha='center',va='center',color='w',transform=ax1.transAxes,fontsize=14)
ax2.text(0.1,0.9,'b',ha='center',va='center',color='k',transform=ax2.transAxes,fontsize=14)
ax3.text(0.05,0.9,'c',ha='center',va='center',color='k',transform=ax3.transAxes,fontsize=14)
ax4.text(0.05,0.9,'d',ha='center',va='center',color='k',transform=ax4.transAxes,fontsize=14)

fig4.subplots_adjust(wspace=0.35)
fig4.savefig('../figures/davids_leaky_contour',dpi=300,bbox_inches='tight')

##################################################
# demo plot showing Rstar values w.r.t. 
#hydrogen peroxide concentration

f5,ax5 = plt.subplots()

Hs = np.linspace(0,75,100)
Pstars = ksp*(dp+kdam*Hs)/(k2p-(dp+kdam*Hs))
Sstar = kss*ds/(k2s-ds)

ax5.axhline(Sstar,label='$Synechococcus$ R$^{*}$',c='k',lw=2)
ax5.plot(Hs,Pstars,label='$Prochlorococcus$ R$^{*}$',c='k',ls='--',lw=2)

ax5.set_xlabel(hlab)
ax5.set_ylabel(r'Minimal nutrient requirement, R$^{*}$ ($\mu$M)')

l5 = ax5.legend()
l5.draw_frame(False)

f5.savefig('../figures/demplot',dpi=300,bbox_inches='tight')

