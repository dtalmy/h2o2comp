"""
Reproduces figures 2-6 of the H2O2 competition paper.

Run from src/ directory:
    python simulations.py
"""

from matplotlib import gridspec
from functions import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.patches import ConnectionPatch
from scipy.integrate import odeint

##################################
# fonts and helpers
##################################

plt.rcParams.update({
    'font.size': 13,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.titlesize': 14,
    'axes.titlesize': 13,
})

# format contour labels as 10^n or c x 10^n
def sci_fmt(x):
    if x == 0:
        return '0'
    exp = int(np.floor(np.log10(abs(x))))
    coef = x / 10.0**exp
    if abs(coef - 1.0) < 1e-3:
        return rf'$10^{{{exp}}}$'
    if abs(coef - round(coef)) < 1e-3:
        coef_str = f'{int(round(coef))}'
    else:
        coef_str = f'{coef:.1f}'.rstrip('0').rstrip('.')
    return rf'${coef_str}{{\times}}10^{{{exp}}}$'

# write S*/P*/N*/H* labels just outside the right edge of an axis
def label_star(ax,y,text,color):
    ax.text(1.01,y,text,transform=ax.get_yaxis_transform(),
            va='center',ha='left',color=color,fontsize=12,fontweight='bold',
            clip_on=False)

# background for panel letters
letterbox = dict(boxstyle='round,pad=0.25',facecolor='white',edgecolor='none',alpha=0.85)

##################################
# parameters and labels
##################################

step = 0.001
ndays = 200
mtimes = np.linspace(0,ndays,int(ndays/step))
inits = (1e5,1e5,1e-4,1)  # P0, S0 (cells mL-1); N0 (mM); H0 (pmol mL-1)

SNs = np.linspace(1.0e-4,5.0e-4,num=10)
Shs = np.linspace(0,500,num=10)
Z = np.zeros((SNs.shape[0],Shs.shape[0]),float)
Nc,Pc,Sc,Hc = Z.copy(),Z.copy(),Z.copy(),Z.copy()

cmin,cmax = 1e+5,1e+8
nmin,nmax = 1e-7,1e-3
hmin,hmax = 1e-2,1e+3

SNlab = r'Nitrogen supply rate (mM day$^{-1}$)'
SHlab = r'H$_2$O$_2$ supply rate (pmol mL$^{-1}$ day$^{-1}$)'
nlab = r'Nitrogen concentration (mM)'
plab = r'$Prochlorococcus$ (cells mL$^{-1}$)'
slab = r'$Synechococcus$ (cells mL$^{-1}$)'
clab = r'Cells (mL$^{-1}$)'
hlab = r'H$_2$O$_2$ concentration (pmol mL$^{-1}$)'
cmlabs = [nlab,plab,slab,hlab]
cbarlab = 'Fraction of cyanobacteria cells \n that are $Prochlorococcus$'

##################################
# figure 2: R* vs H2O2 concentration
##################################

f2,ax_r = plt.subplots()
Hs_axis = np.linspace(0,160,200)
Pstars = ksp*(dp+kdam*Hs_axis)/(mumaxp-(dp+kdam*Hs_axis))
Sstar_init = kss*ds/(mumaxs-ds)

ax_r.axhline(Sstar_init,label='$Synechococcus$ R$^{*}$',c='k',lw=2)
ax_r.plot(Hs_axis,Pstars,label='$Prochlorococcus$ R$^{*}$',c='k',ls='--',lw=2)
ax_r.set_xlabel(hlab)
ax_r.set_ylabel(r'Minimal nutrient requirement, R$^{*}$ (mM)')
l2 = ax_r.legend()
l2.draw_frame(False)
f2.savefig('../figures/figure2',dpi=300,bbox_inches='tight')

##################################
# figure containers for 3, 4, 5
##################################

fig3,axs3 = plt.subplots(1,2,figsize=(13,4.8),dpi=300)
fig3.subplots_adjust(wspace=0.55)

fig4 = plt.figure(figsize=(14,14),dpi=200)
outer_gs = gridspec.GridSpec(2,1,height_ratios=[1,2.4],hspace=0.22)
top_gs = gridspec.GridSpecFromSubplotSpec(1,3,subplot_spec=outer_gs[0],wspace=0.1)
ax_top = fig4.add_subplot(top_gs[1])
ax_top.set_box_aspect(1)
bottom_gs = gridspec.GridSpecFromSubplotSpec(1,3,subplot_spec=outer_gs[1],wspace=0.65)
inner_axes = []
for i in range(3):
    inner_gs = gridspec.GridSpecFromSubplotSpec(3,1,subplot_spec=bottom_gs[i],hspace=0.35)
    inner_axes.append([fig4.add_subplot(inner_gs[j]) for j in range(3)])

fig5,axc = plt.subplots(2,2,figsize=(13,11),dpi=300)
axc = axc.flatten()

##################################
# contour simulations
##################################

# figure 3a: no syn detox
phi = 0.0
Za = np.zeros_like(Z)
for (i,SN) in zip(range(SNs.shape[0]),SNs):
    for (j,Sh) in zip(range(Shs.shape[0]),Shs):
        params = [ksp,kss,mumaxp,mumaxs,dp,ds,kdam,deltah,phi,rho,SN,Sh,Qnp,Qns,0,0]
        leaky = odeint(leak,inits,mtimes,args=(params,))
        Psc,Ssc = leaky[:,0],leaky[:,1]
        Za[i,j] = np.mean(Ssc[-50:])/(np.mean(Ssc[-50:])+np.mean(Psc[-50:]))

grid = axs3[0].pcolormesh(Shs,SNs,Za,vmin=0,vmax=np.max(Za),cmap='summer',shading='auto')
axs3[0].set(xlabel=SHlab); axs3[0].set(ylabel=SNlab)
fig3.colorbar(grid,ax=axs3[0],label=cbarlab)

# figure 3b, figure 4 heatmap, figure 5: with syn detox
phi = 1.7e-6
for (i,SN) in zip(range(SNs.shape[0]),SNs):
    for (j,Sh) in zip(range(Shs.shape[0]),Shs):
        params = [ksp,kss,mumaxp,mumaxs,dp,ds,kdam,deltah,phi,rho,SN,Sh,Qnp,Qns,0,False]
        leaky = odeint(leak,inits,mtimes,args=(params,))
        Psc,Ssc,Nsc,Hsc = leaky[:,0],leaky[:,1],leaky[:,2],leaky[:,3]
        Z[i,j] = np.mean(Ssc[-50:])/(np.mean(Ssc[-50:])+np.mean(Psc[-50:]))
        Nc[i,j],Pc[i,j],Sc[i,j],Hc[i,j] = np.mean(Nsc[-50:]),np.mean(Psc[-50:]),np.mean(Ssc[-50:]),np.mean(Hsc[-50:])

##################################
# figure 3 panel b
##################################

grid = axs3[1].pcolormesh(Shs,SNs,Z,vmin=0,vmax=np.max(Z),cmap='summer',shading='auto')
axs3[1].set(xlabel=SHlab); axs3[1].set(ylabel=SNlab)
fig3.colorbar(grid,ax=axs3[1],label=cbarlab)
for (ax,l) in zip(axs3,'ab'):
    ax.text(0.05,0.9,l,ha='center',va='center',c='w',
            transform=ax.transAxes,fontsize=16,fontweight='bold')
fig3.savefig('../figures/figure3.png',bbox_inches='tight',dpi=300)

##################################
# figure 4: heatmap + 3-column dynamics
##################################

grid = ax_top.pcolormesh(Shs,SNs,Z,vmin=0,vmax=np.max(Z),cmap='summer',shading='auto')
ax_top.set(xlabel=SHlab); ax_top.set(ylabel=SNlab)
fig4.colorbar(grid,ax=ax_top,label=cbarlab)
ax_top.text(0.06,0.92,'a',ha='center',va='center',color='k',
            transform=ax_top.transAxes,fontsize=16,fontweight='bold',
            bbox=letterbox)

def plot_dynamics(col_idx,Sh,SN,star_func,letter,arc):
    params = [ksp,kss,mumaxp,mumaxs,dp,ds,kdam,deltah,phi,rho,SN,Sh,Qnp,Qns,0,False]
    Nstar,Pstar,Sstar,Hstar = star_func(params)
    sol = odeint(leak,inits,mtimes,args=(params,))
    Ps,Ss,Ns,Hs_ = sol[:,0],sol[:,1],sol[:,2],sol[:,3]

    a1,a2,a3 = inner_axes[col_idx]
    nlab_wrapped = 'Nitrogen concentration\n(mM)'
    hlab_wrapped = 'H$_2$O$_2$ concentration\n(pmol mL$^{-1}$)'
    a1.set(ylabel=clab); a2.set(ylabel=nlab_wrapped); a3.set(ylabel=hlab_wrapped)
    a3.set_xlabel('Time (days)')
    a1.tick_params(labelbottom=False); a2.tick_params(labelbottom=False)

    a1.plot(mtimes,Ps,linewidth=3,color='g')
    a1.plot(mtimes,Ss,linewidth=3,color='orange')
    a2.plot(mtimes,Ns,linewidth=3,color='purple')
    a3.plot(mtimes,Hs_,linewidth=3,color='red')

    a1.axhline(Sstar,color='orangered',linestyle='--',linewidth=2)
    a1.axhline(Pstar,color='darkgreen',linestyle='--',linewidth=2)
    a2.axhline(Nstar,color='purple',linestyle='--',linewidth=2)
    a3.axhline(Hstar,color='red',linestyle='--',linewidth=2)
    label_star(a1,Sstar,'S*','orangered')
    label_star(a1,Pstar,'P*','darkgreen')
    label_star(a2,Nstar,'N*','purple')
    label_star(a3,Hstar,'H*','red')

    a1.set_ylim([cmin,cmax]); a2.set_ylim([nmin,nmax]); a3.set_ylim([hmin,hmax])
    a1.semilogy(); a2.semilogy(); a3.semilogy()

    ax_top.plot(Sh,SN,'o',color='blue',markersize=8)
    con = ConnectionPatch((Sh,SN),(0.5,1.12),axesA=ax_top,axesB=a1,
                          coordsA=ax_top.transData,coordsB='axes fraction',
                          color='black',arrowstyle='->',linewidth=2,
                          connectionstyle=arc)
    fig4.add_artist(con)
    a1.text(0.06,0.92,letter,ha='center',va='center',color='k',
            transform=a1.transAxes,fontsize=16,fontweight='bold',
            bbox=letterbox)

plot_dynamics(0,Shs[0], SNs[2],Pwins,  'b','arc3,rad=0.15')
plot_dynamics(1,Shs[-3],SNs[4],Swins,  'c','arc3,rad=-0.15')
plot_dynamics(2,Shs[2], SNs[7],Coexist,'d','arc3,rad=-0.15')

# bbox_inches='tight' breaks the ConnectionPatch arrows on this figure
fig4.savefig('../figures/figure4.png',dpi=200)

##################################
# figure 5: 2x2 concentration heatmaps
##################################

for (ax,C,cm,cmlab) in zip(axc,[Nc,Pc,Sc,Hc],['Purples','Blues','Reds','Greens'],cmlabs):
    grid = ax.pcolormesh(Shs,SNs,C,norm=colors.LogNorm(),cmap=cm,shading='auto')
    ax.set(xlabel=SHlab); ax.set(ylabel=SNlab)
    fig5.colorbar(grid,ax=ax,label=cmlab)

import matplotlib.patheffects as pe
Pc_m = np.ma.masked_where(Pc < 1e+2,Pc)
contour_levels = [
    [1e+5, 1e+6, 5e+6],         # Pro
    [5e+5, 1e+6, 2e+6, 3e+6],   # Syn
    [50, 100, 200, 300],        # H2O2
]
label_fmts = [sci_fmt, sci_fmt, '%d']
label_effects = [pe.withStroke(linewidth=3, foreground='black')]
for (ax,C,levs,fmt) in zip(axc[1:], [Pc_m,Sc,Hc], contour_levels, label_fmts):
    contour = ax.contour(Shs,SNs,C,levels=levs,colors='white',linewidths=1.8)
    contour.set(path_effects=[pe.withStroke(linewidth=3.2, foreground='black')])
    labels = ax.clabel(contour, inline=True, fontsize=11, fmt=fmt, colors='white')
    for t in labels:
        t.set_path_effects(label_effects)

for (ax,l) in zip(axc,'abcd'):
    ax.text(0.08,0.92,l,ha='center',va='center',color='k',
            transform=ax.transAxes,fontsize=16,fontweight='bold',
            bbox=letterbox)

fig5.subplots_adjust(wspace=0.45,hspace=0.3)
fig5.savefig('../figures/figure5.png',bbox_inches='tight',dpi=300)

##################################
# figure 6: EZ55 + dynamics
##################################

fig6 = plt.figure(figsize=(13,13),dpi=200)
gs6 = gridspec.GridSpec(
    4,6,
    height_ratios=[1.5,0.15,1.2,1.2],
    width_ratios=[1,1,1,1,1,1],
    hspace=0.35,wspace=1.0
)
ax6_1 = fig6.add_subplot(gs6[0,0:3])
ax6_2 = fig6.add_subplot(gs6[0,3:6])
ax6_3 = fig6.add_subplot(gs6[2,1:5])
ax6_4 = fig6.add_subplot(gs6[3,1:5])

EZ55s = np.linspace(1e+4,5e+5,num=10)
Shs   = np.linspace(0,500,num=10)
Z6    = np.zeros((EZ55s.shape[0],Shs.shape[0]),float)
HOOH  = np.zeros_like(Z6)
SN    = 0.00025

step6 = 0.001; ndays6 = 800
mtimes6 = np.linspace(0,ndays6,int(ndays6/step6))

for (i,EZ55) in zip(range(EZ55s.shape[0]),EZ55s):
    for (j,Sh) in zip(range(Shs.shape[0]),Shs):
        params = [ksp,kss,mumaxp,mumaxs,dp,ds,kdam,deltah,phi,rho,SN,Sh,Qnp,Qns,EZ55,False]
        leaky = odeint(leak,inits,mtimes6,args=(params,))
        Psc,Ssc,Hsc = leaky[:,0],leaky[:,1],leaky[:,3]
        Ssc_av,Psc_av = np.mean(Ssc[-10:]),np.mean(Psc[-10:])
        Z6[i,j]   = Ssc_av/(Ssc_av+Psc_av)
        HOOH[i,j] = Hsc[-1]

grid = ax6_1.pcolormesh(Shs,EZ55s/1e+5,np.where(Z6==-1,np.nan,Z6),
                        vmin=0,vmax=np.max(Z6),cmap='summer',shading='auto')
ax6_1.set(xlabel=r'Supply H$_2$O$_2$ (pmol mL$^{-1}$ day$^{-1}$)')
ax6_1.set(ylabel='EZ55 cell density (x10$^5$ mL$^{-1}$)')
fig6.colorbar(grid,cmap='summer',label=cbarlab,ax=ax6_1)

grid = ax6_2.pcolormesh(Shs,EZ55s/1e+5,HOOH,vmin=0,cmap='Reds',shading='auto')
ax6_2.set(xlabel=r'Supply H$_2$O$_2$ (pmol mL$^{-1}$ day$^{-1}$)')
ax6_2.set(ylabel='EZ55 cell density (x10$^5$ mL$^{-1}$)')
fig6.colorbar(grid,cmap='summer',label=r'H$_2$O$_2$ concentration (pmol mL$^{-1}$)',ax=ax6_2)

# second y-axis: abiotic decay rate (13e-6 * EZ55)
ax6_1t = ax6_1.twinx()
ax6_1t.spines['left'].set_position(('outward',60))
ax6_1t.spines['left'].set_visible(True)
ax6_1t.yaxis.set_label_position('left')
ax6_1t.yaxis.set_ticks_position('left')
ax6_1t.spines['right'].set_visible(False)
grid = ax6_1t.pcolormesh(Shs,EZ55s*13e-6,np.where(Z6==-1,np.nan,Z6),
                         vmin=0,vmax=np.max(Z6),cmap='summer',shading='auto')
ax6_1t.set_ylabel(r'abiotic decay rate (day$^{-1}$)')

# Pro-wins dynamics
Sh,EZ55 = Shs[4],EZ55s[6]
ax6_1.plot(Sh,EZ55/1e+5,'o',color='blue',markersize=8)
con = ConnectionPatch((Sh,EZ55/1e+5),(-0.2,0.5),axesA=ax6_1,axesB=ax6_3,
                      coordsA=ax6_1.transData,coordsB='axes fraction',
                      color='black',arrowstyle='->',linewidth=2,
                      connectionstyle='arc3,rad=0.3')
fig6.add_artist(con)
params = [ksp,kss,mumaxp,mumaxs,dp,ds,kdam,deltah,phi,rho,SN,Sh,Qnp,Qns,EZ55,False]
leaky = odeint(leak,inits,mtimes6,args=(params,))
Nstar,Pstar,Sstar,Hstar = Pwins(params)
Ps,Ss = leaky[:,0],leaky[:,1]
ax6_3.plot(mtimes6,Ps,linewidth=3,color='g',label='$Prochlorococcus$')
ax6_3.plot(mtimes6,Ss,linewidth=3,color='orange',label='$Synechococcus$')
ax6_3.axhline(Sstar,color='orangered',linestyle='--',linewidth=2)
ax6_3.axhline(Pstar,color='darkgreen',linestyle='--',linewidth=2)
label_star(ax6_3,Sstar,'S*','orangered')
label_star(ax6_3,Pstar,'P*','darkgreen')
l3 = ax6_3.legend(); l3.draw_frame(False)

# Syn-wins dynamics
Sh,EZ55 = Shs[7],EZ55s[0]
ax6_1.plot(Sh,EZ55/1e+5,'o',color='blue',markersize=8)
con = ConnectionPatch((Sh,EZ55/1e+5),(-0.2,0.5),axesA=ax6_1,axesB=ax6_4,
                      coordsA=ax6_1.transData,coordsB='axes fraction',
                      color='black',arrowstyle='->',linewidth=2,
                      connectionstyle='arc3,rad=0.4')
fig6.add_artist(con)
params = [ksp,kss,mumaxp,mumaxs,dp,ds,kdam,deltah,phi,rho,SN,Sh,Qnp,Qns,EZ55,False]
Nstar,Pstar,Sstar,Hstar = Swins(params)
leaky = odeint(leak,inits,mtimes6,args=(params,))
Ps,Ss = leaky[:,0],leaky[:,1]
ax6_4.plot(mtimes6,Ps,linewidth=3,color='g')
ax6_4.plot(mtimes6,Ss,linewidth=3,color='orange')
ax6_4.axhline(Sstar,color='orangered',linestyle='--',linewidth=2)
ax6_4.axhline(Pstar,color='darkgreen',linestyle='--',linewidth=2)
label_star(ax6_4,Sstar,'S*','orangered')
label_star(ax6_4,Pstar,'P*','darkgreen')

ax6_3.set(xlabel='Time (days)',ylabel=clab)
ax6_4.set(xlabel='Time (days)',ylabel=clab)
ax6_3.semilogy(); ax6_4.semilogy()
ax6_3.set_ylim([1e+4,1e+8]); ax6_4.set_ylim([1e+4,1e+8])

ax6_1t.text(0.08,0.92,'a',ha='center',va='center',color='k',
            transform=ax6_1t.transAxes,fontsize=16,fontweight='bold',
            bbox=letterbox)
ax6_2.text(0.08,0.92,'b',ha='center',va='center',color='k',
           transform=ax6_2.transAxes,fontsize=16,fontweight='bold',
           bbox=letterbox)
ax6_3.text(0.04,0.92,'c',ha='center',va='center',color='k',
           transform=ax6_3.transAxes,fontsize=16,fontweight='bold',
           bbox=letterbox)
ax6_4.text(0.04,0.92,'d',ha='center',va='center',color='k',
           transform=ax6_4.transAxes,fontsize=16,fontweight='bold',
           bbox=letterbox)

fig6.subplots_adjust(wspace=0.4)
fig6.savefig('../figures/figure6.png',dpi=200)
