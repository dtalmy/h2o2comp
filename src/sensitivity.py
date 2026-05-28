"""
Sensitivity analysis for figure S4: figure 3b's heatmap rerun under +/- 50%
perturbations of kdam, phi_s, phi_b (via EZ55), and dp.

Run from src/ directory:
    python sensitivity.py
"""

from functions import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

##################################
# fonts
##################################

plt.rcParams.update({
    'font.size': 13,
    'axes.labelsize': 14,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'legend.fontsize': 12,
    'figure.titlesize': 14,
    'axes.titlesize': 13,
})

##################################
# grid
##################################

step = 0.001
ndays = 200
mtimes = np.linspace(0,ndays,int(ndays/step))
inits = (1e5,1e5,1e-4,1)

SNs = np.linspace(1.0e-4,5.0e-4,num=10)
Shs = np.linspace(0,500,num=10)

SHlab = r'H$_2$O$_2$ supply rate (pmol mL$^{-1}$ day$^{-1}$)'
SNlab = r'Nitrogen supply rate (mM day$^{-1}$)'
cbarlab = 'Fraction of cyanobacteria cells \n that are $Prochlorococcus$'

# EZ55 background for the phi_b row (other rows use EZ55 = 0)
EZ55_base = 1e5

##################################
# sweeps
##################################

multipliers = [0.5, 1.5]

# (latex label, code key, EZ55 baseline for that row)
sweeps = [
    (r'$k_{dam}$',  'kdam',  0),
    (r'$\phi_s$',   'phi',   0),
    (r'$\phi_b$',   'phi_b', EZ55_base),
    (r'$\delta_p$', 'dp',    0),
]

def run_grid(kdam_v,phi_v,dp_v,ez55_v):
    Z = np.zeros((SNs.shape[0],Shs.shape[0]),float)
    for i,SN in enumerate(SNs):
        for j,Sh in enumerate(Shs):
            params = [ksp,kss,mumaxp,mumaxs,dp_v,ds,kdam_v,deltah,phi_v,rho,
                      SN,Sh,Qnp,Qns,ez55_v,False]
            sol = odeint(leak,inits,mtimes,args=(params,))
            Psc,Ssc = sol[:,0],sol[:,1]
            Z[i,j] = np.mean(Ssc[-50:])/(np.mean(Ssc[-50:])+np.mean(Psc[-50:]))
    return Z

results = {}
for (latex,code_name,ez55_row) in sweeps:
    for mult in multipliers:
        kdam_v,phi_v,dp_v,ez55_v = kdam,phi,dp,ez55_row
        if code_name == 'kdam':
            kdam_v = kdam * mult
        elif code_name == 'phi':
            phi_v = phi * mult
        elif code_name == 'phi_b':
            # scaling EZ55 is equivalent to scaling the 13e-6 detox coefficient
            ez55_v = EZ55_base * mult
        elif code_name == 'dp':
            dp_v = dp * mult
        results[(code_name,mult)] = run_grid(kdam_v,phi_v,dp_v,ez55_v)

all_vals = np.concatenate([Z.ravel() for Z in results.values()])
vmax = float(np.max(all_vals))
print('shared colorbar vmax:',vmax)

##################################
# plot 4x2 figure
##################################

fig,axs = plt.subplots(4,2,figsize=(11,16),dpi=200)
last_grid = None
for row,(latex,code_name,ez55_row) in enumerate(sweeps):
    for col,mult in enumerate(multipliers):
        ax = axs[row,col]
        Z = results[(code_name,mult)]
        last_grid = ax.pcolormesh(Shs,SNs,Z,vmin=0,vmax=vmax,cmap='summer',shading='auto')
        sign = r'$-50\%$' if mult == 0.5 else r'$+50\%$'
        ax.set_title(latex+' '+sign)
        if row == 3:
            ax.set_xlabel(SHlab)
        if col == 0:
            ax.set_ylabel(SNlab)

fig.subplots_adjust(right=0.86,left=0.10,top=0.96,bottom=0.05,
                    wspace=0.25,hspace=0.45)
cbar_ax = fig.add_axes([0.89,0.10,0.025,0.82])
cbar = fig.colorbar(last_grid,cax=cbar_ax)
cbar.set_label(cbarlab)

fig.savefig('../figures/figureS4.png',dpi=200)
print('saved figureS4.png')
