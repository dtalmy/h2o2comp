import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm

##################################
# load data and tidy
##################################

raw = pd.read_excel('../data/Ma_fig_1_raw_data_compiled.xlsx',sheet_name='raw_data')

# long-form replicates, drop missing or non-positive counts
tidy = raw.melt(id_vars=['strain','temperature','treatment','Time(days)'],
                value_vars=['rep1','rep2','rep3'],
                var_name='replicate',value_name='cells')
tidy = tidy.dropna(subset=['cells','Time(days)','strain','temperature','treatment'])
tidy = tidy[tidy['cells'] > 0].copy()

tidy['Time(days)'] = pd.to_numeric(tidy['Time(days)'],errors='coerce')
tidy['treatment']  = pd.to_numeric(tidy['treatment'],errors='coerce')
tidy['temperature'] = tidy['temperature'].astype(str).str.strip()

# 26D treatment is dropped: protocol variant unknown
tidy = tidy[tidy['temperature'] != '26D'].copy()

tidy['ln_cells'] = np.log(tidy['cells'])

##################################
# helpers
##################################

def fit_log_linear(g):
    X = sm.add_constant(g['Time(days)'].astype(float))
    y = g['ln_cells'].astype(float)
    return sm.OLS(y,X).fit()

def temp_sort_key(s):
    digits = ''.join(c for c in str(s) if (c.isdigit() or c == '.'))
    try:
        return float(digits)
    except ValueError:
        return np.inf

##################################
# fit mu from zero-H2O2 controls
##################################

controls = tidy[tidy['treatment'] == 0]
mu_models = {}
all_models = {}   # keyed by (strain, temp, h2o2)
mu_rows = []

for (strain,temp),g in controls.groupby(['strain','temperature'],sort=False):
    if len(g) < 3 or g['Time(days)'].nunique() < 2:
        continue
    m = fit_log_linear(g)
    mu_models[(strain,temp)] = m
    all_models[(strain,temp,0)] = m
    mu_rows.append({
        'strain': strain,
        'temperature': temp,
        'mu': m.params['Time(days)'],
        'mu_se': m.bse['Time(days)'],
        'mu_ci_low': m.conf_int().loc['Time(days)',0],
        'mu_ci_high': m.conf_int().loc['Time(days)',1],
        'r2': m.rsquared,
        'n_obs': len(g),
    })

mu_table = pd.DataFrame(mu_rows)

##################################
# fit kdam from H2O2-treated curves
##################################

treated = tidy[tidy['treatment'] > 0]
kdam_rows = []

for (strain,temp,h2o2),g in treated.groupby(['strain','temperature','treatment'],sort=False):
    if (strain,temp) not in mu_models:
        continue
    if len(g) < 3 or g['Time(days)'].nunique() < 2:
        continue

    mu_m = mu_models[(strain,temp)]
    tr_m = fit_log_linear(g)
    all_models[(strain,temp,h2o2)] = tr_m

    mu = mu_m.params['Time(days)']
    mu_se = mu_m.bse['Time(days)']
    slope_eff = tr_m.params['Time(days)']
    slope_eff_se = tr_m.bse['Time(days)']

    # kdam = (mu - effective slope) / [H2O2]
    kdam = (mu - slope_eff) / h2o2
    kdam_se = np.sqrt(mu_se**2 + slope_eff_se**2) / h2o2
    kdam_ci_low  = kdam - 1.96*kdam_se
    kdam_ci_high = kdam + 1.96*kdam_se
    kdam_nonzero = not (kdam_ci_low <= 0 <= kdam_ci_high)

    kdam_rows.append({
        'strain': strain,
        'temperature': temp,
        'h2o2': h2o2,
        'mu': mu,
        'mu_se': mu_se,
        'slope_eff': slope_eff,
        'slope_eff_se': slope_eff_se,
        'kdam': kdam,
        'kdam_se': kdam_se,
        'kdam_ci_low': kdam_ci_low,
        'kdam_ci_high': kdam_ci_high,
        'kdam_nonzero': kdam_nonzero,
        'r2': tr_m.rsquared,
        'n_obs': len(g),
    })

coef = pd.DataFrame(kdam_rows)

##################################
# write tabular outputs
##################################

tidy.to_csv('../data/tidy_growth_data.csv',index=False)
mu_table.to_csv('../data/mu_from_zero_h2o2_controls.csv',index=False)
coef.to_csv('../data/kdam_estimates.csv',index=False)

# kdam point estimate with 95% CI, sample size, and R^2
kdam_table = (coef[['strain','temperature','h2o2',
                    'kdam','kdam_ci_low','kdam_ci_high','n_obs','r2']]
              .rename(columns={'h2o2':'h2o2_pmol_per_mL',
                               'kdam':'kdam_mL_per_pmol_per_day',
                               'kdam_ci_low':'ci95_low',
                               'kdam_ci_high':'ci95_high'})
              .sort_values(['strain','temperature','h2o2_pmol_per_mL'])
              .reset_index(drop=True))
kdam_table.to_csv('../data/kdam_table_clean.csv',index=False)

with pd.ExcelWriter('../data/prochlorococcus_kdam_outputs.xlsx',engine='openpyxl') as w:
    tidy.to_excel(w,sheet_name='tidy_data',index=False)
    mu_table.to_excel(w,sheet_name='mu_controls',index=False)
    coef.to_excel(w,sheet_name='kdam_estimates',index=False)
    kdam_table.to_excel(w,sheet_name='kdam_table_clean',index=False)

##################################
# kdam heatmap (MIT9312 and MED4)
##################################

strains = ['MIT9312','MED4']

# strain heatmap matrix, with cells masked where CI crosses zero
def heatmap_matrix(strain):
    sub = coef[coef['strain'] == strain]
    temps = sorted(sub['temperature'].astype(str).unique(),key=temp_sort_key)
    ros = sorted(sub['h2o2'].dropna().unique())
    mat = pd.DataFrame(np.nan,index=ros,columns=temps,dtype=float)
    for _,r in sub.iterrows():
        if r['kdam_nonzero']:
            mat.loc[r['h2o2'],str(r['temperature'])] = r['kdam']
    return mat

mats = [heatmap_matrix(s) for s in strains]
finite_vals = np.concatenate([m.to_numpy()[np.isfinite(m.to_numpy())] for m in mats])
vmax = np.nanmax(finite_vals) if finite_vals.size else 1.0

f,ax = plt.subplots(1,2,figsize=[12,5.2],constrained_layout=True)

for (a,strain,mat) in zip(ax,strains,mats):
    data = np.ma.masked_invalid(mat.to_numpy(dtype=float))
    x = np.arange(data.shape[1]+1)
    y = np.arange(data.shape[0]+1)
    pcm = a.pcolormesh(x,y,data,shading='flat',vmin=0,vmax=vmax)
    a.set_xticks(np.arange(data.shape[1])+0.5)
    a.set_xticklabels(mat.columns.tolist())
    a.set_yticks(np.arange(data.shape[0])+0.5)
    a.set_yticklabels([str(int(v)) if float(v).is_integer() else str(v) for v in mat.index])
    a.set_xlabel('Temperature treatment')
    a.set_ylabel(r'H$_2$O$_2$ treatment (pmol mL$^{-1}$)')
    a.set_title(strain)
    a.invert_yaxis()
    a.set_xticks(np.arange(data.shape[1]+1),minor=True)
    a.set_yticks(np.arange(data.shape[0]+1),minor=True)
    a.grid(which='minor',color='white',linewidth=1.0)
    a.tick_params(which='minor',bottom=False,left=False)

cbar = f.colorbar(pcm,ax=ax,shrink=0.95)
cbar.set_label(r'$k_{dam}$ (mL pmol$^{-1}$ day$^{-1}$)')
f.suptitle(r'Estimated $k_{dam}$ by temperature and H$_2$O$_2$ treatment',y=1.03)

f.savefig('../figures/figureS3',dpi=300,bbox_inches='tight')

##################################
# time-dependent dynamics per strain (SI, portrait)
##################################

from matplotlib.lines import Line2D

letters = 'abcdefghijklmnopqrstuvwxyz'

si_fig = {'MIT9312':'S1','MED4':'S2'}

for strain in strains:
    sdf = tidy[tidy['strain'] == strain]
    temps = sorted(sdf['temperature'].unique(),key=temp_sort_key)
    h2o2s = sorted(sdf['treatment'].dropna().unique())

    cmap = plt.get_cmap('viridis',max(len(h2o2s),2))
    h2o2_to_color = {h: cmap(i) for i,h in enumerate(h2o2s)}

    ncols = 2
    nrows = int(np.ceil(len(temps)/ncols))
    f,ax = plt.subplots(nrows,ncols,figsize=[7,2.4*nrows],sharey=True)
    ax = np.atleast_2d(ax)
    axes_flat = ax.flatten()

    for (a,temp,l) in zip(axes_flat,temps,letters):
        tsub = sdf[sdf['temperature'] == temp]
        for h2o2 in h2o2s:
            hsub = tsub[tsub['treatment'] == h2o2]
            if hsub.empty:
                continue
            c = h2o2_to_color[h2o2]
            a.scatter(hsub['Time(days)'],hsub['cells'],s=18,color=c,alpha=0.85)
            key = (strain,temp,h2o2)
            if key in all_models:
                m = all_models[key]
                tx = np.linspace(hsub['Time(days)'].min(),hsub['Time(days)'].max(),50)
                ty = np.exp(m.params['const'] + m.params['Time(days)']*tx)
                a.plot(tx,ty,color=c,lw=1.5)
        a.semilogy()
        a.text(0.05,0.92,l,transform=a.transAxes,fontsize=12,fontweight='bold')
        a.text(0.95,0.05,'T = '+str(temp),transform=a.transAxes,ha='right',va='bottom',fontsize=10)

    for a in axes_flat[len(temps):]:
        a.axis('off')

    for a in ax[-1,:]:
        a.set_xlabel('Time (days)')
    for a in ax[:,0]:
        a.set_ylabel('Cells (mL$^{-1}$)')

    handles = [Line2D([0],[0],marker='o',color=h2o2_to_color[h],lw=1.5,
                      label=str(int(h))+r' pmol mL$^{-1}$') for h in h2o2s]
    leg = f.legend(handles=handles,title=r'H$_2$O$_2$',loc='upper center',
                   ncol=len(h2o2s),bbox_to_anchor=(0.5,1.0),frameon=False)

    f.suptitle('$'+strain+'$',y=1.02,fontsize=12)
    f.tight_layout()
    f.savefig('../figures/figure'+si_fig[strain],dpi=300,bbox_inches='tight')
    plt.close(f)

