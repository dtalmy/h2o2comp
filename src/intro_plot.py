import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates

##################################
# load data and wrangle
##################################

# cells
md = pd.read_excel('../data/data_MEGA.xlsx',sheet_name='Morris_2016_biomass')
autdf = md.dropna(subset=['Pro/mL']) # just the autotrophs

# mlds
mlddf = pd.read_excel('../data/data_MEGA.xlsx',sheet_name='morris_2016_mlds')
mlddf = mlddf.transpose()
mlddf.columns = mlddf.iloc[0]
mlddf = mlddf.iloc[1:]
mlddf.index.name = 'Time'
mlddf = mlddf.reset_index()
mlddf.index.name = None

# hooh
hmd = pd.read_excel('../data/data_MEGA.xlsx',sheet_name='depth_HOOH_insitu_2016Morris')

##################################
# plot
##################################

f,ax = plt.subplots(1,2,figsize=[9,4.5])
ax[0].set_xlabel('Time')
ax[1].set_xlabel('Time')
ax[0].set_ylabel('H$_2$O$_2$ concentration (nM)')
ax[1].set_ylabel('cells (mL$^{-1}$)')

for (a,l) in zip(ax,'ab'):
    a.text(0.1,0.9,l,transform=a.transAxes)

###########################
# HOOH

# Merge max depth into data and convert all depth data to be floats
merged = hmd.merge(mlddf, on='Time', how='left')
merged['Depth'] = merged['Depth'].astype(float)
merged['Mixed Layer Depth'] = merged['Mixed Layer Depth'].astype(float)

# Filter where depth <= max_depth
hmld = merged[merged['Depth'] <= merged['Mixed Layer Depth']].drop(columns='Mixed Layer Depth')

# Group and compute mean and std (or sem)
grouped = hmld.groupby('Time')['initial  [HOOH]']
mean = grouped.mean()
err = grouped.std()  # or grouped.sem() for standard error

# Plot with log-scaled y-axis
ax[0].errorbar(mean.index, mean.values, yerr=err.values, fmt='-o')
ax[0].set_yscale('log')
ax[0].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
ax[0].tick_params(axis='x', rotation=45)

###########################
#Pro

# Group and compute mean and std (or sem)
grouped = autdf.groupby('Time')['Pro/mL']
mean = grouped.mean()
err = grouped.std()  # or grouped.sem() for standard error

# Plot with log-scaled y-axis
ax[1].errorbar(mean.index, mean.values, yerr=err.values, fmt='-o')
ax[1].set_yscale('log')
ax[1].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
ax[1].tick_params(axis='x', rotation=45)

f.subplots_adjust(bottom=0.2,wspace=0.3)
f.savefig('../figures/intro_plot',bbox_inches='tight',dpi=300)
