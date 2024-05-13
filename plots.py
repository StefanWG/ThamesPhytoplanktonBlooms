import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
from copy import deepcopy
from utils import *

plt.rcParams['text.usetex'] = True

import warnings
warnings.filterwarnings("ignore")

data = readData()

##########################
# Multiple paramter plots
##########################
fig = plt.figure(figsize=(9, 12))

vars = ["diatoms", "n-chloro", "p-chloro", "totCyano"]
labs = ["Diatoms", "Nano-chlorophytes", "Pico-chlorophytes", "Cyanobacteria"]
for var, row, c in zip(vars, [1,2,3,4], ["b","r", "g", "y"]):
    for col in [1,2,3]:
        increases = deepcopy(data)
        increases["diff"] = increases[var].diff()
        diffs = increases[increases["diff"] > 0]

        ax = fig.add_subplot(4,3,(row-1)*3+col)
        if col == 1: # Temp and Flow, Algae type text
            ax.scatter(diffs["temp"], diffs["flow"], diffs["diff"]/1000,
                    edgecolors='black', c=[c]*len(diffs))
            ax.text(17, 250, s=labs[row-1], ha="center", color=c)
        elif col == 2:  # Temp and Sun
            ax.scatter(diffs["temp"], diffs["sun3d"], diffs["diff"]/1000,
                    edgecolors='black', c=[c]*len(diffs))
        elif col == 3: # Flow and Sun, Legend
            scatter = ax.scatter(diffs["flow"], diffs["sun3d"], s=diffs["diff"]/1000,
                    edgecolors='black', c=[c]*len(diffs))
            kw = dict(prop="sizes", num=3, color=c, func=lambda s: s*1000)
            handles, labels = scatter.legend_elements(**kw)
            ax.legend(handles, labels, loc="upper right")    
        if row == 4: # X axis labels
            xlabels = ["Water Temperature (°C)", "Water Temperature (°C)", r'River Flow (m$^3$/s)']
            ax.set_xlabel(xlabels[col-1])          
# Y axis labels           
fig.tight_layout(w_pad=2.0)
fig.text(-0.01, 0.52, r'River Flow (m$^3$/s)', va='center', rotation='vertical')
fig.text(0.34, 0.52, 'Solar Radiation (3 Day Mean)', va='center', rotation='vertical')
fig.text(0.67, 0.52, 'Solar Radiation (3 Day Mean)', va='center', rotation='vertical')


plt.savefig("figures/multiParamCells.jpg", bbox_inches='tight')
# TODO: Size of text, border on legends, scales dont match, add boxes (using thresholds)

##########################
# Time series plots
##########################

weeksSinceStart = [((i - data["date"][0]).days + 1) // 7 + 1 for i in data["date"]]

dates = [data["date"][0] + pd.Timedelta(weeks=i) for i in range(max(weeksSinceStart))]
xvars = ["diatoms", "n-chloro", "p-chloro", "totCyano", "totPhyto"]
# xvars = ["diatoms"]
for xvar in xvars:
    fig = plt.figure(figsize=(12, 8)) 
    vars = [xvar, "temp", "flow", "sun3d", "nutrients"]
    labels = {"diatoms":"Diatom Cells\n(m/l)", 
            "temp":"Water Temperature\n(°C)",
            "flow":r'River Flow' "\n" r'(m$^3$/s)',
            "sun3d":r'Solar Radiation' "\n" r'(W/m$^2$)',
            "n-chloro":"Nano-chlorophytes\n(m/l)", 
            "p-chloro":"Pico-chlorophytes\n(m/l)",
            "totCyano":"Cyanobacteria\n(m/l)",
            "totPhyto":"Total Phytoplankton\n(m/l)"
        }
    axes = []
    smoothed, peaks, mins, maxs  = getBlooms(data["date"], data[xvar], 3)
    years = mdates.YearLocator()
    colors = ["tab:brown", "black", "tab:blue", "tab:olive"]


    for idx, v in enumerate(vars):
        if v =="nutrients":
            ax = fig.add_subplot(len(vars),1,idx+1)
            # Create second y-axis
            ax1 = ax.twinx()
            axes.append(ax)
            axes.append(ax1)
            #Nitrate
            varMap = {a:b for a,b, in zip(weeksSinceStart, data["nitrate"])}
            var = np.array([varMap[i+1] if i+1 in varMap else np.nan for i in range(max(weeksSinceStart))])
            ax.plot(dates, var, label=r'Nitrate (N0$_3$)', color="tab:cyan")
            #Silicon
            varMap = {a:b for a,b, in zip(weeksSinceStart, data["silicon"])}
            var = np.array([varMap[i+1] if i+1 in varMap else np.nan for i in range(max(weeksSinceStart))])
            ax.plot(dates, var, label="Silicon (Si)", color="tab:gray")
            #Phos
            varMap = {a:b for a,b, in zip(weeksSinceStart, data["SRP"])}
            var = np.array([varMap[i+1] if i+1 in varMap else np.nan for i in range(max(weeksSinceStart))])
            ax1.plot(dates, var, label="Soluble Reactive Phosphorus (SRP)", color="tab:purple")
            #Labels
            ax.set_ylabel(r'N0$_3$ and Si' "\n(mg/l)")
            ax1.set_ylabel("SRP\n" r'($\mu$m/l)')

        else:
            # Plot observed time series data
            varMap = {a:b for a,b, in zip(weeksSinceStart, data[v])}
            var = np.array([varMap[i+1] if i+1 in varMap else np.nan for i in range(max(weeksSinceStart))])
            ax = fig.add_subplot(len(vars),1,idx+1)
            ax.plot(dates, var, color = colors[idx], label="Observed")
            ax.set_ylabel(labels[v])

            # Add thresholds
            if v != vars[0]:
                minThresh, maxThresh = getThresholds(data[v], mins, maxs, weeksSinceStart)
                if v != "flow":
                    ax.axhline(minThresh, ls="--", color = colors[idx], label="Min Threshold")
                if v != "sun3d":
                    ax.axhline(maxThresh, ls=":", color = colors[idx], label="Max Threshold")
            axes.append(ax)
        # Add geometry for blooms and cessations
        for i in range(len(peaks)):
            ax.axvspan(dates[mins[i]], dates[maxs[i]], color="#999999")
            ax.axvline(dates[peaks[i]], color="#444444", linestyle=":")
        ax.xaxis.set_major_locator(years)
        


    axes[-1].set_xlabel("Year")
    fig.align_ylabels()
    # Legend
    handles, labels = axes[1].get_legend_handles_labels()
    handles1, labels1 = axes[-2].get_legend_handles_labels()
    handles2, labels2 = axes[-1].get_legend_handles_labels()
    order = [0,3,1,4,2,5]
    h = handles+handles1+handles2
    l = labels+labels1+labels2
    fig.legend([h[idx] for idx in order], [l[idx] for idx in order], loc = "lower center", ncol=3)
    handles = [h.set_color("tab:red") for h in handles]

    plt.savefig(f'figures/{xvar.replace("-","")}.jpg', bbox_inches="tight")

print("Figures have been saved to ./figure directory")

