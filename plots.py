import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
from copy import deepcopy
from utils import *
import sys
import kneed

# TODO: Comment? And improve this strucutre
def plotMethodComparison():
    for loc in ["runnymede", "ock", "ray", "cherwell", "kennet", "thame"]:
        data = readData(loc)
        weeksSinceStart = [((i - data["date"][0]).days + 1) // 7 + 1 for i in data["date"]]

        dates = [data["date"][0] + pd.Timedelta(weeks=i) for i in range(max(weeksSinceStart))]

        fig = plt.figure()
        i = 1
        for xvar in ["diatoms", "p-chloro", "n-chloro","totCyano", "totPhyto", "chl"]:
            ax = fig.add_subplot(6,2,i*2-1)
            # Kneedle
            smoothed, peaks, mins, maxs = getBlooms(data["date"], data[xvar], method="kneedle",sigma=3, kneedleS=1)

            [ax.axvspan(dates[mins[i]], dates[maxs[i]], ymin=0, ymax=.5, color="r", label="kneedle") for i in range(len(mins))]
            # Peacks
            smoothed, peaks, mins, maxs = getBlooms(data["date"], data[xvar], method="find_peaks",sigma=3)
            [ax.axvspan(dates[mins[i]], dates[maxs[i]], ymin=.5, ymax=1, color="b", label="find_peaks") for i in range(len(mins))]

            handles, labels = ax.get_legend_handles_labels()
            ax.set_yticks([])
            ax.set_ylabel(xvar)

            ax = fig.add_subplot(6,2,i*2)
            s = smooth(np.sort(smoothed), 3)
            perc = np.linspace(0,1,len(s))
            ax.plot(s, perc, label=xvar)
            knee = kneed.KneeLocator(s, perc).knee

            ax.vlines(knee, 0, 1,linestyles=":")
            i += 1

        ax.legend(handles[0:1] + handles[-1:], labels[0:1] + labels[-1:], loc="upper center", ncol=2, bbox_to_anchor=(0.5, -0.38))
        fig.suptitle(loc)
        fig.savefig(f"figures/{loc}_metComp.png")
        plt.close()

print(f"Figures have been saved to {output}")

