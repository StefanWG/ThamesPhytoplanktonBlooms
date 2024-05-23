from utils import *
from dataProcessing import readData
from thresholds import *
from copy import deepcopy
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
from PIL import Image

plt.rcParams['text.usetex'] = True


def multiParamPlot(location, outputFolder):
    '''
    Create multi-parameter bubble plots.
    Each row corresponds to a target variable.
    Each column corresponds to a combination of 
    two parameters.

    Inputs:
        - location (str): location
        - outputFolde (str): output folder
    '''
    # TODO: Size of text, border on legends, scales dont match, add boxes (using thresholds)
    # TODO: Size of bubbles
    data = readData(location)
    fig = plt.figure(figsize=(9, 12))

    for idx, var in enumerate(TARGETVARS):
        row = idx+1
        for col in [1,2,3]:
            increases = deepcopy(data)
            increases["diff"] = increases[var].diff()
            diffs = increases[increases["diff"] > 0]

            ax = fig.add_subplot(len(TARGETVARS),3,(row-1)*3+col)
            diffScale = 100/max(diffs["diff"])
            c = COLORS[idx]

            if col == 1: # Temp and Flow, Algae type text
                ax.scatter(diffs["temp"], diffs["flow"], diffs["diff"]*diffScale,
                    edgecolors='black', c=[c]*len(diffs))
                ax.text(max(diffs["temp"])*17/25, max(diffs["flow"])*5/6, 
                        s=LABELS[var], ha="center", color=c)
            elif col == 2:  # Temp and Sun
                ax.scatter(diffs["temp"], diffs["sun3d"], diffs["diff"]*diffScale,
                        edgecolors='black', c=[c]*len(diffs))
            elif col == 3: # Flow and Sun, Legend
                scatter = ax.scatter(diffs["flow"], diffs["sun3d"], diffs["diff"]*diffScale,
                        edgecolors='black', c=[c]*len(diffs))
                kw = dict(prop="sizes", num=3, color=c, func=lambda s: s / diffScale)
                handles, labels = scatter.legend_elements(**kw)
                ax.legend(handles, labels, loc="upper right")    
            if row == len(TARGETVARS): # X axis labels
                xlabels = ["Water Temperature (°C)", "Water Temperature (°C)", r'River Flow (m$^3$/s)']
                ax.set_xlabel(xlabels[col-1])   
    fig.tight_layout(w_pad=3.0)
    fig.text(-0.02, 0.52, r'River Flow (m$^3$/s)', va='center', rotation='vertical')
    fig.text(0.33, 0.52, 'Sunlight Hours (3 Day Mean)', va='center', rotation='vertical')
    fig.text(0.67, 0.52, 'Sunlight Hours (3 Day Mean)', va='center', rotation='vertical')


    plt.savefig(f"{outputFolder}/{location}_multiParamCells.jpg", 
                bbox_inches='tight')
    plt.close()

def timeSeriesPlot(location, outputFolder):
    '''
    Generate time series plot.

    Inputs:
        - location (str): location
        - ouputFolder (str): output folder
    '''
    # TODO: CHECK FILE PATHS
    data = readData(location)

    weeksSinceStart = getWeeksSinceStart(data["date"])
    dates = getGapFilledDates(weeksSinceStart, data["date"][0])

    for tvar in TARGETVARS:
        fig = plt.figure(figsize=(12, 8)) 

        axes = []
        smoothed, peaks, mins, maxs  = getBlooms(data["date"], data[tvar], method="find_peaks",sigma=SIGMA)
        years = mdates.YearLocator()
        vars = [tvar, "temp", "flow", "sun3d", "nutrients"]

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
                ax.plot(dates, var, color = COLORSTAB[idx], label="Observed")
                ax.set_ylabel(LABELSUNITS[v])

                # Add thresholds
                if v != vars[0]:
                    minThresh, maxThresh = getThresholds(data[v], mins, maxs, weeksSinceStart)
                    if v != "flow":
                        ax.axhline(minThresh, ls="--", color = COLORSTAB[idx], label="Min Threshold")
                    if v != "sun3d":
                        ax.axhline(maxThresh, ls=":", color = COLORSTAB[idx], label="Max Threshold")
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

        plt.savefig(f"{outputFolder}/{location}_{tvar}.jpg", 
                    bbox_inches="tight")
        plt.close()

def bloomTimingComparison(location, outputFolder):
    data = readData(location)

    weeksSinceStart = getWeeksSinceStart(data["date"])
    dates = getGapFilledDates(weeksSinceStart, data["date"][0])
    fig = plt.figure()

    for idx, tvar in enumerate(TARGETVARS):
        ax = fig.add_subplot(len(TARGETVARS), 2, idx*2+1)
        for idxm, method in enumerate(BLOOMMETHODS):
            smoothed, peaks, mins, maxs = getBlooms(data["date"], data[tvar], method=method,sigma=3)
            ymin, ymax = idxm / len(BLOOMMETHODS), (idxm + 1) / len(BLOOMMETHODS)
            for i in range(len(peaks)):
                ax.axvspan(dates[mins[i]], dates[maxs[i]], ymin=ymin, ymax=ymax, 
                        color=COLORS[idxm], label=method) 

            ax.set_yticks([])
            ax.set_ylabel(tvar)
        if idx == len(TARGETVARS) - 1:
            handles, labels = ax.get_legend_handles_labels()
            idxl = [np.where(np.array(labels) == m)[0][0] for m in BLOOMMETHODS]
            ax.legend(handles = list(np.array(handles)[idxl]),
                      labels = list(np.array(labels)[idxl]),
                      ncol=3, bbox_to_anchor=(1.25, -0.38))

        ax = fig.add_subplot(len(TARGETVARS),2,(idx+1)*2)
        smoothed = getSmoothedData(data["date"], data[tvar], SIGMA)
        smoothedSorted = np.sort(data[tvar])
        ax.plot(smoothedSorted, np.linspace(0, 1, len(smoothedSorted)))
        knee = getKneedleThreshold(smoothed)
        infl = getInflectionThreshold(smoothed)
        ax.vlines(knee, 0, 1,linestyles=":",label="knee")
        ax.vlines(infl, 0, 1,linestyles="--",label="inflection")

        if idx == len(TARGETVARS) - 1:
            ax.legend(ncol=2, bbox_to_anchor=(1, -0.38))
    fig.savefig(f"{outputFolder}/{location}_methodComp.jpg")
    plt.close()

def thresholds(outputFolder):
    '''
    Create plots illustrating the thresholds produced my different
    method of defining blooms.
    '''
    thresholds = []
    for loc in LOCATIONS:
        data = readData(loc)
        weeksSinceStart = getWeeksSinceStart(data["date"])
        for tvar in TARGETVARS:
            vars = ["temp", "flow", "sun3d"]
            for method in ["kneedle", "find_peaks", "inflection"]:
                smoothed, peaks, mins, maxs = getBlooms(data["date"], data[tvar], method=method,sigma=3)
                for idx, v in enumerate(vars):
                    minThresh, maxThresh = getThresholds(data[v], mins, maxs, weeksSinceStart)
                    thresholds.append({"location":loc,"target":tvar, "var":v, "method":method, "min":minThresh,"max":maxThresh})

    
    df = pd.DataFrame.from_records(thresholds)
    mask = df["var"] == "flow"
    df.loc[mask, "min"] = 0

    mfig = plt.figure(figsize=(12,12))
    mfig.tight_layout()
    gspec = gridspec.GridSpec(3, 3, wspace=-.2, hspace=-.5)

    for i, tvar in enumerate(TARGETVARS):
        w = df[df["target"]==tvar]
        fig = mfig.add_subfigure(gspec[i])
        fig.suptitle(tvar)
        axes = {
            "temp":fig.add_subplot(3,1,1),
            "flow":fig.add_subplot(3,1,2),
            "sun3d":fig.add_subplot(3,1,3)
        }

        scales = {
            "temp":30,
            "flow":100,
            "sun3d":20
        }
        locDict = {l:i*4 for i,l in enumerate(df["location"].unique())}
        offsets = {
            "find_peaks":-1,
            "kneedle":0,
            "inflection":1
        }
        colors = {
            "find_peaks":"b",
            "kneedle":"r",
            "inflection":"g"
        }

        for idx, row in w.iterrows():
            ax = axes[row["var"]]
            loc = row["location"]
            ax.axvline(x=locDict[loc] + offsets[row["method"]],ymin=row["min"]/scales[row["var"]], 
                       ymax=row["max"]/scales[row["var"]], color=colors[row["method"]],lw=9)
            
        for k, ax in axes.items():
            ax.set_ylabel(k)
            ax.set_ylim(0, scales[k])
            ax.set_xlim(0-2, 5*4+2)
            ax.set_xticks(list(locDict.values()), list(locDict.keys()))
        mfig.savefig(f"{outputFolder}/thresholds.jpg")
        plt.close()


def groupPlots(outputFolder):
    '''
    Create with tiling of all images of same type.

    Inputs:
        - outputFolder (str): output folder
    '''
    for var in TARGETVARS + ["linModel", "methodComp"]:
        fig = plt.figure(figsize=(20,20))

        for i, loc in enumerate(["cherwell", "kennet", "ock", "ray", "runnymede","thame"]):
            ax = fig.add_subplot(2,3,i+1)
            img = np.asarray(Image.open(f'figures/{loc}_{var}.jpg'))
            ax.axis('off')
            ax.imshow(img)
            ax.set_title(loc)
        fig.tight_layout()
        fig.subplots_adjust(bottom=0, top=1, hspace= -0.68 if var != "linModel" else -.5)
        fig.savefig(f"{outputFolder}/{var}.jpg", bbox_inches="tight")
        plt.close()


