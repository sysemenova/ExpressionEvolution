import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
wow = pd.read_csv("rel606_data.csv")


# Compare carbon sources, boxplot
yes = False
if yes:
    options = ["glucose", "gluconate", "lactate", "glycerol"]
    dt = wow[wow["Type"].isin(options)]
    #dt = dt[dt["Phase"] == "exponential"]
    #print(dt.head())
    ax = sns.boxplot(x = "Type", y = "R2s", data = dt)
    ax.set_title("One Carbon Sources")
    #ax.show()
    plt.show()

# MGSO4
# High Mg Strangeness Plot
yes = False
if yes:
    fig, ax = plt.subplots()
    dt = wow[wow["Mg Type"] == "highMg"]
    dt = dt[dt["Phase"] == "stationary"]
    
    scatter = ax.scatter(np.log(dt["Concentration"]), dt["Correlation to glucose_24.1"])
    ax.set_ylabel("Correlation to best stationary glucose expression")
    ax.set_xlabel("MgSO4 Concentration, log")
    ax.set_title("High MgSO4 Stress, Stationary Phase")
    plt.show()

# Simple scatter with number as color
yes = False
if yes:
    fig, ax = plt.subplots()
    dt = wow[wow["Type"] == "mgso4"]
    scatter = ax.scatter(dt["Growth Time"], dt["R2s"], c = dt["Concentration"])
    handles, labels = scatter.legend_elements(prop = "colors", alpha = 0.6)
    legend2 = ax.legend(handles, labels, title = "Concentration")
    ax.set_ylabel("R value")
    ax.set_xlabel("Growth Time (h)")
    ax.set_title("MgSO4: Time vs R")
    plt.show()
    
# Concentration vs r with time as color, subplots by phase.
# Change to NACL if necessary
yes = True
if yes:
    data = wow[wow["Type"] == "nacl"]
    grpd = data.groupby("Phase")
    fig, ax = plt.subplots(2, 1, sharex = True, sharey = True)
    x = 0
    for name, dt in grpd:
        scatter = ax[x].scatter(dt["Concentration"], dt["R2s"], c = dt["Growth Time"])
        handles, labels = scatter.legend_elements(prop = "colors")
        legend2 = ax[x].legend(handles, labels, title = "Growth Time")
        
        ax[x].set_ylabel("$R^2$ value")
        ax[x].set_xlabel("Concentration (mM)")
        names = {"exponential": "Exponential Phase", "stationary": "Stationary Phase"}
        title = "NaCl " + names[name] + ": Concentration vs $R^2$"
        ax[x].set_title(title)
        x += 1
    plt.show()


# ALL
colors = {"glucose": "#ef3a38", "glycerol": "#3893d2", "nacl": "#f1ea49", "lactate": "#2d2b2e", "gluconate": "#b240e3",
          "mgso4": "#eb9e16"}

# Time vs r with type as color. Toggle "log" on and off
yes = False
if yes:
    #colors = {"glucose": u"red", "glycerol": u"blue", "nacl": u"yellow", "lactate": u"black", "gluconate": "#b240e3", "mgso4": u"orange"}
    fig, ax = plt.subplots()
    color_map = wow["Type"].map(colors)
    scatter = ax.scatter(wow["Growth Time"], wow["R2s"], c = color_map)
    makers = [plt.Line2D([0,0],[0,0], color = color, marker = 'o', linestyle = '') for color in colors.values()]
    ax.legend(makers, colors.keys(), numpoints = 1)
    ax.set_ylabel("R value")
    ax.set_xlabel("Growth Time (h)")
    #ax.set_xscale("log")
    ax.set_title("Growth Time vs R")
    plt.show()

# Time vs r with type as color. But each gets their own subplot
yes = False
if yes:
    #colors = {"glucose": u"red", "glycerol": u"blue", "nacl": u"yellow", "lactate": u"black", "gluconate": u"purple", "mgso4": u"orange"}
    caps = {"glucose": "Glucose", "glycerol": "Glycerol", "nacl": "NaCl", "lactate": "Lactate", "gluconate": "Gluconate", "mgso4": "MgSO4"}
    fig, axs = plt.subplots(3, 2, sharex = True, sharey = True)
    
    grpd = wow.groupby("Type")
    x = 0
    y = 0
    for name, data in grpd:
        if name == "mgso4" or name == "nacl":
            scatter = axs[x, y].scatter(data["Growth Time"], data["R2s"], c = data["Concentration"], norm = matplotlib.colors.LogNorm())
            handles, labels = scatter.legend_elements(prop = "colors", alpha = 0.6)
            legend2 = axs[x, y].legend(handles, labels, title = "Concentration")
        else:
            axs[x, y].scatter(data["Growth Time"], data["R2s"], c = colors[name])
        axs[x, y].set_ylabel("$R^2$ value")
        axs[x, y].set_xlabel("Growth Time (h)")
        axs[x, y].set_xscale("log")
        title = "Growth time vs $R^2$: " + caps[name]
        axs[x, y].set_title(title)
        x += 1
        if x == 3:
            x = 0
            y = 1
    
    plt.show()
    
# Time vs r with batch as color
yes = False
if yes:
    #colors = {"glucose": u"red", "glycerol": u"blue", "nacl": u"yellow", "lactate": u"black", "gluconate": u"purple", "mgso4": u"orange"}
    fig, axs = plt.subplots(3, 2, sharex = True, sharey = True)
    
    grpd = wow.groupby("Type")
    x = 0
    y = 0
    for name, data in grpd:
        scatter = axs[x, y].scatter(data["Growth Time"], data["R2s"], c = data["Batch"])
        handles, labels = scatter.legend_elements(prop = "colors")
        legend2 = axs[x, y].legend(handles, labels, title = "Batch number")
        axs[x, y].set_ylabel("R value")
        axs[x, y].set_xlabel("Growth Time (h)")
        axs[x, y].set_xscale("log")
        title = "Growth time vs R: " + name
        axs[x, y].set_title(title)
        x += 1
        if x == 3:
            x = 0
            y = 1
    
    plt.show()

# Similarity to best condition
yes = False
if yes:
    fig, ax = plt.subplots()
    #colors = {"glucose":"red", "glycerol":"blue", "nacl":"yellow", "lactate":"black", "gluconate":"purple", "mgso4":"orange"}
    scatter = ax.scatter(wow["Similarity to mgso4_400_7.2"], wow["R2s"], c = wow["Type"].map(colors))
    makers = [plt.Line2D([0,0],[0,0], color = color, marker = 'o', linestyle = '') for color in colors.values()]
    ax.legend(makers, colors.keys(), numpoints = 1)
    ax.set_ylabel("R value")
    ax.set_xlabel("Similarity to MgSO4 at 400 mM at 7 hours (Euclidean distance)")
    ax.set_title("Similarity to Best Condition vs R")
    plt.show()

# Similarity to best condition but phase is color
yes = False
if yes:
    fig, ax = plt.subplots()
    colors = {"stationary":"#ef3a38", "exponential":"#5dd95d", "late_stationary":"#f1ea49"}
    scatter = ax.scatter(wow["Similarity to mgso4_400_7.2"], wow["R2s"], c = wow["Phase"].map(colors))
    makers = [plt.Line2D([0,0],[0,0], color = color, marker = 'o', linestyle = '') for color in colors.values()]
    ax.legend(makers, colors.keys(), numpoints = 1)
    ax.set_ylabel("$R^2$ value")
    ax.set_xlabel("Similarity to MgSO4 at 400 mM at 7 hours (Euclidean distance)")
    ax.set_title("Similarity to Best Condition vs $R^2$, Phase Colored")
    plt.show()

# Similarity to best glucose
yes = False
if yes:
    fig, ax = plt.subplots()
    #colors = {"glucose":"red", "glycerol":"blue", "nacl":"yellow", "lactate":"black", "gluconate":"purple", "mgso4":"orange"}
    scatter = ax.scatter(wow["Similarity to glucose_4"], wow["R2s"], c = wow["Type"].map(colors))
    makers = [plt.Line2D([0,0],[0,0], color = color, marker = 'o', linestyle = '') for color in colors.values()]
    ax.legend(makers, colors.keys(), numpoints = 1)
    ax.set_ylabel("R value")
    ax.set_xlabel("Similarity to Glucose base at 4 hours (Euclidean distance)")
    ax.set_title("Similarity to Best Glucose Condition vs R")
    plt.show()
    




