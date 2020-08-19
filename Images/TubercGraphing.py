import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
os.chdir("C:\\Users\\sysem\\Documents\\ExpressionEvolution\\Images")
wow = pd.read_csv("tuberc_r2.csv")

typecolors = {"Dextrose Exponential Phase": "#FF3333",
              "Dextrose Stationary Phase": "#FF8F33",
              "Dextrose Late Stationary Phase": "#FFF633",
              "Lipid Exponential Phase": "#33AEFF",
              "Lipid Stationary Phase": "#33FDFF",
              "Lipid Late Stationary Phase": "#33FFA2"}

barchart = True
# Simple bar chart
if barchart:
    fig, ax = plt.subplots()
    dt = wow.sort_values("R2s")
    scatter = ax.bar(dt["Conditions"], dt["R2s"])
    justtypes = list(dt["Type"])
    for i in range(0, len(justtypes)):
        scatter[i].set_color(typecolors[justtypes[i]])
    makers = [plt.Line2D([0, 0], [0, 0], color=color, marker='o', linestyle='') for color in typecolors.values()]
    ax.legend(makers, typecolors.keys(), numpoints=1)
    ax.set_title("R^2 Of Conditions")
    ax.set_ylabel("R^2")
    ax.set_xlabel("Condition")
    plt.xticks(rotation=20, horizontalalignment="right")
    plt.show()
