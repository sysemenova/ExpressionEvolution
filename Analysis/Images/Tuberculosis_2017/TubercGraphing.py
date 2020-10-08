import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
wow = pd.read_csv("tuberc_2017_r2.csv")

typecolors = {"Dextrose Exponential Phase": "#ef3a38",
              "Dextrose Stationary Phase": "#eb9e16",
              "Dextrose NRP1 Phase": "#f1ea49",
              "Lipid Exponential Phase": "#5dd95d",
              "Lipid Stationary Phase": "#3893d2",
              "Lipid NRP1 Phase": "#b240e3"}

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
    ax.set_title("$R^2$ Of Conditions")
    ax.set_ylabel("$R^2$")
    ax.set_xlabel("Condition")
    plt.xticks(rotation=20, horizontalalignment="right")
    plt.show()
