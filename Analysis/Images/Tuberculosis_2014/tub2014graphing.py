import pandas as pd
import matplotlib.pyplot as plt
# import numpy as np
wow = pd.read_csv("tub2014.csv")

typecolors = {"Dextrose Exponential Phase": "#eb9e16",
              "Dextrose Stationary Phase": "#f1ea49",
              "Fatty Acids Exponential Phase": "#5dd95d",
              "Fatty Acids Stationary Phase": "#3893d2"}

barchart = True
# Simple bar chart
if barchart:
    fig, ax = plt.subplots()
    dt = wow
    scatter = ax.bar(dt["Conditions"], dt["R2s"])
    justtypes = list(dt["Conditions"])
    for i in range(0, len(justtypes)):
        scatter[i].set_color(typecolors[justtypes[i]])
    makers = [plt.Line2D([0, 0], [0, 0], color=color, marker='o', linestyle='') for color in typecolors.values()]
    ax.legend(makers, typecolors.keys(), numpoints=1)
    ax.set_title("$R^2$ Of Conditions")
    ax.set_ylabel("$R^2$")
    ax.set_xlabel("Condition")
    plt.xticks(rotation=10, horizontalalignment="right", fontsize = 7)
    plt.show()