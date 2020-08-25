import pandas as pd
import matplotlib.pyplot as plt
# import numpy as np
wow = pd.read_csv("tub2014.csv")

typecolors = {"Dextrose_Exp": "#83f556",
              "Dextrose_sta": "#51a630",
              "Fatty_acids_exp": "#6acdf7",
              "Fatty_acids_sta": "#2e8bb3"}

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
    ax.set_title("R^2 Of Conditions")
    ax.set_ylabel("R^2")
    ax.set_xlabel("Condition")
    plt.xticks(rotation=20, horizontalalignment="right")
    plt.show()