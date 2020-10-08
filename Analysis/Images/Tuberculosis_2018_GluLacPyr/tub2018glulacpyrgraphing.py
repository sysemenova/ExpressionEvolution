import os
import pandas as pd
import matplotlib.pyplot as plt
# import numpy as np
wow = pd.read_csv("mtbhr2.csv")

typenames = {"gluc_exp": "Glucose 6h",
             "gluc_sta": "Glucose 24h",
             "lac_exp": "Lactose 6h",
             "lac_sta": "Lactose 24h",
             "pyr_exp": "Pyruvate 6h",
             "pyr_sta": "Pyruvate 24h"}
typecolors = {"gluc_exp": "#83f556",
             "gluc_sta": "#51a630",
             "lac_exp": "#6acdf7",
             "lac_sta": "#2e8bb3",
             "pyr_exp": "#fc3535",
             "pyr_sta": "#b32525"}

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
    ax.legend(makers, typenames.values(), numpoints=1)
    ax.set_title("$R^2$ Of Conditions")
    ax.set_ylabel("$R^2$")
    ax.set_xlabel("Condition")
    plt.xticks(rotation=20, horizontalalignment="right", fontsize = 7)
    plt.show()