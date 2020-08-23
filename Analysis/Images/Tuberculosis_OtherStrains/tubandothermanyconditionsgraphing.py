import os
import pandas as pd
import matplotlib.pyplot as plt
# import numpy as np
os.chdir("C:\\Users\\sysem\\Documents\\ExpressionEvolution\\Analysis\\Images\\Tuberculosis_OtherStrains")
wow = pd.read_csv("mtb_and_others_r2.csv")

typecolors = {"0.1% Byturate": "#83f556",
              "0.2% Glucose + 0.1% Butyrate": "#51a630",
              "0.4% Glucose": "#6acdf7",
              "High Iron": "#2e8bb3",
              "Low Iron, 1 Day": "#fc3535",
              "Low Iron, 1 Week": "#b32525",
              "Tolaxapol, pH 7": "#666666",
              "Tolaxapol, pH 5.5": "#222222"}

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