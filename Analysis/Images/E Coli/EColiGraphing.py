import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
os.chdir("C:\\Users\\visitor\\Documents\\ExLinComb\\Images")
wow = pd.read_csv("ecoli_data.csv")

typenames = {"limite_carbon": "Limited Carbon",
             "one_carbon": "One Carbon Source",
             "rich": "Rich Media",
             "stationary": "Stationary Phase",
             "stress": "Stress Conditions"}
typecolors = {"limite_carbon": "yellow",
              "one_carbon": "blue",
              "rich": "green",
              "stationary": "orange",
              "stress": "red"}

barchart = False
barchartw0 = False
partialsbar = True
simil = False


# Bar chart
if partialsbar:
    fig, ax = plt.subplots()
    dt = wow.sort_values("LB Partials")
    print(dt.head())
    scatter = ax.bar(dt["Conditions"], dt["LB Partials"])
    justtypes = list(dt["Type"])
    for i in range(0, len(justtypes)):
        scatter[i].set_color(typecolors[justtypes[i]])
    makers = [plt.Line2D([0,0],[0,0], color = color, marker = 'o', linestyle = '') for color in typecolors.values()]
    ax.legend(makers, typenames.values(), numpoints = 1)
    ax.set_title("How Much Do These Conditions Add to LB?")
    ax.set_ylabel("Partial Value")
    ax.set_xlabel("Condition")
    plt.xticks(rotation = 20, horizontalalignment = "right")
    plt.show()
    
# Bar chart
if barchart:
    fig, ax = plt.subplots()
    dt = wow.sort_values("R2s")
    print(dt.head())
    scatter = ax.bar(dt["Conditions"], dt["R2s"])
    justtypes = list(dt["Type"])
    for i in range(0, len(justtypes)):
        scatter[i].set_color(typecolors[justtypes[i]])
    makers = [plt.Line2D([0,0],[0,0], color = color, marker = 'o', linestyle = '') for color in typecolors.values()]
    ax.legend(makers, typenames.values(), numpoints = 1)
    ax.set_title("R^2 Of Conditions")
    ax.set_ylabel("R^2")
    ax.set_xlabel("Condition")
    plt.xticks(rotation = 20, horizontalalignment = "right")
    plt.show()

# Bar chart with 0s
if barchartw0:
    fig, ax = plt.subplots()
    dt = wow.sort_values("R2s w 0")
    #print(dt.head())
    scatter = ax.bar(dt["Conditions"], dt["R2s w 0"])
    justtypes = list(dt["Type"])
    for i in range(0, len(justtypes)):
        scatter[i].set_color(typecolors[justtypes[i]])
    makers = [plt.Line2D([0,0],[0,0], color = color, marker = 'o', linestyle = '') for color in typecolors.values()]
    ax.legend(makers, typenames.values(), numpoints = 1)
    ax.set_title("R^2 Of Conditions (with 0 filled)")
    ax.set_ylabel("R^2")
    ax.set_xlabel("Condition")
    plt.xticks(rotation = 20, horizontalalignment = "right")
    plt.show()

# Similarity to LB
if simil:
    fix, ax = plt.subplots()
    groups = wow.groupby("Type")
    for name, dt in groups:
        scatter = ax.scatter(dt["Similarity to LB"], dt["R2s"], label = typenames[name])
    ax.legend(title = "Condition Type")
    ax.set_xlabel("Euclidian distance to LB expression")
    #ax.set_xscale("log")
    ax.set_ylabel("R^2")
    ax.set_title("Similarity to LB vs R^2")
    plt.show()
    
