import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
os.chdir("C:\\Users\\visitor\\Documents\\ExLinComb\\Images")
wow = pd.read_csv("r2.csv")

typenames = {"limite_carbon": "Limited Carbon",
             "one_carbon": "One Carbon Source",
             "rich": "Rich Media",
             "stationary": "Stationary Phase",
             "stress": "Stress Conditions"}

barchart = False
simil = True

# Bar chart
if barchart:
    fig, ax = plt.subplots()
    groups = wow.groupby("Type")
    for name, dt in groups:
        scatter = ax.bar(dt["Conditions"], dt["R2s"], label = typenames[name])
    ax.legend()
    ax.set_title("R^2 Of Conditions")
    ax.set_ylabel("R^2")
    ax.set_xlabel("Condition")
    plt.xticks(rotation = 75)
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
    
