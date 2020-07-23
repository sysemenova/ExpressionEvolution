import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
os.chdir("C:\\Users\\visitor\\Documents\\ExLinComb")
wow = pd.read_csv("rel06r2cs.csv")

# MGSO4
# Time vs r with concentration as color.
yes = False
if yes:
    fig, ax = plt.subplots()
    scatter = ax.scatter(wow["Growth Time"], wow["Mgso4 r2"], c = wow["Concentration"])
    handles, labels = scatter.legend_elements(prop = "colors", alpha = 0.6)
    legend2 = ax.legend(handles, labels, title = "Concentration")
    ax.set_ylabel("R value")
    ax.set_xlabel("Growth Time (h)")
    ax.set_title("MgSO4: Time vs R")
    plt.show()

# Similarity to glucose vs r with concentration as color
yes = False
if yes:
    fig, ax = plt.subplots()
    tempConcentration = np.log(wow["Concentration"])
    scatter = ax.scatter(wow["Similarity to glucose_4"], wow["Mgso4 r2"], c = tempConcentration)
    handles, labels = scatter.legend_elements(prop = "colors")
    legend2 = ax.legend(handles, labels, title = "Concentration")
    ax.set_ylabel("R value")
    ax.set_xlabel("Similarity to glucose_4")
    ax.set_title("MgSO4: Similarity to glucose vs r (with concentration)")
    plt.show()
    
# Similarity to glucose vs concentraton with batch as color (copy of Dinara's plots)
yes = True
if yes:
    fig, ax = plt.subplots()
    temp = wow[wow["Type"] == "mgso4"]
    temp = temp[temp["Concentration"] >= 0.8]
    scatter = ax.scatter(temp["Concentration"], temp["Similarity to glucose_4"], c = temp["Batch"])
    handles, labels = scatter.legend_elements(prop = "colors")
    legend2 = ax.legend(handles, labels, title = "Batch")
    ax.set_xlabel("Concentration")
    ax.set_xscale("log")
    ax.set_ylabel("Similarity to glucose_4")
    ax.set_title("MgSO4: Similarity to glucose vs concentration (with batch number)")
    plt.show()
    
# Concentration vs r with time as color
yes = False
if yes:
    fig, ax = plt.subplots()
    scatter = ax.scatter(wow["Concentration"], wow["Mgso4 r2"], c = wow["Growth Time"])
    handles, labels = scatter.legend_elements(prop = "colors", alpha = 0.6)
    legend2 = ax.legend(handles, labels, title = "Growth Time")
    ax.set_ylabel("R value")
    ax.set_xlabel("Concentration (mM)")
    ax.set_title("MgSO4: Concentration vs R")
    plt.show()

# NACL
yes = False
if yes:
    fig, ax = plt.subplots()
    onlyNacl = wow[wow["Type"] == "nacl"]
    scatter = ax.scatter(onlyNacl["Growth Time"], onlyNacl["R2s"], c = onlyNacl["Concentration"])
    handles, labels = scatter.legend_elements(prop = "colors", alpha = 0.6)
    legend2 = ax.legend(handles, labels, title = "Concentration")
    ax.set_ylabel("R value")
    ax.set_xlabel("Growth Time (h)")
    ax.set_title("NaCl: Time vs R")
    plt.show()

# ALL
# Time vs r with type as color. Toggle "log" on and off
yes = False
if yes:
    colors = {"glucose": u"red", "glycerol": u"blue", "nacl": u"yellow", "lactate": u"black", "gluconate": u"purple", "mgso4": u"orange"}
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
    colors = {"glucose": u"red", "glycerol": u"blue", "nacl": u"yellow", "lactate": u"black", "gluconate": u"purple", "mgso4": u"orange"}
    fig, axs = plt.subplots(3, 2, sharex = True, sharey = True)
    
    grpd = wow.groupby("Type")
    x = 0
    y = 0
    for name, data in grpd:
        axs[x, y].scatter(data["Growth Time"], data["R2s"], c = colors[name])
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
    
# Time vs r with type as color and batch as color
yes = False
if yes:
    colors = {"glucose": u"red", "glycerol": u"blue", "nacl": u"yellow", "lactate": u"black", "gluconate": u"purple", "mgso4": u"orange"}
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
    colors = {"glucose":"red", "glycerol":"blue", "nacl":"yellow", "lactate":"black", "gluconate":"purple", "mgso4":"orange"}
    scatter = ax.scatter(wow["Similarity to mgso4_400_7.2"], wow["R2s"], c = wow["Type"].map(colors))
    makers = [plt.Line2D([0,0],[0,0], color = color, marker = 'o', linestyle = '') for color in colors.values()]
    ax.legend(makers, colors.keys(), numpoints = 1)
    ax.set_ylabel("R value")
    ax.set_xlabel("Similarity to MgSO4 at 400 mM at 7 hours (Euclidean distance)")
    ax.set_title("Similarity to Best Condition vs R")
    plt.show()

# Similarity to best glucose
yes = False
if yes:
    fig, ax = plt.subplots()
    colors = {"glucose":"red", "glycerol":"blue", "nacl":"yellow", "lactate":"black", "gluconate":"purple", "mgso4":"orange"}
    scatter = ax.scatter(wow["Similarity to glucose_4"], wow["R2s"], c = wow["Type"].map(colors))
    makers = [plt.Line2D([0,0],[0,0], color = color, marker = 'o', linestyle = '') for color in colors.values()]
    ax.legend(makers, colors.keys(), numpoints = 1)
    ax.set_ylabel("R value")
    ax.set_xlabel("Similarity to Glucose base at 4 hours (Euclidean distance)")
    ax.set_title("Similarity to Best Glucose Condition vs R")
    plt.show()
    
# Plot only glucose and Mg
yes = False
if yes:
    colors = {"glucose": u"red", "glycerol": u"blue", "nacl": u"yellow", "lactate": u"black", "gluconate": u"purple", "mgso4": u"orange"}
    fig, ax = plt.subplots()
    options = ["glucose", "mgso4"]
    wow = wow[wow["Type"].isin(options)]
    
    color_map = wow["Type"].map(colors)
    scatter = ax.scatter(wow["Growth Time"], wow["R2s"], c = color_map)
    makers = [plt.Line2D([0,0],[0,0], color = color, marker = 'o', linestyle = '') for color in colors.values()]
    ax.legend(makers, colors.keys(), numpoints = 1)
    ax.set_ylabel("R value")
    ax.set_xlabel("Growth Time (h)")
    ax.set_xscale("log")
    ax.set_title("Growth Time vs R of Glucose and MgSO4 (logged)")
    plt.show()










    
