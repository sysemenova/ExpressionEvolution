# Documentation here

## Step 0: Import necessary libraries
print("Downloading libraires...")
import pandas as pd
import numpy as np
#from keras.layers import Input, Dense
#from keras.models import Model
from sklearn import linear_model
#from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
#from mpl_toolkits import mplot3d
import os
print("Done. Beginning program...")

def pr_r2(in_r2, data, which):
    # That - 2 is because Gene Name and Evo Rate are at the beginning. Will likely
    # Have to debug
    print("   ", which, "R^2: %1.4f" %(in_r2[data.columns.get_loc(which) - 2]))

def pr_r2_all(in_r2, data):
    cols = data.columns.copy()
    cols = cols[2:]
    zipped = zip(in_r2.copy(), cols)
    sor = sorted(zipped)
    sor.reverse()
    print("R^2:")
    for i in sor:
        print("   %-23s %1.4f" %(i[1]+":", i[0]))

def pr_coef_r2(coefficients, in_r2, data, cols):
    zipped = zip(coefficients.copy(), cols.copy())
    sor = sorted(zipped)
    print("Coefficients and R^2: ")
    for i in sor:
        print("   %-22s coef: %1.4f; R^2: %1.4f" %(i[1], i[0], in_r2[data.columns.get_loc(i[1])-2]))
    
def pr_coef(coefficients, cols):
    zipped = zip(coefficients.copy(), cols.copy())
    sor = sorted(zipped)
    print("Coefficients: ")
    for i in sor:
        print("   %-23s %1.4f" %(i[1]+":", i[0]))

def highest_coef(num, coefficients, cols):
    zipped = zip(coefficients.copy(), cols.copy())
    sor = sorted(zipped)
    toret = []
    for i in range(num):
        toret.append(sor[i][1])
    return toret

def lin_model(data, columns):
    expression_data = data[columns]
    evo_rates = data["Evo Rate"]
    reg_model = linear_model.LinearRegression()
    reg_model.fit(expression_data, evo_rates)
    print("R^2: ", reg_model.score(expression_data, evo_rates))
    #print("Intercept: ", reg_model.intercept_)
    coef = reg_model.coef_
    pr_coef(coef, columns)
    
    lin_comb = np.multiply(coef[0], expression_data.iloc[:, 0])
    for i in range(1, len(columns)):
        lin_comb = np.add(lin_comb, np.multiply(coef[i], expression_data.iloc[:, i]))

    plt.scatter(lin_comb, evo_rates, color = 'red')
    plt.show()
    
    return (coef, columns)


"""
Program notes:
 - May change this to a loop that goes through all folders in
   this directory. The current structure goes through species
   based on the list right below program notes
"""
# List of all species that we have. May change.
species_names = ["Escherichia_coli"]

# Get the full path before jumping into species
cur_path = os.getcwd()

# Go through all species
for species_name in species_names:
    print()
    print("Current species:", species_name)
    
    ## Step 1: Obtain data
    """
    Step 1 Notes:
     - Will have to calculate the evolution rates in this step
    """

    # Get the current path and create the full path
    full_path = cur_path + "/" + species_name

    # If path does not exist, raise an exception. Else, change the path
    if not os.path.exists(full_path):
        raise Exception("The path", full_path, "does not exist.")
    else:
        os.chdir(full_path)

    # Go through all csv files and import their data. Keep track of num datasets
    data = {}
    i = 0
    for filename in os.listdir():
        if filename.endswith(".csv"):
            data[str(i)] = pd.read_csv(filename)
            i += 1
    num_datasets = i


    ## Step 2: Combine data
    """
    Step 2 Notes:
     - Figure out a way to merge just based on first column, or ensure that
       all first columns are named the same thing.
    """
    # Merge all datasets by their column named "Gene Name"
    full_data = data["0"]
    for i in range(1, num_datasets):
        full_data = pd.merge(full_data, data[str(i)], how = "outer", on = "Gene Name")
        
    ## Step 3.1: Process data
    """
    Step 3 Notes:
    """
    len_bef = len(full_data)
    
    expression_columns = []
    full_data.dropna(inplace = True)
    num_expression_columns = 0
    #ss = StandardScaler()
    
    for i in list(full_data):
        if i != "Gene Name":
            if i != "Evo Rate":
                    expression_columns.append(i)
                    num_expression_columns += 1
            # This section log's it if necessary
            #full_data = full_data[full_data[i] != 0]
            #full_data[i] = np.log(full_data[i])

            # Rank data here
            full_data[i] = full_data[i].rank()
            
    len_aft = len(full_data)
    print("Num rows dropped:", len_bef-len_aft)
    if float(len_bef-len_aft)/float(len_aft) > 0.05:
        print("Warning: around", 100*(len_bef-len_aft)/len_aft, "percent of your data has been dropped.")
        input()
    else:
        print("Dropped rows within acceptable parameters.")
        print()
    print("Processed data head: ")
    print(full_data.head(8))
    print()
    
    ## Step 4: Create a full linear regression model  
    print("Full model:")
    cur_model = lin_model(full_data, expression_columns)
    expression_data = full_data[expression_columns]
    evo_rates = full_data["Evo Rate"]
    indiv_r2 = []
    for i in range(0, num_expression_columns):
        r = linear_model.LinearRegression()
        r.fit(expression_data.iloc[:, [i]], evo_rates)
        indiv_r2.append(r.score(expression_data.iloc[:, [i]], evo_rates))
    print()

    # Edit the commands later
    print("Commands: drop __, __, ...; undrop __, __, ...; only __, __, ...; r2 __; next")
    flag = True
    dropped = []
    while flag:
        user_inp = input("Command: ")
        in_list = user_inp.replace(',','').split()
        valid_commands = ["drop", "r2", "undrop", "only", "next", "pr", "save", "highest"]
        if in_list[0].lower() not in valid_commands:
            print("Not valid command.")
        else:
            if in_list[0].lower() == "next":
                flag = False
            elif in_list[0].lower() == "r2" and in_list[1].lower() == "all":
                pr_r2_all(indiv_r2, full_data)
            elif in_list[0].lower() == "pr":
                if in_list[1].lower() == "coef":
                    pr_coef(cur_model[0], cur_model[1])
                elif in_list[1].lower() == "r2":
                    for i in cur_model[1]:
                        pr_r2(indiv_r2, full_data, i)
                elif in_list[1].lower() == "coef_r2":
                    pr_coef_r2(cur_model[0], indiv_r2, full_data, cur_model[1])
            elif in_list[0].lower() == "highest":
                cols = highest_coef(int(in_list[1]), cur_model[0], cur_model[1])
                cur_model = lin_model(full_data, cols)
            else:
                # Checks if all inputs are okay
                good_input = True
                for i in range(1, len(in_list)):
                    in_list[i] = in_list[i].replace("_", " ")
                    if not in_list[i] in full_data.columns:
                        good_input = False
                # If the input is fine, contiue
                if good_input:
                    com = in_list[0].lower()
                    if com == "drop":
                        for i in range(1, len(in_list)):
                            dropped.append(in_list[i])
                            expression_columns.remove(in_list[i])
                        cur_model = lin_model(full_data, expression_columns)
                    elif com == "r2":
                        for i in range(1, len(in_list)):
                            pr_r2(indiv_r2, full_data, in_list[i])
                    elif com == "undrop":
                        for i in range(1, len(in_list)):
                            expression_columns.append(in_list[i])
                            dropped.remove(in_list[i])
                        cur_model = lin_model(full_data, expression_columns)
                    elif com == "only":
                        cur_model = lin_model(full_data, in_list[1:])
                    
                    else:
                        print("Command", com, "is incomplete.")
                else:
                    print("Not valid inputs.")
        print()
    ## Step 5: Decide what to do from here
    print()
    print()
print("Program is over.")
input("Press enter to exit. ")
