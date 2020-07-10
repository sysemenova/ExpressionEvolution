# Documentation here

## Step 0: Import necessary libraries
print("Downloading libraires...")
print("Pandas... ", end = '')
import pandas as pd
print("NumPy... ", end = '')
import numpy as np
print("SKLearn... ", end = '')
from sklearn import linear_model
from sklearn import metrics
print("SciPy... ", end = '')
from scipy.optimize import lsq_linear
from scipy import stats
print("MatPlotLib... ", end = '')
import matplotlib.pyplot as plt
print("OS... ")
import os

print("Done. Beginning program...")

evo_rate_name = "Evo Rate"

def pr_r2(in_r2, data, which):
    # That - 2 is because Gene Name and Evo Rate are at the beginning. Will likely
    # have to debug
    print(" ", which, "R: %1.4f, p-val: %1.4f" %(in_r2[data.columns.get_loc(which) - 2][0], in_r2[data.columns.get_loc(which) - 2][1]))

def pr_r2_all(in_r2, data):
    cols = data.columns.copy()
    cols = cols[2:]
    zipped = zip(in_r2.copy(), cols)
    sor = sorted(zipped)
    sor.reverse()
    print("R and p-val:")
    for i in sor:
        print("   %-23s %1.4f, %1.4f" %(i[1]+":", i[0][0], i[0][1]))

def pr_coef_r2(coefficients, in_r2, data, cols):
    zipped = zip(coefficients.copy(), cols.copy())
    sor = sorted(zipped)
    print("Coefficients and R: ")
    for i in sor:        
        print("   %-22s coef: %1.4f; R: %1.4f, p-val: %1.4f" %(i[1], i[0], in_r2[data.columns.get_loc(i[1])-2][0],in_r2[data.columns.get_loc(i[1])-2][1]))
    
def pr_coef(coefficients, cols):
    zipped = zip(coefficients.copy(), cols.copy())
    sor = sorted(zipped)
    print("Coefficients: ")
    for i in sor:
        print("   %-23s %1.8f" %(i[1]+":", i[0]))

def highest_coef(num, coefficients, cols):
    zipped = zip(coefficients.copy(), cols.copy())
    sor = sorted(zipped)
    toret = []
    for i in range(num):
        toret.append(sor[i][1])
    return toret

def lin_model(data, columns, alph = 1.0):
    expression_data = data[columns]
    evo_rates = data[evo_rate_name]
    
    reg_model = linear_model.LinearRegression()
    reg_model.fit(expression_data, evo_rates)
    print("---LINEAR REGRESSION---")
    print("R^2: ", reg_model.score(expression_data, evo_rates))
    coef = reg_model.coef_
    pred = reg_model.predict(data[columns])
    pr_coef(coef, columns)

    reg_model = linear_model.Ridge(alpha = alph)
    reg_model.fit(expression_data, evo_rates)
    print("---RIDGE---")
    print("R^2: ", reg_model.score(expression_data, evo_rates))
    coef = reg_model.coef_
    pr_coef(coef, columns)

    # Create the linear combination and plot it
    lin_comb = np.multiply(coef[0], expression_data.iloc[:, 0])
    for i in range(1, len(columns)):
        lin_comb = np.add(lin_comb, np.multiply(coef[i], expression_data.iloc[:, i]))

    #plt.scatter(lin_comb, evo_rates, color = 'red')
    #plt.show()
    
    return (coef, columns, pred)

def partial(data, x1, x2, x3):
    r12 = stats.pearsonr(data[x1], data[x2])[0]
    r13 = stats.pearsonr(data[x1], data[x3])[0]
    r23 = stats.pearsonr(data[x2], data[x3])[0]
    r123 = (r12-r13*r23)/np.sqrt((1-r13*r13)*(1-r23*r23))
    return r123

"""
Program notes:
"""
# List of all species that we have. May change to going through all folders
# marked with something (so Images isn't gone through.)
species_names = ["Escherichia_coli_REL606"]

# Get the full path before jumping into species
cur_path = os.getcwd()

# Go through all species
for species_name in species_names:
    print()
    print("Current species:", species_name)
    
    ## Step 1: Obtain data
    """
    Step 1 Notes:
     - Will have to calculate the evolution rates in this step for some
       species.
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
        print(filename)
        if filename.endswith(".csv") and filename[0:2] != "NO":
            print(filename)
            data[str(i)] = pd.read_csv(filename)
            i += 1
    num_datasets = i


    ## Step 2: Combine data
    """
    Step 2 Notes:
     - All datasets must have a column named "Gene Name"
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
    print(len_bef)
    expression_columns = []
    full_data.dropna(inplace = True)
    num_expression_columns = 0
    for i in list(full_data):
        if i != "Gene Name":
            if i != evo_rate_name:
                    expression_columns.append(i)
                    num_expression_columns += 1
            # This section log's it if necessary
            #full_data = full_data[full_data[i] != 0]
            #full_data[i] = np.log(full_data[i])

            # Rank data here
            full_data[i] = full_data[i].rank()
            
    len_aft = len(full_data)
    print("Num rows dropped:", len_bef-len_aft)
    if float(len_bef-len_aft)/float(len_bef) > 0.05:
        print("Warning: around", 100*(len_bef-len_aft)/len_bef, "percent of your data has been dropped.")
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
    # Create all of the individual r^2
    # May edit to make the indiv_r2 a dictionary instead of a list
    expression_data = full_data[expression_columns]
    evo_rates = full_data["Evo Rate"]
    indiv_r2 = []
    for i in expression_columns:
        r = stats.pearsonr(full_data[i], evo_rates)
        #print(i, r)
        indiv_r2.append(r)
    print()

    # Edit the commands later
    print("Commands: next, drop, undrop, full, only, highest, alpha, partial, r2, pr, save, help.")
    flag = True
    dropped = []
    while flag:
        user_inp = input("Command: ")
        in_list = user_inp.replace(',','').split()
        valid_commands = ["drop","alpha", "partial", "full", "r2", "undrop", "only", "next", "pr", "save", "highest", "help"]
        if in_list[0].lower() not in valid_commands:
            print("Invalid command.")
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
                else:
                    print(in_list[1], "is not a valid print statement.")
            elif in_list[0].lower() == "highest":
                cols = highest_coef(int(in_list[1]), cur_model[0], cur_model[1])
                cur_model = lin_model(full_data, cols)
            elif in_list[0].lower() == "help":
                print("Help currently unavailable.")
            elif in_list[0].lower() == "save":
                if in_list[1].lower() == "r2":
                    r2 = pd.DataFrame()
                    r2s = []
                    for i in indiv_r2:
                        r2s.append(i[0])
                    r2["Conditions"] = expression_columns.copy()
                    r2["R2s"] = r2s
                    user_file = input("Name the file (with .csv): ")
                    r2.to_csv(user_file)
            elif in_list[0].lower() == "alpha":
                if in_list[1].lower() == "tune":
                    y_pred = cur_model[2]
                    se = metrics.mean_squared_error(y_pred, evo_rates) * len(y_pred)
                    print("Squared error:", se)
                    if len(in_list) == 2:
                        alph = 0.01*se
                    else:
                        tuner = float(in_list[2])
                        alph = tuner*se
                    cur_model = lin_model(full_data, cur_model[1], alph)
                else:
                    cur_model = lin_model(full_data, cur_model[1], float(in_list[1]))
            elif in_list[0].lower() == "full":
                cols = expression_columns.copy()
                cur_model = lin_model(full_data, cols)
            else:
                # Checks if all inputs are okay
                good_input = True
                for i in range(1, len(in_list)):
                    
                    if not in_list[i] in expression_columns and in_list[i] != "all":
                        # RIGHT HERE: CHECK CAPITALIZATION
                        # Note: whatever gets here won't go through the commands
                        in_list[i] = in_list[i].replace("_", " ")
                        if not in_list[i] in expression_columns and in_list[i] != "all":
                            good_input = False
                # If the input is fine, contiue
                if good_input:
                    # Note: no danger of key errors
                    com = in_list[0].lower()
                    if com == "drop":
                        cols = cur_model[1].copy()
                        for i in range(1, len(in_list)):
                            cols.remove(in_list[i])
                        cur_model = lin_model(full_data, cols)
                    elif com == "r2":
                        for i in range(1, len(in_list)):
                            pr_r2(indiv_r2, full_data, in_list[i])
                    elif com == "undrop":
                        cols = cur_model[1].copy()
                        for i in range(1, len(in_list)):
                            cols.append(in_list[i])
                        cur_model = lin_model(full_data, cols)
                    elif com == "only":
                        cur_model = lin_model(full_data, in_list[1:])
                    elif com == "partial":
                        if in_list[1].lower() == "all":
                            partls = pd.DataFrame()
                            pd.set_option("display.max_rows", None, "display.max_columns", None)
                            partls["Accounted for"] = expression_columns.copy()
                            for cur in expression_columns:
                                temp = []
                                for acc in expression_columns:
                                    if acc != cur:
                                        temp.append(partial(full_data, evo_rate_name, cur, acc))
                                    else:
                                        temp.append(0)
                                partls[cur] = temp
                            print(partls)
                            #partls.to_csv('partials.csv')
                        elif in_list[2].lower() == "all":
                            print("Partial correlations for", in_list[1])
                            for i in expression_columns:
                                if i != in_list[1]:
                                    print("   %-23s %1.4f" %(i+":", partial(full_data, evo_rate_name, in_list[1], i)))
                        else:
                            print("Partial correlation:", partial(full_data, evo_rate_name, in_list[1], in_list[2]))
                    else:
                        print("Command", com, "is unavailable.")
                else:
                    print("Invalid inputs.")
        print()
    ## Step 5: Decide what to do from here
    print()
    print()
print("Program is over.")
input("Press enter to exit. ")
