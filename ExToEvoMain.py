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

    # THIS CODE IS FOR PLOTS: COMMENT AND UNCOMMENT AS NECESSARY
    #lin_comb = np.multiply(coef[0], expression_data.iloc[:, 0])
    #for i in range(1, len(columns)):
        #lin_comb = np.add(lin_comb, np.multiply(coef[i], expression_data.iloc[:, i]))
    
    #plt.scatter(expression_data[columns[0]], evo_rates, color = 'red', marker = '.')
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
# List of all species to go through -- change as necessary
species_names = ["Escherichia_coli"]

# Get the full path before jumping into species
cur_path = os.getcwd()

# Go through all species
for species_name in species_names:
    print()
    print("Current species:", species_name)
    
    ## Step 1: Obtain data
    
    # Create the full path
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
        if filename.endswith(".csv") and filename[0:2] != "NO":
            print(filename)
            data[str(i)] = pd.read_csv(filename)
            i += 1
    num_datasets = i
    print()

    ## Step 2: Combine data
    
    # Merge all datasets by their column named "Gene Name"
    full_data = data["0"]
    for i in range(1, num_datasets):
        full_data = pd.merge(full_data, data[str(i)], how = "outer", on = "Gene Name")

    if evo_rate_name not in full_data.columns:
        print("This program is about to crash because there is no column titled", evo_rate_name)
    else:
        colsies = [evo_rate_name] + [col for col in full_data if col != evo_rate_name]
        full_data = full_data[colsies]
        
    ## Step 3.1: Process data

    ## IF AVERAGING, SAY TRUE
    av = False

    
    len_bef = len(full_data)    # Keep track of how much was dropped
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

            # Rank data here if necessary
            if not av:
                full_data[i] = full_data[i].rank()
    full_data.dropna(inplace = True)


    # AVERAGING DATA FOR REL606
    if species_name == "Escherichia_coli_REL606" and av:
        conds = {
            "gluc_exp": ["glucose_3", "glucose_3.1", "glucose_3.2",
                       "glucose_4", "glucose_4.1", "glucose_4.2",
                       "glucose_5", "glucose_5.1", "glucose_5.2",
                       "glucose_6", "glucose_6.1", "glucose_6.2",
                       "glucose_8", "glucose_8.1", "glucose_8.2"],
            "gluc_sta": ["glucose_24", "glucose_24.1", "glucose_24.2",
                         "glucose_48", "glucose_48.1", "glucose_48.2"],
            "gluc_latesta": ["glucose_168", "glucose_168.1", "glucose_168.2",
                             "glucose_336", "glucose_336.1", "glucose_336.2"],
            "glyc_exp": ["glycerol_5", "glycerol_7", "glycerol_8", "glycerol_10", "glycerol_14",
                        "glycerol_5.1", "glycerol_7.1", "glycerol_8.1", "glycerol_10.1", "glycerol_14.1",
                        "glycerol_5.2", "glycerol_7.2", "glycerol_8.2", "glycerol_10.2", "glycerol_14.2"],
            "glyc_sta":["glycerol_24", "glycerol_48", "glycerol_168", "glycerol_336",
                       "glycerol_24.1", "glycerol_48.1", "glycerol_168.1",
                       "glycerol_24.2", "glycerol_48.2", "glycerol_168.2"],
            #"nacl_5_exp": ["glucose_5.5", "glucose_6.2", "glucose_5.5.1"],
            #"nacl_5_sta": ["glucose_29", "glucose_28", "glucose_29.1"],
            "nacl_100_exp": ["nacl_100_6.5"],
            #"nacl_200_exp": ["nacl_200_6", "nacl_200_8", "nacl_200_8.1"],
            "nacl_300_exp": ["nacl_300_8", "nacl_300_10", "nacl_300_10.1"],
            "nacl_100_sta": ["nacl_100_28", "nacl_100_29", "nacl_100_29.1"],
            #"nacl_200_sta": ["nacl_200_28", "nacl_200_29", "nacl_200_29.1"],
            "nacl_300_sta": ["nacl_300_28", "nacl_300_29", "nacl_300_29.1"],
            "lact_exp": ["lactate_9", "lactate_8", "lactate_9.1"],
            "lact_sta": ["lactate_29", "lactate_29.1","lactate_29.2"],
            "gluconate_exp": ["gluconate_6.1", "gluconate_6.1"],
            "gluconate_sta": ["gluconate_27", "gluconate_27.1", "gluconate_27.2"],
            "mgso4_0.08_exp": ["mgso4_0.08_5.5.1", "mgso4_0.08_5.5.2", "mgso4_0.08_5.5.3",
                               "mgso4_0.08_5.5.4", "mgso4_0.08_5.5.5"],
            "mgso4_0.08_sta": ["mgso4_0.08_26", "mgso4_0.08_26.1", "mgso4_0.08_26.2", "mgso4_0.08_28",
                               "mgso4_0.08_28.1", "mgso4_0.08_28.2"],
            #"mgso4_0.8_exp": ["glucose_5.5.2", "glucose_5.5.3", "glucose_5.5.4"],
            #"mgso4_0.8_sta": ["glucose_28.3", "glucose_28.2", "glucose_28.1"],
            "mgso4_8_exp": ["mgso4_8_5.5", "mgso4_8_5.5.1", "mgso4_8_5.5.2"],
            #"mgso4_50_exp": ["mgso4_50_5", "mgso4_50_5.1", "mgso4_50_5.2"],
            #"mgso4_200_exp": ["mgso4_200_5", "mgso4_200_5.1", "mgso4_200_5.2"],
            "mgso4_400_exp": ["mgso4_400_7", "mgso4_400_7.1", "mgso4_400_7.2"],
            "mgso4_8_sta": ["mgso4_8_28", "mgso4_8_28.1", "mgso4_8_28.2"],
            #"mgso4_50_sta": ["mgso4_50_28", "mgso4_50_28.1", "mgso4_50_28.2"],
            #"mgso4_200_sta": ["mgso4_200_28", "mgso4_200_28.1", "mgso4_200_28.2"],
            "mgso4_400_sta": ["mgso4_400_28", "mgso4_400_28.1", "mgso4_400_28.2"],
            "mgso4_0.005_exp": ["mgso4_0.005_5.1", "mgso4_0.005_5.2"],
            #"mgso4_0.01_exp": ["mgso4_0.01_5.5", "mgso4_0.01_5.5.1", "mgso4_0.01_5.5.2"],
            #"mgso4_0.02_exp": ["mgso4_0.02_5.5", "mgso4_0.02_5.5.1", "mgso4_0.02_5.5.2"],
            #"mgso4_0.04_exp": ["mgso4_0.04_5.5", "mgso4_0.04_5.5.1", "mgso4_0.04_5.5.2"],
            "mgso4_0.005_sta": ["mgso4_0.005_26", "mgso4_0.005_26.1", "mgso4_0.005_26.2"],
            #"mgso4_0.01_sta": ["mgso4_0.01_26", "mgso4_0.01_26.1", "mgso4_0.01_26.2"],
            #"mgso4_0.02_sta": ["mgso4_0.02_26", "mgso4_0.02_26.1", "mgso4_0.02_26.2"],
            #"mgso4_0.04_sta": ["mgso4_0.04_26", "mgso4_0.04_26.1", "mgso4_0.04_26.2"]
            }
        # Removed:
        # gluconate_6 (MURI_091); mgso4_0.08_5.5 (MURI_130)
        # nacl_300_10 (MURI_072); mgso4_0.005_5 (MURI_142)
        
        avg_exp = pd.DataFrame()
        for k in conds.keys():
            avg_exp[k] = full_data[conds[k]].mean(axis = 1)
        expression_columns = list(avg_exp.columns)
        print(avg_exp.head(8))
        avg_exp["Gene Name"] = full_data["Gene Name"]
        avg_exp["Evo Rate"] = full_data["Evo Rate"]
        full_data = avg_exp
        # Reorder
        colsies = [evo_rate_name, "Gene Name"] + [col for col in full_data if col != evo_rate_name and col != "Gene Name"]
        full_data = full_data[colsies]

        # Comment out if ranking is unecessary
        for i in full_data:
            full_data[i] = full_data[i].rank()
    

    len_aft = len(full_data)
    print("Num rows dropped:", len_bef-len_aft)
    if float(len_bef-len_aft)/float(len_bef) > 0.05:
        print("Warning: around", 100*(len_bef-len_aft)/len_bef, "percent of your data has been dropped.")
        print()
    else:
        print("Dropped rows within acceptable parameters.")
        print()
    print("Processed data: ")
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
        r = (r[0]*r[0], r[1])
        indiv_r2.append(r)
    print()

    # Edit the commands later
    print("Commands: next, drop, undrop, full, only, highest, alpha, partial, r2, pr, save, help.")
    flag = True
    dropped = []
    while flag:
        user_inp = input("Command: ")
        in_list = user_inp.replace(',','').split()
        valid_commands = ["drop","alpha", "correlation_to", "compare_to", "partial", "full", "r2", "undrop", "only", "next", "pr", "save", "highest", "help"]
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
                        in_list[i] = in_list[i].replace("_", " ")
                        if not in_list[i] in expression_columns and in_list[i] != "all":
                            good_input = False
                # If the input is fine, contiue
                if good_input:
                    # Note: no danger of key errors
                    com = in_list[0].lower()
                    if com == "drop":
                        cols = cur_model[1].copy()
                        print(cols)
                        for i in range(1, len(in_list)):
                            cols.remove(in_list[i])
                        cur_model = lin_model(full_data, cols)
                    elif com == "compare_to":
                        comp_cond = in_list[1]
                        # dist = np.linalg.norm(a-b)
                        toret = pd.DataFrame()
                        temp = []
                        for i in expression_columns:
                            temp.append(np.linalg.norm(full_data[comp_cond]-full_data[i]))
                        user_file = input("Name of file (with .csv please): ")
                        toret["Conditions"] = expression_columns
                        name = "Similarity to " + comp_cond
                        toret[name] = temp
                        toret.to_csv(user_file)
                    elif com == "correlation_to":
                        comp_cond = in_list[1]
                        # r = stats.pearsonr(full_data[i], evo_rates)
                        toret = pd.DataFrame()
                        temp = []
                        for i in expression_columns:
                            temp.append(stats.pearsonr(full_data[comp_cond], full_data[i])[0])
                        user_file = input("Name of file (with .csv please): ")
                        toret["Conditions"] = expression_columns
                        name = "Correlation to " + comp_cond
                        toret[name] = temp
                        toret.to_csv(user_file)
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
                            user_name = input("Enter the name of the file (with .csv): ")
                            partls.to_csv(user_name)
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
    print()
    print()
print("Program is over.")
input("Press enter to exit. ")
