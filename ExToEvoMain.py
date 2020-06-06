# Documentation here

## Step 0: Import necessary libraries
print("Downloading libraires...")
import pandas as pd
import numpy as np
#from keras.layers import Input, Dense
#from keras.models import Model
from sklearn import linear_model
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import os
import sys
print("Done. Beginning program...")

"""
Program notes:
 - May change this to a loop that goes through all folders in
   this directory. The current structure goes through species
   based on the list right below program notes
"""
# List of all species that we have. May change.
species_names = ["Testing"]

# Go through all species
for species_name in species_names:
    print("Current species:", species_name)
    ## Step 1: Obtain data
    """
    Step 1 Notes:
     - Will have to calculate the evolution rates in this step
    """

    # Get the current path and create the full path
    cur_path = os.getcwd()
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
    #print("Unprocesed merged data: ")
    #print(full_data.head(8))
    ## Step 3.1: Process data
    """
    Step 3 Notes:
     - Should I drop rows with null values? (data.dropna() if yes)
     - What do I do with 0 values when logging?
         - Add an epsilon ep = lowest expression that isn't 0
     - How to normalize properly
    """
    expression_columns = []
    # Drop NA
    full_data.dropna(inplace = True)
    # Add epsilon to all 0's, log the result, and normalize it
    # Should epsilon be something lower? Higher?
    # This might be unnecessary
    num_expression_columns = 0
    ss = StandardScaler()
    for i in list(full_data):
        if i != "Gene Name":
            if i != "Evo Rate":
                    expression_columns.append(i)
                    num_expression_columns += 1
            # This entire section may be unecessary
            #full_data.loc[full_data[i] == 0, i] = 0.0001
            full_data[i] = np.log(full_data[i])
            #full_data[i] = np.sqrt(np.sqrt(full_data[i]))
            
            #full_data[i] = ss.fit_transform(full_data[[i]])
    print("Processed data: ")
    print(full_data.head(8))
    # Get the number of expression columns
    #num_inputs = num_expression_columns
       
    

    ## Step 3.2: Section out data
    expression_data = full_data[expression_columns]
    evo_rates = full_data["Evo Rate"]
    # X = dataset[['part1']]
    # y = dataset[['part2']]
    # y = train_data.pop('Gene Data').values gets you that

    ## Step 4: Create a linear regression model
    """
    model_input = Input(shape = (num_inputs, ))
    model_output = Dense(1, activation = 'linear')(model_input)
    model = Model(inputs = model_input, outputs = model_output)

    sgd = keras.optimizers.SGD()
    model.compile(optimizer = sgd, loss = 'mse', metrics = ['mse'])
    model.fit(expression_data, evo_rates, epochs = 30, verbose = 2)
    # Have to figure out how to properly show everything

    weights = None
    for layer in model.layers:
        weights = layer.get_weights()
    # Filename for weights
    weights_filename = species_name + ".h5"
    model.save_weights(weights_filename)
    # Use matplotlib to plot these weights
    print(weights)
    """
    
    reg_model = linear_model.LinearRegression()
    reg_model.fit(expression_data, evo_rates)
    print("R^2: ", reg_model.score(expression_data, evo_rates))
    print("Intercept: ", reg_model.intercept_)
    coef = reg_model.coef_
    print("Coefficients: ", coef)

    # Create the linear combination properly right here
    print(coef[0])
    lin_comb = np.multiply(coef[0], expression_data.iloc[:, 0])
    print(len(expression_data))
    for i in range(1, num_expression_columns):
        lin_comb = np.add(lin_comb, np.multiply(coef[i], expression_data.iloc[:, i]))
    # Figure out how to properly display the coefficients
    input()
    fig, axs = plt.subplots(3)
    axs[0].scatter(lin_comb, evo_rates, color = 'red')
    axs[1].scatter(expression_data.iloc[:, 0], evo_rates)
    axs[2].scatter(expression_data.iloc[:, 1], evo_rates)
    #plt.plot(expression_data, reg_model.predict(expression_data), color = 'blue')
    plt.show()

    input("Happy?")

    ## Step 5: Decide what to do from here




