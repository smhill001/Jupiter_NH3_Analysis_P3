def concat_csv(pathin,featuretype='fNH3'):
    import pandas as pd
    import glob
    import os
    
    # Find all CSV files in the directory
    files=os.listdir(pathin)
    filescsv = [item for item in files if ".csv" in item]
    filesblob = [item for item in filescsv if "blob" in item]
    filesfeature = [item for item in filesblob if featuretype in item]
    
    print(filesfeature)    
    # Read each file into a list of DataFrames and concatenate
    df = pd.concat([pd.read_csv(pathin+f) for f in filesfeature], ignore_index=True)
    
    # Save to a new CSV
    df.to_csv(pathin+featuretype+"combined.csv", index=False)