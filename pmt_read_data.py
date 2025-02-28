import pandas as pd
import numpy as np

def pmt_read_data(file_name):
    # Load the data
    data = pd.read_csv(file_name, delim_whitespace=True)
    
    # Extract the header (first 4 rows)
    header = data.iloc[:4]
    
    # Extract the main data (from row 4 onward)
    data_clean = data.iloc[4:].reset_index(drop=True)
    
    # Extract gate time from the header (convert to float)
    gate_time = float(header.iloc[0, 0])
    
    # Rename the data column
    data_clean.columns = ['counts']
    
    # Compute time array
    index_list = np.arange(len(data_clean))  # Instead of using `index.tolist()`
    time = index_list * gate_time / 1000  # Convert to seconds
    
    return time, data_clean['counts']