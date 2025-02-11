import os
import pandas as pd
import sys

#SCRIPT che gira ricorsivamente una directory, prende i file .emapper.annotations e filtra attraverso alcune colonne di interesse. Dopodichè processa i dati rimanenti in un csv apposito.
#SCRIPT that recursively traverses a directory, retrieves .emapper.annotations files, and filters them based on specific columns of interest. Then, it processes the remaining data into a dedicated CSV file.


# Controlla se i parametri sono stati passati
if len(sys.argv) != 3:
    print("Errore: Devi passare il percorso del file di input e il percorso del file di output.")
    sys.exit(1)

# Leggi i parametri dalla riga di comando
input_file = sys.argv[1]  # File di input
output_file = sys.argv[2]  # File di output

df_final = pd.DataFrame()

# This will store the final DataFrame, where each row will be a file and columns are Preferred_name values
all_preferred_names = set()  # To store all unique 'Preferred_name' values

# First loop: Read all files and collect Preferred_names
file_data = {}

for i in os.listdir(input_file):
    if i.endswith("emapper.annotations"):
        # Read the file
        df = pd.read_csv(input_file + i, sep="\t", header=0, skiprows=4)
        
        # Check if necessary columns exist
        if 'Preferred_name' not in df.columns or '#query' not in df.columns or 'score' not in df.columns:
            raise ValueError(f"Le colonne 'Preferred_name' o 'query' o 'score' non sono presenti nel file: {file_path + i}")
        
        # Filter the dataframe: remove rows with gene "-" and keep only score >= 100
        df_sel1 = df[(df['Preferred_name'] != "-") & (df['score'] > 100)]
        
        # Select only relevant columns and remove duplicates
        df_sel2 = df_sel1[['Preferred_name', '#query']].drop_duplicates()
        
        # Group by 'Preferred_name' and aggregate queries into a list
        df_raggruppato = df_sel2.groupby('Preferred_name')['#query'].apply(list).reset_index()
        
        # Add the preferred names of this file to the overall set
        all_preferred_names.update(df_raggruppato['Preferred_name'])
        
        # Store the grouped data by file name
        file_data[i] = df_raggruppato.set_index('Preferred_name')['#query'].to_dict()

# Convert the set of all preferred names into a sorted list
all_preferred_names = sorted(list(all_preferred_names))

# Create an empty DataFrame with files as rows and Preferred_names as columns
df_final = pd.DataFrame(index=file_data.keys(), columns=all_preferred_names)

# Fill the DataFrame with the appropriate #query values
for file_name, preferred_name_queries in file_data.items():
    for preferred_name, queries in preferred_name_queries.items():
            df_final.at[file_name, preferred_name] = ', '.join(map(str, queries))

df_final = df_final.fillna("NaN")
df_final = df_final.sort_index(axis=0).sort_index(axis=1)
df_final.columns = df_final.columns.str.strip()
df_final.index = df_final.index.str.replace('.emapper.annotations', '.fna.faa')

df_final.to_csv(output_file, index=True)

print(df_final)


# In[ ]:




