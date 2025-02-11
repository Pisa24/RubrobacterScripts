import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#SCRIPT: prende in input un file csv formato da sole due colonne: NCBIGenome e Numero di proteine(intercambiabile facilmente con numero di geni). Dopo esegue uno studio per creare un'istogramma
#SCRIPT: Takes as input a CSV file consisting of only two columns: NCBIGenome and Number of proteins (easily interchangeable with the number of genes). Then, it performs an analysis to generate a histogram.


# Caricare il file CSV
file_path = ""   #Inserisci qui il tuo path
df = pd.read_csv(file_path)

# Rinominare le colonne per una gestione più semplice
df.rename(columns={"File": "NCBI_Genome", "Numero di Proteine": "Proteins"}, inplace=True)

# Creare i bin per il numero di proteine a intervalli di 180 partendo da 1000
bins = list(range(900, df["Proteins"].max() + 200, 200))

# Calcolare l'istogramma come percentuale
hist, bin_edges = np.histogram(df["Proteins"], bins=bins)
percent_hist = (hist / len(df)) * 100

# Creare il grafico con etichette più dettagliate sugli intervalli
plt.figure(figsize=(12, 6))
plt.bar(bin_edges[:-1], percent_hist, width=180, align='edge', edgecolor='black')

# Impostare i tick per mostrare tutti gli intervalli di interesse
plt.xticks(bin_edges, rotation=45)

# Etichettare gli assi
plt.xlabel("Number of protein")
plt.ylabel("NCBI Rubrobacter proteomes percentage")
plt.title("Distribution of protein in NCBI Rubrobacter proteomes")

# Mostrare il grafico
plt.show()

