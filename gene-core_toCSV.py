import pandas as pd

#SCRIPT per confrontare la colonna KEGG_ko di un file csv con un file di testo con i geni separati da una virgola. Ideato per estrarre i gene-core dal dataset e creare un'csv di output.
#SCRIPT to compare the KEGG_ko column of a CSV file with a text file containing genes separated by commas. Designed to extract core genes from the dataset and generate an output CSV.

# Funzione per leggere i geni dal file di testo (separati da virgola)
def get_genes_from_txt(file_txt):
    with open(file_txt, 'r') as f:
        # Leggi tutto il file e suddividi per virgola, rimuovendo eventuali spazi extra
        genes = f.read().strip().split(',')
    return [gene.strip() for gene in genes]

# Funzione principale
def extract_genes_from_csv(csv_file, txt_file, output_csv):
    # Carica il CSV
    df = pd.read_csv(csv_file)
    
    # Rimuovi eventuali spazi extra nei nomi delle colonne
    df.columns = df.columns.str.strip()

    # Mostra i nomi delle colonne per verificare il nome corretto della colonna
    print("Colonne del CSV:", df.columns)

    # Ottieni la lista dei geni dal file di testo
    genes_from_txt = get_genes_from_txt(txt_file)
    
    # Supponiamo che la prima colonna contenga i geni
    first_column = df.iloc[:, 0]  # Prima colonna del CSV
    
    # Verifica se esiste la colonna 'KEGG_ko'
    if 'KEGG_ko' not in df.columns:
        raise ValueError("La colonna 'KEGG_ko' non è presente nel file CSV.")

    # Filtra il CSV per prendere solo le righe in cui i geni nella prima colonna corrispondono a quelli nel file di testo
    core_df = df[first_column.isin(genes_from_txt)]
    
    # Estrai la prima colonna (i geni) e la colonna 'KEGG_ko'
    core_df_filtered = core_df[[first_column.name, 'KEGG_ko']]  # Usa il nome della prima colonna e 'KEGG_ko'
    
    # Scrivi il risultato in un nuovo file CSV
    core_df_filtered.to_csv(output_csv, index=False)

    print(f"File CSV filtrato salvato come {output_csv}")

# Esegui la funzione
csv_file = ""  # Percorso del tuo CSV
txt_file = ""  # Percorso del file di testo con i geni
output_csv = ""  # Percorso del file di output

extract_genes_from_csv(csv_file, txt_file, output_csv)

