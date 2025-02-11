
import os
import pandas as pd
import sys
from Bio import SeqIO
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch


#LO_SCRIPT il più doloroso, il primo. Fa decisamente troppe cose: Ricerca nella directory specificata tutti i file .emapper.annotations Estrae e aggrega dati di interesse dalle colonne Preferred_name, query, # score, filtra solo i record con score > 100. Crea un DataFrame con i geni presenti nei vari file.
#Genera un file CSV aggregato (table.csv) che contiene una tabella con i geni identificati e le sequenze proteiche associate. Esclude valori mancanti e ordina alfabeticamente i dati
#Estrae sequenze proteiche dai file .faa. Per ogni gene trovato, cerca nel CSV i corrispondenti file .faa. Se il file .faa esiste, estrae l'ID della sequenza proteica corrispondente. Salva la sequenza in un file di output nella cartella ./out/<nome_gene>
#Classifica i geni in categorie biologiche
#Genera un grafico a torta della distribuzione genica
#Fornisce informazioni dettagliate sui geni. Cerca un gene specifico e indica la sua categoria di appartenenza. Se presente nei file .emapper.annotations, estrae e stampa la sua descrizione o EC number
#Input richiesti: input_directory (directory contenente i file .emapper.annotations), gene_column_name (opzionale, nome di una colonna di interesse nel CSV)
#Output generati: table.csv: file CSV con le informazioni aggregate, ./out/grouped_genes.txt: lista dei geni divisi per categoria, ./out/gene_info.csv: dettagli sui geni estratti dai file .emapper.annotations,
#./out/gene_distribution_piechart.png: grafico della distribuzione genica.

#THE_SCRIPT, the most painful one, the first. It does far too many things:
#Searches the specified directory for all .emapper.annotations files.
#Extracts and aggregates relevant data from the Preferred_name, query, and score columns, filtering only records with score > 100.
#Creates a DataFrame with the genes present in the various files.
#Generates an aggregated CSV file (table.csv), which contains a table with the identified genes and associated protein sequences. It excludes missing values and sorts data alphabetically.
#Extracts protein sequences from .faa files.
#For each identified gene, it looks for the corresponding .faa file in the CSV.
#If the .faa file exists, it extracts the corresponding protein sequence ID.
#Saves the sequence in an output file within ./out/<gene_name>.
#Classifies genes into biological categories.
#Generates a pie chart of gene distribution.
#Provides detailed information about genes.
#Searches for a specific gene and indicates its category.
#If present in the .emapper.annotations files, it extracts and prints its description or EC number.
#Required Input:
#input_directory (directory containing .emapper.annotations files)
#gene_column_name (optional, name of a column of interest in the CSV)
#Generated Output:
#table.csv: aggregated CSV file with extracted gene data.
#./out/grouped_genes.txt: list of genes categorized by type.
#./out/gene_info.csv: detailed gene information extracted from .emapper.annotations files.
#./out/gene_distribution_piechart.png: pie chart of gene distribution.

def cerca_colonna_csv(file_csv, nome_colonna):
    """Carica il CSV e verifica se la colonna esiste."""
    try:
        df = pd.read_csv(file_csv)
        if nome_colonna not in df.columns:
            raise ValueError(f"Gene '{nome_colonna}' not found.")
            sys.exit(1)
    except FileNotFoundError:
        print(f"Error: file {file_csv} not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error while reading the CSV: {e}")
        sys.exit(1)

def estrai_id_e_sequenza(file_faa, sequenza_da_cercare):
    """Estrae l'ID e la sequenza proteica dal file .faa, cercando il nome della sequenza."""
    try:
        for record in SeqIO.parse(file_faa, "fasta"):
            if sequenza_da_cercare == record.id:
                return record.id, str(record.seq)
    except FileNotFoundError:
        print(f"Error: file {file_faa} not found.")
        return None, None
    return None, None

def crea_file_output(sequenza, nome_file_output, sequenza_id, sequenza_proteica):
    """Crea un file di output con l'ID trovato e la sequenza proteica."""
    os.makedirs(os.path.dirname(nome_file_output), exist_ok=True)  # Crea la directory se non esiste
    with open(nome_file_output, 'w') as f:
        f.write(f">{sequenza_id}\n")
        for i in range(0, len(sequenza_proteica), 60):
            f.write(f"{sequenza_proteica[i:i+60]}\n")

def analizza_csv(file_csv, nome_colonna, input_dir, output_dir):
    """Funzione principale per analizzare il CSV e creare i file di output."""
    # Carica il CSV e verifica se la colonna esiste
    df = pd.read_csv(file_csv)
    if nome_colonna not in df.columns:
        print(f"Error: Gene column '{nome_colonna}' not found in the CSV. No output will be generated.")
        return  # Se la colonna non esiste, esce dalla funzione senza fare nulla

    # Se la colonna esiste, crea la cartella di output
    os.makedirs(output_dir, exist_ok=True)  # Crea la directory di output se non esiste

    # Cicla su tutte le righe e per ogni valore non NaN nella colonna
    for index, row in df.iterrows():
        sequenze = str(row[nome_colonna]).split(',')  # Separa le sequenze in caso siano divise da virgole
        
        for sequenza in sequenze:
            sequenza = sequenza.strip()  # Rimuove eventuali spazi prima o dopo la sequenza
            if pd.notna(sequenza) and sequenza != "-":
                
                # Ottieni il nome del file .faa dall'indice della riga (è il nome completo del file)
                nome_file_faa = row.iloc[0]  # Usa l'indice del dataframe come nome del file
                
                # Costruisci il percorso completo del file .faa aggiungendo l'indice alla directory
                file_faa = os.path.join(input_dir, nome_file_faa)  # La path completa con nome del file
                
                # Verifica che il file .faa esista
                if os.path.exists(file_faa):
                    sequenza_id, sequenza_proteica = estrai_id_e_sequenza(file_faa, sequenza)
                    if sequenza_id:
                        # Crea il nome del file di output
                        nome_file_output = os.path.join(output_dir, f"{sequenza}_{nome_file_faa}")
                        
                        # Crea la directory di output solo se la sequenza è stata trovata
                        os.makedirs(os.path.dirname(nome_file_output), exist_ok=True)
                        
                        # Crea il file di output solo se una sequenza è stata trovata
                        crea_file_output(sequenza, nome_file_output, sequenza_id, sequenza_proteica)
                        print(f"Created file: {nome_file_output}")
                    else:
                        print(f"Gene '{sequenza}' not found in file {nome_file_faa}, skipping...")


def leggi_emapper_annotations(input_dir):
    """Legge i file .emapper.annotations e crea un DataFrame."""
    all_preferred_names = set()
    file_data = {}

    for filename in os.listdir(input_dir):
        if filename.endswith("emapper.annotations"):
            df = pd.read_csv(os.path.join(input_dir, filename), sep="\t", header=0, skiprows=4)

            if 'Preferred_name' not in df.columns or '#query' not in df.columns or 'score' not in df.columns:
                raise ValueError(f"Column 'Preferred_name', 'query' or 'score' not found in: {filename}")
            
            # Filtra e aggrega i dati
            df_filtered = df[(df['Preferred_name'] != "-") & (df['score'] > 100)]
            df_grouped = df_filtered[['Preferred_name', '#query']].drop_duplicates().groupby('Preferred_name')['#query'].apply(list).reset_index()
            
            all_preferred_names.update(df_grouped['Preferred_name'])
            file_data[filename] = df_grouped.set_index('Preferred_name')['#query'].to_dict()

    return sorted(all_preferred_names), file_data

def crea_csv(file_data, all_preferred_names, output_file):
    """Crea il CSV finale con le informazioni aggregate se il file non esiste già."""
    if os.path.exists(output_file):
        print(f"CSV file '{output_file}' already exists. Skipping CSV creation.")
        return pd.read_csv(output_file, index_col=0)  # Leggi il CSV esistente senza rifarlo

    df_final = pd.DataFrame(index=file_data.keys(), columns=all_preferred_names)

    for file_name, preferred_name_queries in file_data.items():
        for preferred_name, queries in preferred_name_queries.items():
            df_final.at[file_name, preferred_name] = ', '.join(map(str, queries))

    df_final = df_final.fillna("-")
    df_final = df_final.sort_index(axis=0).sort_index(axis=1)
    df_final.columns = df_final.columns.str.strip()
    
    # Sostituisce l'estensione .emapper.annotations con .fna.faa
    df_final.index = df_final.index.str.replace('.emapper.annotations', '.fna.faa', regex=False)
    
    df_final.to_csv(output_file, index=True)
    print(f"Created CSV file: {output_file}")
    
    return df_final
    

def distribuzione_gene_pie(df_final, file_data, input_dir, output_grouped_genes_file='./out/grouped_genes.txt', output_gene_info_file='./out/gene_info.csv'):
    """
    Conta la percentuale di file in cui ogni gene (colonna) è presente (escludendo "-") 
    e crea un grafico a torta (pie chart) della distribuzione in 4 categorie.
    """
    # Verifica se il grafico esiste già
    pie_chart_file = './out/gene_distribution_piechart.png'
    if os.path.exists(pie_chart_file):
        print(f"Pie chart '{pie_chart_file}' already exists. Skipping pie chart creation.")
        return  # Esci senza rifare il grafico

    # Calcola la percentuale di file (righe) in cui ogni gene (colonna) è presente (escludendo "-")
    max_righe = len(df_final)
    counts = df_final.apply(lambda x: (x != "-").sum(), axis=0)  # Conta le righe non vuote per ciascun gene
    percentuali = (counts / max_righe) * 100  # Calcola la percentuale per ciascun gene

    # Inizializza i gruppi di geni per categoria
    categorie = {"core": [], "soft-core": [], "shell": [], "cloud": []}
    gene_info = {}  # Dizionario per memorizzare informazioni uniche sui geni (evita duplicati)

    # Raggruppa i geni nelle 4 categorie in base alla loro percentuale di presenza nei file
    for gene, percentuale in percentuali.items():
        if percentuale >= 90:
            categorie["core"].append(gene)
        elif 80 <= percentuale < 90:
            categorie["soft-core"].append(gene)
        elif 10 <= percentuale < 80:
            categorie["shell"].append(gene)
        else:
            categorie["cloud"].append(gene)

        # Raccogli informazioni per ogni gene presente in file .emapper.annotations
        for filename, preferred_name_queries in file_data.items():
            if gene in preferred_name_queries:
                # Leggi il file .emapper.annotations
                file_path = os.path.join(input_dir, filename)
                df = pd.read_csv(file_path, sep="\t", header=0, skiprows=4)

                # Se il gene è trovato, prendi la riga e la salva (evita duplicati)
                gene_row = df[df['Preferred_name'] == gene].iloc[0]
                # Salva solo la prima riga incontrata (evita duplicati di gene)
                if gene not in gene_info:
                    gene_info[gene] = gene_row

    # Calcola il numero di geni per categoria e la percentuale rispetto al totale
    total_genes = len(df_final.columns)
    distribuzione_percentuale = {
        key: (len(genes) / total_genes) * 100
        for key, genes in categorie.items()
    }

    # Crea il grafico a torta
    plt.figure(figsize=(8, 8))
    pie_values = [len(categorie[key]) for key in categorie]
    pie_labels = [f"{key}: {len(categorie[key])} Genes ({distribuzione_percentuale[key]:.2f}%)" for key in categorie]
    plt.pie(pie_values, labels=pie_labels, autopct='%1.1f%%', startangle=90, colors=plt.cm.Paired.colors)

    # Aggiungi la leggenda con tutte le categorie
    legend_labels = [
        f"{key}: {len(categorie[key])} Genes ({distribuzione_percentuale[key]:.2f}%)"
        for key in categorie
    ]
    legend_patches = [Patch(color=plt.cm.Paired.colors[i], label=legend_labels[i]) for i in range(len(categorie))]
    plt.legend(legend_patches, legend_labels, title="Gene Categories", bbox_to_anchor=(1.05, 1), loc='best', fontsize=10)

    # Aggiungi la seconda legenda con gli intervalli di percentuale delle categorie
    interval_legend = {
        "core": "90% or more",
        "soft-core": "Between 80% and 90%",
        "shell": "Between 10% and 80%",
        "cloud": "Less than 10%"
    }
    interval_legend_labels = [
        f"{category}: {interval_legend[category]}" for category in categorie
    ]
    interval_legend_patches = [
        Patch(color=plt.cm.Paired.colors[i], label=interval_legend_labels[i]) for i in range(len(categorie))
    ]
    plt.legend(
        interval_legend_patches, interval_legend_labels, title="Category Ranges", 
        bbox_to_anchor=(1.05, 0), loc='best', fontsize=10
    )

    plt.title('Genes Distribution across Categories Based on File Presence', fontsize=16)
    plt.tight_layout()
    plt.savefig(pie_chart_file)
    print(f"Graphic saved in {pie_chart_file}")

    # Se i file di output già esistono, non crearli nuovamente
    if os.path.exists(output_gene_info_file):
        print(f"Gene info CSV already exists at {output_gene_info_file}. Skipping CSV creation.")
    else:
        # Scrivi i nomi dei geni raggruppati in un file
        with open(output_grouped_genes_file, 'w') as f:
            for category, genes in categorie.items():
                f.write(f"{category.capitalize()}: {', '.join(genes)}\n")
        print(f"Gene groups have been written to: {output_grouped_genes_file}")

        # Salva le informazioni dettagliate dei geni in un file CSV
        gene_info_df = pd.DataFrame.from_dict(gene_info, orient='index')
        gene_info_df.to_csv(output_gene_info_file)
        print(f"Gene detailed info saved in: {output_gene_info_file}")

def stampa_info_gene(df_final, gene, file_data, input_dir):
    """
    Stampa la categoria in cui il gene specificato è presente (core, soft-core, shell, cloud),
    la percentuale di presenza e il contenuto della colonna 'Definition' o 'EC' del primo file che contiene il gene.
    La colonna 'Definition' ha priorità sulla colonna 'EC' se entrambe sono presenti in uno o più file.
    """
    # Calcola la percentuale di presenza del gene (escludendo "-")
    max_righe = len(df_final)
    counts = df_final.apply(lambda x: (x != "-").sum(), axis=0)  # Conta le righe non vuote per ciascun gene
    percentuali = (counts / max_righe) * 100  # Calcola la percentuale per ciascun gene

    # Verifica se il gene è presente nella colonna
    if gene in df_final.columns:
        percentuale = percentuali[gene]

        # Determina la categoria del gene in base alla percentuale
        if percentuale >= 90:
            categoria = "core"
        elif 80 <= percentuale < 90:
            categoria = "soft-core"
        elif 10 <= percentuale < 80:
            categoria = "shell"
        else:
            categoria = "cloud"
        
        # Stampa la categoria del gene
        print(f"Gene '{gene}' is categorized as: {categoria} ({percentuale:.2f}% presence in the files)")

        # Flag per verificare se il gene è stato trovato in almeno un file con una colonna 'Definition' o 'EC'
        found_info = False

        # Cerca il gene in tutti i file
        for filename, preferred_name_queries in file_data.items():
            if gene in preferred_name_queries:
                # Leggi il file .emapper.annotations
                file_path = os.path.join(input_dir, filename)
                df = pd.read_csv(file_path, sep="\t", header=0, skiprows=4)

                # Controlla se la colonna 'Description' esiste
                if 'Description' in df.columns:
                    Description = df[df['Preferred_name'] == gene]['Description'].values
                    if Description.size > 0:
                        print(f"Gene '{gene}' description: {Description[0]}")
                        found_info = True
                        break  # Esce dal ciclo se trova la definizione

                # Se 'Definition' non esiste, controlla la colonna 'EC'
                elif 'EC' in df.columns:
                    ec = df[df['Preferred_name'] == gene]['EC'].values
                    if ec:
                        print(f"Gene '{gene}' EC: {ec[0]}")
                        found_info = True
                        break  # Esce dal ciclo se trova l'EC
                
        if not found_info:
            print(f"Gene '{gene}' does not have 'Description' or 'EC' column in any of the files.")
    else:
        print(f"Gene '{gene}' not found in the data.")


def main():
    # Verifica che il numero di argomenti sia corretto
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Error: Expected a directory and a gene name (optional).")
        print("Usage: python script.py <input_directory> <gene_column_name (optional)>")
        sys.exit(1)

    input_dir = sys.argv[1]  # Directory contenente i file .emapper.annotations
    
    # Usa il secondo argomento se presente, altrimenti imposta un valore predefinito
    nome_colonna = sys.argv[2] if len(sys.argv) == 3 else None  # Se non fornito, imposta a None
    
    output_file = './table.csv'  # Nome del file CSV di output



    # Solo se 'nome_colonna' è specificato, crea la directory di output
    if nome_colonna:
        # Leggi il CSV prima di creare la cartella
        df = pd.read_csv(output_file)
        if nome_colonna not in df.columns:
            print(f"Error: Gene column '{nome_colonna}' not found in the CSV. Skipping output creation.")
            sys.exit(1)  # Se la colonna non esiste, esci subito

        output_dir = os.path.join('./out', nome_colonna)  # Directory di output per il gene
        os.makedirs(output_dir, exist_ok=True)  # Crea la directory se non esiste
        print(f"Processing for gene column: {nome_colonna}")
    else:
        output_dir = None  # Non crea una directory se non c'è nome_colonna
        print("No gene column specified, skipping CSV analysis.")

    # Leggi e analizza i file .emapper.annotations
    all_preferred_names, file_data = leggi_emapper_annotations(input_dir)

    # Crea il CSV con i dati aggregati
    df_final = crea_csv(file_data, all_preferred_names, output_file)
    
    # Se il nome della colonna è stato specificato, esegui l'analisi del CSV
    if nome_colonna:
        # Analizza il CSV e crea i file di output per ogni sequenza
        analizza_csv(output_file, nome_colonna, input_dir, output_dir)
        stampa_info_gene(df_final, nome_colonna, file_data, input_dir)

    # Crea la distribuzione dei geni in base alla loro presenza nei file
    distribuzione_gene_pie(df_final, file_data, input_dir)



if __name__ == "__main__":
    main()

