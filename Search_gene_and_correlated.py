
import os
import pandas as pd
from pathlib import Path

#SCRIPT: cerca ricorsivamente un gene all'interno del dataset per poi correlarlo con geni che contengono *stringa. In output avremo una tabella presenza/assenza csv. Con tutti i file che contengono il gene
#cercato e tutti i geni correlati come colonne.
#SCRIPT: Recursively search for a gene within the dataset and correlate it with genes containing *string. The output will be a presence/absence CSV table, listing all files that contain the searched gene and #all correlated genes as columns.


# Funzione per cercare "New Delhi metallo-beta-lactamase" nel file GFF
def find_NDM1_in_gff(gff_file):
    found = False
    other_metallo = {}  # Dizionario per raccogliere i valori product metallo-beta-lactamase senza NDM-1
    with open(gff_file, 'r') as file:
        for line in file:
            if not line.startswith('#'):  # Salta le righe di commento
                columns = line.strip().split('\t')
                if len(columns) >= 9:  # COLONNE NEL GFF
                    feature_type = columns[2]  # Colonna C (indice 2)
                    if feature_type == "CDS":  # Se è un CDS
                        attributes = columns[8]  # Colonna I (indice 8)
                        # Dividi gli attributi usando il punto e virgola
                        attributes_list = attributes.split(';')
                        
                        # Controlla se uno degli attributi è 'product=New Delhi metallo-beta-lactamase'
                        for attr in attributes_list:
                            if "product=" in attr:
                                product_value = attr.split('=')[1]  # Ottieni il valore del product
                                if "metallo-beta-lactamase NDM-1" in product_value:  # MODIFICA(COSA CERCA)
                                    found = True
                                elif "beta-lactamase" in product_value and "NDM-1" not in product_value: #MODIFICA(COSA CORRELA)
                                    if product_value not in other_metallo:
                                        other_metallo[product_value] = 0  # Inizializza a 0
    return found, other_metallo

# Funzione per eseguire il processo nella directory (ricorsivamente)
def search_genes_in_directory(directory):
    results = []
    file_count = 0  # Contatore dei file esaminati
    product_columns = {}  # Dizionario per memorizzare le colonne di prodotto

    # Usa os.walk() per esplorare tutte le cartelle e i file nella directory
    for root, dirs, files in os.walk(directory):
        for filename in files:
            if filename.endswith(".gff"):  # Verifica che sia un file GFF
                file_count += 1
                file_path = os.path.join(root, filename)  # Percorso completo del file
                file_name = Path(file_path)

                # Trova NDM-1 e altri metallo-beta-lactamase
                found, other_metallo = find_NDM1_in_gff(file_path)

                # Creazione della riga per questo file con 0 iniziali
                file_result = {'File': file_name.parent.name, 'NDM-1': 0}   #MODIFICA

                # Aggiungi le colonne per i metallo-beta-lactamase trovati
                for product_value in other_metallo:
                    # Se non abbiamo ancora una colonna per questo prodotto, creiamola
                    if product_value not in product_columns:
                        product_columns[product_value] = 0  # Aggiungi la colonna inizialmente a 0
                    file_result[product_value] = 0  # Imposta inizialmente 0 per ogni prodotto trovato

                # Se NDM-1 è trovato, aggiorniamo la colonna 'NDM-1' con 1
                if found:
                    file_result['NDM-1'] = 1    #MODIFICA 

                # Se altri prodotti sono trovati, metti 1 nella loro colonna
                for product_value in other_metallo:
                    if product_value in file_result:
                        file_result[product_value] = 1

                # Aggiungi la riga alla lista dei risultati
                results.append(file_result)

                # Aggiungi una stampa ogni 100 file analizzati
                if file_count % 200 == 0:
                    print(f"{file_count} file checked...")

    # Creazione del DataFrame e salvataggio in un CSV
    # Assicurati di includere tutte le colonne di prodotto trovate
    if results:
        df = pd.DataFrame(results)

        # Aggiungi tutte le colonne di prodotto che non sono ancora presenti (se non trovate nei file)
        for product_value in product_columns:
            if product_value not in df.columns:
                df[product_value] = 0  # Aggiungi colonna con valore 0 per tutti i file

        # Riorganizza le colonne per avere 'File' e 'NDM-1' all'inizio, seguiti dai prodotti
        columns_order = ['File', 'NDM-1'] + list(product_columns.keys())
        df = df[columns_order]

 	# Aggiungi una riga con la somma delle colonne numeriche
        sum_row = df.sum(numeric_only=True)  # Somma delle colonne numeriche
        sum_row['File'] = 'Sum Row'  # Aggiungi una label per la riga delle somme
        df = df.append(sum_row, ignore_index=True)  # Aggiungi la riga delle somme al DataFrame


        # Salva il CSV
        output_path = ""   #INSERISCI L'OUTPUT PATH
        df.to_csv(output_path, index=False)
        print(f"File CSV created: {output_path}")
    else:
        print("no gene found.")

# Esegui la ricerca nella directory specificata
directory_input = ""  #INPUT DIR
search_genes_in_directory(directory_input)

