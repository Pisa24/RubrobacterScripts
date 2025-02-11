import pandas as pd

#SCRIPT: Filtraggio e estrazione di dati dal CSV
#Legge un file .txt contenente una lista di nomi separati da virgole.
#Cerca questi nomi nella quarta colonna di un file .emapper.annotations (CSV separato da tabulazioni).
#Per ogni corrispondenza trovata, estrae le colonne A (indice 0) e D (indice 3) e le salva in un nuovo CSV chiamato output_file.csv con due colonne:Gene_cloud(modificabile),Query (relativa alla corrispondenza)
#Creazione di un nuovo file .faa filtrato. Legge il file CSV filtrato per estrarre una lista di query. Cerca le query nel file .faa (formato FASTA).
#Se trova una corrispondenza tra la query e l'ID della sequenza nel .faa:
#Crea una nuova intestazione e scrive la sequenza FASTA con righe di 60 caratteri nel nuovo file di output .faa.

#SCRIPT: Filtering and extraction of data from the CSV
#Reads a .txt file containing a list of names separated by commas.
#Searches for these names in the fourth column of a .emapper.annotations file (CSV separated by tabulations).
#For each match found, extracts columns A (index 0) and D (index 3) and saves them in a new CSV called output_file.csv with two columns: Gene_cloud (changeble), Query (corresponding match).
#Creation of a new filtered .faa file. Reads the filtered CSV file to extract a list of queries. Searches for the queries in the .faa file (FASTA format).
#If it finds a match between the query and the sequence ID in the .faa:
#Creates a new header and writes the FASTA sequence with lines of 60 characters in the new output .faa file.

# Funzione per leggere il file txt e restituire una lista di nomi
def read_names_from_txt(txt_file):
    with open(txt_file, 'r') as f:
        names = f.read().strip().split(',')
    return [name.strip() for name in names]

# Funzione per cercare corrispondenze nel file CSV
def find_matches(txt_file, csv_file, output_file):
    # Leggi i nomi dal file txt
    names = read_names_from_txt(txt_file)
    
    # Leggi il file CSV (specificando il separatore tabulazione)
    df = pd.read_csv(csv_file, sep='\t')  # Usa '\t' per separazione tabulazioni
    
    # Lista per i risultati da scrivere nel nuovo CSV
    results = []
    
    # Cicla attraverso ogni nome
    for name in names:
        # Trova tutte le righe nel CSV dove la quarta colonna (indice 3) corrisponde al nome
        matches = df[df.iloc[:, 3].str.contains(name, case=False, na=False)]
        
        # Per ogni corrispondenza, aggiungi le colonne A e D al risultato
        for _, match in matches.iterrows():
            results.append([match.iloc[3], match.iloc[0]])  # Colonna A (index 0) e Colonna D (index 3)
    
    # Crea un DataFrame dai risultati
    result_df = pd.DataFrame(results, columns=['Gene_cloud', 'Query'])
    
    # Scrivi i risultati su un nuovo file CSV
    result_df.to_csv(output_file, index=False)
    print(f"Risultati salvati su {output_file}")

# Esegui la funzione con i tuoi file
txt_file = ''  # Sostituisci con il tuo file di testo
csv_file = ''  # Sostituisci con il tuo file .annotations
output_file = ''  # Il file CSV dove verranno salvati i risultati
find_matches(txt_file, csv_file, output_file)


# In[39]:


from Bio import SeqIO
import textwrap

# Funzione per leggere il file CSV e restituire una lista di query (seconda colonna)
def read_queries_from_csv(csv_file):
    import pandas as pd
    df = pd.read_csv(csv_file)
    return df['Query'].tolist()  # Assumendo che la colonna 'Query' esista nel CSV

# Funzione per cercare le query nel file .faa e scrivere il nuovo file .faa
def search_and_create_faa(csv_file, faa_file, output_faa_file):
    queries = read_queries_from_csv(csv_file)
    
    # Leggi il file .faa (formato FASTA)
    faa_sequences = list(SeqIO.parse(faa_file, "fasta"))
    
    # Lista per raccogliere le sequenze trovate
    output_sequences = []
    
    # Contatore progressivo per gli ID
    progressivo = 1
    
    # Cicla attraverso le query
    for query in queries:
        # Cicla attraverso tutte le sequenze nel file .faa
        for seq_record in faa_sequences:
            if query in seq_record.id:
                # Costruisci una nuova intestazione con il formato richiesto
                new_header = f">{query} #{progressivo}"
                
                # Formatta la sequenza in blocchi di 60 caratteri
                wrapped_sequence = textwrap.fill(str(seq_record.seq), width=60)
                
                # Aggiungi la sequenza e l'intestazione alla lista di output
                output_sequences.append(f"{new_header}\n{wrapped_sequence}")
                progressivo += 1  # Incrementa il contatore progressivo
    
    # Scrivi il risultato nel nuovo file .faa
    with open(output_faa_file, 'w') as output_handle:
        output_handle.write("\n".join(output_sequences))
    
    print(f"File {output_faa_file} creato con successo!")

# Esegui la funzione con i tuoi file
csv_file = ''  # Il tuo file CSV
faa_file = ''  # Il tuo file .faa
output_faa_file = ''  # File di output

search_and_create_faa(csv_file, faa_file, output_faa_file)




