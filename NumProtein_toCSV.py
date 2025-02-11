import os
import csv

#SCRIPT che data una directory entra ricorsivamente in ciascun file FASTA e conta il numero di proteine in ciascun file
#Here is a Bash script that recursively scans a directory, finds all FASTA files, and counts the number of protein sequences in each file. The script assumes that the FASTA files are in standard format, where each protein sequence starts with a ">" character.

def count_proteins_in_faa(directory, output_file):
    protein_counts = []
    
    for filename in os.listdir(directory):
        if filename.endswith(".fna") | filename.endswith(".faa"):  # Controlla che il file sia di tipo .faa
            filepath = os.path.join(directory, filename)
            try:
                with open(filepath, "r", encoding="utf-8") as file:
                    count = sum(1 for line in file if line.startswith(">"))
                protein_counts.append([filename, count])
            except Exception as e:
                print(f"Errore nella lettura del file {filename}: {e}")
    
    with open(output_file, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["File", "Numero di Proteine"])
        writer.writerows(protein_counts)

if __name__ == "__main__":
    directory = input("Inserisci il percorso della directory: ")
    output_file = input("Inserisci il percorso del file di output CSV: ")
    count_proteins_in_faa(directory, output_file)
    print(f"Risultati salvati in {output_file}")

