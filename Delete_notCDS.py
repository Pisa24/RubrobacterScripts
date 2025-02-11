import os
#SCRIPT che cancella tutti i file che non abbiano intestazione cds
#SCRIPT that deletes all files that do not have a "cds" header.

def delete_non_GCA_GCF_files(directory):
    for root, dirs, files in os.walk(directory):
        for file in files:
            if not (file.startswith("cds")):
                file_path = os.path.join(root, file)
                try:
                    os.remove(file_path)
                    print(f"Eliminato: {file_path}")
                except Exception as e:
                    print(f"Errore eliminando {file_path}: {e}")

if __name__ == "__main__":
    directory = input("Inserisci il percorso della directory: ")
    if os.path.isdir(directory):
        delete_non_GCA_GCF_files(directory)
    else:
        print("Il percorso specificato non è una directory valida.")

