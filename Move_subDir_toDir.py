import os
import shutil

#SCRIPT di supporto, apre una directory, ogni cartella sottodirectory ed estrae un solo file da ciascuna. Alla fine li sposta tutti in una directory desiderata.
#Support script that opens a directory, accesses each subdirectory, and extracts a single file from each. Finally, it moves all extracted files to a designated directory.

def move_and_rename_files(input_directory):
    if not os.path.isdir(input_directory):
        print("Errore: la directory di input non esiste.")
        return
    
    for folder in os.listdir(input_directory):
        folder_path = os.path.join(input_directory, folder)
        
        if os.path.isdir(folder_path):  # Controlla che sia una cartella
            files = os.listdir(folder_path)
            
            if len(files) == 1:  # Assicura che ci sia un solo file
                file_name = files[0]
                original_file_path = os.path.join(folder_path, file_name)
                new_file_name = f"{folder}_{file_name}"
                new_file_path = os.path.join(input_directory, new_file_name)
                
                shutil.move(original_file_path, new_file_path)
                print(f"Spostato: {original_file_path} -> {new_file_path}")
            else:
                print(f"Attenzione: La cartella '{folder}' non contiene esattamente un file.")

# Esempio di utilizzo
input_directory = ""
move_and_rename_files(input_directory)

