import ete3
import dendropy
import numpy as np
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import fcluster


#SCRIPT: prende in input i file newick di due alberi e opera diverse operazioni: distanza (Robinson-Foulds), matrice dei cladi, dendogramma delle differenze tra cladi.
#SCRIPT: Takes as input the Newick files of two trees and performs various operations: distance (Robinson-Foulds), clade matrix, dendrogram of clade differences.

# 🔹 Specifica i percorsi dei file Newick direttamente qui
tree1_path = "" 
tree2_path = "" 

def load_tree(file_path):
    """ Carica un albero filogenetico da un file Newick. """
    return ete3.Tree(file_path, format=1)


def calculate_rf_distance(tree1, tree2):
    """ Calcola la distanza di Robinson-Foulds tra due alberi non radicati. """
    result = tree1.robinson_foulds(tree2, unrooted_trees=True)
    rf_distance, max_rf = result[:2]  # Prendiamo solo i primi due valori
    return rf_distance / max_rf  # Normalizziamo la distanza


def create_clade_matrix(tree):
    """ Crea una matrice binaria rappresentante i cladi dell'albero. """
    taxa = tree.get_leaf_names()
    taxa_index = {taxa[i]: i for i in range(len(taxa))}
    num_taxa = len(taxa)
    num_clades = len(tree.get_descendants())

    matrix = np.zeros((num_taxa, num_clades), dtype=int)
    
    for i, node in enumerate(tree.get_descendants()):
        if not node.is_leaf():
            leaves = node.get_leaf_names()
            for leaf in leaves:
                matrix[taxa_index[leaf], i] = 1
    
    return matrix, taxa

def plot_heatmap(matrix, title):
    """ Visualizza una heatmap basata sulla matrice dei cladi. """
    plt.figure(figsize=(10, 6))
    sns.heatmap(matrix, cmap="coolwarm", cbar=True)
    plt.title(title)
    plt.xlabel("Cladi")
    plt.ylabel("Foglie")
    plt.show()


def hierarchical_clustering_colored(matrix, taxa_labels, num_clusters=5):
    """ Esegue il clustering gerarchico e colora i cladi principali. """
    distances = ssd.pdist(matrix, metric="euclidean")
    linkage_matrix = sch.linkage(distances, method="ward")
    
    plt.figure(figsize=(12, 6))
    dendrogram = sch.dendrogram(linkage_matrix, labels=taxa_labels, orientation="top", leaf_rotation=90, leaf_font_size=10,
                                color_threshold=0.7 * max(linkage_matrix[:,2]))  # Soglia per colorare i cluster
    
    plt.title(f"Dendrogramma delle Differenze nei Cladi (Colorato, {num_clusters} gruppi)")
    plt.xlabel("Cladi")
    plt.ylabel("Distanza")
    plt.show()


def compare_trees(tree1_path, tree2_path):
    """ Carica, confronta e visualizza i due alberi filogenetici. """
    
    # Caricamento alberi
    tree1 = load_tree(tree1_path)
    tree2 = load_tree(tree2_path)


    # Calcolo distanza di Robinson-Foulds
    rf_distance = calculate_rf_distance(tree1, tree2)
    print(f"Robinson-Foulds Distance: {rf_distance:.4f}")

    # Creazione delle matrici binarie dei cladi
    matrix1, taxa1 = create_clade_matrix(tree1)
    matrix2, taxa2 = create_clade_matrix(tree2)

    # Visualizzazione heatmap delle matrici
    plot_heatmap(matrix1, "Matrice dei Cladi - Albero 1")
    plot_heatmap(matrix2, "Matrice dei Cladi - Albero 2")

    # Calcolo e visualizzazione differenza tra le due matrici
    diff_matrix = np.abs(matrix1 - matrix2)
    plot_heatmap(diff_matrix, "Differenza tra le Matrici dei Cladi")

    # Clustering gerarchico della differenza tra le matrici
    hierarchical_clustering_colored(diff_matrix, taxa_labels = taxa1)

# 🔹 Esegui direttamente lo script senza passare argomenti
compare_trees(tree1_path, tree2_path)

