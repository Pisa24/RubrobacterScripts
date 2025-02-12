library(TreeDist)
library(ape)
library(ComplexHeatmap)

#SCRIPT: prende in input due alberi in formato newick. Calcola: informazione condivisa, crea la matrice differenza tra le foglie,
# matrici binarie dei cladi, differenza di matrice binaria (heatmap)
#SCRIPT: Takes as input two trees in Newick format. Computes: shared information, creates the difference matrix between leaves, 
# binary clade matrices, binary matrix difference (heatmap).


# Caricare i due alberi
tree1 <- ape::read.tree("")
tree2 <- ape::read.tree("")

#FOGLIE COMUNI ENTROPIA DEL CLUSTERING E INFORMAZIONE CONDIVISA
# Identificare le foglie in comune tra i due alberi
commonTips <- intersect(tree1$tip.label, tree2$tip.label)
print(paste("Numero di foglie in comune:", length(commonTips)))

# Calcolare l'entropia di clustering per ciascun albero
entropy_tree1 <- ClusteringEntropy(tree1)
entropy_tree2 <- ClusteringEntropy(tree2)

print(paste("Entropia di clustering per tree1:", entropy_tree1))
print(paste("Entropia di clustering per tree2:", entropy_tree2))

# Calcolare l'informazione condivisa tra i due alberi
mutual_info <- MutualClusteringInfo(tree1, tree2)
print(paste("Informazione mutua (in bit):", mutual_info))

# Normalizzare rispetto all'informazione massima possibile
normalized_mutual_info <- MutualClusteringInfo(tree1, tree2, normalize = TRUE)
print(paste("Informazione mutua normalizzata:", normalized_mutual_info))



#MATRICE DELLE DIFFERENZE TRA FOGLIE. 
# Assumiamo che le foglie siano le stesse nei due alberi
dist_matrix_tree1 <- cophenetic(tree1)
dist_matrix_tree2 <- cophenetic(tree2)

# Calcolare la differenza assoluta tra le due matrici di distanza
diff_matrix <- abs(dist_matrix_tree1 - dist_matrix_tree2)
corr_matrix <- cor(diff_matrix, method = "pearson")

ComplexHeatmap::Heatmap(as.matrix(corr_matrix), cluster_rows = FALSE,
                        cluster_columns = FALSE, name = "Pearson correlation") 




#MATRICE BINARIA PER VEDERE SE DUE FOGLIE APPARTENGONO ALLO STESSO CLADE

num_tips <- Ntip(tree1)  # Numero di foglie
num_nodes <- Nnode(tree1)  # Numero totale di nodi


# Stampiamo per vedere i limiti corretti
print(paste("Numero di foglie:", num_tips))
print(paste("Numero totale di nodi:", num_nodes))


internal_nodes <- (num_tips + 1):(num_tips + num_nodes - 1)

# Ora possiamo estrarre un clade a partire da un nodo interno valido
clade_tree1 <- extract.clade(tree1, node = internal_nodes[1])  # Primo nodo interno valido
clade_tree2 <- extract.clade(tree2, node = internal_nodes[1])


clade_membership_matrix <- function(tree) {
  num_tips <- Ntip(tree)
  num_nodes <- Nnode(tree)
  internal_nodes <- (num_tips + 1):(num_tips + num_nodes - 1)  # Solo nodi interni
  
  # Creare una matrice binaria con le foglie nei cladi
  matrix <- matrix(0, nrow = num_tips, ncol = length(internal_nodes))
  rownames(matrix) <- tree$tip.label
  colnames(matrix) <- paste0("Clade_", internal_nodes)
  
  for (i in seq_along(internal_nodes)) {
    clade <- extract.clade(tree, node = internal_nodes[i])$tip.label
    matrix[rownames(matrix) %in% clade, i] <- 1
  }
  
  return(matrix)
}


# MATRICI DEI CLADI DEI SINGOLI ALBERI
ComplexHeatmap::Heatmap(
  as.matrix(clade_matrix_tree1), 
  name = "Clade Matrix Tree 1",
  row_names_side = "left",
  column_names_side = "top"
)

ComplexHeatmap::Heatmap(
  as.matrix(clade_matrix_tree2), 
  name = "Clade Matrix Tree 2",
  row_names_side = "left",
  column_names_side = "top"
)

#HEATMAP DIFFERENZA ALBERI DEI CLADI
# Creare le matrici per ciascun albero
clade_matrix_tree1 <- clade_membership_matrix(tree1)
clade_matrix_tree2 <- clade_membership_matrix(tree2)
diff_clade_matrix <- abs(clade_matrix_tree1 - clade_matrix_tree2)

# Visualizza la heatmap delle differenze nei cladi
ComplexHeatmap::Heatmap(
  as.matrix(diff_clade_matrix), 
  name = "HeatMap differenza cladi.",
  row_names_side = "left",
  column_names_side = "top"
)



dist_clades <- dist(diff_clade_matrix, method = "euclidean")

# Clustering gerarchico
hc <- hclust(dist_clades)

# Visualizza il dendrogramma dei cladi
plot(hc, main = "Dendrogramma delle differenze nei Cladi", xlab = "Cladi", ylab = "Distanza")

