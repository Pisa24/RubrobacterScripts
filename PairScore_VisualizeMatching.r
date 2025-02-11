
library(TreeDist)
library(ape)

#SCRIPT: Studio di PairScores e VisualizeMatching, salva su file.txt PairScores
#SCRIPT: study of PairScores and VisualizeMatching, saves PairScores to a .txt file.


# Caricamento degli alberi
tree_genoma_codificante <- read.tree("")
tree_proteine <- read.tree("")

# Calcola solo le informazioni di pair scores
pair_scores <- ClusteringInfoDistance(tree_genoma_codificante, tree_proteine, reportMatching = TRUE)

# Salva l'output in un file di testo
output_file <- ""    # inserisci output file.txt
write.table(pair_scores, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Conferma della scrittura del file
cat("Le informazioni di pair scores sono state salvate in:", output_file, "\n")

SharedPhylogeneticInfo(tree_genoma_codificante, tree_proteine)
attr(SharedPhylogeneticInfo(tree_genoma_codificante, tree_proteine, reportMatching = TRUE), "matching")
VisualizeMatching(SharedPhylogeneticInfo, tree_genoma_codificante, tree_proteine) # Which clades are matched?
