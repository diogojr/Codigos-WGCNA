####################################################
## Prepara��o dos dados de express�o para an�lise ##
####################################################
# Carregar o pacote WGCNA
library(WGCNA);
enableWGCNAThreads();
# O script seguinte � importante, n�o omitir.
options(stringsAsFactors = FALSE);
# Ler o arquivo do data set
DOISm = read.table("Painel_E-MTAB-62_323_amostras_1293_genes.csv", header = TRUE, sep="\t") 
# Remover os dados auxiliares e transportar os dados de express�o para an�lise futura.
datExpr = as.data.frame(t(DOISm[, -c(1)]));
names(datExpr) = DOISm$Hybridization_REF;
rownames(datExpr) = names(DOISm)[-c(1)];
# Fazer a checagem para genes e amostras com muitos dados inexistentes:
gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK
# Se a afirmativa for verdadeira [TRUE], todos os genes passaram pelos cortes. 
# Fazer o agrupamento (clusteriza��o) das amostras para observar a exist�ncia de outliers. 
library(flashClust) 
sampleTree = flashClust(dist(datExpr), method = "average");
# Abrir uma janela de sa�da gr�fica para visualiza��o dos clusters.
sizeGrWindow(12,9)
par(cex = 0.6); 
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Agrupamento da amostra para detec��o de outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
# Plotar uma linha para mostrar o corte
abline(h = 50, col = "red");
# Determinar o grupamento sob a linha
clust = cutreeStatic(sampleTree, cutHeight = 50, minSize = 10) 
table(clust)
## Clust
##  0   1 
## 11 312
# O grupamento 1 cont�m as amostras que ser�o mantidas
keepSamples = (clust==1)
dadosExpr = datExpr[keepSamples, ] 
nGenes = ncol(dadosExpr)
nSamples = nrow(dadosExpr)
# Salvar os dados de express�o relevantes para utiliza��o posterior.
save(dadosExpr, sampleTree, clust, file = "DOISm-dataInput.RData")
