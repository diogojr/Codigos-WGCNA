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

###############################################
## Constru��o da rede e detec��o dos m�dulos ##
###############################################
# Se a sess�o n�o estiver carregada, carregar o arquivo salvo
lnames = load(file = "DOISm-dataInput.RData");
# A vari�vel lnames cont�m os nomes das vari�veis carregadas
lnames
# Escolher um conjunto de poderes do limiar frouxo
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Chamar a fun��o de an�lise de topologia da rede
sft = pickSoftThreshold(dadosExpr, powerVector = powers, verbose = 5)
# Plote os resultados 
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# �ndice de ajuste da topologia scale-free em fun��o do limiar frouxo
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Limiar frouxo(poder)",ylab = "Modelo de ajuste da topologia scale Free, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red");
# Esta linha corresponde a utiliza��o de um ponto de corte do h em R^2
abline(h=0.88,col="red")
# Conectividade media como uma fun��o do limiar frouxo
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Limiar frouxo (poder)",ylab="Conectividade M�dia", type="n", main = paste("Conectividade M�dia"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
# Calcular as adjac�ncias, utilizando o limiar frouxo 7
beta1 = 7;
adjacency = adjacency(dadosExpr, power = beta1);
# Criar um gr�fico da topologia scale-free.
par(mfrow=c(1,1))
scaleFreePlot(adjacency, main=paste("Limiar frouxo, poder=",beta1), truncated=F); 
# Transformar a adjac�ncia em sobreposi��o topol�gica
TOM = TOMsimilarity(adjacency); 
dissTOM = 1-TOM
# Cria��o do dendograma hier�rquico dos genes
library(flashClust)
# Chamar a fun��o de agrupamento hier�rquico
geneTree = flashClust(as.dist(dissTOM), method = "average");
# Plotar a �rvore de agrupamento resultante (dendograma)
sizeGrWindow(12,9)
pdf(file="cluster_genico.pdf",w=12,h=9)
plot(geneTree, xlab="", sub="", main = "Agrupamento g�nico baseado na dissimilaridade TOM", labels = FALSE, hang = 0.04);
# Escolher o tamanho m�nimo dos m�dulos. No caso escolhemos m�dulos relativamente grandes
minModuleSize = 30;
# Identifica��o do m�dulo utilizando o corte din�mico de �rvore
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);
table(dynamicMods)
write.csv(dynamicMods, file="dynamicMods.csv")
# Converter as legendas num�ricas em cores 
dynamicColors = labels2colors(dynamicMods) 
table(dynamicColors)
write.csv(dynamicColors, file="dynamicColors.csv")
# Plotar o dendograma com as cores abaixo
sizeGrWindow(12,9)
plotDendroAndColors(geneTree, dynamicColors, "Corte din�mico de �rvore", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Dendograma g�nico e cores dos m�dulos")
###############################
## C�lculo dos m�dulos eigengene ##
###############################
# Calcular os eigengenes e checar os escores do PCA dos eigengenes
MEList = moduleEigengenes(dadosExpr, colors = dynamicColors)
MEs = MEList$eigengenes
signif(cor(MEs, use="p"), 2)
MEdata = moduleEigengenes(datExpr, moduleColors, softPower = softPower)
varExpl = t(MEList$varExplained)
rownames(varExpl) = colnames(MEList$eigengenes)
write.table(varExpl, file = paste("eigengene.PC1.scores.txt", sep=""), sep = "\t", quote = F, col.names=F)
minVarExpl = 0.5
validMods = rownames(varExpl)[varExpl >= minVarExpl]
# Calcular o escore do module membership para cada gene (correla��o com o m�dulo eigengene)
modMembership = as.matrix(cor(dadosExpr, MEs, use = 'p'))
modMembershipP = as.matrix(corPvalueStudent(modMembership, nrow(dadosExpr)))
names(dynamicColors) = colnames(dadosExpr)
modMembershipGenes = t(sapply(rownames(modMembership), function(g) {
  m = paste("ME", dynamicColors[g], sep="")
  if(!(m %in% c("MENA", "MEgrey")) & m %in% colnames(modMembership)) {
    c(modMembership[g, m], modMembershipP[g, m])
  } else {
    c(NA,NA)
  }
}))
colnames(modMembershipGenes) = c("modMembershipR", "modMembershipP")
write.table(cbind(symbol = as.character(mget(rownames(modMembershipGenes), ifnotfound = NA)), modMembershipGenes, module = dynamicColors),
            file = paste("wgcna.modmembership.txt", sep=""), sep="\t", quote=F, col.names=NA)

# Calcular a dissimilaridade dos m�dulos eigengene
MEDiss = 1-cor(MEs);
# Agrupar os m�dulos eigengene
METree = flashClust(as.dist(MEDiss), method = "average");
# Plotar o resultado
sizeGrWindow(7, 6)
pdf(file="clustering_eigengenes.pdf",w=7,h=6)
plot(METree, main = "Clustering dos m�dulos eigengenes", xlab = "", sub = "")
# O ponto de corte escolhido foi de 0,2, o que corresponde a uma correla��o pairwise correlation de 0,8, para fazer a fus�o
MEDissThres = 0.2
# Plotar a linha de corte no dendograma
abline(h=MEDissThres, col = "red")
# Chamar a fun��o de fus�o autom�tica
merge = mergeCloseModules(dadosExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# As cores dos m�dulos fundidos
mergedColors = merge$colors;
# Nas an�lises subsequentes, ser�o utilizadas as cores dos m�dulos fundidas no mergedColors
# Eigengenes dos novos m�dulos fundidos
mergedMEs = merge$newMEs;
# Plotar o resultado
sizeGrWindow(12, 9)
pdf(file="clustering_merged.pdf",w=12,h=9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Corte din�mico de �rvore", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
# Criar um gr�fico de dispers�o das amostras ao longo dos m�dulos eigengenes
# Colocar os histogramas na diagonal
panel.hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
# Colocar as correla��es (absolutas) nos pain�is superiores, com o tamanho proporcional �s correla��es
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) { 
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  r <- abs(cor(x, y)) 
  txt <- format(c(r, 0.123456789), digits = digits)[1] 
  txt <- paste0(prefix, txt) 
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt) 
  text(0.5, 0.5, txt, cex = cex.cor * r) 
} 
MEordered = MEList[,METree$order]
pairs(MEordered, upper.panel = panel.smooth, lower.panel = panel.cor, diag.panel=panel.hist, main="Rela��o entre m�dulos eigengenes")
# Salvar as vari�veis relevantes:
# Renomear para moduleColors
moduleColors = mergedColors
# Construir r�tulos num�ricos correspondentes �s cores
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Salvar as cores dos m�dulos e os r�tulos para utiliza��o posterior
save(MEs, moduleLabels, moduleColors, geneTree, file = "DOISm-networkConstruction.RData")

######################################################
## Visualiza��o da rede utilizando as fun��es WGCNA ##
######################################################
# Se a sess�o n�o estiver carregada, carregar o arquivo salvo:
lnames = load(file = "DOISm-dataInput.RData");
# A vari�vel lnames cont�m os nomes das vari�veis carregadas
lnames
# CArregar os dados da rede
lnames = load(file = "DOISm-networkConstruction.RData"); 
lnames
nGenes = ncol(dadosExpr) 
nSamples = nrow(dadosExpr)
# Calcular a sobreposi��o topol�gica
beta1 = 7;
dissTOM = 1-TOMsimilarityFromExpr(dadosExpr, power = beta1);
# Transformar a dissTOM com um limiar frouxo para fazer com que conex�es mais fortes sejam mais vis�veis no heatmap
plotTOM = dissTOM^beta1;
# Definir o diagonal para NA para um gr�fico mais agrad�vel
diag(plotTOM) = NA;
# Chamar a fun��o de gr�fico
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Gr�fico tipo heatmap da rede, todos os genes")
# Criar um gr�fico de escala multidimensional (MDS) utilizando duas dimens�es
cmd1 = cmdscale(as.dist(dissTOM), 3)
par(mfrow = c(1, 1))
plot(cmd1, col = moduleColors, main = "Gr�fico MDS", xlab = "Dimens�o 1", ylab = "Dimens�o 2")

#################################################################################
## Produ��o de tabelas de genes para utiliza��o em servi�os e programas online ##
#################################################################################
# Se a sess�o n�o estiver carregada, carregar o arquivo salvo:
lnames = load(file = "DOISm-dataInput.RData");
# A vari�vel lnames cont�m os nomes das vari�veis carregadas
lnames
# Carregar os dados da rede
lnames = load(file = "DOISm-networkConstruction.RData"); 
# A vari�vel lnames cont�m os nomes as vari�veis carregadas
lnames
# Fazer a leitura do arquivo de anota��o
GeneAnnotation = read.csv(file = "E-MTAB-62_e_DOISm-identificadores.csv");
# Corresponder as amostras no data set com os IDs do arquivo de anota��o
probes = names(dadosExpr)
probes2GeneAnnotation = match(probes, GeneAnnotation$Hybridization_REF)
# Fazer a correspond�ncia dos Locus Link IDs
allLLIDs = GeneAnnotation$Entrez_Gene_ID[probes2GeneAnnotation];
# Escolher os m�dulos de interesse
intModules = c("turquoise", "blue", "brown", "yellow", "green")
for (module in intModules){
  # Selecionar as amostras nos m�dulos
  modGenes = (moduleColors==module)
  # Obter os c�digos de identifica��o Entrez
  modLLIDs = allLLIDs[modGenes];
  # Escrever cada m�dulo em um arquivo
  fileName = paste("Entrez_Gene_IDs-", module, ".txt", sep=""); 
  write.table(as.data.frame(modLLIDs), file = fileName, row.names = FALSE, col.names = FALSE)
}
# Escrever todos os genes em um arquivo �nico
fileName = paste("Entrez_Gene_IDs-all.txt", sep=""); 
write.table(as.data.frame(allLLIDs), file = fileName, row.names = FALSE, col.names = FALSE)

################################################
## An�lise de enriquecimento diretamente no R ##
################################################
# Essa � uma tabela contendo os 10 e os 30 melhores termos para cada m�dulo presente no moduleColors. Utilizar um ou outro para a constru��o da tabela de enriquecimento.
# Os 10 primeiros genes com melhor valor-p
GOenr = GOenrichmentAnalysis(moduleColors, allLLIDs, organism = "human", nBestP = 10);
# Os 30 primeiros genes com melhor valor-p
GOenr = GOenrichmentAnalysis(moduleColors, allLLIDs, organism = "human", nBestP = 30);
# A fun��o acima retorna uma longa lista e os componentes mais interessantes est�o agrupados utilizando esta fun��o
tab = GOenr$bestPTerms[[4]]$enrichment
# Os nomes das colunas de cada tabela podem ser acessados por
names(tab)
write.table(tab, file = "GOEnrichmentTable.csv", sep = ",", quote = TRUE, row.names = FALSE)

###########################################################################
## Cria��o de arquivos para visualiza��o utilizando o programa Cytoscape ##
###########################################################################
# Carregar o pacote WGCNA
library(WGCNA);
enableWGCNAThreads();
## Allowing parallel execution with up to 3 working processes.
# O script seguinte � importante, n�o omitir.
options(stringsAsFactors = FALSE);
# Carregar os dados da rede 
lnames = load(file = "DOISm-dataInput.RData");
# A vari�vel lnames cont�m os nomes as vari�veis carregadas
lnames
# Load network data saved in the second part.
lnames = load(file = "DOISm-networkConstruction.RData"); 
lnames
# Fazer a leitura a partir do arquivo de anota��o
GeneAnnotation = read.csv(file = "E-MTAB-62_e_DOISm-identificadores.csv");
# Selecionar os m�dulos de interesse. Primeiramente, foi selecionado o m�dulo verde (m�dulo eigengene), depois todos os m�dulos juntos.
# modules = c("yellow", "turquoise", "green", "brown", "blue")
# modules = c("green")
# Selecionar os genes dos m�dulos
probes = names(dadosExpr);
inModule=is.finite(match(moduleColors,modules));
modProbes=probes[inModule];
match1=match(modProbes,GeneAnnotation$Hybridization_REF);
modGenes=GeneAnnotation$gene_symbol[match1]
# Calcular a adjac�ncia, utilizando o limiar frouxo 7
beta1 = 7;
adjacency = adjacency(dadosExpr, power = beta1);
# Transformar a adjac�ncia em sobreposi��o topol�gica
TOM = TOMsimilarity(adjacency); 
dissTOM = 1-TOM
# Selecionar as sobreposi��es topol�gicas correspondentes
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
# Exportar as redes em listas de bordas (edge) e n�s (node) para o Cytoscape
# O limiar 0.02 foi escolhido por detectar um n�mero satisfat�rio de genes dos m�dulos
cyt = exportNetworkToCytoscape(modTOM, edgeFile=paste("DOISm_Cytoscape_Input-Edge-",paste(modules,collapse="-"),".txt",sep=""), nodeFile=paste("DOISm_Cytoscape_Input-Node-",paste(modules,collapse="-"),".txt",sep=""), weighted = TRUE, threshold = 0.02, nodeNames=modProbes, altNodeNames = modGenes, nodeAttr = moduleColors[inModule])