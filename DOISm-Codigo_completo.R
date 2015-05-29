####################################################
## Preparação dos dados de expressão para análise ##
####################################################
# Carregar o pacote WGCNA
library(WGCNA);
enableWGCNAThreads();
# O script seguinte é importante, não omitir.
options(stringsAsFactors = FALSE);
# Ler o arquivo do data set
DOISm = read.table("Painel_E-MTAB-62_323_amostras_1293_genes.csv", header = TRUE, sep="\t") 
# Remover os dados auxiliares e transportar os dados de expressão para análise futura.
datExpr = as.data.frame(t(DOISm[, -c(1)]));
names(datExpr) = DOISm$Hybridization_REF;
rownames(datExpr) = names(DOISm)[-c(1)];
# Fazer a checagem para genes e amostras com muitos dados inexistentes:
gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK
# Se a afirmativa for verdadeira [TRUE], todos os genes passaram pelos cortes. 
# Fazer o agrupamento (clusterização) das amostras para observar a existência de outliers. 
library(flashClust) 
sampleTree = flashClust(dist(datExpr), method = "average");
# Abrir uma janela de saída gráfica para visualização dos clusters.
sizeGrWindow(12,9)
par(cex = 0.6); 
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Agrupamento da amostra para detecção de outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
# Plotar uma linha para mostrar o corte
abline(h = 50, col = "red");
# Determinar o grupamento sob a linha
clust = cutreeStatic(sampleTree, cutHeight = 50, minSize = 10) 
table(clust)
## Clust
##  0   1 
## 11 312
# O grupamento 1 contém as amostras que serão mantidas
keepSamples = (clust==1)
dadosExpr = datExpr[keepSamples, ] 
nGenes = ncol(dadosExpr)
nSamples = nrow(dadosExpr)
# Salvar os dados de expressão relevantes para utilização posterior.
save(dadosExpr, sampleTree, clust, file = "DOISm-dataInput.RData")

###############################################
## Construção da rede e detecção dos módulos ##
###############################################
# Se a sessão não estiver carregada, carregar o arquivo salvo
lnames = load(file = "DOISm-dataInput.RData");
# A variável lnames contém os nomes das variáveis carregadas
lnames
# Escolher um conjunto de poderes do limiar frouxo
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Chamar a função de análise de topologia da rede
sft = pickSoftThreshold(dadosExpr, powerVector = powers, verbose = 5)
# Plote os resultados 
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Índice de ajuste da topologia scale-free em função do limiar frouxo
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Limiar frouxo(poder)",ylab = "Modelo de ajuste da topologia scale Free, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red");
# Esta linha corresponde a utilização de um ponto de corte do h em R^2
abline(h=0.88,col="red")
# Conectividade media como uma função do limiar frouxo
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Limiar frouxo (poder)",ylab="Conectividade Média", type="n", main = paste("Conectividade Média"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
# Calcular as adjacências, utilizando o limiar frouxo 7
beta1 = 7;
adjacency = adjacency(dadosExpr, power = beta1);
# Criar um gráfico da topologia scale-free.
par(mfrow=c(1,1))
scaleFreePlot(adjacency, main=paste("Limiar frouxo, poder=",beta1), truncated=F); 
# Transformar a adjacência em sobreposição topológica
TOM = TOMsimilarity(adjacency); 
dissTOM = 1-TOM
# Criação do dendograma hierárquico dos genes
library(flashClust)
# Chamar a função de agrupamento hierárquico
geneTree = flashClust(as.dist(dissTOM), method = "average");
# Plotar a árvore de agrupamento resultante (dendograma)
sizeGrWindow(12,9)
pdf(file="cluster_genico.pdf",w=12,h=9)
plot(geneTree, xlab="", sub="", main = "Agrupamento gênico baseado na dissimilaridade TOM", labels = FALSE, hang = 0.04);
# Escolher o tamanho mínimo dos módulos. No caso escolhemos módulos relativamente grandes
minModuleSize = 30;
# Identificação do módulo utilizando o corte dinâmico de árvore
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);
table(dynamicMods)
write.csv(dynamicMods, file="dynamicMods.csv")
# Converter as legendas numéricas em cores 
dynamicColors = labels2colors(dynamicMods) 
table(dynamicColors)
write.csv(dynamicColors, file="dynamicColors.csv")
# Plotar o dendograma com as cores abaixo
sizeGrWindow(12,9)
plotDendroAndColors(geneTree, dynamicColors, "Corte dinâmico de árvore", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Dendograma gênico e cores dos módulos")
###############################
## Cálculo dos módulos eigengene ##
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
# Calcular o escore do module membership para cada gene (correlação com o módulo eigengene)
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

# Calcular a dissimilaridade dos módulos eigengene
MEDiss = 1-cor(MEs);
# Agrupar os módulos eigengene
METree = flashClust(as.dist(MEDiss), method = "average");
# Plotar o resultado
sizeGrWindow(7, 6)
pdf(file="clustering_eigengenes.pdf",w=7,h=6)
plot(METree, main = "Clustering dos módulos eigengenes", xlab = "", sub = "")
# O ponto de corte escolhido foi de 0,2, o que corresponde a uma correlação pairwise correlation de 0,8, para fazer a fusão
MEDissThres = 0.2
# Plotar a linha de corte no dendograma
abline(h=MEDissThres, col = "red")
# Chamar a função de fusão automática
merge = mergeCloseModules(dadosExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# As cores dos módulos fundidos
mergedColors = merge$colors;
# Nas análises subsequentes, serão utilizadas as cores dos módulos fundidas no mergedColors
# Eigengenes dos novos módulos fundidos
mergedMEs = merge$newMEs;
# Plotar o resultado
sizeGrWindow(12, 9)
pdf(file="clustering_merged.pdf",w=12,h=9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Corte dinâmico de árvore", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
# Criar um gráfico de dispersão das amostras ao longo dos módulos eigengenes
# Colocar os histogramas na diagonal
panel.hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
# Colocar as correlações (absolutas) nos painéis superiores, com o tamanho proporcional às correlações
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
pairs(MEordered, upper.panel = panel.smooth, lower.panel = panel.cor, diag.panel=panel.hist, main="Relação entre módulos eigengenes")
# Salvar as variáveis relevantes:
# Renomear para moduleColors
moduleColors = mergedColors
# Construir rótulos numéricos correspondentes às cores
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Salvar as cores dos módulos e os rótulos para utilização posterior
save(MEs, moduleLabels, moduleColors, geneTree, file = "DOISm-networkConstruction.RData")

######################################################
## Visualização da rede utilizando as funções WGCNA ##
######################################################
# Se a sessão não estiver carregada, carregar o arquivo salvo:
lnames = load(file = "DOISm-dataInput.RData");
# A variável lnames contém os nomes das variáveis carregadas
lnames
# CArregar os dados da rede
lnames = load(file = "DOISm-networkConstruction.RData"); 
lnames
nGenes = ncol(dadosExpr) 
nSamples = nrow(dadosExpr)
# Calcular a sobreposição topológica
beta1 = 7;
dissTOM = 1-TOMsimilarityFromExpr(dadosExpr, power = beta1);
# Transformar a dissTOM com um limiar frouxo para fazer com que conexões mais fortes sejam mais visíveis no heatmap
plotTOM = dissTOM^beta1;
# Definir o diagonal para NA para um gráfico mais agradável
diag(plotTOM) = NA;
# Chamar a função de gráfico
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Gráfico tipo heatmap da rede, todos os genes")
# Criar um gráfico de escala multidimensional (MDS) utilizando duas dimensões
cmd1 = cmdscale(as.dist(dissTOM), 3)
par(mfrow = c(1, 1))
plot(cmd1, col = moduleColors, main = "Gráfico MDS", xlab = "Dimensão 1", ylab = "Dimensão 2")

#################################################################################
## Produção de tabelas de genes para utilização em serviços e programas online ##
#################################################################################
# Se a sessão não estiver carregada, carregar o arquivo salvo:
lnames = load(file = "DOISm-dataInput.RData");
# A variável lnames contém os nomes das variáveis carregadas
lnames
# Carregar os dados da rede
lnames = load(file = "DOISm-networkConstruction.RData"); 
# A variável lnames contém os nomes as variáveis carregadas
lnames
# Fazer a leitura do arquivo de anotação
GeneAnnotation = read.csv(file = "E-MTAB-62_e_DOISm-identificadores.csv");
# Corresponder as amostras no data set com os IDs do arquivo de anotação
probes = names(dadosExpr)
probes2GeneAnnotation = match(probes, GeneAnnotation$Hybridization_REF)
# Fazer a correspondência dos Locus Link IDs
allLLIDs = GeneAnnotation$Entrez_Gene_ID[probes2GeneAnnotation];
# Escolher os módulos de interesse
intModules = c("turquoise", "blue", "brown", "yellow", "green")
for (module in intModules){
  # Selecionar as amostras nos módulos
  modGenes = (moduleColors==module)
  # Obter os códigos de identificação Entrez
  modLLIDs = allLLIDs[modGenes];
  # Escrever cada módulo em um arquivo
  fileName = paste("Entrez_Gene_IDs-", module, ".txt", sep=""); 
  write.table(as.data.frame(modLLIDs), file = fileName, row.names = FALSE, col.names = FALSE)
}
# Escrever todos os genes em um arquivo único
fileName = paste("Entrez_Gene_IDs-all.txt", sep=""); 
write.table(as.data.frame(allLLIDs), file = fileName, row.names = FALSE, col.names = FALSE)

################################################
## Análise de enriquecimento diretamente no R ##
################################################
# Essa é uma tabela contendo os 10 e os 30 melhores termos para cada módulo presente no moduleColors. Utilizar um ou outro para a construção da tabela de enriquecimento.
# Os 10 primeiros genes com melhor valor-p
GOenr = GOenrichmentAnalysis(moduleColors, allLLIDs, organism = "human", nBestP = 10);
# Os 30 primeiros genes com melhor valor-p
GOenr = GOenrichmentAnalysis(moduleColors, allLLIDs, organism = "human", nBestP = 30);
# A função acima retorna uma longa lista e os componentes mais interessantes estão agrupados utilizando esta função
tab = GOenr$bestPTerms[[4]]$enrichment
# Os nomes das colunas de cada tabela podem ser acessados por
names(tab)
write.table(tab, file = "GOEnrichmentTable.csv", sep = ",", quote = TRUE, row.names = FALSE)

###########################################################################
## Criação de arquivos para visualização utilizando o programa Cytoscape ##
###########################################################################
# Carregar o pacote WGCNA
library(WGCNA);
enableWGCNAThreads();
## Allowing parallel execution with up to 3 working processes.
# O script seguinte é importante, não omitir.
options(stringsAsFactors = FALSE);
# Carregar os dados da rede 
lnames = load(file = "DOISm-dataInput.RData");
# A variável lnames contém os nomes as variáveis carregadas
lnames
# Load network data saved in the second part.
lnames = load(file = "DOISm-networkConstruction.RData"); 
lnames
# Fazer a leitura a partir do arquivo de anotação
GeneAnnotation = read.csv(file = "E-MTAB-62_e_DOISm-identificadores.csv");
# Selecionar os módulos de interesse. Primeiramente, foi selecionado o módulo verde (módulo eigengene), depois todos os módulos juntos.
# modules = c("yellow", "turquoise", "green", "brown", "blue")
# modules = c("green")
# Selecionar os genes dos módulos
probes = names(dadosExpr);
inModule=is.finite(match(moduleColors,modules));
modProbes=probes[inModule];
match1=match(modProbes,GeneAnnotation$Hybridization_REF);
modGenes=GeneAnnotation$gene_symbol[match1]
# Calcular a adjacência, utilizando o limiar frouxo 7
beta1 = 7;
adjacency = adjacency(dadosExpr, power = beta1);
# Transformar a adjacência em sobreposição topológica
TOM = TOMsimilarity(adjacency); 
dissTOM = 1-TOM
# Selecionar as sobreposições topológicas correspondentes
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
# Exportar as redes em listas de bordas (edge) e nós (node) para o Cytoscape
# O limiar 0.02 foi escolhido por detectar um número satisfatório de genes dos módulos
cyt = exportNetworkToCytoscape(modTOM, edgeFile=paste("DOISm_Cytoscape_Input-Edge-",paste(modules,collapse="-"),".txt",sep=""), nodeFile=paste("DOISm_Cytoscape_Input-Node-",paste(modules,collapse="-"),".txt",sep=""), weighted = TRUE, threshold = 0.02, nodeNames=modProbes, altNodeNames = modGenes, nodeAttr = moduleColors[inModule])