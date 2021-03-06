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