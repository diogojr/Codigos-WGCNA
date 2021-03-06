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