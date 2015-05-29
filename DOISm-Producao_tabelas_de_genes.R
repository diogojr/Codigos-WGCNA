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