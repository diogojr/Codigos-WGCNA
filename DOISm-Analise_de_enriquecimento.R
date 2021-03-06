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