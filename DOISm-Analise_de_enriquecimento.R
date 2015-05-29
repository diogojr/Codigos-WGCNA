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