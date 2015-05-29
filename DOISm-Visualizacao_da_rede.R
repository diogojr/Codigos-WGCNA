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