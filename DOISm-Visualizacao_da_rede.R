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