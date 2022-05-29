library(ontologyIndex)
library(ontologyPlot)
# library(apcluster)

#go <- get_ontology("/Users/eric/Documents/UTFPR/Mestrado/Projeto_do_mestrado/Conjunto_dados/termosGO-GOGO/go.obo", propagate_relationships=c("is_a", "part_of"), extract_tags="everything")
# go <- get_ontology("Datasets/go_termos/go-new.obo", propagate_relationships=c("is_a", "part_of"), extract_tags="everything")
go <- get_ontology("inst/extdata/go-new.obo", propagate_relationships=c("is_a", "part_of"), extract_tags="everything")

#adicionando o número de filhos a propriedades do termos
go$number_of_children <- sapply(go$children, length)

outliers <- function(x) {
  if(any(is.na(x)))
    stop("x is missing values")
  if(!is.numeric(x))
    stop("x is not numeric")
  Q3 <- quantile(x,0.75)
  Q1 <- quantile(x,0.25)
  IQR <- (Q3-Q1)
  left <- (Q1-(1.5*IQR))
  right <- (Q3+(1.5*IQR))
  c(x[(x > left) & (x < right)])
}

levelTree <- function(goTerm){
	goTerm <- goTerm[[1]]
	namespace <- go$namespace[[goTerm]]
	ancestors <- get_ancestors(go, goTerm)
	ancestorsBP<-c()
	for(x in ancestors){
		if(go$namespace[[x]] == namespace){
			ancestorsBP<-c(ancestorsBP, x)
		}
	}
	goCopy<-go
	for(term in ancestorsBP){
		if(length(goCopy$parents[[term]]) < 1){
			goCopy$depth[term] <- 1
		}else{
			parents <- intersect(goCopy$parents[[term]], ancestorsBP)
			parentsMaxDepth <- max(goCopy$depth[parents])
			goCopy$depth[term] <- parentsMaxDepth + 1
		}
	}
	return(goCopy$depth[[goTerm]])
}

go$level<-sapply(go$id, levelTree)

getLvl <- function(goTerm){
	return(go$level[goTerm][[1]])
}


similarity <- function(go1, go2){
	# print(paste("Termo GO 1: ", go1))
	# print(paste("Termo GO 2: ", go2))
	vectorA <- get_ancestors(go, go1)
	vectorB <- get_ancestors(go, go2)
	matrixWeightA <- matrix(0, nrow = length(vectorA), ncol = length(vectorA), dimnames = list(vectorA, vectorA))
	listGO<-intersect(vectorA, vectorB)
	# print(listGO)
	lvlsAncestrais <- sapply(listGO, getLvl)
	# print(lvlsAncestrais)
	if(length(listGO) == 0){
		maiorLvlAncestral <- 0
	}else{
		maiorLvlAncestral <- max(lvlsAncestrais)
	}
	#peso das arestas
	for(nodePai in vectorA){
		for(nodeFilho in vectorA){
			if(match(nodePai, go$part_of[nodeFilho][[1]], nomatch = 0) > 0){
				matrixWeightA[nodePai, nodeFilho] <- go$number_of_children[nodePai][[1]]^(0.6 - 1)
			}else{
				if(match(nodePai, go$is_a[nodeFilho][[1]], nomatch = 0) > 0){
					matrixWeightA[nodePai, nodeFilho] <- go$number_of_children[nodePai][[1]]^(0.8 - 1)
				}
			}	
		}
	}
	matrixWeightA <- matrixWeightA[,c(length(vectorA):1)] #inverte a ordem das colunas da matriz
	sValuesA <- matrix(0, nrow = 1, ncol = length(vectorA), dimnames = list("Values", rev(vectorA))) #cria uma matriz de uma linha para o S-values
	#termoA <- "GO:0005975" #termos em analise, como exemplos
	#termoB <- "GO:1901135"

	#s-values
	for(nodePai in rev(vectorA)){
			if(nodePai == go1){
				sValuesA[,nodePai] <- 1
			}else{
				auxSValues<-sValuesA
				for(x in 1:length(vectorA)){
					auxSValues[,x] <- auxSValues[x] * matrixWeightA[nodePai, x]
				}
				sValuesA[,nodePai] <- max(auxSValues)
			}
	}
	matrixWeightB <- matrix(0, nrow = length(vectorB), ncol = length(vectorB), dimnames = list(vectorB, vectorB))

	#peso das arestas
	for(nodePai in vectorB){
		for(nodeFilho in vectorB){
			if(match(nodePai, go$part_of[nodeFilho][[1]], nomatch = 0) > 0){
				matrixWeightB[nodePai, nodeFilho] <- go$number_of_children[nodePai][[1]]^(0.6 - 1)
			}else{
				if(match(nodePai, go$is_a[nodeFilho][[1]], nomatch = 0) > 0){
					matrixWeightB[nodePai, nodeFilho] <- go$number_of_children[nodePai][[1]]^(0.8 - 1)
				}
			}	
		}
	}


	matrixWeightB <- matrixWeightB[,c(length(vectorB):1)] #inverte a ordem das colunas da matriz
	sValuesB <- matrix(0, nrow = 1, ncol = length(vectorB), dimnames = list("Values", rev(vectorB))) #cria uma matriz de uma linha para o S-values
	#termoA <- "GO:0005975" #termos em analise, como exemplos
	#termoB <- "GO:1901135"
	
	#s-values
	for(nodePai in rev(vectorB)){
			if(nodePai == go2){
				sValuesB[,nodePai] <- 1
			}else{
				auxSValues<-sValuesB
				for(x in 1:length(vectorB)){
					auxSValues[,x] <- auxSValues[x] * matrixWeightB[nodePai, x]
				}
				sValuesB[,nodePai] <- max(auxSValues)
			}
	}
	#Ancestrais em comum
	ancestrais <- ""
	ancestrais <- intersect(get_ancestors(go, c(go1)), get_ancestors(go, c(go2)))
	if(length(ancestrais) != 0){
		svaluesAncestrais <- 0
		for(x in ancestrais){
			svaluesAncestrais <- sValuesA[,x] + sValuesB[,x] + svaluesAncestrais
		}
		# print(svaluesAncestrais)
		# print((sum(sValuesA) + sum(sValuesB)))

		similarity <- svaluesAncestrais/(sum(sValuesA) + sum(sValuesB))
		# print(similarity)
	}else{
		similarity <- 0
	}
	vetorAux<-c()
	for(x in ancestrais){
		vetorAux<-c(vetorAux, levelTree(x))
	}
	# similarity<-similarity*(1-exp(-1)*((1^maiorLvlAncestral)/factorial(maiorLvlAncestral)))


	return(similarity)
}

# levelTree <- function(goTerm){
# 	goTerm <- goTerm[[1]]
# 	ancestors <- get_ancestors(go, goTerm)
# 	ancestorsBP<-c()
# 	for(x in ancestors){
# 		if(go$namespace[[x]] == "biological_process"){
# 			ancestorsBP<-c(ancestorsBP, x)
# 		}
# 	}
# 	goCopy<-go
# 	for(term in ancestorsBP){
# 		if(length(goCopy$parents[[term]]) < 1){
# 			goCopy$depth[term] <- 1
# 		}else{
# 			parents <- intersect(goCopy$parents[[term]], ancestorsBP)
# 			parentsMaxDepth <- max(goCopy$depth[parents])
# 			goCopy$depth[term] <- parentsMaxDepth + 1
# 		}
# 	}
# 	return(goCopy$depth[[goTerm]])
# }

calcSimilarityGene <- function(g1, g2){
	 tryCatch(
        expr = {
            sumBMgene1 <- 0
			# print(paste("Lista Genes: ", g1[1]))
			# print(paste("Lista Genes: ", g2[1]))
			listSimilarities<-c()

			ListGene1Edited <- paste("G1.", g1, sep="")
			ListGene2Edited <- paste("G2.", g2, sep="")

			listGenes <- unique(c(ListGene1Edited, ListGene2Edited))

			squadMatrix <- matrix(nrow = length(listGenes), ncol=length(listGenes), dimnames = list(listGenes, listGenes))

			valoresSimilarity<-c()

			for(gi in listGenes){
				for(gk in listGenes){
					if((gi == gk)){
						squadMatrix[gi,gk] <- NA
					}else{
						if(substr(gi, 1, 2) == substr(gk, 1, 2)){
							squadMatrix[gi,gk] <- NA
						}else{
							resultSimilarity <- similarity(substr(gi, 4, 13), substr(gk, 4, 13))
							squadMatrix[gi,gk] <- resultSimilarity
						}
						
					}
				}
			}

			# print(squadMatrix)

			for(gi in listGenes){
				maxRow<-max(squadMatrix[gi,][!is.na(squadMatrix[gi,])])
				for(gk in listGenes){
					maxColumn<-max(squadMatrix[,gk][!is.na(squadMatrix[,gk])])
					if(!is.na(squadMatrix[gi,gk])){
						if((maxRow != squadMatrix[gi,gk]) && (maxColumn != squadMatrix[gi,gk])){
							squadMatrix[gi,gk]<-0
						}
					}
				}
			}

			print("Max das combinacoes")
			print(squadMatrix)

			squadMatrix[squadMatrix == 0] <- NA

			for(x in 1:length(squadMatrix[,1])){
				for(y in 1:length(squadMatrix[1,])){
					if(!is.na(squadMatrix[x,y])){
						valoresSimilarity<-c(valoresSimilarity, squadMatrix[x,y])
					}
				}
			}

			print("valoresSimilarity")
			print(valoresSimilarity)

			apres <- apcluster(squadMatrix, details=TRUE)

			# print(squadMatrix)

			size <- length(apres)
			# print(paste("size", size))
			print(apres)

			mediaClusters<-c()

			verificadorCluster <- 0

			if((size > 1) && (length(squadMatrix[1,]) > 2)){
				for(x in 1:size){
					termos <- names(apres[[x]])
					print("Termos presentes no cluster")
					print(termos)
					listaTermosGODistinct <- c()
					for(termo in termos){
						listaTermosGODistinct <- c(listaTermosGODistinct, substr(termo, 1, 3))
					}
					
					print("Lista de genes no cluster")
					print(listaTermosGODistinct)

					print("Lista de genes distintos no cluster")
					print(unique(listaTermosGODistinct))

					numTermosGODistinct <- length(unique(listaTermosGODistinct))

					if(numTermosGODistinct < 2){
						verificadorCluster <- verificadorCluster + 1
					}else{
						values <- c()

						sizeCluster <- length(termos)

						if(sizeCluster != 1){
							for(i in 1:(sizeCluster-1)){
								for(k in (1+i):sizeCluster){
									go1<-termos[i]
									go2<-termos[k]
									values<-c(values, squadMatrix[go1,go2])
								}
							}
						}else{
							values <- 0
						}

						values <- values[!is.na(values)]
						# print(values)
						mediaClusters <- c(mediaClusters, mean(values))
					}
				}
			}else{
				mediaClusters <- mean(valoresSimilarity)
			}


			print("Cluster sem termos dos dois genes")
			print(verificadorCluster)

			print("Numero de clusters formados")
			print(size)
			
			if(verificadorCluster == size){
				similarityGene <- mean(valoresSimilarity)
			}else{
				# print(mediaClusters)
				mediaClusters <- mediaClusters[!is.na(mediaClusters)]
				similarityGene <- max(mediaClusters)
			}

			print("Similaridade")
			print(similarityGene)

			return(round(similarityGene, digits=3))
        },
        error = function(e){
            message('Caught an error!')
            print(e)
            valoresSimilarity<-c()
            for(gi in listGenes){
				for(gk in listGenes){
					if((gi == gk)){
						squadMatrix[gi,gk] <- NA
					}else{
						if(substr(gi, 1, 2) == substr(gk, 1, 2)){
							squadMatrix[gi,gk] <- NA
						}else{
							resultSimilarity <- similarity(substr(gi, 4, 13), substr(gk, 4, 13))
							squadMatrix[gi,gk] <- resultSimilarity
							valoresSimilarity <- c(valoresSimilarity, resultSimilarity)
						}
						
					}
					
				}
			}

			for(gi in listGenes){
				maxRow<-max(squadMatrix[gi,][!is.na(squadMatrix[gi,])])
				for(gk in listGenes){
					maxColumn<-max(squadMatrix[,gk][!is.na(squadMatrix[,gk])])
					if(!is.na(squadMatrix[gi,gk])){
						if((maxRow != squadMatrix[gi,gk]) && (maxColumn != squadMatrix[gi,gk])){
							squadMatrix[gi,gk]<-0
						}
					}
				}
			}

			squadMatrix[squadMatrix == 0] <- NA

			# print(round(mean(valoresSimilarity[!is.na(valoresSimilarity)], na.rm=TRUE), digits=3))
            return(round(mean(valoresSimilarity[!is.na(valoresSimilarity)], na.rm=TRUE), digits=3))
        }
    )    
}

similarityGene <- function(g1, g2){
	
	g1Name<-g1[1]
	g2Name<-g2[1]

	print(paste("Gene: ", g1[1]))
	print(paste("Gene: ", g2[1]))

	write(paste(g1Name, g2Name, sep=" - "), file="testString.txt",append=TRUE)
	g1<-g1[2:length(g1)]
	g2<-g2[2:length(g2)]

	g1MF<-c()
	g1CC<-c()
	g1BP<-c()

	for(goTerm in g1){
		if(go$namespace[goTerm] == "biological_process"){
			g1BP<-c(g1BP, goTerm)
		}else{
			if(go$namespace[goTerm] == "molecular_function"){
				g1MF<-c(g1MF, goTerm)
			}else{
				if(go$namespace[goTerm] == "cellular_component"){
					g1CC<-c(g1CC, goTerm)
				}
			}
		}
	}

	g2MF<-c()
	g2CC<-c()
	g2BP<-c()

	for(goTerm in g2){
		if(go$namespace[goTerm] == "biological_process"){
			g2BP<-c(g2BP, goTerm)
		}else{
			if(go$namespace[goTerm] == "molecular_function"){
				g2MF<-c(g2MF, goTerm)
			}else{
				if(go$namespace[goTerm] == "cellular_component"){
					g2CC<-c(g2CC, goTerm)
				}
			}
		}
	}

	if((length(g1MF) > 0 ) && (length(g2MF) > 0)){
		print("molecular_function")
		write("molecular_function", file="testString.txt",append=TRUE)
		simMF <- calcSimilarityGene(g1MF, g2MF)
	}else{
		write("molecular_function", file="testString.txt",append=TRUE)
		simMF <- 0
	}
	if((length(g1CC) > 0 ) && (length(g2CC) > 0)){
		print("cellular_Component")
		write("cellular_Component", file="testString.txt",append=TRUE)
		simCC <- calcSimilarityGene(g1CC, g2CC)
	}else{
		write("cellular_Component", file="testString.txt",append=TRUE)
		simCC <- 0
	}
	if((length(g1BP) > 0 ) && (length(g2BP) > 0)){
		print("biological_process")
		write("biological_process", file="testString.txt",append=TRUE)
		simBP <- calcSimilarityGene(g1BP, g2BP)
	}else{
		write("biological_process", file="testString.txt",append=TRUE)
		simBP <- 0
	}

	#similarityGene <- c(paste(g1Name, " - ", g2Name, "\nMolecular_Function: ", simMF),paste("Cellular_Component: ", simCC),paste("Biological_Process: ", simBP))
	similarityGene <- c(paste(g1Name, g2Name, "BPO", simBP, "CCO", simCC, "MFO", simMF, sep=","))
	return(similarityGene)
	
}

arquivo<-"../dataset/mannose/genes-Go.txt"

numLines <- length(readLines(arquivo));
csv <- read.csv(file=arquivo, header=FALSE, sep=",")
resultados<-c()
for(x in 1:(numLines-1)){
	for(y in seq((x+1),numLines)){
		g1<-scan(arquivo,what="character", sep=",", skip = (x-1), nlines = 1)
		g2<-scan(arquivo,what="character", sep=",", skip = (y-1), nlines = 1)
		print(paste("Numero:",x,"de", numLines))
		print(paste("Numero:",y,"de", numLines))
		resultados<-c(resultados, similarityGene(g1,g2))
	}
}

fileConn<-file("../results/mannose.txt")
writeLines(resultados, fileConn)
close(fileConn)

