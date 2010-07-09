random.polychor.pa <- function(nvar="NULL", n.ss="NULL", nrep, nstep="NULL", data.matrix, q.eigen, r.seed="NULL") {
	start.t<- Sys.time()
	cat("computation starts at:", format(start.t, "%X"), "\n", "\n")
	flush.console()
	require(psych)
	require(nFactors)
	cat("\n")

	for (z in 1:ncol(data.matrix)) {
		if(is.numeric(data.matrix[,z])==FALSE) {
			data.matrix[,z]<-as.numeric(data.matrix[,z]) 
		}
	}

	data.matrix <- as.matrix(data.matrix)
	if (!is.null(dimnames(data.matrix))) {
		dimnames(data.matrix) <- list(NULL, NULL)
		data.matrix
		}
	if (is.null(q.eigen)) {
		stop("the quantile to be extracted is not declared. Please provide a valid quantile")
		}
	if (r.seed=="NULL") {  set.seed(1335031435)
		}
		else { set.seed(r.seed)
		}
	data.matrix.0<-na.omit(data.matrix)
	n.ss<-nrow(data.matrix)
	cat("number of units (rows) in data.matrix:", n.ss, "\n")
	if (nrow(data.matrix.0) < nrow(data.matrix)) {
		cat("########################################################################", "\n")
		cat("######################## Missing Values (NA) detected", "\n")
		cat("######################## LISTWISE deletion needed. The new sample size is:", nrow(data.matrix.0), "\n")
		cat("########################################################################", "\n")
		data.matrix<-data.matrix.0
		n.ss <- nrow(data.matrix.0)
		}
	nvar<-ncol(data.matrix)
	cat("number of variables (cols) in data.matrix:", nvar, "\n", "\n")
	flush.console()

	### loop che effettua un confronto tra item.tab fornita e quella ricavata da data.matrix per capire che se i dati di data.matrix e item.tab coincidono
	item.tab.ex<-matrix(,ncol(data.matrix), 2)
	for (h in 1:ncol(data.matrix)) {
		item.tab.ex[h,2]<-nrow(as.data.frame(table(data.matrix[,h])))
		item.tab.ex[h,1]<-h
		}
	item.tab<-matrix(,nrow(table(item.tab.ex[,2])), 2)
	item.tab.ex<-(as.data.frame(table(item.tab.ex[,2])))
	item.tab[,1]<-item.tab.ex[,2]
	item.tab[,2]<-as.numeric(levels(item.tab.ex[,1]))
	colnames(item.tab)<-c("Items", "Categories")
	row.labels<-paste(c(1:nrow(item.tab)), "GROUP")
	rownames(item.tab)<-c(row.labels)
	cat("the following table shows the groups of items with diffent number of categories found in your data.matrix:", "\n")
	print(item.tab)
	flush.console()
	cat("\n")

	item.tab <- as.matrix(item.tab)
	if (!is.null(dimnames(item.tab))) {
		dimnames(item.tab) <- list(NULL, NULL)
		item.tab
		}

	############ inizio definizione delle matrici che saranno utilizzate
	eigen.data <- matrix(0, nvar, nrep)
	eigen.data1 <- matrix(0, nvar, nrep)
	eigen.data.pca <- matrix(0, nvar, nrep)
	eigen.data1.pca <- matrix(0, nvar, nrep)
	f1.poly.cor<- matrix(0, nvar, )
	f1.cor<- matrix(0, nvar, )
	st.matrix<- matrix(0, nvar, 15)
	st.matrix.pca<- matrix(0, nvar, 15)
	############ fine

	############ INIZIO loop che genera nrep matrici di dati casuali
	for (j in 1:nrep) {
		matrix.3 <- matrix(0, n.ss, nvar)

		######################## INIZIO loop che genera una matrice di dati discreti casuali
		u<-1
		t<-1
		while(u <= nrow(item.tab)) {
			for (z in 1:(item.tab[u,1])) {
				matrix.3[,t]<-round(runif(n.ss, min=1, max=item.tab[u,2]))
				if(t < sum(item.tab[,1])) t<-t+1
			}
			u<-u+1
		}
		######################## FINE loop che genera una matrice di dati discreti casuali

		f1.poly.cor<-poly.mat(matrix.3, short = TRUE, std.err = FALSE, ML = FALSE)
		f1.cor<-cor(matrix.3)
		eigen.data[,j] <-eigen(corFA(f1.poly.cor))$values
		eigen.data1[,j] <-eigen(corFA(f1.cor))$values
		eigen.data.pca[,j] <-eigen(cor(f1.poly.cor))$values
		eigen.data1.pca[,j] <-eigen(cor(f1.cor))$values
		if (j==1) {
			end.pt<-Sys.time()
			estimated.t<-difftime(end.pt, start.t)
			estimated.total<-estimated.t*nrep
			estimated.t<-as.difftime(estimated.t, format = "%H-%M-%S", units = "mins")
			estimated.total<-round((as.difftime(estimated.total, format = "%H %M %S", units = "mins")), digits=2)
			cat("The whole simulation will take no less than", estimated.total/60, "min.", "to terminate", "\n", "\n")
			flush.console()
		}
	}
	############ FINE loop che genera nrep matrici di dati casuali

	for (col in 1:nvar) {
		eigen.data.t <- t(eigen.data)
		st.matrix[col,1]<- mean(eigen.data.t[,col])
		st.matrix[col,2]<- sd(eigen.data.t[,col])
		st.matrix[col,3]<- quantile(eigen.data.t[,col], .95)
		st.matrix[col,4]<- quantile(eigen.data.t[,col], q.eigen)
		st.matrix[col,5]<- (col)
		eigen.data1.t <- t(eigen.data1)
		st.matrix[col,6]<- mean(eigen.data1.t[,col])
		st.matrix[col,7]<- sd(eigen.data1.t[,col])
		st.matrix[col,8]<- quantile(eigen.data1.t[,col], .95)
		st.matrix[col,9]<- quantile(eigen.data1.t[,col], q.eigen)

		eigen.data.pca.t <- t(eigen.data.pca)
		st.matrix.pca[col,1]<- mean(eigen.data.pca.t[,col])
		st.matrix.pca[col,2]<- sd(eigen.data.pca.t[,col])
		st.matrix.pca[col,3]<- quantile(eigen.data.pca.t[,col], .95)
		st.matrix.pca[col,4]<- quantile(eigen.data.pca.t[,col], q.eigen)
		st.matrix.pca[col,5]<- (col)
		eigen.data1.pca.t <- t(eigen.data1.pca)
		st.matrix.pca[col,6]<- mean(eigen.data1.pca.t[,col])
		st.matrix.pca[col,7]<- sd(eigen.data1.pca.t[,col])
		st.matrix.pca[col,8]<- quantile(eigen.data1.pca.t[,col], .95)
		st.matrix.pca[col,9]<- quantile(eigen.data1.pca.t[,col], q.eigen)
	}
	colnames(st.matrix) <- c("P.SimMeanEigen", "P.SimSDEigen",  "P.Sim95perc", "P.SimQuant", "Factor", "C.SimMeanEigen", "C.SimSDEigen",  "C.Sim95perc", "C.SimQuant", "Emp.Polyc.Eigen", "Emp.Pears.Eigen", "Eigen.diff.Polyc", "Eigen.diff.Pears", "Nr.poly.pa.fact", "Nr.pear.pa.fact")
	st.matrix
	colnames(st.matrix.pca) <- c("P.SimMeanEigen.PCA", "P.SimSDEigen.PCA",  "P.Sim95perc.PCA", "P.SimQuant.PCA", "Factor.PCA", "C.SimMeanEigen.PCA", "C.SimSDEigen.PCA",  "C.Sim95perc.PCA", "C.SimQuant.PCA", "Emp.Polyc.Eigen.PCA", "Emp.Pears.Eigen.PCA", "Eigen.diff.Polyc.PCA", "Eigen.diff.Pears.PCA", "Nr.poly.pa.fact.PCA", "Nr.pear.pa.fact.PCA")
	st.matrix.pca
	matrix.cor1 <- poly.mat(data.matrix, short = TRUE, std.err = FALSE, ML = FALSE)
	matrix.cor2 <- cor(data.matrix)

	### begin of polychor.map function
	######################## inizio funzione che calcola il Velicer MAP .
	map.polychor <- function(poly.matrix, corr.matrix) {
		loadings.map<- matrix(0, nvar, nvar)
		loadings.map1<- matrix(0, nvar, nvar)
		fm <- matrix(0, nvar, 5)
		eigen.map <- eigen(poly.matrix)
		eigen.map1 <- eigen(corr.matrix)
		loadings.map <- (eigen.map$vectors %*% (sqrt(diag(eigen.map$values))))
		loadings.map1 <- (eigen.map1$vectors %*% (sqrt(diag(eigen.map1$values))))
		fm[1,2]<- (sum(poly.matrix^2)-nvar) / (nvar*(nvar-1))
		fm[1,3]<- (sum(poly.matrix^4)-nvar) / (nvar*(nvar-1))
		fm[1,4]<- (sum(corr.matrix^2)-nvar) / (nvar*(nvar-1))
		fm[1,5]<- (sum(corr.matrix^4)-nvar) / (nvar*(nvar-1))
		for (m in 1:(nvar-1)) {
			A <- loadings.map[,1:m]
			partcov <- poly.matrix - (A %*% t(A))
			d <- diag(1 / sqrt(diag(partcov)))
			pr <- d %*% partcov %*% d
			fm[m+1, 2] <- (sum(pr^2)-nvar) / (nvar*(nvar-1))
			fm[m+1, 3] <- (sum(pr^4)-nvar) / (nvar*(nvar-1))
					### ripeto il tutto per la matrice di corr  ###
			A1 <- loadings.map1[,1:m]
			partcov1 <- corr.matrix - (A1 %*% t(A1))
			d1 <- diag(1 / sqrt(diag(partcov1)))
			pr1 <- d1 %*% partcov1 %*% d1
			fm[m+1, 4] <- (sum(pr1^2)-nvar) / (nvar*(nvar-1))
			fm[m+1, 5] <- (sum(pr1^4)-nvar) / (nvar*(nvar-1))
		}
		minfm.map <- fm[1,2]
		minfm4.map <- fm[1,3]
		minfm.map1 <- fm[1,4]
		minfm4.map1 <- fm[1,5]
		nfactors.map <- 0
		nfactors4.map <- 0
		nfactors.map1 <- 0
		nfactors4.map1 <- 0
		for (s in 1:nrow(fm)) {
			fm[s,1] <- s - 1
			if ( fm[s,2] < minfm.map ) {
				minfm.map = fm[s,2]
				nfactors.map = s - 1
			}
		}
		for (s in 1:nrow(fm)) {
			fm[s,1] <- s - 1
			if ( fm[s,3] < minfm4.map ) {
				minfm4.map = fm[s,3]
				nfactors4.map = s - 1
			}
		}
						### ripeto il tutto per la matrice di corr  ###
		for (s in 1:nrow(fm)) {
			fm[s,1] <- s - 1
			if ( fm[s,4] < minfm.map1 ) {
				minfm.map1 = fm[s,4]
				nfactors.map1 = s - 1
			}
		}
		for (s in 1:nrow(fm)) {
			fm[s,1] <- s - 1
			if ( fm[s,5] < minfm4.map1 ) {
				minfm4.map1 = fm[s,5]
				nfactors4.map1 = s - 1
			}
		}
		cat("******* RESULTS:", "\n")
		cat("# of factors (PCA) for Velicer MAP criterium (Polychoric corr): ", nfactors.map, "\n")
		cat("# of factors (PCA) for Velicer MAP(4th power)(Polychoric corr): ", nfactors4.map, "\n")
		cat("# of factors (PCA) for Velicer MAP criterium (Pearson corr)...: ", nfactors.map1, "\n")
		cat("# of factors (PCA) for Velicer MAP(4th power)(Pearson corr)...: ", nfactors4.map1, "\n")
		colnames(fm) <- c("Factor", "POLY.MAP.squared", "POLY.MAP.4th", "CORR.MAP.squared", "CORR.MAP.4th")
		max.nr <- rbind(nfactors.map, nfactors4.map, nfactors.map1, nfactors4.map1)
		return(fm[(1:(max(max.nr)+2)),])
	}
	### end of polychor.map function
	map.result<-map.polychor(poly.matrix=matrix.cor1, corr.matrix=matrix.cor2)

	st.matrix[,10]<- eigen(corFA(matrix.cor1))$values
	st.matrix[,11]<- eigen(corFA(matrix.cor2))$values
	for (riga in 1:nvar) {
		st.matrix[riga,12]<-(st.matrix[riga,10] - st.matrix[riga,4])
		st.matrix[riga,14]<-(st.matrix[riga,12]>0)
		st.matrix[riga,13]<-(st.matrix[riga,11] - st.matrix[riga,9])
		st.matrix[riga,15]<-(st.matrix[riga,13]>0)
        }
        st.matrix.pca[,10]<- eigen(cor(matrix.cor1))$values
	st.matrix.pca[,11]<- eigen(cor(matrix.cor2))$values
	for (riga in 1:nvar) {
		st.matrix.pca[riga,12]<-(st.matrix.pca[riga,10] - st.matrix.pca[riga,4])
		st.matrix.pca[riga,14]<-(st.matrix.pca[riga,12]>0)
		st.matrix.pca[riga,13]<-(st.matrix.pca[riga,11] - st.matrix.pca[riga,9])
		st.matrix.pca[riga,15]<-(st.matrix.pca[riga,13]>0)
	}
	### stabilisce il nr di fattori da ritenere per matrici polycoriche
	coin<-1
	nr_fact<-0
	righetta<-1
	while(coin==1) {
		nr_fact<-righetta
		st.matrix[righetta,14]<-(st.matrix[righetta,12]>0)
		coin<-st.matrix[righetta,14]
		righetta<-righetta+1
	}
	### stabilisce il nr di fattori da ritenere per matrici correlazioni pearson
	coin<-1
	nr_fact1<-0
	righetta<-1
	while(coin==1) {
		nr_fact1<-righetta
		st.matrix[righetta,15]<-(st.matrix[righetta,13]>0)
		coin<-st.matrix[righetta,15]
		righetta<-righetta+1
	}
	coin.pca<-1
	nr_fact.pca<-0
	righetta.pca<-1
	while(coin.pca==1) {
		nr_fact.pca<-righetta.pca
		st.matrix.pca[righetta.pca,14]<-(st.matrix.pca[righetta.pca,12]>0)
		coin.pca<-st.matrix.pca[righetta.pca,14]
		righetta.pca<-righetta.pca+1
	}

	coin.pca<-1
	nr_fact1.pca<-0
	righetta.pca<-1
	while(coin.pca==1) {
		nr_fact1.pca<-righetta.pca
		st.matrix.pca[righetta.pca,15]<-(st.matrix.pca[righetta.pca,13]>0)
		coin.pca<-st.matrix.pca[righetta.pca,15]
		righetta.pca<-righetta.pca+1
	}
#results in the plot
        plot(st.matrix[,5], st.matrix[,10], type="b", xlim=c(1, max(st.matrix[,5])), ylim=c((min(st.matrix)), (max(st.matrix[,c(4, 9, 10, 11)]))), xlab="# factors", ylab="eigenvalues", main="Parallel Analysis")
	points(st.matrix[,5], st.matrix[,4], type="b", pch=8)
	points(st.matrix[,5], st.matrix[,11], type="b", pch=20, col="red")
	points(st.matrix[,5], st.matrix[,9], type="b", pch=2, col="red")
	perc <- paste(q.eigen*100, "° perc. Polychoric corr. Sim. FA", sep="")
	perc1 <- paste(q.eigen*100,"° perc. Pearson corr. Sim. FA", sep="")
        res<-(nr_fact-1)
	res1<-(nr_fact1-1)
        res.pca<-(nr_fact.pca-1)
	res1.pca<-(nr_fact1.pca-1)
	ris <- paste("# factors with Polyc.PA: ", res, sep="")
	ris1 <- paste("# factors with Pear.PA: ", res1, sep="")
	legend(x="topright", c("Polychoric corr. Empirical FA","Pearson corr. Empirical FA",perc,perc1,ris,ris1), col=c(1,2,1,2,1,2), pch=c(1, 20, 8, 2))
	abline(h=1)
	abline(h=0)
# results in the text
	cat("# of factors (PCA) for PA method (Random Polychoric Corr.)....: ", (nr_fact.pca-1), "\n")
	cat("# of factors (PCA) for PA method (Random Pearson Corr.).......: ", (nr_fact1.pca-1), "\n")
	cat("# of factors for PA method (Random Polychoric Corr.)..........: ", (nr_fact-1), "\n")
	cat("# of factors for PA method (Random Pearson Corr.).............: ", (nr_fact1-1), "\n", "\n")
	end.t<-Sys.time()
	elapsed.t<-difftime(end.t, start.t)
	cat("Elapsed Time:", elapsed.t, "\n")
	table.result.poly<- cbind(st.matrix[(1:(res+1)), c(5, 10, 1, 2, 4)])
	table.result.corr<- cbind(st.matrix[(1:(res1+1)), c(5, 11, 6, 7, 9)])
	table.result.pca.poly<- cbind(st.matrix.pca[(1:(res.pca+1)), c(5, 10, 1, 2, 4)])
	table.result.pca.corr<- cbind(st.matrix.pca[(1:(res1.pca+1)), c(5, 11, 6, 7, 9)])
	list(MAP.selection=map.result, POLYCHORIC=table.result.poly, PEARSON=table.result.corr, PLYCHORIC.PCA=table.result.pca.poly, PEARSON.PCA=table.result.pca.corr)
}
