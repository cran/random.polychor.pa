random.polychor.pa <- function(nvar, n.ss, nrep, nstep, data.matrix, q.eigen, r.seed="NULL") {   # definisce il nome e il numero dei parametri per iniziare la simulazione
	start.t<- Sys.time()   # prendo nota del tempo di inizio della simulazione
	data.matrix <- as.matrix(data.matrix)
	if (!is.null(dimnames(data.matrix))) {
		dimnames(data.matrix) <- list(NULL, NULL)
		data.matrix
		}
#	if (na.fail(data.matrix))) {
#		stop("Some missing values (NA) were found in the provided dataset. Please remove NA values form the dataset before running the function.")
#		}
	if (is.null(nvar)) {
		stop("number of variables is missing. Please provide the number of variables")
		}
	if (is.null(n.ss)) {
		stop("number of subjects is missing. Please provide the number of subjects")
		}
	if (is.null(nstep)) {
		stop("number of replication is missing. Please provide the number of replication")
		}
	if (is.null(q.eigen)) {
		stop("the quantile to be extracted is not declared. Please provide a valid quantile")
		}
	if (r.seed=="NULL") {  set.seed(1335031435) 	# necessario per permettere la riproducibilità della soluzione
		}
		else { set.seed(r.seed)
		}
	check <- (nrow(data.matrix)==n.ss)
	check1 <- (ncol(data.matrix)==nvar)
	if (check==0 | check1==0)
		stop("number of variables and/or subjects not equal to the entered data matrix")
	data.matrix.0<-na.omit(data.matrix)
	if (nrow(data.matrix.0) < nrow(data.matrix)) {
		cat("########################################################################", "\n")
		cat("######################## Missing Values (NA) detected", "\n")
		cat("######################## LISTWISE deletion needed. The new sample size is:", nrow(data.matrix.0), "\n")
		cat("########################################################################", "\n")
		data.matrix<-data.matrix.0
		n.ss <- nrow(data.matrix.0)
		}
	require(psych)   # carico la libreria "psych", serve per la funzione poly.mat()
	require(nFactors)   # carico la libreria "nFactors", serve per la funzione hetcor() e eigen()

	############ inizio definizione delle matrici che saranno utilizzate
	eigen.data <- matrix(0, nvar, nrep)
	eigen.data1 <- matrix(0, nvar, nrep)
	f1.poly.cor<- matrix(0, nvar, )
	f1.cor<- matrix(0, nvar, )
	st.matrix<- matrix(0, nvar, 15)
	treshold<- matrix(0, (nstep-1), 1)
	############ fine

	############ INIZIO loop che genera nrep matrici di dati casuali
	for (j in 1:nrep) {
		matrix.1 <- matrix(0, n.ss, nvar)
		matrix.2 <- matrix(0, n.ss, nvar)
		matrix.3 <- matrix(0, n.ss, nvar)

		############	INIZIO loop che genera una matrice di dati categoriali casuali
		for (i in 1:nvar) {
			matrix.1[,i]<-rnorm(n.ss)   # determina n valori random (per ogni nvar) distribuiti normalmente
			matrix.2[,i]<-pnorm(matrix.1[,i])   # determina la p degli n valori

			############	INIZIO programma che calcola soglie e step degli n*nvar valori e li trasforma in numeri discreti (i.e., 1, 2 ..., nstep)
			for (t in 1:(nstep-1)) treshold[t,]<- (t/nstep) # determina il valore(t/nstep) e il numero (nstep-1) delle soglie della scala di risposta dell'item in funzione degli step
				treshold
			for (b in 1:nvar) {   # loop che ha l'obiettivo di riempire due matrici (eigen.data, eigen.data1) con gli autovalori delle matrici casuali generate
				for (a in 1:n.ss) {   # per ogni riga di matrix.2
					tampone<-matrix(0,(nstep-1),1)  # definisce una matrice che serve a registrare momentaneamente se il valore di matrix.2 è sopra o sotto la soglia
					trick=2   # imposta un valore di default per la variabile
					for (s in 1:(nstep-1)) {   # loop che serve a stabilire per ogni valore di matrix.2 se è sopra o sotto soglia
		                               # funzionamento: il valore di matrix.2 viene confrontato con ogni soglia per capire
		                               # a quale soglia si avvicina di più
						tampone[s,1]<- (treshold[s,1]-matrix.2[a,b])  # confronta la soglia con il valore contenuto in matrix.2
#						cat(tampone[s,1], "\n")   # stampa i valori di tampone per capire se ci sono errori nel calcolo (da eliminare)
					}
#				cat(which.min(abs(tampone[,1])), "\n")   # stampa a video il il numero di riga che contiene il valore (in assoluto) più vicino a 0
				trick<-(tampone[which.min(abs(tampone[,1])),1]<0)   # valuta se il valore (all'interno di tampone) più vicino a 0 è negativo
				if(trick==1) {   # se il valore è negativo allora significa che per determinare lo step corretto è necessario aggiungere 1
					matrix.3[a,b]<-(which.min(abs(tampone[,1]))+1)   # aggiunge uno allo step
					}
					else { matrix.3[a,b]<-(which.min(abs(tampone[,1])))
					}   # altrimenti lascia lo step così
				}
			}
			############	FINE programma che calcola soglie e step
		}
		############	FINE loop che genera una matrice di dati categoriali casuali

		f1.poly.cor<-poly.mat(matrix.3, short = TRUE, std.err = FALSE, ML = FALSE) # matrice di correlazion di Polycorica
		f1.cor<-cor(matrix.3) # matrice di correlazione di Pearson
		eigen.data[,j] <-eigen(corFA(f1.poly.cor))$values  # estrae gli eigenvalues dalla matrice di correlazion Polycorica
		eigen.data1[,j] <-eigen(corFA(f1.cor))$values  # estrae gli eigenvalues dalla matrice di correlazion di Pearson
		if (j==1) {
			end.pt<-Sys.time()   # fine del tempo necessario ai calcoli
			estimated.t<-difftime(end.pt, start.t)   # calcola il tempo impiegato
			estimated.total<-estimated.t*nrep   # calcola il tempo impiegato
			estimated.t<-as.difftime(estimated.t, format = "%H-%M-%S", units = "secs")
			estimated.total<-round((as.difftime(estimated.total, format = "%H %M %S", units = "mins")/60), digits=0)
			cat("The whole simulation will take no less than", estimated.total, "min.", "to terminate", "\n")   # stampa il tempo impiegato
#			cat("The whole simulation will terminate not before:", format(start.t, "%X"), "+", estimated.total/60, "min.", "\n")   # stampa il tempo impiegato
		}
	}
	############ FINE loop che genera nrep matrici di dati casuali

	for (col in 1:nvar) {   # loop che riempie la matrice
		eigen.data.t <- t(eigen.data) # trasposta della matrice
		st.matrix[col,1]<- mean(eigen.data.t[,col])   # matrice polycorica: dato il primo autovalore estratto, calcola la media tra tutti gli nrep autovalori estratti come primi (e così via per il II autovalore estratto ...)
		st.matrix[col,2]<- sd(eigen.data.t[,col])   # matrice polycorica: stesso ma calcola la deviazione standard per il I autovalore estratto, ecc.
		st.matrix[col,3]<- quantile(eigen.data.t[,col], .95)   # matrice polycorica: calcola il 95° percentile del I autovalore estratto
		st.matrix[col,4]<- quantile(eigen.data.t[,col], q.eigen)   # matrice polycorica: calcola il 99° percentile del I autovalore estratto
		st.matrix[col,5]<- (col)
		eigen.data1.t <- t(eigen.data1) # trasposta della matrice
		st.matrix[col,6]<- mean(eigen.data1.t[,col])   # matrice correlazione: dato il primo autovalore estratto, calcola la media tra tutti gli nrep autovalori estratti come primi (e così via per il II autovalore estratto ...)
		st.matrix[col,7]<- sd(eigen.data1.t[,col])   # matrice correlazione: stesso ma calcola la deviazione standard per il I autovalore estratto, ecc.
		st.matrix[col,8]<- quantile(eigen.data1.t[,col], .95)   # matrice correlazione: calcola il 95° percentile del I autovalore estratto
		st.matrix[col,9]<- quantile(eigen.data1.t[,col], q.eigen)   # matrice correlazione: calcola il 99° percentile del I autovalore estratto
	}
	colnames(st.matrix) <- c("P.SimMeanEigen", "P.SimSDEigen",  "P.Sim95perc", "P.SimQuant", "Factor", "C.SimMeanEigen", "C.SimSDEigen",  "C.Sim95perc", "C.SimQuant", "Emp.Polyc.Eigen", "Emp.Pears.Eigen", "Eigen.diff.Polyc", "Eigen.diff.Pears", "Nr.poly.pa.fact", "Nr.pear.pa.fact")
	st.matrix   # mostra avideo la matrice della parallel analysis
	matrix.cor1 <- poly.mat(data.matrix, short = TRUE, std.err = FALSE, ML = FALSE)
	matrix.cor2 <- cor(data.matrix)

	### begin of polychor.map function
	######################## inizio funzione che calcola il Velicer MAP .
	map.polychor <- function(poly.matrix, corr.matrix) {   # definisce il nome e il numero dei parametri per iniziare la simulazione
		loadings.map<- matrix(0, nvar, nvar)
		loadings.map1<- matrix(0, nvar, nvar)
		fm <- matrix(0, nvar, 5)
		eigen.map <- eigen(poly.matrix)  # estrae gli eigenvalues e gli eigenvectors dalla matrice di correlazion Polycorica
		eigen.map1 <- eigen(corr.matrix)  # estrae gli eigenvalues e gli eigenvectors dalla matrice di correlazion Polycorica
		loadings.map <- (eigen.map$vectors %*% (sqrt(diag(eigen.map$values))))
		loadings.map1 <- (eigen.map1$vectors %*% (sqrt(diag(eigen.map1$values))))
		fm[1,2]<- (sum(poly.matrix^2)-nvar) / (nvar*(nvar-1))
		fm[1,3]<- (sum(poly.matrix^4)-nvar) / (nvar*(nvar-1))
		fm[1,4]<- (sum(corr.matrix^2)-nvar) / (nvar*(nvar-1))
		fm[1,5]<- (sum(corr.matrix^4)-nvar) / (nvar*(nvar-1))
		for (m in 1:(nvar-1)) {   # loop che riempie la matrice
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
		for (s in 1:nrow(fm)) {   # loop che riempie la matrice
			fm[s,1] <- s - 1
			if ( fm[s,2] < minfm.map ) {
				minfm.map = fm[s,2]
				nfactors.map = s - 1
			}
		}
		for (s in 1:nrow(fm)) {   # loop che riempie la matrice
			fm[s,1] <- s - 1
			if ( fm[s,3] < minfm4.map ) {
				minfm4.map = fm[s,3]
				nfactors4.map = s - 1
			}
		}
						### ripeto il tutto per la matrice di corr  ###
		for (s in 1:nrow(fm)) {   # loop che riempie la matrice
			fm[s,1] <- s - 1
			if ( fm[s,4] < minfm.map1 ) {
				minfm.map1 = fm[s,4]
				nfactors.map1 = s - 1
			}
		}
		for (s in 1:nrow(fm)) {   # loop che riempie la matrice
			fm[s,1] <- s - 1
			if ( fm[s,5] < minfm4.map1 ) {
				minfm4.map1 = fm[s,5]
				nfactors4.map1 = s - 1
			}
		}
		cat("# of factors for Velicer MAP criterium (Polychoric corr): ", nfactors.map, "\n")
		cat("# of factors for Velicer MAP(4th power)(Polychoric corr): ", nfactors4.map, "\n")
		cat("# of factors for Velicer MAP criterium (Pearson corr)   : ", nfactors.map1, "\n")
		cat("# of factors for Velicer MAP(4th power)(Pearson corr)   : ", nfactors4.map1, "\n")
		colnames(fm) <- c("Factor", "POLY.MAP.squared", "POLY.MAP.4th", "CORR.MAP.squared", "CORR.MAP.4th")
		fm[,1] <- 1+fm[,1]
		max.nr <- rbind(nfactors.map, nfactors4.map, nfactors.map1, nfactors4.map1)
		return(fm[(1:(max(max.nr)+1)),])
	}
	### end of polychor.map function
	map.result<-map.polychor(poly.matrix=matrix.cor1, corr.matrix=matrix.cor2)

	st.matrix[,10]<- eigen(corFA(matrix.cor1))$values   # aggiunge gli autovalori della matrice policorica dei dati empirici
	st.matrix[,11]<- eigen(corFA(matrix.cor2))$values   # aggiunge gli autovalori della matrice correlazioni Pearson dei dati empirici
	for (riga in 1:nvar) {   # loop che riempie la matrice
		st.matrix[riga,12]<-(st.matrix[riga,10] - st.matrix[riga,4])
		st.matrix[riga,14]<-(st.matrix[riga,12]>0)
		st.matrix[riga,13]<-(st.matrix[riga,11] - st.matrix[riga,9])
		st.matrix[riga,15]<-(st.matrix[riga,13]>0)
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
	plot(st.matrix[,5], st.matrix[,10], type="b", xlim=c(1, max(st.matrix[,5])), ylim=c((min(st.matrix)), (max(st.matrix[,c(4, 9, 10, 11)]))), xlab="# factors", ylab="eigenvalues", main="Parallel Analysis")
	points(st.matrix[,5], st.matrix[,4], type="b", pch=8)
	points(st.matrix[,5], st.matrix[,11], type="b", pch=20, col="red")
	points(st.matrix[,5], st.matrix[,9], type="b", pch=2, col="red")
	perc <- paste(q.eigen*100, "° perc. Polychoric corr. Sim. FA", sep="")
	perc1 <- paste(q.eigen*100,"° perc. Pearson corr. Sim. FA", sep="")
	res<-(nr_fact-1)
	res1<-(nr_fact1-1)
	ris <- paste("# factors with Polyc.PA: ", res, sep="")
	ris1 <- paste("# factors with Pear.PA: ", res1, sep="")
	legend(x="topright", c("Polychoric corr. Empirical FA","Pearson corr. Empirical FA",perc,perc1,ris,ris1), col=c(1,2,1,2,1,2), pch=c(1, 20, 8, 2))
	abline(h=1)
	abline(h=0)
	cat("# of factors for PA method (Random Polychoric Corr.)    : ", (nr_fact-1), "\n")   # stampa il nr di fattori da ritenere
	cat("# of factors for PA method (Random Pearson Corr.)       : ", (nr_fact1-1), "\n")   # stampa il nr di fattori da ritenere
	end.t<-Sys.time()   # fine del tempo necessario ai calcoli
	elapsed.t<-difftime(end.t, start.t)   # calcola il tempo impiegato
	cat("Elapsed Time:", elapsed.t, "\n")   # stampa il tempo impiegato
	super.matrix <- merge(st.matrix, map.result, by="Factor")
	table.result.poly<- cbind(st.matrix[(1:(res+1)), c(5, 10, 1, 2, 4)])
	table.result.corr<- cbind(st.matrix[(1:(res1+1)), c(5, 11, 6, 7, 9)])
	list(MAP.selection=map.result, POLYCHORIC=table.result.poly, PEARSON=table.result.corr)
}

