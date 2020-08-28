###### make_annotation 
## make annotation will return a set of GRanges accordingly to the 
## parameters given to the function. Each GRange represents a single gene. 
## In case more than one transcript is annotated for a gene, all transcripts are
## merged in order to obtain the largest region possible associated to 
## that gene. In case a gene has non overlapping transcripts associated, the 
## gene is removed.
## Arguments:
##   - "annotation" can be either "singleChr" or "genome"
##   - "singleChrName" must specify chromosome name (in case "singleChr" is selected)
##   - "annotationFeature" a character can be either "TSS", "GB", "TES"
##       > "TSS" returns the TSS region defined by the coordinates defined in 
##			promoter_limits
##       > "TES" returns the TES region defined by the coordinates defined in 
##			promoter_limits
##       > "GB" returns the region between the TSS and the TES
##       > "tx" returns the whole transcript, extended in the promoter and termination
##         regions, according to promoter_limits and termination_limits definitions
##   - "promoter_limits" a numeric of length 2, defines the upstream and
## 		downstream distances from the annotated transcription-start-site
##   - "termination_limits" a numeric of length 2, defines the upstream and
## 		downstream distances from the annotated transcription-end-site
##   - "txdb" a TxDb object
make_annotation <- function(annotation, singleChrName, annotationFeature, 
	promoter_limits=c(50, 300), termination_limits=c(1000, 4000),
	txdb=TxDb.Mmusculus.UCSC.mm9.knownGene) {

	promoter_start <- promoter_limits[1]
	promoter_end   <- promoter_limits[2]

	termination_start <- termination_limits[1]
	termination_end   <- termination_limits[2]

	message('Creating annotation...')
	reduceLevels <- function(grange) {
		grange@seqnames <- droplevels(grange@seqnames)
		grange@seqinfo <- grange@seqinfo[levels(grange@seqnames)]
		grange }
	# library(TxDb.Mmusculus.UCSC.mm9.knownGene)
	# txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
	tx <- transcripts(txdb)
	names(tx) <- select(txdb, keys = tx$tx_name
		, keytype = "TXNAME", columns = "GENEID")$GENEID
	if( annotation == 'singleChr' )
	{######## single chromosome
		tx <- reduceLevels(tx[seqnames(tx)==singleChrName])
	}
	# remove NA named genes
	tx <- tx[!is.na(names(tx))]
	######## select a unique feature per each gene
	tx <- reduce(split(tx, names(tx)))
	# remove genes which still have more than one transcript after reducing...
	tx <- unlist(tx[elementNROWS(tx)==1])
	annotationGR <- tx
	######## create either gene body or promoter annotation		
	annotationGR <- switch(annotationFeature
		# for GB, reduce and unlist the GRList split by the name of the gene
		, 'GB'={
			# remove the promoter region from the gene body
			#### PLUS
			start(annotationGR[strand(tx)=='+']) <- # min(TTS+prom,TES)
				apply(cbind(
					start(tx[strand(tx)=='+'])+(promoter_end+1)
					, end(tx[strand(tx)=='+'])+1
					), 1, min)
			end(annotationGR[strand(tx)=='+']) <- # max(TTS+prom,TES-term)
				apply(cbind(
					start(tx[strand(tx)=='+'])+(promoter_end+1)
					, end(tx[strand(tx)=='+'])-termination_start
					), 1, max)
			#### MINUS
			start(annotationGR[strand(tx)=='-']) <- 
				apply(cbind(
					end(tx[strand(tx)=='-'])-promoter_end-1
					, start(tx[strand(tx)=='-'])+termination_start
					), 1, min)
			end(annotationGR[strand(tx)=='-']) <- 
				apply(cbind(
					end(tx[strand(tx)=='-'])-(promoter_end+1)
					, start(tx[strand(tx)=='-'])-1
					), 1, max)
			trim(annotationGR)
		}
		, 'TSS'={
			# keep only the promoter region (plus strand)
			start(annotationGR[strand(tx)=='+']) <- 
				start(tx[strand(tx)=='+'])-promoter_start
			end(annotationGR[strand(tx)=='+']) <- 
				apply(cbind(
					start(tx[strand(tx)=='+'])+promoter_end
					, end(tx[strand(tx)=='+'])
					), 1, min)
			# keep only the promoter region (minus strand)
			end(annotationGR[strand(tx)=='-']) <- 
				end(tx[strand(tx)=='-'])+promoter_start
			start(annotationGR[strand(tx)=='-']) <- 
				apply(cbind(
					end(tx[strand(tx)=='-'])-promoter_end
					, start(tx[strand(tx)=='-'])
					), 1, max)
			trim(annotationGR)
		}
		, 'TES'={
			# keep only the promoter region (plus strand)
			end(annotationGR[strand(tx)=='+']) <- 
				end(tx[strand(tx)=='+'])+termination_end
			start(annotationGR[strand(tx)=='+']) <- 
				apply(cbind(
					start(tx[strand(tx)=='+'])+promoter_end+2
					, end(tx[strand(tx)=='+'])-(termination_start-1)
					), 1, max)
			# keep only the promoter region (minus strand)
			start(annotationGR[strand(tx)=='-']) <- 
				start(tx[strand(tx)=='-'])-termination_end
			end(annotationGR[strand(tx)=='-']) <- 
				apply(cbind(
					end(tx[strand(tx)=='-'])-promoter_end-2
					, start(tx[strand(tx)=='-'])+(termination_start-1)
					), 1, min)
			trim(annotationGR)
		}, 'tx'={
			# keep only the promoter region (plus strand)
			start(annotationGR[strand(tx)=='+']) <- 
				start(tx[strand(tx)=='+'])-promoter_start
			end(annotationGR[strand(tx)=='+']) <- 
				end(tx[strand(tx)=='+'])+termination_end
			# keep only the promoter region (minus strand)
			end(annotationGR[strand(tx)=='-']) <- 
				end(tx[strand(tx)=='-'])+promoter_start
			start(annotationGR[strand(tx)=='-']) <- 
				start(tx[strand(tx)=='-'])-termination_end
			trim(annotationGR)
		}
		)

	return(annotationGR)

}

###### make_transcripts_annotation (modify)
## make annotation will return a set of GRanges accordingly to the 
## parameters given to the function. Each GRange represents a single gene. 
## In case more than one transcript is annotated for a gene, all transcripts are
## merged in order to obtain the largest region possible associated to 
## that gene. In case a gene has non overlapping transcripts associated, the 
## gene is removed.
## Arguments:
##   - "annotation" can be either "singleChr" or "genome"
##   - "singleChrName" must specify chromosome name (in case "singleChr" is selected)
##   - "annotationFeature" a character can be either "TSS", "GB", "TES"
##       > "TSS" returns the TSS region defined by the coordinates defined in 
##			promoter_limits
##       > "TES" returns the TES region defined by the coordinates defined in 
##			promoter_limits
##       > "GB" return the region between the TSS and the TES
##   - "promoter_limits" a numeric of length 2, defines the upstream and
## 		downstream distances from the annotated transcription-start-site
##   - "termination_limits" a numeric of length 2, defines the upstream and
## 		downstream distances from the annotated transcription-end-site
##   - "txdb" a TxDb object
make_transcripts_annotation <- function(txdb) {
	print('Creating the annotation...')
	################
	## define 5'UTR, exons, introns, 3'UTR
	#################
	myGenetx <- reduce(transcriptsBy(txdb,'gene'))
	myGenetx <- unlist(myGenetx[elementNROWS(myGenetx)==1])
	myGenetxPlus <- myGenetx[strand(myGenetx)=='+']
	myGenetxMinus <- myGenetx[strand(myGenetx)=='-']
	unambiguousTx <- names(myGenetx)
	####################
	## 5'UTR - exons ######
	###################
	print("1/6 : 5'UTR/exons...")
	fiveUtrGR <- fiveUTRsByTranscript(txdb, use.names=TRUE)
	fiveUtrMatch <- select(txdb, keys = names(fiveUtrGR), keytype = "TXNAME", columns = "GENEID")
	#### annotate the 5utr db (plus strand)
	myFiveUtrDB <- fiveUtrMatch[fiveUtrMatch$GENEID%in%names(myGenetxPlus),]
	myFiveUtrDB$start <- sapply(start(fiveUtrGR[myFiveUtrDB$TXNAME]),min)
	myFiveUtrDB$exons <- sapply(start(fiveUtrGR[myFiveUtrDB$TXNAME]),length)
	myFiveUtrDB$txstart <- start(myGenetxPlus[myFiveUtrDB$GENEID])
	# select the 5'UTR that match the start of the annotation
	myFiveUtrDB <- myFiveUtrDB[myFiveUtrDB$start == myFiveUtrDB$txstart,]
	# remove the identical 5'UTR annotations
	myFiveUtrDB <- myFiveUtrDB[!duplicated(myFiveUtrDB[,-1]),]
	# keep only the most complete annotation
	myFiveUtrDB <- do.call('rbind',lapply(split(myFiveUtrDB,myFiveUtrDB$GENEID), function(x) x[which.max(x$exons),]))
	# create the GRanges
	myFiveUtrGRplus <- fiveUtrGR[myFiveUtrDB$TXNAME]
	names(myFiveUtrGRplus) <- myFiveUtrDB$GENEID
	#### annotate the 5utr db (minus strand)
	myFiveUtrDB <- fiveUtrMatch[fiveUtrMatch$GENEID%in%names(myGenetxMinus),]
	myFiveUtrDB$start <- sapply(end(fiveUtrGR[myFiveUtrDB$TXNAME]),max)
	myFiveUtrDB$exons <- sapply(end(fiveUtrGR[myFiveUtrDB$TXNAME]),length)
	myFiveUtrDB$txstart <- end(myGenetxMinus[myFiveUtrDB$GENEID])
	# select the 5'UTR that match the start of the annotation
	myFiveUtrDB <- myFiveUtrDB[myFiveUtrDB$start == myFiveUtrDB$txstart,]
	# remove the identical 5'UTR annotations
	myFiveUtrDB <- myFiveUtrDB[!duplicated(myFiveUtrDB[,-1]),]
	# keep only the most complete annotation
	myFiveUtrDB <- do.call('rbind',lapply(split(myFiveUtrDB,myFiveUtrDB$GENEID), function(x) x[which.max(x$exons),]))
	# create the GRanges
	myFiveUtrGRminus <- fiveUtrGR[myFiveUtrDB$TXNAME]
	names(myFiveUtrGRminus) <- myFiveUtrDB$GENEID
	#### concatenate
	myFiveUtrGRExons <- c(myFiveUtrGRplus, myFiveUtrGRminus)
	#######################
	## 5'UTR - introns ######
	#####################
	print("2/6 : 5'UTR/introns...")
	myFiveUtrGRIntrons <- psetdiff(unlist(range(myFiveUtrGRExons)),myFiveUtrGRExons)
	myFiveUtrGRIntrons <- myFiveUtrGRIntrons[elementNROWS(myFiveUtrGRIntrons)>0]
	####################
	## 3'UTR - exons ######
	###################
	print("3/6 : 3'UTR/exons...")
	threeUtrGR <- threeUTRsByTranscript(txdb, use.names=TRUE)
	threeUtrMatch <- select(txdb, keys = names(threeUtrGR), keytype = "TXNAME", columns = "GENEID")
	#### annotate the 5utr db (plus strand)
	myThreeUtrDB <- threeUtrMatch[threeUtrMatch$GENEID%in%names(myGenetxPlus),]
	myThreeUtrDB$end <- sapply(end(threeUtrGR[myThreeUtrDB$TXNAME]),max)
	myThreeUtrDB$exons <- sapply(end(threeUtrGR[myThreeUtrDB$TXNAME]),length)
	myThreeUtrDB$txend <- end(myGenetxPlus[myThreeUtrDB$GENEID])
	# select the 5'UTR that match the end of the annotation
	myThreeUtrDB <- myThreeUtrDB[myThreeUtrDB$end == myThreeUtrDB$txend,]
	# remove the identical 5'UTR annotations
	myThreeUtrDB <- myThreeUtrDB[!duplicated(myThreeUtrDB[,-1]),]
	# keep only the most complete annotation
	myThreeUtrDB <- do.call('rbind',lapply(split(myThreeUtrDB,myThreeUtrDB$GENEID), function(x) x[which.max(x$exons),]))
	# create the GRanges
	myThreeUtrGRplus <- threeUtrGR[myThreeUtrDB$TXNAME]
	names(myThreeUtrGRplus) <- myThreeUtrDB$GENEID
	#### annotate the 5utr db (minus strand)
	myThreeUtrDB <- threeUtrMatch[threeUtrMatch$GENEID%in%names(myGenetxMinus),]
	myThreeUtrDB$start <- sapply(start(threeUtrGR[myThreeUtrDB$TXNAME]),min)
	myThreeUtrDB$exons <- sapply(start(threeUtrGR[myThreeUtrDB$TXNAME]),length)
	myThreeUtrDB$txend <- start(myGenetxMinus[myThreeUtrDB$GENEID])
	# select the 5'UTR that match the start of the annotation
	myThreeUtrDB <- myThreeUtrDB[myThreeUtrDB$start == myThreeUtrDB$txend,]
	# remove the identical 5'UTR annotations
	myThreeUtrDB <- myThreeUtrDB[!duplicated(myThreeUtrDB[,-1]),]
	# keep only the most complete annotation
	myThreeUtrDB <- do.call('rbind',lapply(split(myThreeUtrDB,myThreeUtrDB$GENEID), function(x) x[which.max(x$exons),]))
	# create the GRanges
	myThreeUtrGRminus <- threeUtrGR[myThreeUtrDB$TXNAME]
	names(myThreeUtrGRminus) <- myThreeUtrDB$GENEID
	#### concatenate
	myThreeUtrGRExons <- c(myThreeUtrGRplus, myThreeUtrGRminus)
	#######################
	## 3'UTR - introns ######
	#####################
	print("4/6 : 3'UTR/introns...")
	myThreeUtrGRIntrons <- psetdiff(unlist(range(myThreeUtrGRExons)),myThreeUtrGRExons)
	myThreeUtrGRIntrons <- myThreeUtrGRIntrons[elementNROWS(myThreeUtrGRIntrons)>0]
	######################
	## Coding - regions ######
	#######################
	myCodingRegions <- myGenetx
	myFiveUtrGenes <- names(myFiveUtrGRExons)
	myCodingRegions[myFiveUtrGenes] <- psetdiff(myCodingRegions[myFiveUtrGenes],unlist(range(myFiveUtrGRExons)))
	myThreeUtrGenes <- names(myThreeUtrGRExons)
	myCodingRegions[myThreeUtrGenes] <- psetdiff(myCodingRegions[myThreeUtrGenes],unlist(range(myThreeUtrGRExons)))
	####################
	## Coding - exons ######
	###################
	print('5/6 : Creating Coding/exons...')
	myExonsGR <- disjoin(exonsBy(txdb,'gene'))[names(myCodingRegions)]
	myCDSExonsGR <- pintersect(myExonsGR,myCodingRegions,drop.nohit.ranges=TRUE)
	######################
	## Coding - introns ######
	#######################
	print('6/6 : Creating Coding/introns...')
	myCDSIntronsGR <- psetdiff(unlist(range(myCDSExonsGR)),myCDSExonsGR)
	myCDSIntronsGR <- myCDSIntronsGR[elementNROWS(myCDSIntronsGR)>0]
	############
	## Gather ######
	#############
	print('Creating empty GRangesList... (This step could be long, up to 2-3 minutes)')
	emptyGRList <- GRangesList(sapply(unambiguousTx, function(x) GRanges()))
	#    user  system elapsed 
	# 139.860   0.000 139.869
	#+ 31.696   0.420  32.126 
	threeUTRintrons <- threeUTRexons <- CDSintrons <- CDSexons <- fiveUTRintrons <- fiveUTRexons <- emptyGRList
	fiveUTRexons[names(myFiveUtrGRExons)] <- myFiveUtrGRExons
	fiveUTRintrons[names(myFiveUtrGRIntrons)] <- myFiveUtrGRIntrons
	CDSexons[names(myCDSExonsGR)] <- myCDSExonsGR
	CDSintrons[names(myCDSIntronsGR)] <- myCDSIntronsGR
	threeUTRexons[names(myThreeUtrGRExons)] <- myThreeUtrGRExons
	threeUTRintrons[names(myThreeUtrGRIntrons)] <- myThreeUtrGRIntrons
	txAnnotation <- list(
		"fiveUTRexons"=fiveUTRexons,
		"fiveUTRintrons"=fiveUTRintrons,
		"CDSexons"=CDSexons,
		"CDSintrons"=CDSintrons,
		"threeUTRexons"=threeUTRexons,
		"threeUTRintrons"=threeUTRintrons
		)
	print('Done.')
	return(txAnnotation)
}

###### gene length
## returns the length of the genomic region region associated to each gene, 
## as defined in "make_annotation" function
## Arguments:
##   - "annotation" can be either "singleChr" or "genome"
##   - "singleChrName" must specify chromosome name (in case "singleChr" is selected)
gene_length <- function(annotation, singleChrName) {

	message('Creating annotation...')
	reduceLevels <- function(grange) {
		grange@seqnames <- droplevels(grange@seqnames)
		grange@seqinfo <- grange@seqinfo[levels(grange@seqnames)]
		grange }
	library(TxDb.Mmusculus.UCSC.mm9.knownGene)
	txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
	tx <- transcripts(txdb)
	names(tx) <- select(txdb, keys = tx$tx_name
		, keytype = "TXNAME", columns = "GENEID")$GENEID
	if( annotation == 'singleChr' )
	{######## single chromosome
		tx <- reduceLevels(tx[seqnames(tx)==singleChrName])
	}
	# remove NA named genes
	tx <- tx[!is.na(names(tx))]
	######## select a unique feature per each gene
	tx <- reduce(split(tx, names(tx)))
	# remove genes which still have more than one transcript after reducing...
	tx <- unlist(tx[elementNROWS(tx)==1])
	gl <- width(tx)
	names(gl) <- names(tx)
	return(gl)
}