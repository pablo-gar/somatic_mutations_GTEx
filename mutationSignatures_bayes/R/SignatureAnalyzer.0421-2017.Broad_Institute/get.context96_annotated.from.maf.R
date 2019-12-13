get.complementary <- function(x) {
	if (x == "A") {
		y <- "T"
	} else if (x == "T") {
		y <- "A"
	} else if (x == "C") {
		y <- "G"
	} else if (x == "G") {
		y <- "C"
	}
	return(y)
}
get.reverse <- function(x) {
	n.x <- length(x)
	z <- rep(NA,4)
	w <- x
	for (i in 1:n.x) {
		y <- as.vector(unlist(strsplit(as.character(x[i]),"")))
		z[1] <- get.complementary(y[1])
		z[2] <- get.complementary(y[2])
		z[3] <- get.complementary(y[4])
		z[4] <- get.complementary(y[3])
		w[i] <- paste(z[1],z[2],z[3],z[4],sep="")
	}
	return(w)
}

context65 <- read.delim(paste(paste("context65.txt",sep=""),sep=""),header=T,sep="\t")
x <- context65[,2]
x1 <- sapply(x,function(x) strsplit(as.character(x)," in ")[[1]][1]) 
x2 <- sapply(x,function(x) strsplit(as.character(x)," in ")[[1]][2])
context65.flanking <- paste(substring(x2,1,1),substring(x2,3,3),sep="")

context96 <- read.delim(paste(paste("context96.txt",sep=""),sep=""),header=F,sep="\t")
context96 <- as.vector(unlist(context96))
context96.reverse <- get.reverse(context96)
context96.label <- c(context96.reverse[1:48],context96[49:96])
context96.label.reverse <- get.reverse(context96.label)
context96.plus <- paste(context96.label,"+",sep="")
context96.minus <- paste(context96.label,"-",sep="")
context96.reverse.plus <- paste(context96.label.reverse,"+",sep="")
context96.reverse.minus <- paste(context96.label.reverse,"-",sep="")

get.spectrum96.from.maf <- function(maf) {
	maf[,"sample"] <- maf$Tumor_Sample_Barcode
	if ("Variant_Type" %in% colnames(maf)) {
                maf <- maf[maf$Variant_Type %in% "SNP",]
        }
	ref <- toupper(maf[,colnames(maf) %in% c("Reference_Allele")])
	alt <- toupper(maf[,colnames(maf) %in% c("Tumor_Seq_Allele2")])
	context <- toupper(maf[,colnames(maf) %in% c("ref_context")])
	n.context <- length(unlist(strsplit(as.character(context[1]),"")))
	mid <- trunc(n.context/2)+1
	contig <- paste(ref,alt,substring(context,mid-1,mid-1),substring(context,mid+1,mid+1),sep="")
	context96.num <- rep(0,length(contig))
	for (i in 1:96) {
		context96.num[contig %in% c(context96[i],context96.reverse[i])] <- i
	}
	maf[,"context96.num"] <- context96.num
	maf[,"context96.word"] <- contig
	x <- as.data.frame.matrix(table(maf$context96.num,maf$sample))
	if (nrow(x)<96) {
                existing <- which(c(1:96)%in%as.numeric(rownames(x)))
                missing <- which(!c(1:96)%in%as.numeric(rownames(x)))
                n.missing <- length(missing)
                x.missing <- array(0,dim=c(n.missing,ncol(x)))
                rownames(x.missing) <- context96.label[missing]
                colnames(x.missing) <- colnames(x)
                rownames(x) <- context96.label[existing]
                x <- rbind(x,x.missing)
                x <- x[match(context96.label,rownames(x),nomatch=0),]
	} else if (nrow(x) > 96) {
		stop('Unusal context - double-check context information')
        }
	rownames(x) <- context96.label
	return(list(maf,x))
}
