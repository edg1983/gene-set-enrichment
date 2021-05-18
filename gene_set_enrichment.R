##Load results and eventually select desired genes based on significance level
mydata <- read.table(input.file, sep=",", header=T, stringsAsFactors=F)
genes <- mydata$geneName[mydata$F.qvalue < 0.05]

##1.	Load pathways (suppose .gmt file from MSigDB)
inputfile <- "/data0/DATA/MSigDB/c5.bp.v6.0.symbols.gmt"
pathways <- list()
con  <- file(inputfile, open = "r")
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
myvector <- (strsplit(oneLine, "\t"))
pathways[[myvector[[1]][1]]] <- myvector[[1]][3:length(myvector[[1]])][myvector[[1]][3:length(myvector[[1]])] != ""]
}
close(con)

##3.	Load background genes (usually are defined as whole gene tested in results)
backgenes <- unique(mydata$geneName)

##4.	Compute hypergeometric test enrichment
results <- data.frame(category=character(), pvalue=numeric(), results_in_path=numeric(), backgenes_in_path=numeric(), tot_results=numeric(), tot_backgenes=numeric(), stringsAsFactors=FALSE)
	for (h in 1:length(pathways)) {
		successes <- length(which(genes %in% pathways[[h]]))
		tot_path <- length(which(backgenes %in% pathways[[h]]))
		not_path <- length(backgenes) - tot_path
		results[nrow(results)+1,] <- c(names(pathways[h]), phyper(successes-1, tot_path, not_path, length(genes), lower.tail=FALSE), successes, tot_path, length(genes), length(backgenes) )
	}

results$backgenes_in_path <- as.numeric(results$backgenes_in_path)
results$results_in_path <- as.numeric(results$results_in_path)

#filter out categories with less than 10 genes present in background genes
results.filt <- results[results$backgenes_in_path >= 10,] 

#Multiple test correction
results.filt$qvalue <- p.adjust(results.filt$pvalue, method="fdr")

#Save results order by adj P
write.table(results.filt[order(results.filt$qvalue),], file=output.file, sep="\t", row.names=F, quote=F)
