#Modul: Analisis Vaksin Dengue 
#Dataset: GSE18090 (Dengue)
#Platform: Microarray (Affymetrix GPL96)
#Tujuan: Mengidentifikasi Differentially Expressed Genes (DEG) 

#PART A. PENGANTAR KONSEP 

#Analisis ekspresi gen bertujuan untuk membandingkan tingkat ekspresi gen 
#antara dua kondisi biologis (misalnya penderita dengue vs normal) 
#Pada modul ini kita menggunakan pendekatan statistik limma (Linear Models
#for Microarray Data), yang merupakan standar emas untuk data microarray. 

#PART B. PERSIAPAN LINGKUNGAN KERJA (INSTALL & LOAD PACKAGE) 

#Apa itu package? 
#Package adalah kumpulan fungsi siap pakai di R
#Bioinformatika di R sangat bergantung pada package dari CRAN dan Bioconductor 

#1. Install BiocManager (manajer paket Bioconductor) 
#IF adalah struktur logika : “jika kondisi terpenuhi, lakukan aksi”

if (!require("BiocManager", quietly = TRUE))  {
  install.packages("BiocManager") 
}

# 2. Install paket Bioconductor (GEOquery & limma) 
#GEOquery: mengambil data dari database GEO 
#limma: analisis statistik ekspresi gen 

BiocManager::install(c("GEOquery", "limma"), ask = FALSE, update = FALSE) 

#Install annotation package sesuai platform
#GPL96 = Affymetrix Human Genome U133A 
BiocManager::install("hgu133a.db", ask = FALSE, update = FALSE)

#3. Install paket CRAN untuk visualisasi dan manipulasi data 
#phetmap: heatmap ekspresi gen 
#ggplot2: grafik (volcano plot)
#dplyr: manipulasi tabel data 

install.packages(c("pheatmap", "ggplot2", "dplyr"))

#umap: grafik (plot UMAP) 
if (!requireNamespace("umap", quietly = TRUE)) {
  install.packages("umap")
}

#4. Memanggil library 
#library() digunakan agar fungsi di dalam package bisa digunakan
#   Differential expression analysis dengan limma
library(GEOquery)
library(limma)
library(umap)

#PART C. PENGAMBILAN DATA DARI GEO 

#GEO (Gene Expression Omnibus) adalah database publik milik NCBI
#getGEO(): fungsi untuk mengunduh dataset berdasarkan ID GEO
#GSEMatrix = TRUE -> data diambil dalam format ExpressionSet
#AnnotGPL  = TRUE -> anotasi gen (Gene Symbol) ikut diunduh
gset <- getGEO("GSE18090", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# Ikuti name kolom agar sesuai dengan variabel label top tabel
fvarLabels(gset) <- make.names(fvarLabels(gset))

gsms <- "00000000000000000011111111"
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
#ExpressionSet berisi:
# - exprs() : matriks ekspresi gen
# - pData() : metadata sampel
# - fData() : metadata fitur (probe / gen)

#PART D. PRE-PROCESSING DATA EKSPRESI 

# exprs(): mengambil matriks ekspresi gen
# Baris  = probe/gen
# Kolom  = sampel
# group membership for all samples
ex <- exprs(gset)
#Mengapa perlu log2 transformasi?
#Data microarray mentah memiliki rentang nilai sangat besar.
#Log2 digunakan untuk:
#1. Menstabilkan varians
#2. Mendekati asumsi model linear
#3. Memudahkan interpretasi log fold change

#quantile(): menghitung nilai kuantil (persentil)
#as.numeric(): mengubah hasil quantile (yang berupa named vector)
#menjadi vektor numerik biasa agar mudah dibandingkan
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))

#LogTransform adalah variabel logika (TRUE / FALSE)
#Operator logika:
#>  : lebih besar dari
#| | : OR (atau)
#&& : AND (dan)
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }

#PART E. DEFINISI KELOMPOK SAMPEL 
# Masukkan kelompok sampel dan definisikan (pembatasan)
gs <- factor(sml)
groups <- make.names(c("Patient","Normal"))
levels(gs) <- groups
gset$group <- gs

#PART F. DESIGN MATRIX (KERANGKA STATISTIK) 

#model.matrix():
#Membuat matriks desain untuk model linear
#~0 berarti TANPA intercept (best practice limma)
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # hiraukan missing values

# Menenukan precision weights dan menunjukkan plot dari mean-variance trend
v <- vooma(gset, design, plot=T)
# OR weights by group
# v <- voomaByGroup(gset, group=groups, design, plot=T, cex=0.1, pch=".", col=1:nlevels(gs))
v$genes <- fData(gset) # Memberi anotasi gen

#PART G. ANALISIS DIFFERENTIAL EXPRESSION (LIMMA)

# fit linear model
fit  <- lmFit(v)

# set up contrasts of interest dan mengkalkulasi kembali model coefficients
cts <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","GB_ACC","SPOT_ID","Gene.Symbol","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

# Visualisasi dan quality control hasil tes.
# Membangun histogram daari P-values untuk semua genes. Normal test
# Asumsi tidak semua gen merupakan ekspresi diferensial
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
  ylab = "Number of genes", main = "P-adj value distribution")

# Merangkum hasil sebagai "up", "down" atau "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=0)

# Venn diagram
vennDiagram(dT, circle.col=palette())

# Membuat Q-Q plot untuk t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # plih Contrat of Interest
 
# volcano plot menggunakan limma function:
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
  highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

################################################################
# Analisis ekspresi data umum
ex <- exprs(gset)

#PART H. Visualisasi Data
# box-and-whisker plot
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE18090", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE18090", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")

# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 11, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=11", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
col=1:nlevels(gs), title="Group", pt.cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

