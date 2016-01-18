# Festival of Genomics, London, 18th January 2016

## Edinburgh Genomics - Edinburgh Genomics: Oxford Nanopore MinION sequencing - Data handling, analysis and applications

#### View individual files using hdfvew
```sh
hdfview MAP006-1/MAP006-1_downloads/pass/LomanLabz_PC_Ecoli_K12_MG1655_20150924_MAP006_1_5005_1_ch150_file24_strand.fast5 &
hdfview /data/basecall/pc6_groupc/pass/NB02/nanopore2_PoreCamp_GrpC_0829_1_ch501_file24_strand.fast5 &
```

### Looking at Nick's SQK-MAP-006 data
```R
library(poRe)
setwd("C:/Users/Mick/Dropbox (Edinburgh Genomics)/Presentations/FestivalOfGenomics")

# load in pass data
pass <- read.table("pass.meta.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# set standard/expected column names
colnames(pass) <- c("filename","channel_num","read_num","read_start_time",
                    "status","tlen","clen","len2d",
                    "run_id","read_id","barcode","exp_start")

# load in the fail data
fail <- read.table("fail.meta.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# set standard/expected column names
colnames(fail) <- c("filename","channel_num","read_num","read_start_time",
                    "status","tlen","clen","len2d",
                    "run_id","read_id","barcode","exp_start")

# take a look at the data
head(pass)
head(fail)

# count pass and fail reads
nrow(pass)
nrow(fail)

# view pass and fail reads
split.screen(c(1,2))
screen(1)
plot(pass$tlen, pass$clen, xlab="Template length (bp)", ylab="Complement length (bp)", main="Pass reads", pch=16)
screen(2)
plot(fail$tlen, fail$clen, xlab="Template length (bp)", ylab="Complement length (bp)", main="Fail reads", pch=16)
close.screen(all=TRUE)


# find the longest pass 2D read
max(pass$len2d)
pass[pass$len2d==max(pass$len2d),]
fasta <- get_fasta("LomanLabz_PC_Ecoli_K12_MG1655_20150924_MAP006_1_5005_1_ch153_file57_strand.fast5")$'2D'
cat(fasta, file = "longest.fa", sep = "\n", fill = FALSE)

# navigate to http://blast.ncbi.nlm.nih.gov/Blast.cgi

# length histograms
plot.length.histogram(pass)

# plot cumulative yield    
plot.cumulative.yield(pass)

# quality and alignment
https://github.com/arq5x/nanopore-scripts


# Mapping
map <- read.table("MAP006-1.pass.2D.poRe.profile.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

boxplot(map$align_len/map$read_len, main="SQK-MAP-006 alignment length", ylab="% of read length", xlab="Pass", col="skyblue3")

boxplot(100 * map$matches/(map$matches + map$deletions + map$insertions + map$mismatches), main="SQK-MAP-006 accuracy", ylab="% accuracy", xlab="Pass", col="purple2")
                  
nrow(map)
nrow(pass)


# Squiggles
r1 <- "LomanLabz_PC_Ecoli_K12_MG1655_20150924_MAP006_1_5005_1_ch440_file31_strand.fast5"
r2 <- "LomanLabz_PC_Ecoli_K12_MG1655_20150924_MAP006_1_5005_1_ch446_file30_strand.fast5"

t1 <- get.events(r1, path.t = "/Analyses/Basecall_2D_000/", path.c = "/Analyses/Basecall_2D_000/")
t2 <- get.events(r2, path.t = "/Analyses/Basecall_2D_000/", path.c = "/Analyses/Basecall_2D_000/")

split.screen(c(2,1))
screen(1)
plot.squiggle(t1$template, maxseconds=20)
screen(2)
plot.squiggle(t2$template, maxseconds=20)
close.screen(all=TRUE)
```




