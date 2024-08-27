#!/usr/bin/env Rscript

require(reticulate) # need to read python pickles
source('KeyBeanplotFunctions.R')
use_python("~/miniconda3/bin/python", required=T)

setwd("~/afdb_domain/tools/dompair-eval-afdb")
source_python("pkl2dict.py") # imports python function load_domdict(fname:str)

cath <- c('C', 'A', 'T','H')

regexp_match <- function(x, regex){
  #given a regex, return just the first matching part of the string(s) in x.
  regmatches(x = x, m = regexpr(pattern = regex, text = x))
  }

mag <- function(vec) {
  sqrt(sum(vec**2))
  }

unit_vec <- function(vec){
  vec / mag(vec)
  }

vecs2cio <- function(vecs) {
  # remove nans?
  vecs <- matrix(vecs[!is.na(vecs)], ncol=3)

  n <- nrow(vecs)
  uv <- apply(vecs, 1, unit_vec)
  s <- apply(uv, 1, sum, na.rm=F)
  1.0 - mag(s)/ n

  }

map2color_linear <- function(x,pal,limits=NULL){
  if(is.null(limits))
    limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

map2color_log <- function(x,pal,limits=NULL){
  if(is.null(limits))
    limits=range(x)

  col_lims <- axisTicks(usr=limits, log=T, axp=c(limits, 1), nint=length(pal)+1)

  pal[findInterval(x,col_lims, all.inside=TRUE)]
}

nvec <- function(x){
  if (is.vector(x$vectors)){
    stopifnot(length(x$vectors) == 3)
    return(1)
  }
  return(nrow(x$vectors))
}

df_check <- function(x) {
  nums <- c(length(x$aligned_domain_pairs), length(x$choppings), nrow(x$vectors))
  nums
  }

get_cath_name <- function(sfam, tbl, ascend_if_na = FALSE, depth_str='') {
  # first try matching the whole sfam

  n <- tbl[sfam, 'name']
  if (is.na(n) && ascend_if_na){
      sfam_split <- strsplit(sfam,split = '.', fixed = T)[[1]]
      depth <- length(sfam_split)
      n <- get_cath_name(paste(sfam_split[1:(depth-1)], collapse = '.'), tbl, ascend_if_na, depth_str=paste0(" (",cath[depth-1],")"))
  }
  paste0(n, depth_str)
}

# retrieve descriptive names from CATH for the H-families in an ISP.
# If there is no name at the H-level, `ascend_if_na` controls whether we should ascend the CATH hierarchy until a non-empty name is found.
isp_to_names <- function(isp, tbl, sep='-', ascend_if_na = FALSE){

  sfams <- strsplit(x = isp, split = sep, fixed = T)[[1]]
  result <- vector(mode = 'character', length = length(sfams))
  for (i in 1:length(result)){
    result[i] <- get_cath_name(sfams[i], tbl, ascend_if_na)
  }
  result
}

isp_to_pdbs <- function(isp, d, sep=':') {

  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  get_ted_domnum <- function(x){substrRight(x, 5)}

  i <- d[[isp]]

  # all_pairs <- c(i$ref_domain_pair, i$aligned_domain_pairs)
  all_pairs <- i$aligned_domain_pairs

  # ASSUMES that pairs are always domains from the same single chain
  afdb_models <- regexp_match(x = all_pairs, regex = '^AF-[a-zA-Z0-9]*-F1-model_v4')

  splitted <- strsplit(x = all_pairs, split = sep, fixed = T)
  ted_dom_suffixes <- lapply(X = splitted, FUN = function(x) regexp_match(x = x, regex = '_TED[0-9]{2}$'))
  tbl <- cbind(afdb_models, Reduce(f = rbind, x = ted_dom_suffixes, init=data.frame()))
  colnames(tbl) <- c("model", "d1", "d2")

  tbl
}


good_pae_idx <- function (x, thresh){
  return(which(x$pae_score <= thresh))
}


pae_filter_loss <- function (d, cutoff){
  n_isp <- length(d)
  # lengths of unfiltered lists
  l <- sapply(d, function(x){length(x$pae_score)})
  filt_l <- sapply(d, FUN = function(x, thresh){length(good_pae_idx(x, thresh))}, thresh=cutoff)

  still_good_isps <- which(filt_l > 0)
  n_still_good_isps <- length(still_good_isps)
  isp_retention <- n_still_good_isps / n_isp

  sum_l <- sum(l)
  sum_filt_l <- sum(filt_l)
  retention <- sum_filt_l/sum_l
  writeLines(paste0("pAE Cutoff: ", cutoff, " Total: ", sum_l, ", retained: ", sum_filt_l, ", Instance Retention: ", retention, ', ISP retention: ', isp_retention))
  return(retention)
}

pae_filter_one <- function(d_subset, idx_to_retain){
  if (length(idx_to_retain) == 0){
    NA
    }
  out <- list()

  for (k in c("aligned_domain_pairs", "choppings", "pae_score")){
    out[[k]] <- d_subset[[k]][idx_to_retain]
  }
  out$vectors <- d_subset$vectors[idx_to_retain,]
  out
}

# a function to determine how many PDB files we'd need to load to make a figure for a given ISP for the paper
# used interactively, not in script below.
find_total_pdbs_for_hub <- function(sfam, dcath, dafdb){
    idx <- grep(pattern=sfam,x=names(dcath),fixed=T)
    d_subset <- dcath[idx]
    nums  <- sapply(d_subset, nvec)
    total_cath <- sum(nums)
    writeLines('Counts for CATH:')
    print(nums)
    writeLines(paste('Total:', total_cath))

    idx <- grep(pattern=sfam,x=names(dafdb),fixed=T)
    d_subset <- dafdb[idx]
    nums  <- sapply(d_subset, nvec)
    total_afdb <- sum(nums)
    writeLines('Counts for AFDB:')
    print(nums)
    writeLines(paste('Total:', total_afdb))

    writeLines(paste0("Grand total: ", total_cath+total_afdb, " pairs (", (total_cath+total_afdb)*2, " ref + tag PDBs)."))

}


isp_is_homotypic <- function(isp, sep='-'){
    sp <- strsplit(isp, sep, fixed=T) # a list of vectors (should all be of length 2)
    sapply(sp, FUN=function(x){x[1]==x[2]})
}


##########################################
# MAIN

pae_score_cutoff <- 4.0
# suffix for output files
suffix <- paste0("_filt_medianpae_",pae_score_cutoff)

def.par <- par(no.readonly = T)
make_pdf <- T

HOME  <- Sys.getenv("HOME")

afdb_pkl_f <- paste0(HOME,'/afdb_domain/tools/dompair-eval-afdb/all_interaction_vectors_raw_v2.pkl')

writeLines(paste('Reading', afdb_pkl_f, '...'))
d_afdb <- load_domdict(afdb_pkl_f)
d_cath <- load_domdict(paste0(HOME,'/afdb_domain/tools/dompair-eval-cath4.3/dompair_interaction_vector_data.pkl'))

gc() # apparently a lot of memory can be freed at this point!

# names for each node in the CATH hierarchy
cath_names <- read.table(file = "~/afdb_domain/tools/dompair-eval-afdb/cath-names-sane.tsv", header = F, sep = '\t', quote = '', row.names = 1, na.strings = ':')
names(cath_names) <- c("example", "name")
cath_names$name <- gsub(pattern = '^:', replacement = '', x = cath_names$name)

writeLines(paste("Apply median interdomain pAE filter threshold of", pae_score_cutoff,"A..."))

counts_afdb <- sapply(X = d_afdb, FUN = nvec)
writeLines(paste("Original AFDB list has", sum(counts_afdb), "instances across", length(d_afdb), "ISPs."))
filt_idx <- lapply(d_afdb, good_pae_idx, thresh=pae_score_cutoff)


stopifnot(all(names(d_afdb) == names(filt_idx)))
d_afdb_filt <- list()
for (n in names(d_afdb)){
  res <- pae_filter_one(d_afdb[[n]], filt_idx[[n]])
  d_afdb_filt[[n]] <- res
}
gc()
counts_afdb_filt <- sapply(X = d_afdb_filt, FUN = nvec)
bad <- which(counts_afdb_filt == 0)
d_afdb <- d_afdb_filt[-bad]

rm(d_afdb_filt)
gc()
writeLines(paste("After filtering, there are",sum(counts_afdb_filt),"instances across",length(d_afdb),"ISPs."))

writeLines('compute...')
n_afdb <- names(d_afdb)
n_cath <- names(d_cath)

# find the number of occurrences of each ISP
counts_afdb <- sapply(X = d_afdb, FUN = nvec)
counts_cath <- sapply(X = d_cath, FUN = nvec)

pairs_afdb <- strsplit(n_afdb, split = '-', fixed = T)
pairs_cath <- strsplit(n_cath, split = '-', fixed = T)

interactions_counts_afdb <- data.frame(from=sapply(pairs_afdb, `[`, 1), to=sapply(pairs_afdb, `[`, 2), value=counts_afdb)
interactions_counts_cath <- data.frame(from=sapply(pairs_cath, `[`, 1), to=sapply(pairs_cath, `[`, 2), value=counts_cath)

write.table(x = interactions_counts_afdb, file = paste0("isps_counts_afdb",suffix,".tbl"), append = F, quote = F, row.names = F, col.names = F)
write.table(x = interactions_counts_cath, file = "isps_counts_cath.tbl", append = F, quote = F, row.names = F, col.names = F)


# find the interacting superfamily pairs (ISPs) common to the two sets
common <- intersect(n_cath, n_afdb)
unique_to_afdb <- setdiff(n_afdb, common)
unique_to_cath <- setdiff(n_cath, common)

l_afdb <- length(n_afdb)
l_cath <- length(n_cath)

# extract the vector sets for those ISPs
d_comm_afdb <- d_afdb[common]
d_comm_cath <- d_cath[common]

# find the number of occurrences of each ISP
l_comm_afdb <- sapply(X = d_comm_afdb, FUN = nvec)
l_comm_cath <- sapply(X = d_comm_cath, FUN = nvec)

# Sum them to get the grand total number of examples of common ISPs in each set
total_int_afdb <- sum(l_comm_afdb)
total_int_cath <- sum(l_comm_cath)

# Fraction of the total that is contributed by each ISP
frac_comm_afdb <- l_comm_afdb / total_int_afdb
frac_comm_cath <- l_comm_cath / total_int_cath

stopifnot(all(names(l_comm_afdb) == names(l_comm_cath)))

# 'loget' was proposed in one paper as a shorthand for log_2(fold change).
# Needless to say, it didn't catch on.
logets <- log2(l_comm_afdb/l_comm_cath)
logets <- logets[order(logets)]  # increasing order

write(logets, file = 'log2foldchange.tbl', ncolumns = 1, append = F)

lower_i <- which(logets <0)
n_lower <- length(lower_i)
higher_i <- which(logets>0)
n_higher <- length(higher_i)
unchanged_i <- which(logets==0)
n_unchanged <- length(unchanged_i)

writeLines('Plot logets_horiz...')
#stop('Normal stop.')
if (make_pdf)
    pdf(file=paste0("logets_horiz_sorted",suffix,".pdf"), height=7, width=12)
plot(logets, ylab = expression(log[2]("fold change")), pch=16, col = 'deepskyblue', xlim=c(1,4000), las=1, main = "Fold change of raw ISP counts (AFDB vs CATH 4.3)")
abline(h=0, lty=2)
abline(v = range(unchanged_i), lty=3)
text(x = mean(lower_i), y=10, labels=paste(n_lower, 'interactions with higher count in CATH 4.3'), srt=90)
text(x = mean(unchanged_i), y = 10, labels = paste(n_unchanged, 'interactions with equal counts'), srt=90)
text(x = mean(higher_i), y=2.5, labels=paste(n_higher, 'interactions with higher count in AFDB'))
if (make_pdf)
    dev.off()
## AND/OR ...

writeLines('Plot isp_counts_comparison...')
if (make_pdf)
    pdf(file=paste0("isp_counts_comparison",suffix,".pdf"), height=10, width=10)
plot(l_comm_cath, l_comm_afdb,
     log='xy',
     xlab = 'Number of instances in CATH 4.3',
     ylab = 'Number of instances in AFDB',
     main = 'Interacting superfamily pairs common to CATH 4.3 and AFDB',
     pch = 16,
     col = rgb(0,0.7,1,0.5),
     xlim=c(1,1e6), ylim = c(1,1e6)
     )
abline(0,1, lty=2)
if (make_pdf)
    dev.off()

## AND/OR...
writeLines('Plot logets_hist...')
if (make_pdf)
    pdf(file=paste0("logets_hist",suffix,".pdf"), height=5, width=7)
par(mar=c(5,5,1,1))
hist(logets,
     breaks=seq(-10, 20, 1),
     col = 'deepskyblue',
     border = NA,
     xlim = c(-10,20),
     ylim=c(0,500),
     xlab= expression(log[2]("fold change")),
     main=NA
     )
if (make_pdf)
    dev.off()
par(def.par)

writeLines('compute CIO for common ISPs...')
cio_afdb <- sapply(X = d_comm_afdb, FUN = function(x){vecs2cio(x$vectors)})
cio_cath <- sapply(X = d_comm_cath, FUN = function(x){vecs2cio(x$vectors)})

df <- data.frame(counts_afdb = l_comm_afdb, counts_cath = l_comm_cath, cio_afdb, cio_cath)
write.table(x=df, file=paste0("common_isps_counts_cios",suffix,".tbl"), col.names=T, row.names=T, quote=F)

# find the top-n most enriched ISPs
top_n <- 100
most_enriched <- rev(tail(logets, n = top_n))
df_most_enriched <- df[names(most_enriched),]
o <- order(df_most_enriched$counts_cath, decreasing = F)
write.table(x = df_most_enriched[o, ], file = paste0('top', top_n, '_most_enriched_ISPs.tsv'), append = F, quote = F, sep = '\t', row.names = T, col.names = T)


writeLines('plot cio_comparison...')
if (make_pdf)
    pdf(file=paste0("cio_comparison",suffix,".pdf"), height=7, width=7)

# beanplot_params
bw <- 0.05
cut <- 0
what <- c(0,1,1,0)
beanlines <- 'quantiles'
beanline.length <- 0.5
border=NA
# bean fill, ?, ?, beanlines
col <- c('deepskyblue', NA, NA, 'black')


l <- layout(mat = matrix(1:4, nrow=2, byrow=T))
beanplot(df$cio_afdb, border = border, what=what, col = col, beanlines=beanlines, bw=bw, ylim=c(0,1), cut=cut, beanline.length = beanline.length)
plot(df$cio_cath, df$cio_afdb,
     xlab = 'CIO, CATH 4.3',
     ylab = 'CIO, AFDB',
     main = 'CIO comparison on common ISPs\n(0 = complete conservation)',
     pch = 16,
     col = rgb(0,0.7,1,0.5))
plot(0,0, type='n', xlab=NA, ylab=NA, bty='n', xaxt='n', yaxt='n')
beanplot(df$cio_cath, border= border, what=what, col = col, beanlines=beanlines, bw=bw, ylim=c(0,1), cut=cut, beanline.length = beanline.length, horizontal = T)


if (make_pdf)
  dev.off()

par(def.par)


d_cio <- data.frame(afdb=cio_afdb, cath=cio_cath)
d_cio$delta <- d_cio$afdb - d_cio$cath
d_cio_homo <- d_cio[isp_is_homotypic(row.names(d_cio)),]
d_cio_hetero <- d_cio[!isp_is_homotypic(row.names(d_cio)),]


if(make_pdf){
    side=640
    jpeg(file='cio_hist_homo_hetero.jpg', height=side, width=side, quality=100)
}
cex <- 0.5
cex.lab=2
cex.axis=1.5
cex.main=1.5
xaxs='i'
yaxs='i'
par( mar=c(5,12,3,1), las=1)

beanplot(
  d_cio$delta, d_cio_homo$delta, d_cio_hetero$delta,
  border= border, 
  what=what, 
  col = col, 
  beanlines=beanlines, 
  bw=bw, 
  ylim=c(-1,1), 
  xlab=expression(Delta~'CIO'),
  cut=cut, 
  beanline.length = beanline.length, 
  horizontal = T,
  main='Distributions of CIO differences (TED - CATH)',
  names=c('All ISPs', 'Homo-type ISPs', 'Hetero-type ISPs'),
  las=1, 
  cex.axis=1.5,
  cex.lab=1.5
)

if (make_pdf) {
    par(def.par)
    dev.off()
}


# if(make_pdf){
#     side=320
#     jpeg(file='cio_subsets.jpg', height=side, width=3*side, quality=100)
# }
# cex <- 0.5
# cex.lab=2
# cex.axis=1.5
# cex.main=1.5
# xaxs='i'
# yaxs='i'
# par(mfcol=c(1, 3), mar=c(5,5,3,1), las=1)
# plot(cio_cath, cio_afdb, 
#   pch=16, 
#   cex=cex, 
#   xlim=c(0,1), 
#   ylim=c(0,1), 
#   col='deepskyblue', 
#   xlab='CIO, CATH', 
#   ylab='CIO, TED', 
#   main='All common ISPs (n = 3070)', 
#   cex.lab=cex.lab, 
#   cex.axis=cex.axis,
#   cex.main=cex.main,
#   xaxs=xaxs,
#   yaxs=yaxs
#   )
# abline(0,1,lty=2)
# plot(d_cio_homo$cath, d_cio_homo$afdb,
#   pch=16, 
#   cex=cex, 
#   xlim=c(0,1), 
#   ylim=c(0,1), 
#   col='deepskyblue', 
#   xlab='CIO, CATH', 
#   ylab='CIO, TED', 
#   main='Homo-type ISPs (n = 293)', 
#   cex.lab=cex.lab, 
#   cex.axis=cex.axis,
#   cex.main=cex.main,
#   xaxs=xaxs,
#   yaxs=yaxs
#   )
# abline(0,1,lty=2)
# plot(d_cio_hetero$cath, d_cio_hetero$afdb, 
#   pch=16,
#   cex=cex, 
#   xlim=c(0,1), 
#   ylim=c(0,1), 
#   col='deepskyblue', 
#   xlab='CIO, CATH',
#   ylab='CIO, TED', 
#   main='Hetero-type ISPs (n = 2777)', 
#   cex.lab=cex.lab, 
#   cex.axis=cex.axis,
#   cex.main=cex.main,
#   xaxs=xaxs,
#   yaxs=yaxs
#   )
# abline(0,1,lty=2)

# if (make_pdf) {
#     par(def.par)
#     dev.off()
# }

# writeLines('Compute CIO for ALL ISPs, including ones not common to both sets...')
# cio_afdb_all <- sapply(X = d_afdb, FUN = function(x){vecs2cio(x$vectors)})
# cio_cath_all <- sapply(X = d_cath, FUN = function(x){vecs2cio(x$vectors)})

# which_n_cath_homo <- isp_is_homotypic(names(cio_cath_all))
# which_n_cath_hetero <- !which_n_cath_homo

# which_n_afdb_homo <- isp_is_homotypic(names(cio_afdb_all))
# which_n_afdb_hetero <- !which_n_afdb_homo

# cio_cath_all_homo <- cio_cath_all[which_n_cath_homo]
# cio_cath_all_hetero <- cio_cath_all[which_n_cath_hetero]

# cio_afdb_all_homo <- cio_afdb_all[which_n_afdb_homo]
# cio_afdb_all_hetero <- cio_afdb_all[which_n_afdb_hetero]

# if(make_pdf){
#     side=640
#     jpeg(file='cio_subsets_all.jpg', height=side, width=side, quality = 100)
#     par(mar=c(5,8,1,1))
# }
# beanplot(cio_cath_all_homo, cio_cath_all_hetero, cio_afdb_all_homo, cio_afdb_all_hetero, 
#   border= border, 
#   what=what, 
#   col = col, 
#   beanlines=beanlines, 
#   bw=bw, 
#   ylim=c(0,1), 
#   xlab='CIO',
#   cut=cut, 
#   beanline.length = beanline.length, 
#   horizontal = T,
#   #main='Distribution of CIO values',
#   names=c('CATH,\nhomo-type', 'CATH,\nhetero-type', 'TED,\nhomo-type', 'TED,\nhetero-type'),
#   las=1, 
#   cex.axis=1.5,
#   cex.lab=1.5
#   )
# if (make_pdf) {
#     par(def.par)
#     dev.off()
# }

writeLines('Done.')
