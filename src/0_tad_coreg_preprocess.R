library(rtracklayer)
library(txtplot)
library(data.table)
library(tidyverse)
library(magrittr)
library(Matrix)
library(matrixStats)
library(here)
select<-dplyr::select

sum_similiar=T


# fread('ext_data/mm10.cage_peak_phase1and2combined_fair_counts_ann.osc.txt.gz.extract.tsv.gz')
if(!exists('tpms')) tpms = fread('ext_data/mm10.cage_peak_phase1and2combined_fair_tpm_ann_fix1.osc.txt.gz.extract.tsv.gz')
if(!is(counts,'data.table')){
	counts = fread('ext_data/mm10.cage_peak_phase1and2combined_fair_counts_ann.osc.txt.gz.extract.tsv.gz')
}



excludecolpatts = list(
	"Clontech%20Mouse%20Universal%20Reference",
	"unclassified.",
	"whole%20body%",
	"feeder",
	"Universal",
	"SABiosciences",
	"testis",
	"epididymis"
)

for(patt in excludecolpatts){
	counts %<>% select(-matches(patt))
}

tpmmat = counts%>%select(matches('^counts'))%>%mutate_all(as.numeric)%>%as.matrix
tpmmat%<>%set_rownames(counts$short_description)


smalllib = tpmmat%>%colSums%>%`<`(1e6)
tpmmat = tpmmat[,!smalllib]
counts%<>%select(-one_of(names(smalllib)[smalllib]))

nodata = tpmmat%>%rowSums%>%add(1)%>%log2%>%`<`(5)
tpmmat = tpmmat[!nodata,]
counts = counts[!nodata,]

peakexprs = tpmmat%>%rowMedians
rownames(tpmmat) = counts$short_description

if(!file.exists(here('data/peakgr.rds'))){
	peakgr <- {

		peakgr = counts[[1]]%>%str_match('(chr[0-9MXY]+):(\\d+)\\.\\.(\\d+),([+-])')%>%.[,-1]%>%as.data.frame%>%set_colnames(c('seqnames','start','end','strand'))%>%GRanges
		peakgr$peakname = counts$short_description

		library(liftOver)
		chainpath = here('../cortexomics/ext_data/mm9ToMm10.over.chain')
		ch = import.chain(chainpath)
		seqlevelsStyle(cur) = "UCSC"  # necessary
		peakgr = liftOver(peakgr, ch)%>%unlist


		annogr = fread('ext_data/gencode.vM23.annotation.gtf.gz')

		ext_exons = annogr%>%subset(type=='exon')%>%
			subset(start>200)%>%
			split(.,.$transcript_id)%>%
			resize_grl(sum(width(.))+200,fix='end')

		#verify we get all the fat1 peaks
		missingfat1peaks = peakgr%>%subset(peakname%>%str_detect('Fat1'))%>%subsetByOverlaps(inv=T,ext_exons)%>%length
		stopifnot(missingfat1peaks==0)

		fmcol = function(gr,col){
			mcols(gr@unlistData)[[col]][gr@partitioning@end]
		}
		pov = mergeByOverlaps(peakgr[,'peakname'],ext_exons)
		genepeaksdf = tibble(
			peakname=pov$peakname,
			gene_id = fmcol(pov$ext_exons,'gene_id'),
			gene_name = fmcol(pov$ext_exons,'gene_name')
			)%>%
			group_by(peakname)%>%
			distinct%>%
			filter(n_distinct(gene_id)==1)

		peakgr$gene_id = peakgr$peakname%>%tibble(.)%>%left_join(genepeaksdf,by=c(.='peakname'))%>%.$gene_id
		peakgr$gene_name = peakgr$peakname%>%tibble(.)%>%left_join(genepeaksdf,by=c(.='peakname'))%>%.$gene_name
		peakgr$gene_id <- ifelse(is.na(peakgr$gene_id),peakgr$peakname,peakgr$gene_id)
		peakgr$gene_name <- ifelse(is.na(peakgr$gene_name),peakgr$peakname,peakgr$gene_name)
		# peakgr$entrezgene_id <- ifelse(is.na(counts$entrezgene_id),peakgr$peakname,counts$entrezgene_id)
		peakgr$transcript_id = peakgr$gene_id
		# peakgr$gene_id = peakgr$entrezgene_id
		peakgr$type = 'exon'
		# peakgr$gene_name = peakgr$peakname%>%str_extract('[^@]+$')
		peakgr$name = peakgr$peakname
		peakgr%>%subset(seqnames=='chr8')%>%export('data/peak.gtf')
		peakgr%>%subset(seqnames=='chr8')%>%export('data/peak.bed')

		peakgr

	}
	saveRDS(peakgr,here('data/peakgr.rds'))
}else{
	peakgr<-readRDS(here('data/peakgr.rds'))
	stopifnot(peakgr%>%length==	(counts%>%nrow))
}

fmcols <- function(grl,...){
  with(grl@unlistData@elementMetadata,...)[start(grl@partitioning)]
}

peakgrl=peakgr%>%split(.,.$gene_id)

allgenegr = GRanges(
	seqnames(peakgrl)%>%as("CharacterList")%>%.[as(rep(1,length(.)),'IntegerList')]%>%unlist,
	IRanges(
		peakgrl%>%start%>%min,
		peakgrl%>%start%>%max
	),
	strand = strand(peakgrl)%>%as("CharacterList")%>%.[as(rep(1,length(.)),'IntegerList')]%>%unlist
)

allgenegr$peakname<-fmcols(peakgrl,peakname)
allgenegr$gene_id<-fmcols(peakgrl,gene_id)
allgenegr$gene_name<-fmcols(peakgrl,gene_name)
allgenegr = allgenegr[allgenegr$gene_id %>%str_detect("ENSMUSG")]
allgenegr%<>%sort(ignore.strand=TRUE)

stopifnot(allgenegr$peakname%>%str_subset('Fat1')%>%length%>%`>`(0))
stopifnot(allgenegr$peakname%>%str_subset('Zfp42$')%>%length%>%`>`(0))

# if(sumpeaks){
if(TRUE){
	tpmmat =tpmmat[ match(peakgr$peakname,rownames(tpmmat)) , ]
	allgenetpmmat = rowsum(tpmmat,peakgr$gene_id)
	allgenetpmmat %<>% .[allgenegr$gene_id,]
	#rows of the expr mat, and the genegr, now named after the top peak
	allgenegr$peakname%<>%str_replace('.*?@','')
	rownames(allgenetpmmat) = allgenegr$gene_name
	stopifnot('Fat1'%in%rownames(allgenetpmmat))
	stopifnot('Zfp42'%in%rownames(allgenetpmmat))
}else{
	allgenetpmmat<-tpmmat[genepeaksdf$peak,]
	allgenetpmmat%<>%set_rownames(genepeaksdf$short_description)
}
allgenetpmmat_ocols = allgenetpmmat

sim_lim = 0.05
if(sum_similiar){
	libcors = cor(allgenetpmmat[,])
	libhclust = hclust(as.dist(1 - libcors))
	colgrps = cutree(h=sim_lim,libhclust)
	allgenetpmmat = t(rowsum(t(allgenetpmmat),colgrps))
}
stopifnot(ncol(allgenetpmmat) < ncol(allgenetpmmat_ocols))

colnames(allgenetpmmat) <- allgenetpmmat%>%colnames%>%
	str_replace('counts.','')%>%
	str_replace('.mm10.nobar.*$','')%>%
	str_replace_all('%20','_')

file.remove('gene_count_table.tsv')
file.remove('gene_count_table.tsv.gz')

allgenetpmmat[,]%>%as_tibble(rownames='gene')%>%write_tsv('gene_count_table.tsv')
system('gzip gene_count_table.tsv')
normalizePath('gene_count_table.tsv.gz')%>%message

nallgenetpmmat = DESeq2::vst(allgenetpmmat)

rownorm <-function(x) x %>%sweep(.,MARGIN=2,F='/',STAT=)
libsizes = allgenetpmmat%>%{DESeq2::estimateSizeFactorsForMatrix(.)}
allgenetpmmat%<>% sweep(.,MARGIN=1,F='/',STAT=libsizes)
#
rowmeds = nallgenetpmmat%>%rowMedians
rowmaxs = nallgenetpmmat%>%rowMaxs

#now plot
plotfile<- here(paste0('plots/','hkdensplot','.pdf'))
pdf(plotfile)
	dplot = qplot(
		rowmeds,
		rowmaxs,
		geom='density_2d'
	)+
	stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
	scale_x_continuous(paste0('Median Expr'),limits=c(3,10))+
	scale_y_continuous(paste0('Max Expr'),limits=c(3,10))+
	geom_vline(xintercept=5.75,linetype='dashed')+
	ggtitle(paste0('Median Vs Max Expr Fantom Genes'))+
	theme_bw()
	print(dplot)
dev.off()
message(normalizePath(plotfile))

rowmeds = nallgenetpmmat%>%rowMedians
is_ubiq = (rowmeds > 5.75)%>%setNames(rownames(nallgenetpmmat))
allgenegr$gene_name%in%names(is_ubiq)

dir.create('tables')
allgenegr$is_ubiq = is_ubiq[allgenegr$gene_name]
allgenetpmmat_ocols%>%as.data.frame%>%rownames_to_column('gene_name')%>%
	mutate(is_ubiq=is_ubiq[gene_name])%>%
	write_tsv('tables/allgenetpm.tsv')



save.image('data/tad_coreg_preprocess.Rdata')