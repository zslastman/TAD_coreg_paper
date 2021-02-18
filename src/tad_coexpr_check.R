{
library(magrittr)
library(txtplot)
library(tidyverse)
dim = 37
nfact = 3
npoints = 101
factors = sample(1:24,nfact*dim,rep=T)%>%matrix(nrow=dim)
mixes = sample(1:10,nfact*npoints,rep=T)%>%matrix(nrow=nfact)
colnorm <-function(x) x %>%sweep(.,MARGIN=2,F='/',STAT=colSums(.))
mixes%<>%colnorm
# d,m		d,f         f,m

fakedata = factors %*% mixes

library(NMF)

nmfres = nmf(fakedata,nfact)
dim(nmfres@fit@W)==dim(factors)
dim(nmfres@fit@H)==dim(mixes)
nmffacts = nmfres@fit@W
factors = factors[,order(factors[1,],factors[2,])]
nmffacts = nmffacts[,order(nmffacts[1,],nmffacts[2,])]
cor(factors,nmffacts)
txtplot(factors[,1],nmffacts[,1])

pcares = princomp(t(fakedata))

#d,dred
pcares$loadings
#m,dred
pcares$scores

recondata = t(t(pcares$scores %*% t(pcares$loadings))+pcares$center)

txtplot(recondata[1,],t(fakedata)[1,])
}

#TODOS
#rank tads by correlation of genes inside - strip plot
#Are most coregulated genes nearby?
	#co-regged genes not in same tad that often?
#distinction - house o	
#maybe eliminate 
#co-regged tads vs non - difference?

#dotplot of the pair plot - include between neibhoring tad vs in same tad vs all pairs
#density plot of coreg for above 3 categories.
#binarisation of coreg with 4 barplots - in same tad or not, coreg or not
#in randomly shuffled genome how often would coregged genes fall in same tad.
#way of extracting most co-regulated tads.
#Heatmap of fat1 locus coreg vs tad membership.

#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6054162/


library(rtracklayer)
library(txtplot)
library(data.table)
library(tidyverse)
library(magrittr)
library(Matrix)
library(matrixStats)
library(here)
select<-dplyr::select


sumpeaks=TRUE
# filterlow=FALSE
use_vsd = FALSE
# unlog = TRUE
# use_pca=TRUE
# pcadim=829
CODLIM = 0.1
# discretize=FALSE
# plot_cosine=FALSE
# corfun = purrr::partial(cor,method='spearman')
corfun = purrr::partial(cor,method='pearson')

sum_similiar = TRUE
sim_lim = 0.05

ubiqfilts = list(
	allgns = NULL,
	ubiq = function(x) x,
	nonubiq = function(x) ! x
)
ubiq_grp = 'ubiq'
# ubiq_grp = 'highvar'
ubiqfilt = ubiqfilts[[ubiq_grp]]

tadsetnms = c('all','1mb')
tadsetnm=tadsetnms[2]

for(ubiq_grp in names(ubiqfilts)){
	for(tadsetnms in tadsetnm){
	
	ubiqfilt=ubiqfilts[[ubiq_grp]]
vsdstatus = ifelse(use_vsd,'withVSD_','noVSD_')
runfolder = paste0('plots/',
	# ifelse(sumpeaks,'sumpeaks_','bestpeak_'),
	vsdstatus,
	ubiq_grp,'_',
	tadsetnm,'/'
)
dir.create(runfolder,showWarn=F,recursive=T)
#Notes heatmap correlation gene pairing.
{


}

{
neurotads = import('ext_data/Neuron_Cortical.Bonev_2017-raw.domains.bed')
esctads = import('ext_data/mESC.Bonev_2017-raw.domains.bed')
seqlevels(neurotads) %<>% paste0('chr',.)
seqlevels(esctads) %<>% paste0('chr',.)

alltads = neurotads

if(tadsetnm=='1mb'){
	tads2use = alltads%>%subset(width>1e6)
}else{
	tads2use = alltads
}


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

		# pov%>%as.data.frame%>%
			# transmute(peakname = peakgr$peakname[queryHits],fmcols(ext_exons,c(gene_id,gene_name))[subjectHits])

		# genepeaksdf = counts%>%
		# 	select(short_description,entrezgene_id)%>%
		# 	mutate(median_expr = peakexprs)%>%
		# 	mutate(peak=1:n())%>%
		# 	group_by(entrezgene_id)%>%
		# 	slice(which.max(median_expr))
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

}

#get gene level expr, process, filter
{

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

if(sumpeaks){
	tpmmat =tpmmat[ match(peakgr$peakname,rownames(tpmmat)) , ]
	allgenetpmmat = rowsum(tpmmat,peakgr$gene_id)
	allgenetpmmat %<>% .[allgenegr$gene_id,]
	#rows of the expr mat, and the genegr, now named after the top peak
	allgenegr$peakname%<>%str_replace('.*?@','')
	rownames(allgenetpmmat) = allgenegr$peakname
	stopifnot('Fat1'%in%rownames(allgenetpmmat))
	stopifnot('Zfp42'%in%rownames(allgenetpmmat))
}else{
	allgenetpmmat<-tpmmat[genepeaksdf$peak,]
	allgenetpmmat%<>%set_rownames(genepeaksdf$short_description)
}

colnames(allgenetpmmat) <- allgenetpmmat%>%colnames%>%
	str_replace('counts.','')%>%
	str_replace('.mm10.nobar.*$','')%>%
	str_replace_all('%20','_')

if(ubiq_grp=='all'){

	file.remove('gene_count_table.tsv')
	file.remove('gene_count_table.tsv.gz')

	allgenetpmmat[,]%>%as_tibble(rownames='gene')%>%write_tsv('gene_count_table.tsv')
	system('gzip gene_count_table.tsv')
	normalizePath('gene_count_table.tsv.gz')%>%message

}

# if(filterlow){
# 	rmeds = allgenetpmmat%>%{setNames(rowMedians(.),rownames(.))}
# 	# %>%{setNames(.>0,names(.))}
# 	rmedhigh
# 	rmeds%>%log10%>%txtdensity
# 	(rmeds['Zfp42'])
# 	(rmeds['Fat1'])
# 	stopifnot(rmedhigh['Zfp42'])
# 	stopifnot(rmedhigh['Fat1'])
# 	allgenetpmmat%<>%.[rmedhigh,]
# 	allgenegr = allgenegr[rmedhigh]
# }

if(sum_similiar){

	libcors = cor(tpmmat[,])
	libhclust = hclust(as.dist(1 - libcors))
	colgrps = cutree(h=sim_lim,libhclust)
	tpmmat = t(rowsum(t(tpmmat),colgrps))

}


if(use_vsd){
	message('using vsd')
	allgenetpmmat = DESeq2::vst(allgenetpmmat)
}

nallgenetpmmat = if(!use_vsd) DESeq2::vst(allgenetpmmat) else allgenetpmmat

if(!use_vsd){
	rownorm <-function(x) x %>%sweep(.,MARGIN=2,F='/',STAT=)
	libsizes = allgenetpmmat%>%{DESeq2::estimateSizeFactorsForMatrix(.)}
	allgenetpmmat%<>% sweep(.,MARGIN=1,F='/',STAT=libsizes)
}
# if(unlog) allgenetpmmat = 2^allgenetpmmat

# if(discretize){
# 		message('discretize Expression')
# 	library(CoRegNet)
# 	 allgenetpmmat=discretizeExpressionData(allgenetpmmat)
# 	# allgenetpmmat = as.numeric(apply(allgenetpmmat,2,function(x) x>0))
# 	allgenetpmmat = allgenetpmmat[allgenetpmmat%>%rowSds%>%`!=`(0),]
# 	# allgenetpmmat['Zfp42',]%>%table
# 	# allgenetpmmat['Fat1',]%>%table
# 	allgenegr = allgenegr[allgenegr$peakname%in%rownames(allgenetpmmat)]
# }



# if(use_pca){
# 	message('use PCA')
# 	# pcadim = floor(ncol(allgenetpmmat)/2)
	
# 	pca = allgenetpmmat[,]%>%princomp
# 	# recondata = t(t(pca$scores[,1:pcadim] %*% t(pca$loadings[,1:pcadim]) ) +pca$center)

# 	allgenetpmmat = pca$scores[,1:pcadim]
# }

rowmeds = nallgenetpmmat%>%rowMedians
rowmaxs = nallgenetpmmat%>%rowMaxs



#now plot
plotfile<- here(paste0(runfolder,'hkdensplot','.pdf'))
pdf(plotfile)
	dplot = qplot(
		rowmeds,
		rowmaxs,
		geom='density_2d'
	)+
	stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
	scale_x_continuous(paste0('Median Expr'),limits=c(3,10))+
	scale_y_continuous(paste0('Max Expr'),limits=c(3,10))+
	geom_vline(xintercept=4.25,linetype='dashed')+
	ggtitle(paste0('Median Vs Max Expr Fantom Genes'))+
	theme_bw()
	print(dplot)
dev.off()
message(normalizePath(plotfile))


if(!(ubiq_grp=='all')){
	message('High var genes only')
	# coefvar = allgenetpmmat%>%apply(1,function(x) sd(x,na.rm=T)/mean(x,na.rm=T))
	# coefvar%>%txtdensity
	# coefvar%>%median
	rowmeds = nallgenetpmmat%>%rowMedians
	is_ubiq = rowmeds > 4.25
	if(ubiq_grp=='ubiq') keep = is_ubiq
	if(ubiq_grp=='nonubiq') keep = !is_ubiq
	stopifnot(is_ubiq%>%mean%>%between(.58,.59))
	allgenetpmmat = allgenetpmmat[is_ubiq,]
	allgenegr = allgenegr[allgenegr$peakname%in%rownames(allgenetpmmat)]
}

# allgenetpmmat['Zfp42',]%>%table
# allgenetpmmat['Fat1',]%>%table


f1tad = allgenegr%>%
	subset(peakname%>%str_detect('Fat1'))%>%
	findOverlaps(alltads)%>%subjectHits%>%unique
zpftad = allgenegr%>%
	subset(peakname%>%str_detect('Zfp42$'))%>%
	findOverlaps(alltads)%>%subjectHits%>%unique

if(ubiq_grp=='all') stopifnot(identical(f1tad,zpftad))

# stopifnot('Fat1'%in%rownames(allgenetpmmat))


# peakgr%>%{.$name<-.$peakname;.}%>%export('data/tpmpeaks.bed')
# allgenegr%>%{.$name<-.$peakname;.}%>%export('data/allgenegr.bed')

}

################################################################################
########Look at all pairwise combs
################################################################################

highcorlim = 0.15

#
tadovdf = allgenegr%>%findOverlaps(alltads,type='within')%>%as.data.frame%>%set_colnames(c('gene','tad'))
tadovdf_dub = tadovdf%>%left_join(tadovdf,by='tad')
allcorfile = file.exists(here(paste0('data/',vsdstatus,'allcors.rds')))
if(!allcorfile){
	message('getting all cors')
	allcors <- allgenetpmmat%>%t%>%cor
	allcors%<>%Matrix
	saveRDS(allcors,here('data/allcors.rds'))
}else{
	allcors<-readRDS(here('data/allcors.rds'))
	stopifnot(cor(allgenetpmmat[1,],allgenetpmmat[2,])==allcors[2,1])
}


tadovdf = allgenegr%>%findOverlaps(tads2use,type='within')%>%as.data.frame%>%set_colnames(c('gene','tad'))
tadovdfpairs = tadovdf%>%left_join(tadovdf,suffix=c('','2'),by='tad')
tad_coincidenmat = Matrix(FALSE,ncol=length(allgenegr),nrow=length(allgenegr),sparse=TRUE)
#only mark upper triangle, so count once
tad_coincidenmat[as.matrix(tadovdfpairs%>%select(gene,gene2)%>%filter(gene<gene2))] <- TRUE
tad_coincidenmat = tad_coincidenmat | t(tad_coincidenmat)
diag(tad_coincidenmat) = TRUE

highcormat = Matrix(sparse=TRUE,allcors[allgenegr$peakname,allgenegr$peakname]>0.5)

cortadtable = matrix(ncol=2,c(
	sum((!triu(highcormat)) & (!triu(tad_coincidenmat))),
	sum(triu(highcormat) & (!triu(tad_coincidenmat))),
	sum((!triu(highcormat)) & triu(tad_coincidenmat)),
	sum(triu(highcormat) & triu(tad_coincidenmat))
))

sampair = sum(diag(tad_coincidenmat)&  diag(highcormat))
cortadtable[4] = cortadtable[4] - sampair

fisher.test(cortadtable)
cortadtable

tad_coincidenmat <- tad_coincidenmat%>%
	set_colnames(allgenegr$peakname)%>%
	set_rownames(allgenegr$peakname)

ichr='chr1'
allchrs = unique(seqnames(allgenegr))

tad_dist_df = lapply(allchrs,function(chrs){

	chrgenegr = allgenegr%>%subset(seqnames==ichr)
	chrgnms = chrgenegr%>%.$peakname
	chrstarts = chrgenegr%>%resize(1,'center')%>%start
	chrdists = chrstarts%>%dist(method='manhattan')

	chrhighcormat = highcormat[chrgnms,chrgnms]
	chrtad_coincidenmat = tad_coincidenmat[chrgnms,chrgnms]
	chrcormat = allcors[chrgnms,chrgnms]
	chrcormat = Matrix(as.matrix(chrcormat),sparse=TRUE)

	tad_cor_df = summary(Matrix(as.matrix(chrdists),sparse=TRUE))%>%
		as.data.frame%>%
		select(i,j,dist=x)%>%
		left_join(summary(chrhighcormat)%>%rename('highcor':=x))%>%
		left_join(summary(chrtad_coincidenmat)%>%rename('tad':=x))%>%
		left_join(summary(chrcormat)%>%rename('cor':=x))

	tad_cor_df%<>%filter(i!=j)

	tad_cor_df$tad%<>%replace_na(FALSE)
	tad_cor_df$highcor%<>%replace_na(FALSE)

	tad_cor_df

})
tad_dist_df%<>%bind_rows
ntadpairs = tad_dist_df%>%filter(tad)%>%nrow
	}
}

stop()

{

{
################################################################################
########USE CODfinder
################################################################################
#first optimize window size
# chrs='chr8'
# ichr='chr8'
# source('Applications/CODfinder/CODfinder.R')
# chrs%<>%setdiff(c('chrY','chrM'))
# allgenegrsort = allgenegr%>%sort(ignore.strand=TRUE)

# wsizes = 4:8

# testcoddf = mclapply(mc.cores=4,wsizes,function(wsize){
# 	ichr='chr1'
# 	ipeaks = allgenegrsort%>%subset(seqnames==ichr)%>%.$peakname
# 	if(length(ipeaks)<10)return(NULL)
# 	chrexpr = allgenetpmmat[ipeaks,]
# 	chrcor = chrexpr%>%t%>%cor
# 	chrexprfile = str_interp('data/${ichr}_tpm_gene.tsv')
# 	chrexpr%>%t%>%as.data.frame%>%write_tsv((chrexprfile))
# 	codfinder(
# 		chrexprfile, wsize,
# 		outfile = str_interp('data/${ichr}_tpm_gene_codout'),
# 		CORLIM = 0.07
# 	)
# 	message(ichr)
	
# 	coddf = fread(str_interp('data/${ichr}_tpm_gene_codout.COD.csv'))%>%{IRanges(.[[2]],.[[3]])}%>%
# 		as('IntegerList')%>%as.list%>%enframe%>%unnest(value)%>%
# 		set_colnames(c('cod','value'))%>%
# 		as.data.frame%>%
# 		mutate(value = rownames(chrexpr)[.$value])%>%
# 		mutate(cod = paste0(ichr,'_',cod))

# 	chrcor = chrexpr%>%t%>%cor

# 	coddf%>%left_join()

# })



chrs='chr8'
ichr='chr8'
source('Applications/CODfinder/CODfinder.R')
chrs%<>%setdiff(c('chrY','chrM'))
allgenegrsort = allgenegr%>%sort(ignore.strand=TRUE)

coddf = mclapply(mc.cores=1,chrs,function(ichr){
	ipeaks = allgenegrsort%>%subset(seqnames==ichr)%>%.$peakname
	if(length(ipeaks)<10)return(NULL)
	chrexpr = allgenetpmmat[ipeaks,]
	chrcor = chrexpr%>%t%>%cor
	chrexprfile = str_interp('data/${ichr}_tpm_gene.tsv')
	chrexpr%>%t%>%as.data.frame%>%write_tsv((chrexprfile))
	codfinder(
		chrexprfile, 4,
		outfile = str_interp('data/${ichr}_tpm_gene_codout'),
		CORLIM = CODLIM,
		corfun = corfun
	)
	message(ichr)
	
	fread(str_interp('data/${ichr}_tpm_gene_codout.COD.csv'))%>%{IRanges(.[[2]],.[[3]])}%>%
		as('IntegerList')%>%as.list%>%enframe%>%unnest(value)%>%
		set_colnames(c('cod','value'))%>%
		as.data.frame%>%
		mutate(value = rownames(chrexpr)[.$value])%>%
		mutate(cod = paste0(ichr,'_',cod))
})
coddf%<>%bind_rows
coddf%<>%set_colnames(c('cod','peakname'))

# randomcdf = lapply(chrs,function(ichr){
# 	ipeaks = allgenegr%>%subset(seqnames==ichr)%>%.$peakname
# 	if(length(ipeaks)<10)return(NULL)
# 	chrexpr = allgenetpmmat[ipeaks,]
# 	chrexpr = chrexpr[sample(1:nrow(chrexpr),replace=FALSE),]
# 	chrcor = chrexpr%>%t%>%cor
# 	chrexpr%>%t%>%as.data.frame%>%write_tsv(str_interp('data/${ichr}_rand_tpm_gene.tsv'))
# 	codfinder(
# 		str_interp('data/${ichr}_rand_tpm_gene.tsv'), 4,
# 		outfile = str_interp('data/${ichr}_rand_tpm_gene_codout')
# 	)
# 	message(ichr)
	
# 	fread(str_interp('data/${ichr}_rand_tpm_gene_codout.COD.csv'))%>%{IRanges(.[[2]],.[[3]])}%>%
# 		as('IntegerList')%>%as.list%>%enframe%>%unnest(value)%>%
# 		set_colnames(c('cod','value'))%>%
# 		as.data.frame%>%
# 		mutate(value = rownames(chrexpr)[.$value])%>%
# 		mutate(cod = paste0(ichr,'_',cod))
# })

coddf %>% saveRDS(here('data/coddf.rds'))
coddf <- readRDS(here('data/coddf.rds'))
# randomcoddf %>% saveRDS(here('data/randomcoddf.rds'))



# codmat = Matrix(0,length(allgenegr),length(allgenegr),sparse=TRUE)%>%
#  	set_rownames(allgenegr$peakname)%>%
#  	set_colnames(allgenegr$peakname)
# codindmat = coddf%>%bind_rows%>%inner_join(.,.,by='cod')%>%select(2:3)%>%
# 	set_colnames(c('g1','g2'))%>%
# 	mutate(g1=match(g1,allgenegr$peakname))%>%
# 	mutate(g2=match(g2,allgenegr$peakname))%>%
# 	as.matrix
# codmat[codindmat]=1
}

{

fat1ind = allgenegrsort$peakname%>%str_detect('Fat1')%>%which
zfp42ind = allgenegrsort$peakname%>%str_detect('Zfp42$')%>%which


allgenegrsort$peakname%>%str_detect('Zfp42$')%>%which
allgenegrsort$peakname%>%str_detect('Fat1')%>%which

zfp1gr = allgenegrsort[allgenegrsort$peakname%>%str_detect('Zfp42$')%>%which]
fat1gr = allgenegrsort[allgenegrsort$peakname%>%str_detect('Fat1')%>%which]

localgenegr = allgenegrsort[(min(zfp42ind,fat1ind)-100):(max(zfp42ind,fat1ind)+100)]%>%subset(seqnames=='chr8')
localgenegr$order = 1:length(localgenegr)

localcoddf = coddf%>%filter(peakname%in%localgenegr$peakname)%>%bind_rows%>%inner_join(.,.,by='cod')%>%
	select(-cod)%>%mutate(grptype='cod')
localcoddf$peakname.x %<>% match(.,localgenegr$peakname)
localcoddf$peakname.y %<>%  match(.,localgenegr$peakname)
localcoddf%<>%filter(peakname.x>=peakname.y)

tadovdf = localgenegr%>%findOverlaps(tads2use,type='within')%>%as.data.frame%>%set_colnames(c('peakname','tad'))%>%
	mutate(peakname = localgenegr$peakname[.$peakname])
localtaddf = tadovdf%>%filter(peakname%in%localgenegr$peakname)%>%bind_rows%>%inner_join(.,.,by='tad')%>%
	select(-tad)%>%mutate(grptype='tad')
localtaddf$peakname.x %<>% match(.,localgenegr$peakname)
localtaddf$peakname.y %<>%  match(.,localgenegr$peakname)


localtaddf%<>%filter(peakname.x>=peakname.y)


grpsdf=bind_rows(localcoddf%>%mutate(grptype='cod'),localtaddf%>%mutate(grptype='tad'))

gbrks = setNames(
	localgenegr$peakname%>%str_detect(.,'Zfp42|Fat1')%>%which,
	localgenegr$peakname%>%str_subset(.,'Zfp42|Fat1')
)


chrexpr = allgenetpmmat[localgenegr$peakname,]


# if(plot_cosine){
	# chrcor = chrexpr%>%cosdist%>%as.matrix
# }else{
	chrcor = chrexpr%>%t%>%corfun
# }

chrdf = chrcor%>%as.data.frame%>%rownames_to_column('peakname.x')%>%
	gather(peakname.y,cor,-peakname.x)
chrdf$peakname.x %<>% match(.,localgenegr$peakname)
chrdf$peakname.y %<>%  match(.,localgenegr$peakname)
chrdf%<>%filter(peakname.x<=peakname.y)
}

mcorgene = chrdf%>%group_by(peakname.x)%>%filter(peakname.x>85)%>%mutate(nhigh=sum(cor>.9))%>%arrange(desc(cor))%>%ungroup%>%filter(peakname.x==peakname.x[1])%>%.$peakname.x
othmcorgenes = chrdf%>%filter(cor>.9,peakname.x==mcorgene)%>%.$peakname.y
mcorgene%<>%unique
othmcorgenes%<>%unique
mcorgene=145

{
	#now plot
	plotfile<- here(paste0(runfolder,'tadcodgrid','.pdf'))
	pdf(plotfile)
	print(
		ggplot(chrdf, aes(peakname.x, peakname.y)) + 
		geom_raster(data=chrdf,aes(fill=cor))+
		geom_raster(data=grpsdf%>%filter(grptype=='tad'),aes(fill=-1*as.numeric(grptype=='tad')),alpha=I(0.5))+
		geom_raster(data=grpsdf%>%filter(grptype!='tad'),aes(fill=1*as.numeric(grptype=='cod')),alpha=I(0.5))+
		# scale_color_discrete(name='colorname',colorvals)+
		# scale_x_discrete('gene 1',breaks=gbrks,labels=names(gbrks))+
		scale_x_discrete('gene 1')+
		scale_y_discrete('gene 2',breaks=gbrks)+
		ggtitle(paste0('TADs vs CODs'))+
		# scale_fill_gradient2()+
		scale_fill_gradient2(low = "blue", mid='white',high = "red")+
		geom_vline(xintercept=which(localgenegr$peakname%>%str_detect('Zfp42$')),linetype=2)+
		geom_vline(xintercept=which(localgenegr$peakname%>%str_detect('Fat1')),color=I('red'),linetype=2)+
		geom_hline(yintercept=which(localgenegr$peakname%>%str_detect('Fat1')),color=I('red'),linetype=2)+
		geom_hline(yintercept=which(localgenegr$peakname%>%str_detect('Zfp42$')),linetype=2)+
		geom_vline(xintercept=mcorgene,color=I('green'),linetype=2)+
		theme_bw())
	dev.off()
	message(normalizePath(plotfile))
}
chrs = allgenegr%>%seqnames%>%unique
map(chrs,possibly(otherwise=NULL,function(ichr) fread(str_interp('data/${ichr}_tpm_gene_codout.binsignal.csv'))%>%.$mean.cf))%>%unlist%>%mean
# 	mutate(seqnames=ichr)%>%
chrexpr[mcorgene,]%>%txtplot
chrexpr[othmcorgenes[1],]%>%txtplot


mcorgene = chrdf%>%group_by(peakname.x)%>%filter(peakname.x>120)%>%mutate(nhigh=sum(cor>.9))%>%arrange(desc(nhigh))%>%ungroup%>%filter(peakname.x==peakname.x[1])%>%.$peakname.x%>%unique
othmcorgenes = chrdf%>%filter(cor>.9,peakname.x==mcorgene)%>%.$peakname.y%>%unique
chrexpr[mcorgene,]%>%{names(.)[.>4]}

chrexpr[mcorgene,]%>%{names(.)[.>0]}
chrexpr[mcorgene,]%>%.[.==0]
chrexpr[mcorgene,]%>%names
dsetnames = chrexpr%>%colnames
dsetnames[chrexpr[mcorgene,]>4]
dsetnames[chrexpr[mcorgene+1,]>0]
dsetnames[chrexpr[mcorgene+2,]>0]
dsetnames[chrexpr[othmcorgenes,]>0]
dsetnames[chrexpr[othmcorgenes%>%sample(1),]>0]%>%str_detect('testes')

lapply(1:length(othmcorgenes),function(i){dsetnames[chrexpr[othmcorgenes[i],]>0]})%>%unlist%>%unique

(chrexpr[mcorgene,]>0)%>%table


# stop()

################################################################################
########Correlation across boundaries
################################################################################
#shoul probs use only tads with 2 or more genes
{
setstrand<-function(gr,strand){strand(gr)<-strand;gr}
{
lbounds = tads2use%>%resize(1)
rbounds = tads2use%>%resize(1,'end')

allgenegrns = allgenegr%>%setstrand('*')


tadedgeindr = follow(rbounds,allgenegrns)
outtadedgeindr = precede(rbounds,allgenegrns)
okayedges = !((tadedgeindr%>%is.na) | (outtadedgeindr%>%is.na))
tadedgeindr = tadedgeindr[okayedges]
outtadedgeindr = outtadedgeindr[okayedges]
intadedgeindr = follow(allgenegrns[tadedgeindr],allgenegrns)

tadedgeindl = precede(lbounds,allgenegrns)
outtadedgeindl = follow(lbounds,allgenegrns)
okayedges = !((tadedgeindl%>%is.na) | (outtadedgeindl%>%is.na))
tadedgeindl = tadedgeindl[okayedges]
outtadedgeindl = outtadedgeindl[okayedges]
intadedgeindl = precede(allgenegrns[tadedgeindl],allgenegrns)
tadedgeindl = tadedgeindl[(!is.na(tadedgeindl))&(!is.na(outtadedgeindl))&(!is.na(intadedgeindl)) ] 
outtadedgeindl = outtadedgeindl[(!is.na(tadedgeindl))&(!is.na(outtadedgeindl))&(!is.na(intadedgeindl)) ] 
intadedgeindl = intadedgeindl[(!is.na(tadedgeindl))&(!is.na(outtadedgeindl))&(!is.na(intadedgeindl)) ] 

i=1
all(start(allgenegrns[tadedgeindr]) > start(allgenegrns[intadedgeindr]))
all(start(allgenegrns[tadedgeindr]) < start(allgenegrns[outtadedgeindr]))

all(start(allgenegrns[tadedgeindl]) < start(allgenegrns[intadedgeindl]))
all(start(allgenegrns[tadedgeindl]) > start(allgenegrns[outtadedgeindl]))

is2genetad = intadedgeindr %in% tadedgeindl
is1genetad = intadedgeindr %in% outtadedgeindl

tadedgeindl%<>% .[(!is2genetad)&(!is1genetad)]
outtadedgeindl %<>% .[(!is2genetad)&(!is1genetad)]
intadedgeindl %<>% .[(!is2genetad)&(!is1genetad)]



is1genetad = tadedgeindr %in% intadedgeindr
table(is1genetad)

tadedgeindr%<>% .[!is1genetad]
outtadedgeindr %<>% .[!is1genetad]
intadedgeindr %<>% .[!is1genetad]

indexdf = data.frame(
	edge = c(tadedgeindr,tadedgeindl),
	inside = c(intadedgeindr,intadedgeindl),
	outside = c(outtadedgeindr,outtadedgeindl)
)

pcordf = map_df((1:nrow(indexdf)),function(i){
	map_df(.id='tad_boundary_between',c(TRUE,FALSE),function(tadbound){
		sel = if(tadbound) 'outside' else 'inside'
		 c( cor = corfun(
				allgenetpmmat[allgenegr[indexdf$edge[i]]$peakname,],
				allgenetpmmat[allgenegr[indexdf[[sel]][i]]$peakname,],
  			), 
  			width = distance(allgenegrns[indexdf$edge[i]],allgenegrns[indexdf[[sel]][i]])
  		  )
	})
})
indexdf$edge==indexdf$inside

pcordf$tad_boundary_between%<>%`==`(.,1)

allgenegr[indexdf$edge]%>%export('data/edgegenes.bed')
allgenegr[indexdf$inside]%>%export('data/insidegenes.bed')
allgenegr[indexdf$outside]%>%export('data/outsidegenes.bed')
rbounds%>%export('data/rbounds.bed')
lbounds%>%export('data/lbounds.bed')
peakgr%>%{.$name<-.$peakname;.}%>%export('data/tpmpeaks.bed')
allgenegr%>%{.$name<-.$peakname;.}%>%export('data/allgenegr.bed')


}



{

library(ggpubr)
bggdf = pcordf
library(txtplot)
#now plot
plotfile<- here(paste0(runfolder,'boundary_cors','.pdf'))
pdf(plotfile)
ribbonplot = bggdf%>%
	filter(!(tad_boundary_between&(log10(width)%>%{.<4})))%>%
	ggplot(.,aes(y=cor,color=tad_boundary_between,fill=tad_boundary_between,x=log10(width)))+
	geom_smooth()+
	geom_point(alpha=I(0.1))+
	# scale_color_discrete(name='colorname',colorvals)+
	scale_x_continuous(paste0('Log10(distance between Adjacent Cage Peaks)'))+
	scale_y_continuous(paste0('Correlation across all phantom data'))+
	ggtitle(paste0('Correlations Between Neighboring Genes'))+
	theme_bw()
histplot  = bggdf%>%
	ggplot(.,aes(color=tad_boundary_between,fill=tad_boundary_between,x=log10(width)))+
	geom_density()+
	facet_grid(tad_boundary_between~.)+
	# scale_color_discrete(name='colorname',colorvals)+
	scale_x_continuous(paste0('Log10(distance between Adjacent Cage Peaks)'))+
	ggtitle(paste0('Distribution of distance between genes'))+
	theme_bw()
print(ggarrange(plotlist=list(ribbonplot,histplot),ncol=1))
dev.off()
message(normalizePath(plotfile))

}

}

################################################################################
########Now get within tad correlations
################################################################################

{
tadovdf = allgenegr%>%findOverlaps(tads2use,type='within')%>%as.data.frame%>%set_colnames(c('gene','tad'))

tadofinterest = tadovdf%>%filter(allgenegr$peakname[gene]%>%str_detect('Fat1|Zfp42$'))%>%.$tad%>%unique
# tadovdf%>%filter(gene==Zfp42ind)

x=tadovdf%>%{split(.$gene,.$tad)}%>%.[[1]]
allgenetpmmat%<>%set_colnames(NULL)

intratadcors = tadovdf%>%
	{split(.$gene,.$tad)}%>%
	.[map_dbl(.,length)>1]%>%
	map(function(x)cor(t(allgenetpmmat[allgenegr$peakname[x],])))%>%
	map(~ .[lower.tri(.)])%>%
	map(unlist)%>%
	enframe('tad','cor')%>%
	unnest(cor)
#now plot


tadofinterestmean = intratadcors%>%filter(tad==tadofinterest)%>%.$cor%>%mean


plotfile<- here(paste0(runfolder,'intratadcors','.pdf'))
pdf(plotfile,w=5,h=5)
print(
	intratadcors%>%group_by(tad)%>%
	# mutate(mcor=mean(cor>quantile(cor,0.75,na.rm=T)))%>%
	summarise(mcor=mean(cor),n=n())%>%
	ungroup%>%arrange(mcor)%>%mutate(tad=as_factor(tad))%>%
	# group_by(tad)%>%summarise()%>%
	# ggplot(.,aes(x=as.numeric(tad),y=cor))+
	ggplot(.,aes(x=mcor))+
	# geom_linerange()%>%
	# geom_point()+
	geom_histogram()+
	scale_y_continuous(paste0('Frequency'))+
	scale_x_continuous(paste0('Mean(Gene-Gene Correlation)'))+
	ggtitle(paste0('Distribution of average Gene-Gene cor\n
		dashed line is Fat1 Tad'))+
	geom_vline(xintercept=tadofinterestmean,linetype=2)+
	theme_bw()
)
dev.off()
message(normalizePath(plotfile))
}

plotfile<- here(paste0(runfolder,'intratadcors_wsize_2d','.pdf'))
pdf(plotfile,w=5,h=5)
print(
	intratadcors%>%group_by(tad)%>%
	# mutate(mcor=mean(cor>quantile(cor,0.75,na.rm=T)))%>%
	summarise(mcor=mean(cor),n=n())%>%
	ungroup%>%arrange(mcor)%>%mutate(tad=as_factor(tad))%>%
	# group_by(tad)%>%summarise()%>%
	# ggplot(.,aes(x=as.numeric(tad),y=cor))+
	ggplot(.,aes(x=mcor,y=log2(n)))+
	# geom_linerange()%>%
	# geom_point()+
	geom_bin2d()+
	# geom_jitter(width=I(.2))+
	scale_fill_continuous(paste0('Frequency'))+
	scale_x_continuous(paste0('Mean(Gene-Gene Correlation)'))+
	scale_y_continuous('Number of Genes in TAD')+
	ggtitle(paste0('Distribution of average Gene-Gene cor\n
		dashed line is Fat1 Tad'))+
	geom_vline(xintercept=tadofinterestmean,linetype=2)+
	theme_bw()
)
dev.off()
message(normalizePath(plotfile))






}
# gamob = gam(data=tad_dist_df,cor ~ s(log10(dist),bs='cs'))
# tad_dist_df$cor - predict(new.data=tad_dist_df,gamob)

# #now plot
# plotfile<- here(paste0(runfolder,'withintadcordist','.pdf'))
# pdf(plotfile)
# sample(Matrix(allcors)[tad_coincidenmat==1],100))%>%
# 	ggplot(.,aes())+
# 	scale_color_discrete(name='colorname',colorvals)+
# 	scale_x_continuous(paste0('xname'))+
# 	scale_y_continuous(paste0('yname'))+
# 	ggtitle(paste0('title'))+
# 	theme_bw()
# dev.off()
# message(normalizePath(plotfile))

# t_randcorgenes = randcorgenes%>%
# 	left_join(tadovdf,by=c('peakname.x'='gene'))%>%
# 	left_join(tadovdf,by=c('peakname.y'='gene'))%>%
# 	mutate(sametad = isTRUE(tad.x==tad.y) )
# t_randcorgenes%>%.$sametad%>%mean

# t_randncorgenes = randncorgenes%>%
# 	left_join(tadovdf,by=c('peakname.x'='gene'))%>%
# 	left_join(tadovdf,by=c('peakname.y'='gene'))%>%
# 	mutate(sametad = isTRUE(tad.x==tad.y) )

# t_randncorgenes%>%.$sametad%>%mean

# randcorgenes%>%.$sametad%>%mean

# #isintad = overlapsAny(genegr,tads2use)
# # 
# {
# tadgenes = allgenegr%>%subsetByOverlaps(tads2use)
# genegr = allgenegr%>%
# 	# subset(seqnames%in%c('chr1','chr2','chr3','chr4','chr5'))
# 	# subsetByOverlaps(invert=T,alltads)%>%
# 	# sample(length(tadgenes))%>%
# 	identity


# # genegr = genegr[isintad]
# genetpmmat = allgenetpmmat[genegr$id,]

# tadovdf = allgenegr%>%findOverlaps(tads2use,type='within')%>%as.data.frame%>%set_colnames(c('gene','tad'))
# tadovdfpairs = ovdf%>%left_join(ovdf,suffix=c('','2'),by='tad')
# tad_coincidenmat = Matrix(0,ncol=length(genegr),nrow=length(genegr),sparse=TRUE)
# #only mark upper triangle, so count once
# tad_coincidenmat[as.matrix(tadovdfpairs%>%select(gene,gene2)%>%filter(gene<gene2))] <- 1
# tad_coincidenmat
# tadpairs = tad_coincidenmat%>%reshape::melt.matrix(.)%>%filter(X1<X2,value==1)
# nontadpairs = tad_coincidenmat%>%reshape::melt.matrix(.)%>%filter(X1<X2,value==0)

# cormat = cor(t(genetpmmat))
# genetpmmat%>%dim
# cormat%>%dim
# ggdf = bind_rows(
# 	tibble(set='nontad',cor=cormat[as.matrix(nontadpairs[,1:2])],g1=nontadpairs[[1]],g2=nontadpairs[[2]]),
# 	tibble(set='tad',cor=cormat[as.matrix(tadpairs[,1:2])],g1=tadpairs[[1]],g2=tadpairs[[2]])
# )
# ggdf$dist = distance(genegr[ggdf$g1,],genegr[ggdf$g2,])
# ggdf%<>%filter(!is.na(dist))
# ggdf%<>%filter(is.na(dist)<1e6)
# }

# #now plot
# plotfile<- (paste0(runfolder,'tmp','.pdf'))
# pdf(plotfile)
# ggplot(data=ggdf%>%filter(dist<1e7),aes(fill=set,color=set,x=log10(dist),y=cor))+geom_smooth(alpha=I(0.5))+
# 	# geom_point(size=I(0.1))+
# 	theme_bw()
# dev.off()
# message(normalizePath(plotfile))

# library(ape)

# stop()

# genesintads = genegr[ovdf$gene%>%unique,]



# neurotads%>%width%>%log10%>%hist


# randomtads 


# probind = genegr[,3]%>%as.numeric%>%is.na%>%which%>%head(1)
