#This script will run the plots involving HiC signal
#for now may move back intot he 2 tad _v _dist loop later
{
library(rtracklayer)
library(txtplot)
library(data.table)
library(tidyverse)
library(magrittr)
library(Matrix)
library(matrixStats)
library(here)
library(tidyverse)
}

if(!exists('matad_dist_df'))load(here('data/2_tad_v_dist_plot.Rdata'))

library(strawr)
hicchrs = as.character(1:19)
hicchr=hicchrs[1]
hicfile = 'hicdata/CN_mapq30.hic'
strawr::readHicChroms(hicfile)
manum = 2e3


ubiq_setnm = 'ubiq'
ubiq_set=ubiq_sets[[ubiq_setnm]]
tad_dist_df = tad_dist_df_allgenes%>%filter(i%in%ubiq_set,j%in%ubiq_set)
matad_dist_df <- make_matchdistdf(tad_dist_df,tadsetnm,tadgrpcols,distfilt = 1e6)%>%
				filter_bytadsize('1mb')%>%
				identity

hicpairdata <- lapply(hicchrs,function(hicchr){
	message(hicchr)
	fallgenegr_bins = resize(fallgenegr,1,'center')
	hic_res = 10e3
	newpos = start(fallgenegr_bins) %>% `/`(hic_res) %>% round%>%`*`(hic_res)
	fallgenegr_bins = GenomicRanges::shift(fallgenegr_bins,newpos - start(fallgenegr_bins))
	hicchrnm = paste0('chr',hicchr)
	chrgbins <- fallgenegr_bins%>%subset(seqnames%in%hicchrnm)
	chrgbins <- chrgbins%>%subset(gene_name %in% c(matad_dist_df$i,matad_dist_df$j))
	#hic info for that chr
	chrhic <- straw('KR',hicfile,chr1loc=paste0(hicchr,''),
		chr2loc=paste0(hicchr,''),'BP',hic_res,matrix='oe')
	#only bins with genes in them 
	chrhic <- chrhic[(chrhic$x %in% start(chrgbins)) & chrhic$y %in% start(chrgbins),]
	genebin <- tibble(x=start(chrgbins),gene_name=chrgbins$gene_name)
	chrhic <- chrhic%>%filter((y - x )<1e7)
	# browser()
	chr_pairs_hicdf = matad_dist_df%>%
		# sample_n(10)%>%
		select(-weight,-highcor,-tad,-btad)%>%
		filter(i %in% chrgbins$gene_name)%>%
		filter(dist<1e7)%>%
		# head(1000)%>%
		inner_join(genebin,by=c(i='gene_name'))%>%
		inner_join(genebin%>%select(y=x,gene_name),by=c(j='gene_name'))%>%
		inner_join(chrhic,by=c('x','y'))%>%
		select(-x,-y)
	gc()
	message('.')
	chr_pairs_hicdf
})
hicpairdata %<>% bind_rows
hicpairdata%>%saveRDS('data/hicpairdata.rds')
# save.image('data/3_tadsigplot.Rdata')

# c("i", "j", "dist", "highcor", "tad", "btad", "cor", "tadgrp", 
# "weight", "x", "y", "counts")
# chr_pairs_hicdf%>%{cor.test(.$dist,log(.$counts))}
# chr_pairs_hicdf%>%{cor.test(.$cor,log(.$dist))}
# chr_pairs_hicdf%>%{cor.test(.$cor,log(.$counts))}
# chr_pairs_hicdf%>%{txtplot(.$cor,log(.$counts))}
# lm(data=chr_pairs_hicdf,cor ~ log2(dist))%>%anova
# lm(data=chr_pairs_hicdf,cor ~ log2(dist)+log2(counts))%>%anova
# lm(data=chr_pairs_hicdf%>%filter(dist>10e3,dist<10e7),cor ~ log2(dist)+log2(counts))%>%anova
# lm(data=chr_pairs_hicdf,cor ~ log2(dist))%>%anova
# hicpairdata%>%colnames%>%dput
# stopifnot(colnames(chr_pairs_hicdf)==c('tad','btad','dist','tadgrp','tadsrc','tadsig'))



hicpairdata$counts%>%is.finite%>%table
# hicpairdata$dist%>%is.na%>%table
# hicpairdata$dist%>%log10%>%is.finite%>%table

# hicpairdata$counts%<>%pmax()
runfolder = paste0('plots/',
				'ugrp_',ubiq_setnm,'_',
				'tadsize_',tadsetnm,'_',
				'tadsrc_',tadsourcenm,'/'
			)
			dir.create(runfolder,showWarn=F,recursive=T)
{
chr_pairs_hicdf=hicpairdata
make_pair_tadsig_plot<-function(chr_pairs_hicdf,runfolder,tadgrpcols){
				if(!'tadsrc'%in%colnames(chr_pairs_hicdf))chr_pairs_hicdf$tadsrc='.'
				library(zoo)
				plotfile<-here(paste0(runfolder,'tad_sig_dist_',ubiq_setnm,'.pdf'))
				mymanum=manum
				pdf(plotfile)
				plot = chr_pairs_hicdf%>%
					filter(log10(dist)<7)%>%
					group_by(tadgrp)%>%
					# filter((tadgrp=='Within TAD')|tadgrp==tadgrp_match)%>%
					filter(is.finite(counts))%>%
					arrange(dist,by_group=TRUE)%>%
					mutate(tadsig = rollmean(counts,mymanum,na.pad=TRUE))%>%
					# filter(log10(dist)>5.75-0.01,log10(dist)<5.75+0.01)
					# mutate(tadsig = 2^rollmean(log2(counts),manum,na.pad=TRUE))%>%
					ggplot(.,aes(x=log10(dist),y = tadsig,color=tadgrp))+
					# xlab('geomean o/e')+
					ylab(str_interp('rollmean KR norm o/e'))+
					geom_line()+
					facet_grid(~tadsrc)+
					scale_color_manual(values = tadgrpcols)+
					coord_cartesian(xlim=c(4.5,6))+
					# scale_y_continuous('mean count count value')+
					ggtitle(paste0('Tad signal as a function of distance'))+
					theme_bw()
				print(plot)
				dev.off()
				message(plotfile)
			}
make_pair_tadsig_plot(hicpairdata,runfolder,tadgrpcols)
# make_pair_tadsig_plot(chr_pairs_hicdf,runfolder,tadgrpcols)
}
#do genes that are far from one another have a lower than expected 