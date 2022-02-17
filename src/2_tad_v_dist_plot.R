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

select<-dplyr::select



lowcorcol = (c(217,217,217)/256)%>%{rgb(.[1],.[2],.[3])}
highcorcol = (c(96,142,202)/256)%>%{rgb(.[1],.[2],.[3])}
tadcolor = (c(96,142,202)/256)%>%{rgb(.[1],.[2],.[3])}
matchpaircolor = (c(217,217,217)/256)%>%{rgb(.[1],.[2],.[3])}
gpairsamechr = (c(44,105,151)/256)%>%{rgb(.[1],.[2],.[3])}
gpairallcol = (c(0,0,0)/256)%>%{rgb(.[1],.[2],.[3])}


#this is just the .9 quantile
highcorlim = 0.38

if(!exists('allgenetpmmat')) load('data/tad_coreg_preprocess.Rdata')



ubiq_sets = list(
	all=names(is_ubiq),
	ubiq=names(is_ubiq)[is_ubiq],
	non_ubiq=names(is_ubiq)[!is_ubiq]
)

ubiq_setnm='all'
use_vsd= FALSE
tadsetnms = c('all','1mb')
tadsetnm='1mb'
tadsourcenm='Neuron_Cortical'

tadfiles = Sys.glob('ext_data/tadfiles/*.bed')%>%setNames(.,str_extract(basename(.),'[^.]+'))
names(tadfiles)=c('mESC',"mm10_CN_50K", "mm10_E114limb_Q30_50K", "mm10_ES_50K", "Neuron_Cortical")
tadfiles = tadfiles%>%map(import)%>%lapply(function(gr){
	if(is.null(gr$name)) gr else gr[!gr$name %in% c('gap','boundary')]
})
tadovdfs <- readRDS(here('data/tadovdfs.rds'))


message('loading correlation data')
vsdstatus = ifelse(use_vsd,'withVSD_','noVSD_')
allcorfile=here(paste0('data/',vsdstatus,'allcors.rds'))
existsallcorfile = file.exists(allcorfile)
if(!exists('oallcors')){
	if(!existsallcorfile){
		message('getting all cors')
		oallcors <- allgenetpmmat2use%>%t%>%cor
		oallcors%<>%Matrix
		saveRDS(oallcors,(allcorfile))
	}else{
		oallcors<-readRDS((allcorfile))
	}
}


#
for(i in seq_along(tadfiles)) seqlevelsStyle(tadfiles[[i]])<-'UCSC'
for(i in seq_along(tadfiles)) tadfiles[[i]]$big <- width(tadfiles[[i]])>1e6
#
tadpairlist=list()
for(tadsourcenm in names(tadfiles)){
	tadovdf=tadovdfs[[tadsourcenm]]
	tadgr = tadfiles[[tadsourcenm]]


	ubiq_set=ubiq_sets$all

	fallgenegr = allgenegr%>%subset(gene_name %in% ubiq_set)

	allgenetpmmat2use = if(use_vsd)nallgenetpmmat else allgenetpmmat


	g1 = fallgenegr$gene_name[1];g2=fallgenegr$gene_name[2]
	allcors = oallcors[fallgenegr$gene_name,fallgenegr$gene_name]
	stopifnot(cor(allgenetpmmat2use[g1,],allgenetpmmat2use[g2,])==allcors[g2,g1])


	{
	#deffine the gene tad overlaps in a big df
	gstarts = fallgenegr%>%resize(1,'start')
	tadovdf  = gstarts%>%
		findOverlaps(tadgr)%>%
		as.data.frame%>%
		set_colnames(c('gene','tad'))
	tadovdf$gene_name = fallgenegr$gene_name[tadovdf$gene]
	tadovdf$is_ubiq = is_ubiq[tadovdf$gene_name]
	tadovdf$big = tadgr$big[tadovdf$tad]

	#work out distance from the left and right tad boundary for each of these
	tadovdf$ldist = distance(gstarts[tadovdf$gene],resize(tadgr,1,'start')[tadovdf$tad])
	tadovdf$rdist = distance(gstarts[tadovdf$gene],resize(tadgr,1,'end')[tadovdf$tad])
	#define each gene as being a boundary gene
	tadovdf%<>%group_by(tad)%>%mutate(isredge = rdist ==min(rdist))
	tadovdf$rdist <- NULL
	tadovdf%<>%group_by(tad)%>%mutate(isledge = ldist ==min(ldist))
	tadovdf$ldist <- NULL
	
	#Then we can define the pairs as being double boundary genes
	tadovdfpairs = tadovdf%>%
		left_join(tadovdf%>%select(tad,gene,isledge,isredge,gene_name),suffix=c('','2'),by='tad')
	tadovdfpairs%<>%mutate(edgepair = (isledge & isredge2) | (isredge & isledge2))
	tadovdfpairs%<>%select(-isredge,-isledge,-isledge2,-isredge2)
	}

	# gfat1 = fallgenegr$gene_name%>%str_detect('Fat1')%>%which
	# zfp42 = fallgenegr$gene_name%>%str_detect('Zfp42$')%>%which
	# message(gfat1)
	# message(zfp42)
	# tadovdf%>%filter(gene_name=='Fat1')

	# tadovdf%>%filter(gene==gfat1)
	# tadovdf%>%filter(gene==zfp42)
	# fallgenegr%>%subset(gene_name=="Fat1")%>%subsetByOverlaps(tadgr)
	# fallgenegr%>%subset(gene_name=="Zfp42")%>%subsetByOverlaps(tadgr)

	#matrix denoting coincidence in all tads
	make_tad_coincidenmat<-function(fallgenegr,tadovdfpairs){
		tad_coincidenmat = Matrix(FALSE,ncol=length(fallgenegr),nrow=length(fallgenegr),sparse=TRUE)
		#only mark upper triangle, so count once
		tad_coincidenmat[as.matrix(tadovdfpairs%>%select(gene,gene2)%>%filter(gene<gene2))] <- TRUE
		tad_coincidenmat = tad_coincidenmat | t(tad_coincidenmat)
		diag(tad_coincidenmat) = TRUE
		#
		tad_coincidenmat <- tad_coincidenmat%>%
			set_colnames(fallgenegr$gene_name)%>%
			set_rownames(fallgenegr$gene_name)
		tad_coincidenmat
	}
	tad_coincidenmat <- make_tad_coincidenmat(fallgenegr,tadovdfpairs)
	#
	tadovdfpairs$big <- tadgr$big[tadovdfpairs$tad]
	btad_coincidenmat <- make_tad_coincidenmat(fallgenegr,tadovdfpairs%>%subset(big))
	# tad_coincidenmat['Fat1','Zfp42']
	# btad_coincidenmat['Fat1','Zfp42']
	# tadovdf%>%filter(gene_name=='Fat1')
	#
	highcormat = Matrix(sparse=TRUE,allcors[fallgenegr$gene_name,fallgenegr$gene_name]>highcorlim)


	ichr='chr8'
	allchrs = as.character(unique(seqnames(fallgenegr)))
	message('create distance - cor data frame')
	tad_dist_df_allgenes = lapply(X=allchrs,
	#	tad_coincidenmat,btad_coincidenmat,highcormat,fallgenegr,
		FUN=function(ichr){
		message('.')
		#
		chrgenegr = fallgenegr%>%subset(seqnames==ichr)
		if(!(length(chrgenegr)>1)) return(NULL)
		chrgnms = chrgenegr%>%.$gene_name
		chrstarts = chrgenegr%>%resize(1,'center')%>%start
		chrdists = chrstarts%>%dist(method='manhattan')
		# cgfat1 = chrgenegr$gene_name%>%str_detect('Fat1')%>%which
		# czfp42 = chrgenegr$gene_name%>%str_detect('Zfp42$')%>%which
		#
		chrhighcormat = highcormat[chrgnms,chrgnms]
		chrtad_coincidenmat = tad_coincidenmat[chrgnms,chrgnms]
		bigtad_coincidenmat = btad_coincidenmat[chrgnms,chrgnms]
		chrcormat = allcors[chrgnms,chrgnms]
		#
		chrcormat = Matrix(as.matrix(chrcormat),sparse=TRUE)
		#
		#
		tad_cor_df = summary(Matrix(as.matrix(chrdists),sparse=TRUE))%>%
			as.data.frame%>%
			select(i,j,dist=x)%>%
			left_join(summary(chrhighcormat)%>%rename('highcor':=x),c('i','j'))%>%
			left_join(summary(chrtad_coincidenmat)%>%rename('tad':=x),c('i','j'))%>%
			left_join(summary(bigtad_coincidenmat)%>%rename('btad':=x),c('i','j'))%>%
			left_join(summary(chrcormat)%>%rename('cor':=x),c('i','j'))
		#
		tad_cor_df%<>%filter(i!=j)
		tad_cor_df$i = chrgenegr$gene_name[tad_cor_df$i]
		tad_cor_df$j = chrgenegr$gene_name[tad_cor_df$j]
		#
		tad_cor_df$tad%<>%replace_na(FALSE)
		tad_cor_df$highcor%<>%replace_na(FALSE)
		tad_cor_df$btad%<>%replace_na(FALSE)
		#
		tad_cor_df
	})
	tad_dist_df_allgenes%<>%bind_rows
	#
	edgepairdf <- tadovdfpairs%>%ungroup%>%filter(edgepair)%>%select(gene_name,gene_name2,edgepair)
	tad_dist_df_allgenes%<>%left_join(edgepairdf,by=c('i'='gene_name','j'='gene_name2'))
	tad_dist_df_allgenes$edgepair%<>%replace_na(FALSE)

	tadsetnm = tadsetnms[2]
	ubiq_setnm = 'non_ubiq'

	for(ubiq_setnm in names(ubiq_sets)){

		ubiq_set=ubiq_sets[[ubiq_setnm]]

		tad_dist_df = tad_dist_df_allgenes%>%filter(i%in%ubiq_set,j%in%ubiq_set)

		for(tadsetnm in tadsetnms){
			{#DEBUG
			if(is.null(tadpairlist[[tadsourcenm]]))tadpairlist[[tadsourcenm]]=list()
			vsdstatus = ifelse(use_vsd,'withVSD_','noVSD_')
			runfolder = paste0('plots/',
				'ugrp_',ubiq_setnm,'_',
				'tadsize_',tadsetnm,'_',
				'tadsrc_',tadsourcenm,'/'
			)
			dir.create(runfolder,showWarn=F,recursive=T)

			tadgrp_pos='Within TAD'
			tadgrp_match='Not Within (any) TAD\n(distance matched)'
			tadgrp_nt='Non TAD'
			tadgrp_all='All'
			tadgrpcols = c(tadcolor,matchpaircolor,gpairsamechr,gpairallcol)%>%
				setNames(c(tadgrp_pos,tadgrp_match,tadgrp_nt,tadgrp_all))
			#
			make_matchdistdf<-function(tad_dist_df,tadsetnm,tadgrpcols,distfilt= Inf){
				tadgrps <- names(tadgrpcols)
				#
				if(tadsetnm=='1mb')tadsizecol=sym('btad') else tadsizecol=sym('tad')
				#
				distdensfun = tad_dist_df%>%filter(!!tadsizecol)%>%.$dist%>%log10%>%density(.,adjust=0.5)%>%approxfun
				distdensfun2 = function(x)replace_na(distdensfun(x),0)				
				tadmaxdist =  tad_dist_df%>%filter(!!tadsizecol)%>%.$dist%>%max
				matchedpairs = tad_dist_df%>%filter(!!tadsizecol)%>%filter(log10(dist)<log10(tadmaxdist))
				matchedpairs$weight = log10(matchedpairs$dist)%>%distdensfun2
				#
				tadnum = tad_dist_df$tad%>%sum
				tadnum = tadnum*10
				#	
				#
				matcheddistdf = bind_rows(
					tad_dist_df%>%filter(!!tadsizecol)%>%mutate(tadgrp=tadgrps[1]),
					matchedpairs%>%sample_n(min(n(),tadnum),weight=weight,replace=TRUE)%>%mutate(tadgrp=tadgrps[2]),
					tad_dist_df%>%filter(dist<distfilt)%>%filter(!tad)%>%ungroup%>%sample_n(replace=TRUE,min(n(),tadnum))%>%mutate(tadgrp=tadgrps[3]),
					tad_dist_df%>%filter(dist<distfilt)%>%ungroup%>%sample_n(replace=TRUE,min(n(),tadnum))%>%mutate(tadgrp=tadgrps[4])
				)
			}
			#
			manum = 2e3
			#
			make_tadfreq_dplot<-function(matad_dist_df,runfolder,highcorlim,tadsizecol,manum,lowcorcol,highcorcol,pwidth=7){
				if(!'tadsrc'%in%colnames(matad_dist_df))matad_dist_df$tadsrc='.'
				library(zoo)
				plotfile<-here(paste0(runfolder,'allpair_flip_cutoff_ma','.pdf'))
				pdf(plotfile,width=pwidth)
				plot = matad_dist_df%>%
					ungroup%>%
					filter(log10(dist)<6)%>%
					group_by(highcor,tadsrc)%>%
					arrange(dist,by_group=TRUE)%>%
					mutate(tprop = rollmean(tad,manum,na.pad=TRUE))%>%
					ggplot(.,aes(x=log10(dist),y = tprop,color=highcor))+
					scale_color_manual(values=c('FALSE'=lowcorcol,'TRUE'=highcorcol))+
					geom_line()+
					facet_grid(.~tadsrc)+
					coord_cartesian(xlim=c(4.5,6))+
					scale_y_continuous('proportion of genes sharing TADs\n(moving average of 10k genes)')+
					ggtitle(paste0('Tad comembership prob as a function of distance,\nw vs without highcor\nmoving average ',manum,' points'))+
					theme_bw()
				print(plot)
				dev.off()
				message(plotfile)
			}
			filter_bytadsize<- function(matcheddistdf,tadsetnm){
				if(tadsetnm=='1mb'){
					matad_dist_df = matcheddistdf%>%filter((btad) | (!tad))
				}else{
					matad_dist_df = matcheddistdf
				}
			}
			matad_dist_df <- make_matchdistdf(tad_dist_df,tadsetnm,tadgrpcols,distfilt = 1e6)%>%
				filter_bytadsize(tadsetnm)%>%
				identity
			# matad_dist_df <- filter_bytadsize(matcheddistdf,tadsetnm)
			# matad_dist_df <- matcheddistdf
			matad_dist_df %>% make_tadfreq_dplot(runfolder,highcorlim,tadsizecol,manum,lowcorcol,highcorcol)

			#
			if(ubiq_setnm=='all'){tadpairlist[[tadsourcenm]][[tadsetnm]]=matad_dist_df%>%filter(tadgrp!='All')}
			#	
			make_highcorfreq_plot<-function(matad_dist_df,runfolder,tadgrpcols){
				if(!'tadsrc'%in%colnames(matad_dist_df))matad_dist_df$tadsrc='.'
				library(zoo)
				plotfile<-here(paste0(runfolder,'allpair_ma','.pdf'))
				pdf(plotfile)
				plot = matad_dist_df%>%
					filter(log10(dist)<6)%>%
					group_by(tadgrp)%>%
					arrange(dist,by_group=TRUE)%>%
					mutate(cprop = rollmean(cor>highcorlim,manum,na.pad=TRUE))%>%
					ggplot(.,aes(x=log10(dist),y = cprop,color=tadgrp))+
					geom_line()+
					facet_grid(~tadsrc)+
					scale_color_manual(values = tadgrpcols)+
					coord_cartesian(xlim=c(4.5,6))+
					scale_y_continuous('proportion of coexpressed genes\n(moving average of 10k genes)')+
					ggtitle(paste0('Coexpression as a function of distance,\nmoving average 5k points'))+
					theme_bw()
				print(plot)
				dev.off()
				message(plotfile)
			}
			# if((tadsetnm=='1mb') & (ubiq_setnm=='non_ubiq') & (tadsourcenm=='mm10_ES_50K'))browser()
			matad_dist_df%>%filter(!edgepair)%>%make_highcorfreq_plot(runfolder,tadgrpcols)
			}#DEBUG
			# make_matchdistdf<-make_matchdistdf
			matcheddistdf_alldist <- make_matchdistdf(tad_dist_df,tadsetnm,tadgrpcols)
			make_coreg_densplots <- function(matcheddistdf_alldist,thirdgrp=tadgrp_all,filtgrp=tadgrp_nt,highcorlim,runfolder,zoom=F){
				zoomfun = if(!zoom) NULL else coord_cartesian(xlim=c(highcorlim,1),ylim=c(0,1))
				zoomnm = if(!zoom) 'NULL' else '_zoom_'
				#
				grpsizes = matcheddistdf_alldist%>%filter(!tadgrp==filtgrp)%>%group_by(tadgrp)%>%tally
				grpsizes= grpsizes%>%apply(1,paste,collapse=':')
				cordensplot = matcheddistdf_alldist%>%
					filter(!tadgrp==filtgrp)%>%
					ggplot(.,aes(x=cor,color=tadgrp))+
					geom_density(alpha=I(0.8))+
					scale_color_discrete(name='')+
					ggtitle('Distribution of Co-regulation amongst all gene pairs')+
					zoomfun+
					geom_vline(xintercept=highcorlim,linetype='dashed')+
					theme_bw()+
					geom_text(y=Inf,x=Inf,data=tibble(label=paste0(collapse='\n',grpsizes)),aes(label=label),color=I('black'),vjust='top',hjust='right')
				ddensplot = matcheddistdf_alldist%>%
					filter(!tadgrp==filtgrp)%>%
					ggplot(.,aes(x=log10(dist),color=tadgrp))+
					geom_density(alpha=I(0.8))+
					scale_color_discrete(name='')+
					ggtitle('Distribution of distance amongst all gene pairs')+
					geom_vline(xintercept=highcorlim,linetype='dashed')+
					theme_bw()
					#
				plotfile<-here(paste0(runfolder,'cor_distribution_w',zoomnm,str_replace(thirdgrp,' ','_'),'.pdf'))
				pdf(plotfile)
				print(ggpubr::ggarrange(nrow=2,plotlist=list(cordensplot,ddensplot)))
				dev.off()
				message(plotfile)
			}
			make_coreg_densplots(matcheddistdf_alldist,tadgrp_all,tadgrp_nt,highcorlim,runfolder)
			make_coreg_densplots(matcheddistdf_alldist,tadgrp_all,tadgrp_nt,highcorlim,runfolder,zoom=TRUE)
			make_coreg_densplots(matcheddistdf_alldist,tadgrp_nt,tadgrp_all,highcorlim,runfolder)
			make_coreg_densplots(matcheddistdf_alldist,tadgrp_nt,tadgrp_all,highcorlim,runfolder,zoom=TRUE)
			#
			make_tad_fracbarplot <- function(tadcoregtbl,plotname,pcnum=100,ptitle='Fraction of TAD-sharing genes, given shared chromosome'){
				#
				nco_nt = tadcoregtbl[1,2] / sum(tadcoregtbl[1,])
				nco_t = tadcoregtbl[2,2] / sum(tadcoregtbl[2,])
				#	
				tadfracdf = tibble(
					coregulation = c('not co-expressed','co-expressed')%>%as_factor,
					within_tad = c(nco_nt,nco_t)
				)
				# if(pcnum==1000)pcstr = paste0('%','0.') else pcstr='%'
				brks = tadfracdf$within_tad%>%multiply_by(pcnum)%>%max%>%ceiling%>%{0:.}%>%divide_by(pcnum)
				labs = tadfracdf$within_tad%>%multiply_by(pcnum)%>%max%>%ceiling%>%{0:.}%>%paste0(.,'%')
				#now plot
				plotfile<-here(paste0(runfolder,plotname,'.pdf'))
				pdf(plotfile)
				print(
					tadfracdf%>%ggplot(aes(x=coregulation,y=within_tad,fill=coregulation))+
					stat_identity(geom='bar')+
					scale_fill_manual(values=c(lowcorcol,highcorcol))+
					# scale_y_continuous(name='Fraction of Pairs Sharing a TAD',breaks=brks,labels=labs)+
					scale_y_continuous(name='Fraction of Pairs Sharing a TAD')+
					ggtitle(ptitle)+
					theme_bw()
					)
				dev.off()
				message(normalizePath(plotfile))
			}
			tadcoregtbl <- table(cor = tad_dist_df$highcor,tad = tad_dist_df$tad)
			tadcoregtbl %>% make_tad_fracbarplot('coexpr_barplot')
			if(tadsetnm=='1mb'){
				tadmat = btad_coincidenmat
			}else{
				tadmat = tad_coincidenmat
			}
			#tad on columns of this matrix
			sampgenes = rownames(highcormat)%>%sample(5e3)
			tadcoregtbl2 = c(
				sum((!highcormat[sampgenes,sampgenes])&(!tadmat[sampgenes,sampgenes])),
				sum((highcormat[sampgenes,sampgenes])&(!tadmat[sampgenes,sampgenes])),
				sum((!highcormat[sampgenes,sampgenes])&(tadmat[sampgenes,sampgenes])),
				sum((highcormat[sampgenes,sampgenes])&(tadmat[sampgenes,sampgenes]))
			)%>%matrix(ncol=2)
			tadcoregtbl2%>%sweep(1,STAT=rowSums(.),FUN='/')
			tadcoregtbl2 %>% make_tad_fracbarplot('allpairs_coexpr_barplot',ptitle='Fraction of TAD-sharing genes')
	}
}
}

stop()
#

tadsourcecols = c(
	"mESC"='lightpink',
	"Neuron_Cortical"='red',
	"mm10_CN_50K"='green',
	"mm10_E114limb_Q30_50K"='blue',
	"mm10_ES_50K"='darkblue')           
#
# tadsourcecols = c("mm10_CN_50K_Within TAD"='#FF0000',
# 	"mm10_CN_50K_Not Within (any) TAD\n(distance matched)" = '#880000',
# "mm10_CN_50K_Non TAD" = '#220000', 
# "mm10_E114limb_Q30_50K_Within TAD" = '#00FF00', 
# "mm10_E114limb_Q30_50K_Not Within (any) TAD\n(distance matched)" = '#008800',
# "mm10_E114limb_Q30_50K_Non TAD" = '#002200', 
# "mm10_ES_50K_Within TAD" = '#0000FF', 
# "mm10_ES_50K_Not Within (any) TAD\n(distance matched)" = '#000088',
# "mm10_ES_50K_Non TAD" = '#000022'
# )
make_highcorfreq_comp_plot <- function(matad_dist_df,runfolder,tadgrpcols){
                # if(!'tadsrc'%in%colnames(matad_dist_df))matad_dist_df$tadsrc='.'
                library(zoo)
                plotfile<-here(paste0(runfolder,'allpair_ma','.pdf'))
                pdf(plotfile,w=21,h=7)
                ggdf = matad_dist_df%>%
                    filter(log10(dist)<6)%>%
                    group_by(tadgrp,tadsrc)%>%
                    arrange(dist,by_group=TRUE)%>%
                    mutate(cprop = rollmean(cor>highcorlim,manum,na.pad=TRUE))
                plot = ggdf%>%
                    ggplot(.,aes(x=log10(dist),y = cprop,color=tadgrp))+
                    geom_line()+
                   # ggplot(.,aes(x=log10(dist),y = as.numeric(cor>highcorlim),color=tadgrp))+
                    # geom_smooth()+
                    facet_grid(~tadsrc)+
                    scale_color_manual(values = tadgrpcols)+
                    coord_cartesian(xlim=c(4.5,6))+
                    scale_y_continuous('proportion of coexpressed genes\n(moving average of 10k genes)',limits=c(0,0.6))+
                    ggtitle(paste0('Tad comembership prob as a function of distance,\nw vs without highcor\nmoving average 5k points'))+
                    theme_bw()
                print(plot)
                dev.off()
                message(plotfile)
            }

################################################################################
########Now print pannelled plots for comparing tad sources
################################################################################
manum = 2e3	
'plots/tadsourcecomp/'%>%dir.create
#
for(ubiq_setnm in names(ubiq_sets)){
	ubiq_set = ubiq_sets[[ubiq_setnm]]
	for(tadsetnm in tadsetnms){
		runfolder = paste0('plots/tadsourcecomp/',
			'ugrp_',ubiq_setnm,'_',
			'tadsize_',tadsetnm
		)
		dir.create(runfolder,showWarn=F,recursive=T)
		tadsourcecomplist = tadpairlist%>%map(tadsetnm)%>%bind_rows(.id='tadsrc')
		tadsourcecomplist2 = tadsourcecomplist%>%filter(i%in%ubiq_set,j%in%ubiq_set)	
		if(tadsetnm=='1mb'){
				tadsourcecomplist2 = tadsourcecomplist2%>%filter((btad) | (!tad))
   		}
		tadsourcecomplist2%>%	
			make_highcorfreq_comp_plot(runfolder,tadgrpcols)
		tadsourcecomplist2%>%
			make_tadfreq_dplot(runfolder,highcorlim,tadsizecol,manum,lowcorcol,highcorcol,pwidth=21)
		}
}





#oka let's pull info from the hic file
# Definition	RGB colour
# cor<x?	96,142,202
# cor>x?	
# all TADs or >1Mb TADs	96,142,202
# pairs not in same tad that have been chosen by a weighted sampling of pairs using the TAD gene pair separation distribution as a guide	217,217,217
# all gene pairs on same chromosome	44,105,151
# all gene pairs on same chromosome	0,0,0