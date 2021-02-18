
################################################################################
########Table of info about TADS - print this is independent of ubiq subsetting
########and of tad size
################################################################################

neurotads = import('ext_data/Neuron_Cortical.Bonev_2017-raw.domains.bed')
esctads = import('ext_data/mESC.Bonev_2017-raw.domains.bed')
seqlevels(neurotads) %<>% paste0('chr',.)
seqlevels(esctads) %<>% paste0('chr',.)

alltads = neurotads


ubiq_sets = list(
	all=names(is_ubiq),
	ubiq=names(is_ubiq)[is_ubiq],
	non_ubiq=names(is_ubiq)[!is_ubiq]
)	
tadfiles = Sys.glob('ext_data/tadfiles/*.bed')%>%setNames(.,str_extract(basename(.),'[^.]+'))
names(tadfiles)=c('mESC',"mm10_CN_50K", "mm10_E114limb_Q30_50K", "mm10_ES_50K", "Neuron_Cortical")
tadfiles = tadfiles%>%map(import)%>%lapply(function(gr){
	if(is.null(gr$name)) gr else gr[!gr$name %in% c('gap','boundary')]
})
# tadfiles%<>%lapply(function(tads){seqlevels(tads)})
#
for(i in seq_along(tadfiles)) seqlevelsStyle(tadfiles[[i]])<-'UCSC'
for(i in seq_along(tadfiles)) tadfiles[[i]]$big <- width(tadfiles[[i]])>1e6
#
tadovdfs = list()
tadsourcenm=names(tadfiles)[1]
for(tadsourcenm in names(tadfiles)){
	#
	tadgr = tadfiles[[tadsourcenm]]
	#
	tadovdf = allgenegr%>%
		findOverlaps(tadgr,type='within')%>%
		as.data.frame%>%
		set_colnames(c('gene','tad'))
	#
	tadovdf$gene_name = allgenegr$gene_name[tadovdf$gene]
	tadovdf$is_ubiq = is_ubiq[tadovdf$gene_name]
	tadovdf$big = tadgr$big[tadovdf$tad]
	#store this variable
	tadovdfs[[tadsourcenm]]=tadovdf
	#
	tadgnms = tadovdf%>%group_by(tad,is_ubiq)%>%summarise(gnms=paste0(collapse=',',gene_name))%>%
		mutate(is_ubiq=ifelse(is_ubiq,'ubiquitous','non_ubiquitious'))%>%
		spread(is_ubiq,gnms)
	tadovdfpairs = tadovdf%>%left_join(tadovdf,suffix=c('','2'),by='tad')
	#
	tadgnms$chr = as.character(seqnames(tadgr)[tadgnms$tad])
	tadgnms$start = as.character(start(tadgr)[tadgnms$tad])
	tadgnms$end = as.character(end(tadgr)[tadgnms$tad])
	tadgnms$size = as.character(as.numeric(tadgnms$end) - as.numeric(tadgnms$start)+1)
	tadgnms$n_u = tadgnms$ubiquitous%>%str_count(',')%>%add(1)%>%replace_na(0)
	tadgnms$n_nu = tadgnms$non_ubiquitious%>%str_count(',')%>%add(1)%>%replace_na(0)
	tadgnms$n_tot = tadgnms$n_nu+tadgnms$n_u
	#
	tadofinterest = tadovdf%>%filter(allgenegr$gene_name[gene]%>%str_detect('Fat1|Zfp42$'))%>%.$tad%>%unique
	allgenetpmmat%<>%set_colnames(NULL)
	#	
	gset=ubiq_sets[1]
	gsetnm=ubiq_sets%>%names%>%.[1]
	mcorsets = imap(ubiq_sets%>%setNames(c('mcor','mcor_ubi','mcor_nubi')),function(gset,gsetnm){
		intratadcors = tadovdf%>%
			filter(gene_name%in%gset)%>%
			{split(.$gene,.$tad)}%>%
			.[map_dbl(.,length)>1]%>%
			map(function(x)cor(t(allgenetpmmat[allgenegr$gene_name[x],])))%>%
			map(~ .[lower.tri(.)])%>%
			map(unlist)%>%
			enframe('tad','cor')%>%
			unnest(cor)
		intratadcors$tad%<>%as.numeric
		mcorvals = intratadcors%>%group_by(tad)%>%summarise(mcor=mean(cor))
		mcorvals%<>%rename(!!gsetnm:=mcor)
		list(intratadcors,mcorvals)
	})
	mcorvals = mcorsets%>%map(2)%>%reduce(full_join,by='tad')
	intratadcors = mcorsets%>%map(1)
	tadinfodf = tadgnms%>%left_join(mcorvals)
	#
	tadinfofile = paste0('tables/',tadsourcenm,'.tadinfo.tsv')
	#
	tadinfodf%>%
		select(tad,non_ubiquitious,n_nu,ubiquitous,n_u,n_tot,
				chr,start,end,size,mcor_ubi,mcor_nubi,mcor)%>%
		write_tsv(tadinfofile)
	tadinfofile%>%normalizePath%>%message
	################################################################################
	########Now plot distribution within the tad of interest
	################################################################################
	myrgb<-function(r,g,b){
		rgb(r/256,g/256,b/256)
	}
	goicols = list(
		Fat1=myrgb(244,133,32),
		Triml2=myrgb(0,205,109),
		Zfp42=myrgb(0,205,109)
	)
	tadofinterestmeancol = myrgb(96,142,202)


	tadofinterest = tadovdf%>%filter(gene_name%>%str_detect('Fat1|Zfp42$'))%>%.$tad%>%unique
	if(length(tadofinterest)==0){
		message('no tad of interest here - the genes in question dont overlap these tads')
	}else{
		for(i_ubiq_set in seq_along(ubiq_sets)) {
			gsetname=names(ubiq_sets)[i_ubiq_set]

			#	
			intratadcors = mcorsets[[i_ubiq_set]][[1]]
			#now plot
			tadofinterestmean = intratadcors%>%filter(tad==tadofinterest)%>%.$cor%>%mean
			#Fat1 (244,133,32); Triml2 (0,205,109); Zfp42 (0,205,109) - colour the distribution according to these colours	add line on each plot for other two genes using indicated colours. 
			#also add line for average FAT1 TAD coregulation (96,142,202)
			plotfile<- here(paste0('plots/intratadcors_',gsetname,'_',tadsourcenm,'.pdf'))
			pdf(plotfile,w=5,h=5)
			print(
				intratadcors%>%
				group_by(tad)%>%
				summarise(mcor=mean(cor),n=n())%>%
				ungroup%>%
				arrange(mcor)%>%
				mutate(tad=as_factor(tad))%>%
				ggplot(.,aes(x=mcor))+
				geom_histogram()+
				scale_y_continuous(paste0('Frequency'))+
				scale_x_continuous(paste0('Mean(Gene-Gene Correlation)'))+
				ggtitle(paste0('Distribution of average Gene-Gene cor\n
					dashed line is Fat1 Tad'))+
				geom_vline(xintercept=tadofinterestmean,color = I(tadofinterestmeancol),linetype=2)+
				theme_bw()
			)
			dev.off()
			message(normalizePath(plotfile))
		}

		################################################################################
		########Now individual gene plots
		################################################################################
		gois = c('Fat1','Zfp42','Triml2')

		use_vsd=FALSE
		vsdstatus = ifelse(use_vsd,'withVSD_','noVSD_')
		allcorfile=here(paste0('data/',vsdstatus,'allcors.rds'))
			
		allcors<-readRDS((allcorfile))
		igene = gois[1]

		tadnumber = tadovdf%>%filter(gene_name=='Fat1')%>%.$tad
		
		myrgb<-function(r,g,b){
			rgb(r/256,g/256,b/256)
		}
		goicols = list(
			Fat1=myrgb(244,133,32),
			Triml2=myrgb(0,205,109),
			Zfp42=myrgb(0,205,109)
		)
		tadofinterestmeancol = myrgb(96,142,202)

		# Fat1 (244,133,32); Triml2 (0,205,109); Zfp42 (0,205,109) - colour the distribution according to 
		# these colours	add line on each plot for other two genes using indicated colours. 
		# also add line for average FAT1 TAD coregulation (96,142,202)

		for(igene in gois){
			#
			kgenes = gois%>%setdiff(igene)
			jcors = allcors[igene,kgenes]
			kgenes = kgenes[order(jcors)]
			jcors = jcors[order(jcors)]
			#
			plotfile<- here(paste0('plots/cordist_',igene,'_',tadsourcenm,'.pdf'))
			pdf(plotfile,w=5,h=5)
			print(
				allcors[igene,]%>%
				enframe('gene_name','cor')%>%
				ggplot(.,aes(x=cor))+
				geom_histogram(fill=I(goicols[igene]))+
				scale_y_continuous(paste0('Frequency'))+
				scale_x_continuous(paste0('Mean(Gene-Gene Correlation)'))+
				# ggtitle(paste0('Distribution of cors for ',igene,' cor\nred line is mean cor for the tad\n black dasheded lines are ',paste0(kgenes,collapse=','),' '))+
				ggtitle(paste0('Distribution of cors for ',igene))+
				geom_vline(xintercept=jcors[1],linetype=2,color=I(goicols[names(jcors)[1]]))+
				geom_vline(xintercept=jcors[2],linetype=2,color=I(goicols[names(jcors)[2]]))+
				geom_vline(xintercept=tadofinterestmean,linetype=3,color=I(tadofinterestmeancol))+
				theme_bw()
			)
			dev.off()
			message(normalizePath(plotfile))
			#
		}

	}
}
tadovdfs %>% saveRDS(here('data/tadovdfs.rds'))
#TADs Table	tadinfo.tsv	Columns - TAD; non_ubi_gene_names; ubi_gene_names; number_non_ubi_genes; number_ubi_genes; total_genes; chr; start; end; TAD_size; mcorr_non_ubi; mcorr_ubi; mcorr_all)
