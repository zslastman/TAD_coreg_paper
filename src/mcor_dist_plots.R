use_vsd=FALSE


for(ubiq_setnm in names(ubiq_sets)[1]){
	for(tadsetnm in tadsetnms[2]){

	vsdstatus = ifelse(use_vsd,'withVSD_','noVSD_')
	runfolder = paste0('plots/',
		vsdstatus,
		ubiq_setnm,'_',
		tadsetnm,'/'
	)
	dir.create(runfolder,showWarn=F,recursive=T)

	allgenetpmmat2use = if(use_vsd)nallgenetpmmat else allgenetpmmat


	ubiq_set=ubiq_sets[[ubiq_setnm]]
	fallgenegr = allgenegr%>%subset(gene_name %in% ubiq_set)

	allgenetpmmat2use = allgenetpmmat[fallgenegr$gene_name,]

	if(tadsetnm=='1mb'){
		tads2use = alltads%>%subset(width>1e6)
	}else{
		tads2use = alltads
	}
	################################################################################
	########
	################################################################################
		
	tadovdf = fallgenegr%>%findOverlaps(tads2use,type='within')%>%as.data.frame%>%set_colnames(c('gene','tad'))

	tadofinterest = tadovdf%>%filter(fallgenegr$gene_name[gene]%>%str_detect('Fat1|Zfp42$'))%>%.$tad%>%unique
	
	netpmmat2use%<>%set_colnames(NULL)

	tadovdf$gnm = fallgenegr$gene_name[tadovdf$gene]

	intratadcors = map_df(unique(tadovdf$tad)%>%head,function(itad){
		gns = tadovdf$gnm[tadovdf$tad==itad]
		if(length(gns)<2)return(NULL)
		cors = oallcors[gns,gns]%>%.[lower.tri(.)]
		tibble(tad=itad,cor=cors)
	})

	#now plot
	tadofinterestmean = intratadcors%>%filter(tad==tadofinterest)%>%.$cor%>%mean
	#Fat1 (244,133,32); Triml2 (0,205,109); Zfp42 (0,205,109) - colour the distribution according to these colours	add line on each plot for other two genes using indicated colours. 
	#also add line for average FAT1 TAD coregulation (96,142,202)
	
	plotfile<- here(paste0(runfolder,'intratadcors','.pdf'))
	pdf(plotfile,w=5,h=5)
	print(
		intratadcors%>%
		group_by(tad)%>%
		# mutate(mcor=mean(cor>quantile(cor,0.75,na.rm=T)))%>%
		summarise(mcor=mean(cor),n=n())%>%
		ungroup%>%arrange(mcor)%>%mutate(tad=as_factor(tad))%>%
		# group_by(tad)%>%summarise()%>%
		# ggplot(.,aes(x=as.numeric(tad),y=cor))+
		ggplot(.,aes(x=mcor))+
		# geom_linerange()%>%
		# geom_point()+
		geom_histogram(fill=)+
		scale_y_continuous(paste0('Frequency'))+
		scale_x_continuous(paste0('Mean(Gene-Gene Correlation)'))+
		ggtitle(paste0('Distribution of average Gene-Gene cor\n
			dashed line is Fat1 Tad'))+
		geom_vline(xintercept=tadofinterestmean,color = I(96,142,202),linetype=2)+
		theme_bw()
	)
	dev.off()
	message(normalizePath(plotfile))

}}