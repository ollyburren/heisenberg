genes<-c('GENE') ## change
log.dir.stub='' ## change
pdf.file='' ## change 

add.label.plot<-function(b,df,y.name,col.name,region.name){
		text(x=b,y=df[[y.name]],labels=df[[region.name]],
			cex=1,adj=c(-0.5,0.5))
}

pdf(file=pdf.file,paper="a4r",width = 0, height = 0)
for (gene in genes){
	log.dir<-paste(log.dir.stub,gene,'/log/',sep="")
	files<-list.files(path=log.dir,pattern="*.log",full.names=TRUE)

	stats<-do.call("rbind",lapply(files,function(x){
		t<-read.table(file=x,header=FALSE,sep="\t")
	}))
	
	names(stats)<-c('plate','prod','type','length','count')	
	stats<-do.call("rbind",lapply(split(stats,stats$length),function(x){
		x[!duplicated(x$length),]$count<-sum(x$count)
		x[1,]
	}))
	names(stats)<-c('plate','prod','type','length','count')
	r.df<-stats[stats$count>500 | stats$type =='PASSED',]
	r.df<-r.df[order(r.df$count,decreasing=TRUE),]
	title=paste("Insert size distribution of processed reads for",gene)	
	b<-barplot(r.df$count,col=ifelse(r.df$type=='PASSED','green','grey'),
		ylim=c(0,max(r.df$count)*1.1),ylab="Read counts",main=title)
	par(srt=90)
	add.label.plot(b,r.df,'count','col','length')
}
dev.off();
