## CHECK script tomorrow is it dumping out correct figures ?

genes<-c() ## change
pdf.file= ## change
log.dir.stub= ## change
num.per.page<-20
pdf(file=pdf.file,paper="a4r",width = 0, height = 0)
for (gene in genes){
	log.dir<-paste(log.dir.stub,gene,'/log/',sep="")
	files<-list.files(path=log.dir,pattern="*.count",full.names=TRUE)
	
	stats<-do.call("rbind",lapply(files,function(x){
		t<-read.table(file=x,header=FALSE,sep="\t",fill=TRUE,)
    t$V5<-NULL
    plate<-unique(as.character(t$V1))
    t<-as.data.frame(t(t[,c(3,4)]))
    names(t)<-as.character(unname(unlist(t[1,])))
    t<-t[-1,]
    rownames(t)<-plate
    t
		}
		))
  stats[,names(stats)]<-as.numeric(as.character(unlist(stats[,names(stats)])))
  ylimit<-c(0,max(rowSums(stats))*1.05)
  rem<-nrow(stats) %% num.per.page
  whole<-nrow(stats) %/% num.per.page
  if(whole<1){
  	stats$group<-1
  }else{
  	stats$group<- c(unlist(lapply(seq(1:whole),rep,times=num.per.page)),rep(whole+1,rem))
  }
 
  sapply(split(stats,stats$group),function(x){
  	title=paste("Read count distributions through pipeline checkpoints",gene)
  	x$group<-NULL
  	pm<-t(as.matrix(x))
  	barplot(pm,names=row.names(x),las=2,
  	ylim=ylimit,
  	col=rainbow(length(names(x))),legend = rownames(pm),main=title)
  })
}
dev.off();
