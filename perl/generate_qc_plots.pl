use strict;
use File::Temp qw/tempfile/;

my $gf = shift;
my $analysis_dir=shift;
my $out_dir=shift;

my $R_INSERT_SIZE=<<EOR;
genes<-scan('$gf',character())
log.dir.stub='$analysis_dir' ## change
pdf.file='$out_dir/insert_qc.pdf' ## change 

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
	stats<-do.call("rbind",lapply(split(stats,stats\$length),function(x){
		x[!duplicated(x\$length),]\$count<-sum(x\$count)
		x[1,]
	}))
	names(stats)<-c('plate','prod','type','length','count')
	r.df<-stats[stats\$count>500 | stats\$type =='PASSED',]
	r.df<-r.df[order(r.df\$count,decreasing=TRUE),]
	title=paste("Insert size distribution of processed reads for",gene)	
	b<-barplot(r.df\$count,col=ifelse(r.df\$type=='PASSED','green','grey'),
		ylim=c(0,max(r.df\$count)*1.1),ylab="Read counts",main=title)
	par(srt=90)
	add.label.plot(b,r.df,'count','col','length')
}
dev.off();

EOR

my $R_QC_WELL=<<EOR;
genes<-scan('$gf',character())
log.dir.stub='$analysis_dir'
pdf.file='$out_dir/qc_well.pdf'
num.per.page<-20
pdf(file=pdf.file,paper="a4r",width = 0, height = 0)
for (gene in genes){
	log.dir<-paste(log.dir.stub,gene,'/log/',sep="")
	files<-list.files(path=log.dir,pattern="*.count",full.names=TRUE)
	
	stats<-do.call("rbind",lapply(files,function(x){
		t<-read.table(file=x,header=FALSE,sep="\t",fill=TRUE,)
    t\$V5<-NULL
    plate<-unique(as.character(t\$V1))
    t<-as.data.frame(t(t[,c(3,4)]))
    names(t)<-as.character(unname(unlist(t[1,])))
    t<-t[-1,]
    rownames(t)<-plate
    t
		}
		))
  stats[,names(stats)]<-as.numeric(as.character(unlist(stats[,names(stats)])))
  ylimit<-c(0,max(rowSums(stats))*1.05)
  rem<-nrow(stats) \%\% num.per.page
  whole<-nrow(stats) \%/\% num.per.page
  if(whole<1){
  	stats\$group<-1
  }else{
  	stats\$group<- c(unlist(lapply(seq(1:whole),rep,times=num.per.page)),rep(whole+1,rem))
  }
 
  sapply(split(stats,stats\$group),function(x){
  	title=paste("Read count distributions through pipeline checkpoints",gene)
  	x\$group<-NULL
  	pm<-t(as.matrix(x))
  	barplot(pm,names=row.names(x),las=2,
  	ylim=ylimit,
  	col=rainbow(length(names(x))),legend = rownames(pm),main=title)
  })
}
dev.off();

EOR



if(! -e $gf || !-e $analysis_dir || !-e $out_dir){
	print "USAGE $0 genelist analysis_dir output_dir\n";
	exit(1);
}

(my $RSCRIPT=`which Rscript`)=~s/\n//;
if(!$RSCRIPT){
	print "[ERROR] require Rscript to run\n";
	exit(1);
}

my ($fh,$fname) = &get_tmp_file();
print $fh $R_QC_WELL;
print $fh $R_INSERT_SIZE;
close($fh);
my $cmd = "$RSCRIPT --vanilla $fname > $out_dir/generate_qc_plot_error.txt 2>&1";
`$cmd`;
unlink($fname);

sub get_tmp_file{
	my $dir = shift;
	my $template = "GEN_METH_PLOTS_XXXXX";
	return tempfile($template, DIR => $dir) if $dir;
	return tempfile($template, TMPDIR => 1);
}




