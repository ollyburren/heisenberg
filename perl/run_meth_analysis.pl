#!/usr/bin/perl

use strict;
use File::Basename;
use File::Find;
use File::Temp qw/tempfile/;
use POSIX;
use Config::IniFiles;
use Data::Dumper;
use Getopt::Long;
use Cwd 'abs_path';

my $CM_SCRIPT=dirname(abs_path($0))."/call_meth.pl";

if(! -e $CM_SCRIPT){
	die "Cannot find $CM_SCRIPT at $CM_SCRIPT, please check\n";
}

my $USAGE=<<EOL;
$0: Process demultiplexed data and call methylation sites.

$0 -out_dir output_dir --ini inifile.ini --gene_file gene.txt --log_dir log_dir 

	MANDATORY PARAMETERS
		out_dir|o:	path to a location to store intermediate and final analysis files
		ini|i: path to ini file containing configuration options.
		gene_file|g: file containing list of genes to process.
		log_dir|l path to location to store logging information
		data_dir|d absolute path to location of gzipped fastq files
		
	OPTIONAL
		help|h: print this message
		test|t: run one example to check configuration
		prefix|p: plate prefix (DEFAULT: '')
		
	Author: Olly Burren

Last updated 17/10/2013
EOL

####################
#OPTIONS PROCESSING#
####################

my ($out_dir,$ini,$genes,$log_dir,$test,$prefix,$data_dir,$help);

GetOptions(	
	"data_dir|d=s" => \$data_dir,
	"prefix|p=s" => \$prefix,
	"out_dir|b=s" => \$out_dir,
	"ini|i=s" => \$ini,
	"genes_file|g=s" => \$genes,
	"log_dir|l=s"=>\$log_dir,
	"test|t"=>\$test,
	"help|h"=>\$help);

my $ABORT_FLAG= 0 || $help;

unless($ABORT_FLAG){
	if(!$out_dir){
		print "[ERROR] Require --out_dir/-o option: Location for analysis files to be written\n";
		$ABORT_FLAG++;
	}elsif(!$genes){
		print "[ERROR] Require --genes/-g option: Path to list of genes to be analysed\n";
		$ABORT_FLAG++;
	}elsif(!$log_dir){
		print "[ERROR] Require --log_dir/-l option: Path to dir to store q logging info\n";
		$ABORT_FLAG++;
	}elsif(!$ini){
		print "[ERROR] Require --ini option: Path to ini\n";
		$ABORT_FLAG++;
	}elsif(!$data_dir){
		print "[ERROR] Require --data_dir/-d option: Path to data dir containing fq.gz files\n";
		$ABORT_FLAG++;
	}
}
	

unless($ABORT_FLAG){
	if(! -e $ini){
		print "[ERROR] Cannot find $ini ini file\n";
		$ABORT_FLAG++;
	}elsif(! -d $log_dir){
		print "[ERROR] Cannot find log dir $log_dir \n";
		$ABORT_FLAG++;
	}elsif(! -e $genes){
		print "[ERROR] Cannot find list of genes file $genes \n";
		$ABORT_FLAG++;
	}elsif(! -d $data_dir){
		print "[ERROR] Cannot find data_dir $data_dir \n";
		$ABORT_FLAG++;
	}
}

if($ABORT_FLAG){
	print $USAGE;
	exit(1);
}


my $QSUB_CMD="qsub -q all.q -v PERL5LIB=$ENV{PERL5LIB}";
my $cfg = new Config::IniFiles( -file => $ini);
my $TMP_DIR='/tmp/';
## Software only works on gzipped FASTQ
my $SUFFIX='.1.fq.gz';

open(GENES,$genes) || die "Cannot open gene file $genes\n";

while(<GENES>){
	chomp;
	next if /^#/;
  ## in case the file is DOS
  s/\r\n/\n/;
	my $GENE_NAME=$_;
	my @REGEXP=$cfg->val("$GENE_NAME",'wells');
	my $log_dir = $log_dir."/$GENE_NAME/";
	if(! -e $log_dir){
		mkdir($log_dir);
	}else{
		opendir(DH,$log_dir);
		my @delfiles = map{"$log_dir/$_"}grep{/\.log$/}readdir(DH);
		unlink(@delfiles);
		closedir(DH);
	}
	my @files;
	foreach my $re(@REGEXP){
		my $regexp = "$prefix$re$SUFFIX";
		find(sub{push @files,$File::Find::name if /$regexp/;},$data_dir);
	}

	foreach my $f1(@files){
		(my $alt_suffix=$SUFFIX)=~s/1/2/;
		my $rname = basename($f1,$SUFFIX);
		(my $f2=$f1) =~ s/$SUFFIX$/$alt_suffix/;
		my $cmd = "env perl $CM_SCRIPT -b $out_dir -ini $ini -gene $GENE_NAME -f1 $f1 -f2 $f2";
		## uncomment below to run locally rather than on q
		#print "$cmd > $log_dir$rname.log 2>&1";
		#`$cmd > $log_dir$rname.log 2>&1`;
		#next;
		my ($fh,$fname) = &get_tmp_file();
		print($fh $cmd);
		close($fh);
		my $cmd = "$QSUB_CMD -o $log_dir$rname.log -j y $fname";
		debug("$cmd");
		`$cmd`;
		unlink("$fname");
		last if $test;
	}
	last if $test;
}

close(GENES);

sub get_tmp_file{
	my $dir = shift;
	my $template = "CALL_METH_XXXXX";
	return tempfile($template, DIR => $dir) if $dir;
	return tempfile($template, TMPDIR => 1);
}

sub debug{
	my $msg=shift;
	my $ts = POSIX::strftime("[%m/%d/%Y %H:%M:%S]", localtime);
	print join("\t",$ts,$msg)."\n";
}
