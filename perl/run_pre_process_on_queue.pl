use strict;
use File::Basename;
use File::Find;
use File::Temp qw/tempfile/;
use POSIX;


my $TEST=0;
my $PRE_PROC_BIN = '<LOC OF REPO>/perlpre_process_meth_run.pl';
my $ROOT_DIR = "FILL_ME_IN";
my $DATA_DIR = "$ROOT_DIR/DATA/NGS6/";
my $OUT_DIR = "$ROOT_DIR/DATA/NGS6/";
my $LOG_DIR = "$ROOT_DIR/log/NGS6/PRE_PROCESS/";
my $QSUB_CMD="qsub -q all.q -v PERL5LIB=$ENV{PERL5LIB}";
my $TMP_DIR='/tmp';

my $odir = "${OUT_DIR}PRE_PROCESS/";
my $tdir = "${OUT_DIR}INDEX_TRIM/";
my $unpdir = "${OUT_DIR}UNPAIRED/";

for my $d(($LOG_DIR,$OUT_DIR,$odir,$tdir,$unpdir)){
	`mkdir $d` unless -d $d;
}

opendir(DIR, $DATA_DIR) or die $!;

while (my $file = readdir(DIR)) {
	# We only want files
	next unless (-f "$DATA_DIR/$file");
	# Use a regular expression to find files ending in .txt
	next unless ($file =~ m/\.fq\.gz$/);
	print "###PROCESSING $DATA_DIR/$file\n";
	my $stub = basename($file,'.1.fq.gz');
	my $cmd = "env perl $PRE_PROC_BIN $DATA_DIR/$file $tdir $odir $unpdir $LOG_DIR";
	my ($fh,$fname) = &get_tmp_file();
	print($fh $cmd);
	close($fh);
	my $cmd = "$QSUB_CMD -o ${LOG_DIR}$stub.log -j y $fname";
	debug("$cmd");
	`$cmd`;
	unlink("$fname");
	last if $TEST;
}
closedir(DIR);


sub get_tmp_file{
	my $dir = shift;
	my $template = "PREP_METH_XXXXX";
	return tempfile($template, DIR => $dir) if $dir;
	return tempfile($template, TMPDIR => 1);
}

sub debug{
	my $msg=shift;
	my $ts = POSIX::strftime("[%m/%d/%Y %H:%M:%S]", localtime);
	print join("\t",$ts,$msg)."\n";
}
