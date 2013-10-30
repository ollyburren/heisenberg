use strict;
use File::Basename;
use File::Find;
use File::Temp qw/tempfile/;
use POSIX;
use Getopt::Long;
use Cwd 'abs_path';

my $EXTENSION = '.1.fq.gz';
my $TMP_DIR='/tmp';
my $QSUB_CMD="qsub -q all.q -v PERL5LIB=$ENV{PERL5LIB}";
my $TRIMMO_URL="http://www.usadellab.org/cms/?page=trimmomatic";

my $PSCRIPT = dirname(abs_path($0))."/pre_process_meth_run.pl";

my $USAGE=<<EOL;
$0: Pre-process demultiplexed reads and trim to gene-spec adaptor and insert

$0 --data_dir data_dir --log_dir log_dir

	MANDATORY PARAMETERS
		data_dir|d: path to location of demultiplexed data
		log_dir|l: path to location to store loging information
		trimmomatic|tj: path to trimmoatic jar file

	OPTIONAL
		help|h: print this message
		test|t: run one example to check configuration
		
	NOTE:
		Only works on files with a $EXTENSION extension - this can be changed within the script
		Output is written to data_dir/PRE_PROCESS
		Unpaired outut is written to data_dir/UNPAIRED
		
		Requires trimmomatic (tested version 0.30), this can be downloaded from $TRIMMO_URL
		
	Author: Olly Burren
	
Last updated 30/10/2013
EOL

my ($data_dir,$trimmo,$log_dir,$test,$help);

GetOptions(
	"data_dir|d=s" => \$data_dir,
	"log_dir|l=s" => \$log_dir,
	"trimmomatic|tm=s" => \$trimmo,
	"test|t" => \$test,
	"help|h" => \$help);

my $ABORT_FLAG= 0 || $help;

unless($ABORT_FLAG){
	if(!$data_dir){
		print "[ERROR] Require --data_dir/-d option: location of demultiplexed reads with extension $EXTENSION\n";
		$ABORT_FLAG++;
	}elsif(!$log_dir){
		print "[ERROR] Require --log_dir/-l option: path to dir to store q logging info\n";
		$ABORT_FLAG++;
	}elsif(!$trimmo){
		print "[ERROR] Require --trimmomatic/-tm option: path to trimmomatic jar file not found\n";
	}
}

unless($ABORT_FLAG){
	if(! -d $data_dir){
		print "[ERROR] Cannot find $data_dir data_dir\n";
		$ABORT_FLAG++;
	}elsif(! -e $PSCRIPT){
		print "[ERROR] Cannot find preprocess_meth_run script expecting to find at $PSCRIPT\n";
	}elsif(! -e $trimmo){
		print "[ERROR] Cannot find trimmomatic jar file. Is this installed ?\n";
	}
}

if($ABORT_FLAG){
	print $USAGE;
	exit(1);
}


my $odir = "${data_dir}PRE_PROCESS/";
my $unpdir = "${data_dir}UNPAIRED/";

for my $d(($log_dir,$odir,$unpdir)){
	`mkdir $d` unless -d $d;
}

opendir(DIR, $data_dir) or die $!;

while (my $file = readdir(DIR)) {
	# We only want files
	next unless (-f "$data_dir/$file");
	# Use a regular expression to find files ending in .txt
	next unless ($file =~ m/$EXTENSION$/);
	print "###PROCESSING $data_dir/$file\n";
	my $stub = basename($file,$EXTENSION);
	## $trimmo, out_dir, unpaired_dir, log_dir
	my $cmd = "env perl $PSCRIPT $trimmo $data_dir/$file $odir $unpdir $log_dir";
	my ($fh,$fname) = &get_tmp_file();
	debug($cmd);
	print($fh $cmd);
	close($fh);
	my $cmd = "$QSUB_CMD -o ${log_dir}$stub.log -j y $fname";
	debug("$cmd");
	`$cmd`;
	unlink("$fname");
	last if $test;
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
