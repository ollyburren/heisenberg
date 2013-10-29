#!/usr/bin/env perl

use strict;
use File::Basename;
use POSIX;

## PRE PROCESS METH FILES ON THE BASIS OF
## XIN YANG'S PIPELINE
my $JAVA_BIN = '/usr/bin/java';
my $TRIMMO_BIN = "$ENV{HOME}/Trimmomatic-0.30/trimmomatic-0.30.jar";
my $FASTX_TRIM_BIN = "$ENV{HOME}/src/bin/fastx_trimmer";
## -f is number of bases to trim assume -Q only does this if above qscore of 33
## -z outputs file as zip
my $FASTX_PARAM = '-f 10 -Q 33 -z';
my $LINKER_FA = '<FILL IN>' ## location is adaptors dir/linkers.fa
my $TRIMMO_PARAM = " ILLUMINACLIP:$LINKER_FA:1:1:1";
$TRIMMO_PARAM .= ' HEADCROP:18 TRAILING:20 SLIDINGWINDOW:4:20';


## expect e.g. P1_G12.1.fq.gz
my $forward_file = shift;

## expect e.g. P1_G12.2.fq.gz
(my $reverse_file = $forward_file)=~s/1(\.fq\.gz)$/2$1/;

die "Reverse file not found\n" if $reverse_file eq $forward_file;

## note that this script is actually called by 
my $trim_index_dir = shift;
my $out_dir = shift;
my $unpaired_dir = shift;
my $log_dir = shift;

## step 1 after de multiplexing remove first 10 bases if very good quality from 
## forward read. 

my $stub = basename($forward_file,'.1.fq.gz');

my $tfile = $trim_index_dir.basename($forward_file);

my $trim_log = "$log_dir$stub.trim.out";

my $cmd = "zcat $forward_file | $FASTX_TRIM_BIN  $FASTX_PARAM 2> $trim_log > $tfile";

## for NGS6 index barcodes have already been removed - is this part of Xin's script ?

#print "$cmd\n";
#`$cmd`;


## next we remove 18 bp at front of each read (the linker) 
## any runs of bp at end of read below qscore 20 and finally..
## cut as soon as a window  4bp in size with an avg qscore less than 20
## is encountered.

## scoring for thresh is 0.6 for a match, mismatch = -(Q/20). i.e 12bp perfect match == 7.2, our adaptors are 32bp therefore expect score of 19 for perfect match, if we set to 15 then we allow for some mismatching.

my $trimmo_log = "$log_dir$stub.trimmomatic.log";

my $fwd_stub = "$out_dir$stub.1";
my $rev_stub = "$out_dir$stub.2";
my $unp_fwd = "$unpaired_dir$stub.1.unpaired.fq.gz";
my $unp_rev = "$unpaired_dir$stub.2.unpaired.fq.gz";


my $cmd = "$JAVA_BIN -jar $TRIMMO_BIN PE -phred33 -trimlog $trimmo_log $forward_file $reverse_file $fwd_stub.fq.gz  $unp_fwd $rev_stub.fq.gz  $unp_rev  $TRIMMO_PARAM";

print "$cmd\n";
`$cmd`;

## finally gzip log out 

`gzip -f $trimmo_log`;
