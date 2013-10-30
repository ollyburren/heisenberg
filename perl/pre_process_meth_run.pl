#!/usr/bin/env perl

use strict;
use File::Basename;
use POSIX;
use Cwd 'abs_path';

## setup adaptor file
my $running_dir = dirname(abs_path($0));
(my $LINKER_FA = $running_dir)=~s!/[^/]+$!/adaptors/linkers.fa!;

## PRE PROCESS METH FILES ON THE BASIS OF
## XIN YANG'S PIPELINE
(my $JAVA_BIN = `which java`)=~s/\n//;
my $TRIMMO_PARAM = " ILLUMINACLIP:$LINKER_FA:1:1:1";
$TRIMMO_PARAM .= ' HEADCROP:18 TRAILING:20 SLIDINGWINDOW:4:20';

my $trimmo_bin = shift;
## expect e.g. P1_G12.1.fq.gz
my $forward_file = shift;

## expect e.g. P1_G12.2.fq.gz
(my $reverse_file = $forward_file)=~s/1(\.fq\.gz)$/2$1/;

die "Reverse file not found\n" if $reverse_file eq $forward_file;

## note that this script is actually called by run_pre_process_on_queue.pl therefore no arg handling I'm afraid

my $out_dir = shift;
my $unpaired_dir = shift;
my $log_dir = shift;


my $stub = basename($forward_file,'.1.fq.gz');


## We remove 18 bp at front of each read (the linker) 
## any runs of bp at end of read below qscore 20 and finally..
## cut as soon as a window  4bp in size with an avg qscore less than 20
## is encountered.

## scoring for thresh is 0.6 for a match, mismatch = -(Q/20). i.e 12bp perfect match == 7.2, our adaptors are 32bp therefore expect score of 19 for perfect match, if we set to 15 then we allow for some mismatching.

my $trimmo_log = "$log_dir$stub.trimmomatic.log";

my $fwd_stub = "$out_dir$stub.1";
my $rev_stub = "$out_dir$stub.2";
my $unp_fwd = "$unpaired_dir$stub.1.unpaired.fq.gz";
my $unp_rev = "$unpaired_dir$stub.2.unpaired.fq.gz";


my $cmd = "$JAVA_BIN -jar $trimmo_bin PE -phred33 -trimlog $trimmo_log $forward_file $reverse_file $fwd_stub.fq.gz  $unp_fwd $rev_stub.fq.gz  $unp_rev  $TRIMMO_PARAM";

print "$cmd\n";
`$cmd`;

## finally gzip log out 

`gzip -f $trimmo_log`;
