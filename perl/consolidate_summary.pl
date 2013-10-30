#!/usr/bin/perl

use strict;
use File::Find;
use Data::Dumper;

my $dir = shift;
my $odir = shift;

if(! -d $dir || ! -d $odir){
	print "USAGE: $0 analysis_dir output_dir\n";
}

my %results;
find(sub{
	return unless /summary$/;
	print "Processing ".$File::Find::name."\n";
	open(IN,$File::Find::name);
	my $heading;
	while(<IN>){
		chomp;
		my @val=split("\t",$_);
		if($val[1] eq 'gene'){
			$heading = $_;
		}else{
			$results{$val[1]}{$val[0]}=$_;
			$results{$val[1]}{heading}=$heading unless $results{$val[1]}{heading};
		}
		
	}
	close(IN);
	
	},$dir);

foreach my $gene(keys %results){
	open(OUT,">$odir/$gene.tab") || die "Cannot open $gene.tab for writing\n";
	print "#############$gene to $odir/$gene.tab\n";
	print  OUT $results{$gene}{heading}."\n";
	foreach my $well(sort{$a->[1] cmp $b->[1] || $a->[2] cmp $b->[2] || $a->[3] <=> $b->[3]} map{/(.*)_([A-Z])([0-9]+)$/;[$_,$1,$2,$3]}keys %{$results{$gene}}){
		next if $well->[0] eq 'heading';
		print OUT $results{$gene}{$well->[0]}."\n";
	}
	close(OUT);
}
