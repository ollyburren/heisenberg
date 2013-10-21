use XML::LibXML;
use XML::LibXSLT;
use strict;
use Data::Dumper;

## this software attempts to rescue meth status from PE
## reads that FLASH software was not able to combine.
## It takes the contents of nonCombined files from FLASH
## output and attempts to call meth status on each paired 
## read separately. These calls are then combined to yield
## an integrated set of results. 

my $QSCORE=20;

my $MAX_EXTENSION=20;

my $design_file = shift;
my $fq1=shift;
my $fq2=shift;
my $gene=shift || die "No gene name given";

my $XSLT=<<EOXSLT;
<xsl:stylesheet
xmlns:xsl="http://www.w3.org/1999/XSL/Transform" 
version="1.0" 
xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">

  <xsl:output method="html" />
  <xsl:param name="limit2" select="1"/>
  <xsl:strip-space elements="*" />
  
  <xsl:template match="/">
			<xsl:apply-templates select="//w:body" />
  </xsl:template>

  <xsl:template match="w:body"> 
      <xsl:apply-templates select="w:p"/>
		  <xsl:apply-templates select="w:tbl"/>
  </xsl:template>

	<xsl:template match="w:p">
    <xsl:choose>
    	<xsl:when test=".//w:t[starts-with(text(),'&gt;')]">
    		<xsl:element name="header">
    			 <xsl:value-of select="." />
    		</xsl:element>
    	</xsl:when>
    		<xsl:when test=".//\@w:val &#61; &quot;center&quot;">
    		</xsl:when>
    		<xsl:when test="not(normalize-space(w:r))">
    	</xsl:when>
    	<xsl:otherwise>
    		<xsl:element name="sequence">
    			<xsl:apply-templates select="w:r" />
    		</xsl:element>
    		<xsl:text>
    		</xsl:text>
    	</xsl:otherwise>
    </xsl:choose> 
	</xsl:template>

  <xsl:template match="w:r">
    <xsl:apply-templates select="w:t" />
  </xsl:template>

  <xsl:template match="w:t">
      <xsl:apply-templates select="../w:rPr" />
  </xsl:template>

  <xsl:template match="w:rPr">
    <xsl:apply-templates select="w:shd" />
    <xsl:if test="not(w:shd)">
    	<xsl:element name="normal">
      <xsl:value-of select="../w:t" />
      </xsl:element>
    </xsl:if>
  </xsl:template>

  <xsl:template match="w:shd">
  	<xsl:if test="\@w:fill &#61; &quot;00FF00&quot;">
      <adaptor> 
      <xsl:value-of select="../../w:t" />
      </adaptor> 
  	</xsl:if>
  	<xsl:if test="\@w:fill &#61; &quot;00FFFF&quot;">
      <index> 
      <xsl:value-of select="../../w:t" />
      </index> 
  	</xsl:if>
  	<xsl:if test="\@w:fill &#61; &quot;FF0000&quot;">
      <meth> 
      <xsl:value-of select="../../w:t" />
      </meth> 
  	</xsl:if>
  	<xsl:if test="\@w:fill &#61; &quot;FFFF00&quot;">
      <snp> 
      <xsl:value-of select="../../w:t" />
      </snp> 
  	</xsl:if>
  	<xsl:if test="\@w:fill &#61; &quot;FF00FF&quot;">
      <meth> 
      <xsl:value-of select="../../w:t" />
      </meth> 
  	</xsl:if>
	</xsl:template>
</xsl:stylesheet>
EOXSLT




my $string = `unzip -p $design_file word/document.xml`;
my $xslt = XML::LibXSLT->new();
my $dom = XML::LibXML->load_xml(string=>$string);
my $style_doc = XML::LibXML->load_xml(string=>$XSLT, no_cdata=>1);
my $stylesheet = $xslt->parse_stylesheet($style_doc);
my $results = $stylesheet->transform($dom);
my $txt = $stylesheet->output_as_bytes($results);
$txt=~s/(<header>)/<fasta>$1/g;
$txt=~s/(<\/sequence>)/$1<\/fasta>/g;
$dom = XML::LibXML->load_xml(string=>"<run>$txt</run>");
my (%tlu,%met_pos);
foreach my $fasta ($dom->findnodes('/run/fasta')) {
	 my($header) = uc $fasta->findnodes('./header')->to_literal;
	 $header=~s/^>//;
	 $header=~s/(\s+|\-)demeth$//i;
	 my $sequence = $fasta->findnodes('./sequence');
	 my $position=0;
	 foreach my $mu($sequence->get_nodelist){
	 	 my $sequence=$mu->to_literal();
	 	 my $bar = $mu->getChildrenByTagName('*');
	 	 
	 	 my (@met_pos,@snp_pos);
	 	 my $adaptor;
	 	 while(my $foo = $bar->shift){
	 	 	
	 	 	 if($foo->nodeName eq 'meth'){
	 	 	 	 push @met_pos,$position;
	 	 	 	 
	 	 	 }elsif($foo->nodeName eq 'snp'){
	 	 	 	 push @snp_pos,$position;
	 	 	 }
	 	 	 $position += length($foo->to_literal);
	 	 	 $adaptor=$position if $foo->nodeName eq 'adaptor';
	 	 }
	 	 $header=~s/(\s+|\-)meth//i;
	 	 #die Dumper(\@met_pos);
	 	 foreach my $mpos(@met_pos){
	 	 	 my $set=0;
	 	 	 for (my $x=1;$x<$MAX_EXTENSION;$x++){
	 	 	 	 my $tu = substr($sequence,($mpos-$x),2+$x);
	 	 	 	 my $mcount =()= $sequence =~ /$tu/gi;
	 	 	 	 #print "mcount $mcount\n";
	 	 	 	 unless($tlu{$header}{$tu} || $mcount >1){
	 	 	 	 	 ## 0:match_string start pos, 
	 	 	 	 	 ## 1:meth base (C|T)
	 	 	 	 	 ## 2:offset (i.e. where meth base lies in match_string
	 	 	 	 	 ## 3:position relative to first meth position
	 	 	 	 	 $tlu{$header}{$tu}=[$mpos-$x,substr($sequence,$mpos,1),$x,($mpos-$met_pos[0])];
	 	 	 	 	 $met_pos{$header}{$mpos-$met_pos[0]}++;
	 	 	 	 	 $set++;
	 	 	 	 	 last;
	 	 	 	 }
	 	 	 }
	 	 	 if(!$set){
	 	 	 	 print "For $header $mpos cannot find unique tuple\n";
	 	 	 	 exit;
	 	 	 }
	 	 }
	 }
}

## parse out the meth positions

my @rel_pos;
foreach my $rp(sort{$a <=> $b}keys %{$met_pos{$gene}}){
	if($met_pos{$gene}{$rp} == 2){
		push @rel_pos,$rp;
	}else{
		print "Something wrong with design file demeth and meth positions disagree at $rp\n";
	}
}


## read in test fastq files

my @skeys = sort{$tlu{$gene}{$a}->[0] <=> $tlu{$gene}{$b}->[0]}keys %{$tlu{$gene}};


open my $for,"zcat $fq1 |" || die "Cannot open $fq1\n";
open my $rev,"zcat $fq2 |" || die "Cannot open $fq2\n";


## call meth positions for forward ...
my $methf = call_meth($for,$tlu{$gene});
## .. and reverse
my $methr = call_meth($rev,$tlu{$gene},1);

## merge calls 
foreach my $seq(keys %$methf){
	(my $rseq = $seq)=~s!(\s)1\:!${1}2\:!;
	my @out;
	foreach my $rp(@rel_pos){
		my $fm = $methf->{$seq}->{$rp} || 'N';
		my $rm = $methr->{$rseq}->{$rp}|| 'N';
		if($fm eq $rm){
			push @out,$fm;
		}elsif($fm eq 'N'){
			push @out,$rm;
		}elsif($rm eq 'N'){
			push @out,$fm;
		}else{
			## disagreement between forward and reverse strands
			push @out,'M'; 
		}
	}
	print join('',@out)."\n" unless grep{/N|M/}@out;
}
		
sub call_meth{
	my ($fh,$tupleref,$rev)=@_;
	my %tuple = %$tupleref;
	my @found;
	my @len;
	my @meth;
	my %results;
	while(<$fh>){
		chomp;
		if($. % 4==1){
			@len=$_;
		}elsif($. % 4 ==2 ){
			push @len,$_;
		}elsif($. % 4 ==0){
			push @len,$_;
			my @q;
			if($rev){
				$len[1] = join("",map{tr/AGCT/TCGA/;$_;}reverse(split(//,$len[1])));
				@q = reverse map{ord($_)-33}split(//,$len[1]);
			}else{
				@q =  map{ord($_)-33}split(//,$len[2]);
			}
			my $seq_length = length($len[1]);
			for(my $x=0;$x<@skeys;$x++){
				## has the anchor point already been found ?
				## if not we use a regexp to search for it.
				my $rel_pos = @found?$tuple{$skeys[$x]}->[0] - $tuple{$skeys[0]}->[0] + $found[0]:undef;
				last if $rel_pos > $seq_length; ## stop here if we have run out of sequence;
				if(my $m = &find_match(\$skeys[$x],\$len[1],\$rel_pos)){
					push @found,$m;
					my $qual = $q[$tuple{$skeys[$x]}->[2] + $m];
					my $rel_pos =  $tuple{$skeys[$x]}->[3];
					$results{$len[0]}{$rel_pos}=$tuple{$skeys[$x]}->[1] if $qual > $QSCORE;
				}
				## one might consider the following
				## this will be slower but will allow small inserts ??
				## not tested. Find match will be called w/o position so
				## will use regexp to find a match. Slightly dangerous as
				## could cause a majority of matches to not be properly anchored.
				##}else{
				#	@found=();
			  #}
			}
			##reset found after each fastq record
			@found=();
		}
	}
	close($fh);
	return \%results;
}
				
		
sub find_match{
	my ($string,$sequence,$position)=@_;
	if(defined($$position)){
		return $$position if substr($$sequence,$$position,length($$string)) eq $$string;
	}elsif($$sequence=~/($$string)/){
			return $-[0];
	}else{
		return 0;
	}
}
