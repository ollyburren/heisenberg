use XML::LibXML;
use XML::LibXSLT;
use strict;
use Data::Dumper;

my $DELIM = ","; ## CHANGE THIS IF YOU WANT COMMA's


open(GENES,">genes.txt") || die "Cannot open genes.txt\n";

##NOTE THIS XSLT ONLY WORKS WITH LibreOffice Word docs atm

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


my $HEADER=<<EOL;
[GENERAL]
qscore_cutoff=20
cutadapt_bin=<FILL ME IN>
flash_bin=<FILL ME IN>
sickle_bin=<FILL ME IN>


[CUTADAPT ilmn_left]
a=GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG
q=20
O=6
[CUTADAPT ilmn_right]
a=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
q=20
O=6

EOL

my $file = shift;
my $sample_file = shift;
my $string = `unzip -p $file word/document.xml`;
my $xslt = XML::LibXSLT->new();
my $dom = XML::LibXML->load_xml(string=>$string);
my $style_doc = XML::LibXML->load_xml(string=>$XSLT, no_cdata=>1);
my $stylesheet = $xslt->parse_stylesheet($style_doc);
my $results = $stylesheet->transform($dom);;
my $txt = $stylesheet->output_as_bytes($results);
$txt=~s/(<header>)/<fasta>$1/g;
$txt=~s/(<\/sequence>)/$1<\/fasta>/g;
my %CONF;
$dom = XML::LibXML->load_xml(string=>"<run>$txt</run>");
foreach my $fasta ($dom->findnodes('/run/fasta')) {
	 my($header) = $fasta->findnodes('./header')->to_literal;
	 $header=~s/^>//;
	 my $sequence = $fasta->findnodes('./sequence');
	 my $position=0;
	 my ($met_first_pos);
	 foreach my $mu($sequence->get_nodelist){
	 	 my $sequence=$mu->to_literal();
	 	 my $bar = $mu->getChildrenByTagName('*');
	 	
	 	 my (@met_pos,@snp_pos);
	 	 my $adaptor;
	 	 while(my $foo = $bar->shift){
	 	 	 if($foo->nodeName eq 'meth'){
	 	 	 	 $met_first_pos=$position unless $met_first_pos;
	 	 	 	 push @met_pos,$position-$met_first_pos;
	 	 	 	 
	 	 	 }elsif($foo->nodeName eq 'snp'){
	 	 	 	 push @snp_pos,$position-$met_first_pos;
	 	 	 }
	 	 	 $position += length($foo->to_literal);
	 	 	 $adaptor=$position if $foo->nodeName eq 'adaptor';
	 	 }
	 	 next if $header=~m/(\s+|\-)demeth$/i;
	 	 $header=~s/(\s+|\-)meth//i;
	 	 $CONF{uc $header}{met_pos}=\@met_pos;
	 	 $CONF{uc $header}{left}=substr($sequence,$met_first_pos-20,20);
	 	 $CONF{uc $header}{right}=substr($sequence,$met_pos[-1]+$met_first_pos+3,20);
	 	 $CONF{uc $header}{adaptor}=$adaptor;
	 	 $CONF{uc $header}{m}=roundoff($adaptor-$met_first_pos-1,10); ## this reproduces what Chris P has used.
	 	 $CONF{uc $header}{insert_size}=$met_pos[-1]+3;
	 	 $CONF{uc $header}{snp}=\@snp_pos;
	 }
	 
}
 
##next we parse the sample file to work out which genes are assayed in which well

open(IN,$sample_file);

my %lu;
while(<DATA>){
	chomp;
	my ($k,$v)=split("\t",$_);
	$lu{$k}=$v;
}

my %RESULTS;
my @header;
my @gcols;
while(<IN>){
	chomp;
	unless(@header){
		@header=map{uc $_}split($DELIM,$_);
		@gcols=grep{$header[$_]=~/^GENE\s+[0-9]$/i}0..$#header;
	}else{
		my @val = split($DELIM,$_);
		foreach my $i(@gcols){
			next unless $val[$i];
			## if values are different use lookup defined in __DATA__
      if(my $nval=$lu{uc $val[$i]}){
      	push @{$CONF{$nval}{wells}},"$val[0]_$val[1]";
      }else{
      	## otherwise just use raw values and assume that they
      	## are the same between design files and sample file
      	push @{$CONF{uc $val[$i]}{wells}},"$val[0]_$val[1]";
      }
		}
	}
}


print $HEADER;



foreach my $g(keys %CONF){
	print GENES $g."\n";
	my $MP = join(",",@{$CONF{$g}{met_pos}});
	my $WELLS = join("\n",@{$CONF{$g}{wells}});
	my $SNPS = join(",",@{$CONF{$g}{snp}});
	my $ini_sect=<<EOINI;
[$g]
met_pos=$MP
snp_pos=$SNPS
insert_size=$CONF{$g}{insert_size}
wells=<<EOL
$WELLS
EOL
trim=sickle

[CUTADAPT ${g}_left]
g=$CONF{$g}{left}
m=$CONF{$g}{m}
[CUTADAPT ${g}_right]
a=$CONF{$g}{right}
EOINI
	print "$ini_sect\n";
}

sub roundoff{
  my $num = shift;
  my $roundto = shift || 1;
  return int($num/$roundto+0.5)*$roundto;
}

## THIS DATA SECTION IS USED TO MAP FASTA HEADERS IN DESIGN FILE
## TO IDENTIFIERS IN SAMPLE FILE. TAB DELIMITED

## e.g if gene identifier is FOO in sample file BAR in design FASTA HEADER then
## FOO	BAR

close(GENES);

__DATA__
