heisenberg
==========

pipeline to call methylation status for immunophenotyping 

under development, use at your own risk.


basic steps are

1. Open your sample file (xls) and save well to gene designations in a tab delim file (no quotes).
2. Check your design file (docx) and make a note of FASTA headers.
3. open up create_meth_ini.pl and make sure that you fill in __DATA__ section (I need to change this so it works from a file). This maps identifiers used in sample file with those used in the design file.
4. run create_meth_ini.pl <design.file.docx> <sample.file.csv> > my.ini
5. check ini file and fill in any missing bits see call_meth.ini for more details.
6. run run_meth_analysis.pl (see --help for more details). This will run the analysis on the q.
7. run consolidate_summary.pl <analysis base dir>. This will create a tab delim file for each pcr product run of the total number of meth and demethylated sites across all reads within a well (i.e. for an individual).
8. run the R scripts in R directory to graph qc statistics to troubleshoot possible problems.

todo
====

Check for insert size bug when snp in design file.
Make sample file to FASTA header mapping file rather than using DATA tag.
Rscript-fy Rscripts with proper argument handling.
Make xslt more configurable i.e stop hardcoded colours.
Look at trimming based on pcr adaptor sequences prior to PE read assembly to improve yield.
Heuristics for guessing poly-T issues and implementing meth offsets dynamically.
SNP calling ?
Add a serial mode for people w/o access to a SGE q - perhaps provide support for other grid software (see sandman libs).
