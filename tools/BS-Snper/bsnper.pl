use strict;
use Getopt::Long;
use POSIX;
use FindBin '$Bin';
use File::Path;  ## function " mkpath" and "rmtree" deal with directory
use File::Basename qw(basename dirname);


# my ($Help,$fasta,$bam,$mapvalue,$minhetfreq,$minhomfreq,$minquali,$minread2,$mincover,$maxcover,$errorate,$pvalue,$interval,$output,$methcg,$methchg,$methchh);
# my ($regFile,$varOnly);
# GetOptions(
#     "fa:s"=>\$fasta,
# 	"input:s"=>\$bam,
# 	"regions-file:s"=>\$regFile,
# 	"variants-only"=>\$varOnly,
# 	"output:s"=>\$output,
# 	"methcg:s"=>\$methcg,
# 	"methchg:s"=>\$methchg,
# 	"methchh:s"=>\$methchh,
# 	"minhetfreq:f"=>\$minhetfreq,
# 	"minhomfreq:f"=>\$minhomfreq,
# 	"minquali:i"=>\$minquali,
# 	"mincover:i"=>\$mincover,
# 	"maxcover:i"=>\$maxcover,
# 	"minread2:i"=>\$minread2,
# 	"errorate:f"=>\$errorate,
# 	"mapvalue:i"=>\$mapvalue,
# 	"help"=>\$Help
# );

# print "\ndupa\n\n\n";

# print @ARGV + "\n\n";
# if (@ARGV==0) {
# 	print "kkkk\n"
# }

# die `pod2text $0` if (@ARGV==0 || $Help);

# print "\ndupa\n\n\n";

print "\n$Bin\n";

my $bam='/home/thedam/Desktop/Doktorat/CytoMeth/results/clipped/small_FAKE03.clipped.bam';

open INTV,"$Bin/samtools-0.1.19/samtools view -H $bam|" or die $!;
#@SQ     SN:chr3 LN:198295559
while(<INTV>){
	chomp;
	my @a=split;
	my ($chr,$length);
	if(/SN\:(\S+)/){
		$chr=$1;
	}
	if(/LN\:(\d+)/){
		$length=$1;
	}

	if ($chr and $length){
	print "\n\m@a\n";
	print "##contig=<ID=$chr,length=$length>\n";
	}
}
close INTV;