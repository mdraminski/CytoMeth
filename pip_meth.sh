###########
########### moje notataki (Damian)


## BS-Snper
#git clone https://github.com/hellbelly/BS-Snper

# sudo cpan
# install pod2text

bssnper=/home/thedam/Desktop/Doktorat/CytoMeth/tools/BS-Snper/BS-Snper.pl
ref=/home/thedam/Desktop/Doktorat/CytoMeth/referenceData/hg38.fa
in=/home/thedam/Desktop/Doktorat/CytoMeth/results/clipped/small_FAKE03.clipped.bam

perl BS-Snper.pl --fa /home/thedam/Desktop/Doktorat/CytoMeth/referenceData/hg38.fa --input /home/thedam/Desktop/Doktorat/CytoMeth/results/clipped/small_FAKE03.clipped.bam --output snp.candidate.out --methcg meth.cg --methchg meth.chg --methchh meth.chh --minhetfreq 0.1 --minhomfreq 0.85 --minquali 15 --mincover 10 --maxcover 1000 --minread2 2 --errorate 0.02 --mapvalue 20 >SNP.out 2>ERR.log


perl $bssnper --fa $ref --input $in --output snp.candidate.out --methcg meth.cg --methchg meth.chg --methchh meth.chh --minhetfreq 0.1 --minhomfreq 0.85 --minquali 15 --mincover 10 --maxcover 1000 --minread2 2 --errorate 0.02 --mapvalue 20 >SNP.out 2>ERR.log 


perl $bssnper 
--fa $ref 
--input $in 
--output snp.candidate.out 
--methcg meth.cg 
--methchg meth.chg 
--methchh meth.chh 
--minhetfreq 0.1 
--minhomfreq 0.85 
--minquali 15 
--mincover 10 
--maxcover 1000 
--minread2 2 
--errorate 0.02 
--mapvalue 20 
>SNP.out 2>ERR.log 





config$ref_data_intervals_file

  methyl_result_file <- paste0(file.path(config$results_path, config_tools[config_tools$proces=="methratio","temp_results_dirs"], sample_basename),".methylation_results.bed")




# die `pod2text $0` if (@ARGV==0 || $Help);

	if ($chr and $length){ #### Damian added this IF
		print "##contig=<ID=$chr,length=$length>\n";
	}





