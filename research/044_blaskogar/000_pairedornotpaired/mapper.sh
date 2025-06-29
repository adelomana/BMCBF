time /Users/adrian/software/kallisto_old/kallisto/kallisto quant -i /Users/adrian/software/kallisto_old/homo_sapiens/transcriptome.idx -o r1 --single -l 180 -s 20 --verbose --pseudobam -t 8 PE-041_Pre_551_S44_L003_R1_001.fastq.gz
time /Users/adrian/software/kallisto_old/kallisto/kallisto quant -i /Users/adrian/software/kallisto_old/homo_sapiens/transcriptome.idx -o r2 --single -l 180 -s 20 --verbose --pseudobam -t 8 PE-041_Pre_551_S44_L003_R2_001.fastq.gz

cd r1
time samtools view -h pseudoalignments.bam > output1.sam
cd ..

cd r2
time samtools view -h pseudoalignments.bam > output2.sam
cd ..