#!/bin/bash

#SBATCH --job-name=MITF_M_Untreated_IgG_1
#SBATCH --partition=mimir
#SBATCH --nodes=1            
#SBATCH --ntasks-per-node=16    
#SBATCH --hint=nomultithread       
#SBATCH --output=messages/messages.MITF_M_Untreated_IgG_1.out.txt   
#SBATCH --error=messages/messages.MITF_M_Untreated_IgG_1.err.txt 

export LC_ALL=C.UTF-8
export LANG=C.UTF-8

date
pwd
conda info --envs

#
# 1. Create a temporary directory with a unique identifier associated with your jobid
# 
scratchlocation=/scratch/users
if [ ! -d $scratchlocation/$USER ]; then
mkdir -p $scratchlocation/$USER
fi
tdir=$(mktemp -d $scratchlocation/$USER/$SLURM_JOB_ID-XXXX)
echo "the scratch dir is:"
echo $tdir

#
# 2. copy files into scratch
#
echo ""
echo "about to copy files into scratch"
date
cp -rf /hpcdata/Mimir/adrian/door/MITF_M_Untreated_IgG_1 $tdir/.
date
cd $tdir
cd MITF_M_Untreated_IgG_1

#
# 3. call fastp
#
echo ""
echo "about to call fastp"

date
time fastp -i novaseqxplus4_20240513_LH00203_0044_A22F577LT3.s_1_1705_1505.001.R1.fastq.gz -I novaseqxplus4_20240513_LH00203_0044_A22F577LT3.s_1_1705_1505.001.R2.fastq.gz -o clean.R1.fq.gz -O clean.R2.fq.gz --thread 32 --trim_poly_g --detect_adapter_for_pe --allow_gap_overlap_trimming --length_required 25 -h report.html -j report.json
date

#
# 4. call bowtie2
#
echo ""
echo "about to call bowtie2"

date
time bowtie2 -x /users/home/adrian/software/bowtie2/GRCh38_noalt_as/GRCh38_noalt_as -1 clean.R1.fq.gz -2 clean.R2.fq.gz --end-to-end --very-sensitive --no-mixed --no-discordant -p 32 -I 10 -X 700 | samtools view -bS -@ 32 - | samtools sort -@ 32 -o human.bam
date

#
# 5. add read groups
#
echo ""
echo "about to add read groups"
date
time /users/home/adrian/software/java/jdk-21.0.8/bin/java -jar /users/home/adrian/software/picard/picard.jar AddOrReplaceReadGroups I=human.bam O=human.RG.bam RGID=MITF_M_Untreated_IgG_1 RGLB=lib_MITF_M_Untreated_IgG_1 RGPL=ILLUMINA RGPU=unit1 RGSM=MITF_M_Untreated_IgG_1 VERBOSITY=WARNING
date

#
# 6. mark or remove duplicates
#
echo ""
echo "about to mark or remove duplicates"
date
time /users/home/adrian/software/java/jdk-21.0.8/bin/java -jar /users/home/adrian/software/picard/picard.jar MarkDuplicates I=human.RG.bam O=human.RG.marked_dup.bam REMOVE_DUPLICATES=true M=MITF_M_Untreated_IgG_1.marked_dup.info.txt VERBOSITY=WARNING
date

#
# 7. MAPQ filtering
#
echo ""
echo "about to run MAPQ filtering"
date
time samtools view --threads 32 -b -q 30 human.RG.marked_dup.bam > human.RG.marked_dup.q30.bam
date

#
# 8. compute fragment length
#
echo ""
echo "about to compute fragment length"
date
time samtools view -F 0x04 human.RG.marked_dup.q30.bam | awk 'function abs(x){return (x<0?-x:x)}{if($9!=0) sizes[abs($9)]++}END{for(s in sizes) print s, sizes[s]/2}' OFS="\t" | sort -n > fragment_length_info.txt
date

# 
# 9. convert to BEDPE fragments
#
echo ""
echo "about to convert to BEDPE fragments"
date
time samtools sort -n -@ 8 -o human.sorted.bam human.RG.marked_dup.q30.bam
time bedtools bamtobed -bedpe -i human.sorted.bam > human.bedpe
time awk '$1==$4 && $6-$2 < 1000 {print $0}' human.bedpe > human.clean.bed
time cut -f 1,2,6 human.clean.bed | sort -k1,1 -k2,2n -k3,3n > human.fragments.bed
time awk -v w=500 '{print $1, int(($2 + $3)/(2*w))*w + w/2}' human.fragments.bed | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' | sort -k1,1V -k2,2n > human.fragments.500.bed
date

#
# 10. convert to bedGraph
#
echo ""
echo "about to convert to bedGraph"
date
time bedtools genomecov -i human.fragments.bed -g /users/home/adrian/software/bowtie2/GRCh38_noalt_as/GRCh38_noalt_as.genome.sizes -bg > human.bedgraph
date

#
# 11. copy results back to my folders
#
echo ""
echo "about to copy results out of scratch to my dirs"
ls -ltrsah
mkdir /hpcdata/Mimir/adrian/research/keilir/results/MITF_M_Untreated_IgG_1

date
cp fragment_length_info.txt /hpcdata/Mimir/adrian/research/keilir/results/MITF_M_Untreated_IgG_1/.
cp human.fragments.500.bed /hpcdata/Mimir/adrian/research/keilir/results/MITF_M_Untreated_IgG_1/.
cp human.bedgraph /hpcdata/Mimir/adrian/research/keilir/results/MITF_M_Untreated_IgG_1/.
date

#
# 12. clean scratch
#
echo ""
echo "about to remove scratch dir"
rm -rf $tdir
echo "all done."
date
    