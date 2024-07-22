# NAPB_wheat_tauschii_introgression


Ae tauschii Introgressions in Varieties Overley and TAM112
A. Retrieve SRR raw sequence from NCBI
The first step is to retrieve raw sequence data from SRA database-NCBI

#!/bin/bash
#SBATCH --mem=4G
#SBATCH --time=48:00:00
##SBATCH --mem-per-cpu 
##SBATCH --ntasks=2
##SBATCH --ntasks-per-cpu=2
#SBATCH -J=DwldSRA 
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err



sra="/homes/yemaneg/software/SRA/sratoolkit.2.10.8-centos_linux64/bin/fastq-dump-orig.2.10.8"

# write the SRR####### number after the sbatch *.sh command
$sra --split-files --origfmt --gzip $1 
B. Evaluate quality parameter of the reads
#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --job-name=fastpcPE
#SBATCH --out=slurm-%j.out
#SBATCH --error=slurm-%j.err


prefix_in=$1  #The code assumes raw sequence with _1.fastq.gz and _2.fastq.gz
prefix_ou=$2  #The output 


/homes/yemaneg/.conda/envs/fastp/bin/fastp --in1  "$prefix_in"_1.fastq.gz  --in2  \
"$prefix_in"_2.fastq.gz  --detect_adapter_for_pe --out1 "$prefix_ou"_1.fq.gz \
--out2 "$prefix_ou"_2.fq.gz  --unpaired1  "$prefix_ou"_unpaired_1.fq.gz  \
--unpaired2  "$prefix_ou"_unpaired_2.fq.gz  -h  "$prefix_ou"_fastp.html
C. Index the reference genome to be mapped
#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --j index
#SBATCH --out slurm-%j.out
#SBATCH --error slurm-%j.err


module load BWA/0.7.17-GCC-8.3.0


suffixref=$1   
input=$2

cd $(dirname $input) #points the path to where the input came from


bwa index -a bwtsw  -p $suffixref   $input

D. Map filtered reads to iwgsc_2.1 reference genome
#!/bin/bash

#SBATCH --time=160:00:00
#SBATCH --mem=128G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --job-name=align


module load BWA/0.7.17-GCC-8.3.0  
module load SAMtools/1.13-GCC-8.3.0

cd "/bulk/akf/yemaneg/overley/"

REF="/bulk/akf/yemaneg/Ref/iwgsc_refseqv2.1_assembly" #without .fa ... just the prefix of the reference sequence.


file=$(cat list.txt) 

for i in $file
do

echo "preparing bwa for $i" 

bwa mem -t 24 -a -M -R "@RG\tID:$i\tSM:$i" $REF "$i"_1.fq.gz  "$i"_2.fq.gz > "$i".sam

echo "bwa for $i is done; continuing to SAMtools"

samtools sort -@ 20 -m 4G "$i".sam > "$i"_sorted.bam

echo "samtools for $i is done; continuing to the next group!!!"

done

E. Quality filter for high quality reads.
The first step is to filter duplicated reads using Picard and the second is to filter flagged duplicate reads, low quality maps and reads with secondary alignment.

i- Marking duplicates
#!/bin/bash

#SBATCH --time=36:00:00
#SBATCH --ntasks=12
#SBATCH --mem=32GB
#SBATCH --job-name=dupPic


module load GATK/4.1.4.1-GCCcore-8.3.0-Java-11
module load picard/2.21.6-Java-11


input=$1
base=$(basename -s .bam $input)



cd $(dirname $input)   #output where input obtained
bamf=$(basename $input)


java -Xmx24g -jar $EBROOTPICARD/picard.jar MarkDuplicates\
      I=$input\
      O="$base"_marked.bam\
      M="$base"_dup_metrics.txt

echo "duplicate marking completed" 

java -Xmx24g -jar $EBROOTPICARD/picard.jar SortSam \
      I="$base"_marked.bam\
      O="$base"_marked_sorted2.bam \
      SORT_ORDER=coordinate

echo "diplicate marked bam file sorted" 
ii. Filtering for reads with high quality alignments
The following code includes only reads where the mates are properly aligned while excluding secondary and duplicate reads.

#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --job-name=Quality




input_bam=$1
output_bam=$2

echo " Input $input_bam ; output $output_bam; prefix $prefix"

module load SAMtools/1.13-GCC-8.3.0




samtools view -f 0x2 -F 0x808  -F 0x400  -q 60  -@ 12  $input_bam -o $output_bam 
