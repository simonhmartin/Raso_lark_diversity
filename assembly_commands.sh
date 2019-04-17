### Trimmomatic

# Run Trimmomatic

sbatch trim_fastq.sh Skylark_R1_Apollo.fastq.gz Skylark_trimmed_Apollo.fastq.gz
sbatch trim_fastq.sh Skylark_R1_3kb.fastq.gz Skylark_trimmed_3kb.fastq.gz
sbatch trim_fastq.sh Skylark_2nd_R1_Apollo.fastq.gz Skylark_2nd_trimmed_Apollo.fastq.gz
sbatch trim_fastq.sh Skylark_2nd_R1_3kb.fastq.gz Skylark_2nd_trimmed_3kb.fastq.gz

### FastQC

sbatch run_fastqc.sh 

### run allpaths assembly

sbatch prepare_assembly.sh

sbatch run_assembly.sh
