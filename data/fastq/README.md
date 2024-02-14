Small fastq files of 2000 reads each matching the shortened references genomes.

Created via:

    apptainer build sandy.sif docker://galantelab/sandy

    apptainer run sandy.sif genome -v -j 10 -o . -s 42 --coverage 5 data/ref/Pf3D7_01_v3.fa
    # seqkit split2 -1 <(seqkit head -n 10000 out_R1_001.fastq.gz) -2 <(seqkit head -n 10000 out_R2_001.fastq.gz) -p 5 -O data/fastq/pf -e .gz
    seqkit head -n 10000 out_R1_001.fastq.gz > pf_R1_001.fastq.gz
    seqkit head -n 10000 out_R2_001.fastq.gz > pf_R2_001.fastq.gz
    seqkit split2 -1 pf_R1_001.fastq.gz -2 pf_R2_001.fastq.gz -p 5 -O data/fastq/pf -e .gz

    for f in data/fastq/pf/*; do
        if [[ "$f" =~ ^(.*\/)(.*)(R[0-9]_[0-9]+)(.part_[0-9]+)(.fastq.gz)$ ]]; then
            echo mv ${f} ${BASH_REMATCH[1]}${BASH_REMATCH[2]}${BASH_REMATCH[4]#.}_${BASH_REMATCH[3]}${BASH_REMATCH[5]}
        fi
    done

    apptainer run sandy.sif genome -v -j 10 -o . -s 42 --coverage 1 data/ref/GRCh38.chr21.fa.gz
    seqkit head -n 10000 out_R1_001.fastq.gz > human_R1_001.fastq.gz
    seqkit head -n 10000 out_R2_001.fastq.gz > human_R2_001.fastq.gz
    seqkit split2 -1 human_R1_001.fastq.gz -2 human_R2_001.fastq.gz -p 5 -O data/fastq/human -e .gz

    for f in data/fastq/human/*; do
        if [[ "$f" =~ ^(.*\/)(.*)(R[0-9]_[0-9]+)(.part_[0-9]+)(.fastq.gz)$ ]]; then
            echo mv ${f} ${BASH_REMATCH[1]}${BASH_REMATCH[2]}${BASH_REMATCH[4]#.}_${BASH_REMATCH[3]}${BASH_REMATCH[5]}
        fi
    done

    rm out_R*.fastq.gz pf_R*.fastq.gz human_R*.fastq.gz
