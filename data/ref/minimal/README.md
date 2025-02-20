Minimal versions of reference genomes for human (chromosome 21) and Plasmodium falciparum (chromosome 1), as well as a concatenated version.

Created via:

    # samtools faidx PlasmoDB-66_Pfalciparum3D7_Genome.fasta
    samtools faidx PlasmoDB-66_Pfalciparum3D7_Genome.fasta "Pf3D7_01_v3" | bgzip > Pf3D7_01_v3.fa.gz

    # samtools faidx GRCh38.primary_assembly.genome.fa.bgz
    samtools faidx GRCh38.primary_assembly.genome.fa.bgz chr21 | bgzip > GRCh38.chr21.fa.gz

    zcat GRCh38.chr21.fa.gz Pf3D7_01_v3.fa.gz | bgzip > concat.fa.gz

    # variant for uncompressed file, `-` tells cat to read from stdin (i.e. the pipe)
    # zcat GRCh38.chr21.fa.gz | cat Pf3D7_01_v3.fa - | bgzip > concat.fa.gz

Original reference genomes were retrieved via:

    # https://plasmodb.org/plasmo/app/downloads/Current_Release/Pfalciparum3D7/
    wget https://plasmodb.org/common/downloads/release-66/Pfalciparum3D7/fasta/data/PlasmoDB-66_Pfalciparum3D7_Genome.fasta

    # https://www.gencodegenes.org/human/
    # Release 45 (GRCh38.p14)
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh38.primary_assembly.genome.fa.gz
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/MD5SUMS
    md5sum -c MD5SUMS
    # Some tools require bgzip compressed files rather than regular gzip (e.g., samtools faidx), so these were created via:
    gunzip --keep --stdout GRCh38.primary_assembly.genome.fa.gz | bgzip -@ 8 > GRCh38.primary_assembly.genome.fa.bgz
