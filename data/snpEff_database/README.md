# snpEff database structure

SnpEff expects that all the files it needs to be stored in a specific location and with specific filenames, matching the names used in `./config/snpEff.config`. Briefly,

1. A copy of the reference genome should be placed in `./data/snpEff_database/<reference-name>/sequences.fa` (a symlink also works)
2. A gene annotation file should be placed in `./data/snpEff_database/<reference-name>/genes.gff`.
3. Add a CDS fasta file: `./data/snpEff_database/<reference-name>/cds.fa`.
4. Optionally a protein fasta file can also be added: `./data/snpEff_database/<reference-name>/protein.fa`.

The code below can be used to create these files for e.g. the Pv PAM genome. It should be run from inside the `./data/snpEff_database/PlasmoDB-68_PvivaxPAM/` folder.

```
wget wget https://plasmodb.org/common/downloads/release-68/PvivaxPAM/fasta/data/PlasmoDB-68_PvivaxPAM_Genome.fasta
mv PlasmoDB-68_PvivaxPAM_Genome.fasta sequences.fa

wget https://plasmodb.org/common/downloads/release-68/PvivaxPAM/gff/data/PlasmoDB-68_PvivaxPAM.gff
mv PlasmoDB-68_PvivaxPAM.gff genes.gff

wget https://plasmodb.org/common/downloads/release-68/PvivaxPAM/fasta/data/PlasmoDB-68_PvivaxPAM_AnnotatedCDSs.fasta
mv PlasmoDB-68_PvivaxPAM_AnnotatedCDSs.fasta cds.fa

wget https://plasmodb.org/common/downloads/release-68/PvivaxPAM/fasta/data/PlasmoDB-68_PvivaxPAM_AnnotatedProteins.fasta
mv PlasmoDB-68_PvivaxPAM_AnnotatedProteins.fasta protein.fa
```
