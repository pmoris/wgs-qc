# Simple variant calling pipeline for SUMMIT samples

**Note that this script does not yet implement human read removal!**

## Pre-requisities:

- Install conda environment using `conda env create -f env.yml`
- Move or symlink fastq files to `./data/fastq`, e.g. `ln -s /data/antwerpen/grp/ap_itg_mu/projects/2024_SUMMIT_WGS/240404_procomcure/fastq/* . && ln -s /data/antwerpen/grp/ap_itg_mu/projects/2024_SUMMIT_WGS/240404_procomcure_reseq/fastq/* .`
- Move or symlink reference files to `./data/ref`, e.g. `ln -s /data/antwerpen/grp/ap_itg_mu/public_data/reference_genomes/* ./data/ref/`
- Adjust the file paths to the reference fasta files in `config/fastq-screen.conf` (requires absolute paths).

## Running the code:

Each script is wrapped by an associated Slurm script. These can also be bundled into a single Slurm wrapper that calls each component in turn.
