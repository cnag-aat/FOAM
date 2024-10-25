# FOAM
FOAM (Flye ONT Assembly of Mitogenomes). Is a pipeline we have been using to assemble the complete circular mitogenome of different organisms. Although previous versions used a combination of Illumina and ONT reads, due to the increased accuracy of the latest we are only using ONT Q20+ reads. Keep in mind that as each genome and sample poses different challenges is in continuous development.

# Current release
Version 0.5 is available and the different upgrades are shown in the corresponding tag.

# Taxonomic range 
FOAM has been used to assemble the mitogenomes of several ERGA-BGE species. It has been successful across different taxa including insects,sponges, archnids and vertebrates. Previous versions were used in various genomes see the publications below.

# To make it run
1. Clone this repository.
2. Download the refseq81 databases from https://zenodo.org/records/2672835
3. Copy them to your local utils/refseq_dir.
4. Uncompress the refseq databases locally  (example:  tar jvxf refseq81m.tar.bz2). This will allow you to run mitos stand-alone without using the web interface.
5. Create config files for your genome. For that (check examples/config forlder config.yamls and also the foam.cnag_cluster.yaml that works for our slurm cluster in example/configs).
6. Run the pipeline with snakemake pointing to the snakefile in bin.

# Config and key parameters

  - **Reference Length:** set the right length for the reference genome used, example reference_length: 15025 for NC_012901 (i.e. Blattella germanica used as reference for Loboptera ibLobCana)

  - Chose the appropriate genetic code: 2 for Vertebrate 5 for Invertebrate (check code in Genetic_Code_Mitos.txt placed in the utils/ folder)
  - As v0.5 expects ONT Q20+ reads or hifi reads, it will no longer use Illumina data for polishing. Therefore mapping and assembly options below should work for most species:

	  parameters:
		max_read_length: 18000
  		ontfilter_minmatch: 5000
  		ontfilter_minq: 15
  		minimap2_opts: " -ax map-hifi"
  		flye_opts: " --meta --scaffold -i 1"
  		flye_reads: " --nano-hq"
  		polishing_flye: "-i 1"

Flye polishing iterations. 


# Contributors
Fernando Cruz, Jèssica Gómez-Garrido and Tyler Alioto. Genome Assembly and Annotation Team (CNAG).

Visit our website for more info about the team: https://denovo.cnag.cat/about

# References
A chromosome-level reference genome for the common octopus, _Octopus vulgaris_ (Cuvier, 1797)
Dalila Destanović, ProfileDarrin T. Schultz, Ruth Styfhals, Fernando Cruz, Jèssica Gómez-Garrido, Marta Gut, Ivo Gut, Graziano Fiorito, Oleg Simakov, Tyler S. Alioto, Giovanna Ponte, Eve Seuntjens.
doi: https://doi.org/10.1101/2023.05.16.540928 
