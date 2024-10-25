from datetime import datetime
import re
import os
import subprocess # important to call pre-existing modules

#########################################################################
# This version v0.5 works well and has been tested on the rMalMon genome
#########################################################################

#1- allows to call subworkflows (modules)
module polish_workflow:
  snakefile: "../modules/polish.rules.smk"
#evaluate_assemblies.rules.smk
module evaluate_assemblies_workflow:
  snakefile: "../modules/evaluate_assemblies.rules.smk"


#2- main snakemake rules
date = datetime.now().strftime('%Y%m%d.%H%M%S')
#get the bin path for the snakemake file
bin_path = sys.path[0]
#make scripts and utils relative to bin
scripts_dir = bin_path + "/../scripts/"
utils_dir = bin_path + "/../utils/"
keepfiles = False
work_dir = os.getcwd() + "/"
logs_dir = work_dir + "logs/"
if not os.path.exists(logs_dir):
 os.makedirs(logs_dir)


rule format_mito_ref:
  input:
    reference = config["data"]["reference_fasta"],
  output:
    mitoref = work_dir + "s01.1_p00_ref_MT/MT_REF.scaffolds.fa",
  params:
    scripts = scripts_dir,
  log:
    logs_dir+ str(date) + ".j%j.format_mito_ref.out",
    logs_dir+ str(date) + ".j%j.format_mito_ref.err",
  threads: 1
  shell:
    # Just rename it and format the fasta as multiline
    "mkdir -p s01.1_p00_ref_MT; "
    "cd s01.1_p00_ref_MT/; "
    "ln -s {input.reference} . ;"
    "{params.scripts}rename_fasta_seq.pl -f  {input.reference} -n MT_REF | {params.scripts}FastaToTbl | {params.scripts}TblToFasta > {output.mitoref} ;"
    "echo 'reference mitogenome has been formatted into a multiline fasta'; "

rule orient_reference_MT:
  input:
    genome = rules.format_mito_ref.output.mitoref,
  output:
    oriented = work_dir + "s01.1_p00_ref_MT/out/MT_REF.scaffolds.oriented.fa", 
    annotation = work_dir + "s01.1_p00_ref_MT/out/annotation_2/result.gff",
    formatted = work_dir + "s01.1_p00_ref_MT/out/MT_REF.formatted.fa"
  params:
    scripts = scripts_dir,
    refseq_dir = utils_dir + "refseq_dir/",
    genetic_code = config["annotation"]["genetic_code"],
    mitos_options = config["annotation"]["mitos_options"],
    refseq_db = config["annotation"]["refseq_database"],
  conda:
    "../envs/mitos.yaml"
  threads: 8
  log:
    logs_dir+ str(date) + ".j%j.orient_reference_MT.out",
    logs_dir+ str(date) + ".j%j.orient_reference_MT.err",
  shell:
    "echo workdir is ; "
    "/usr/bin/pwd ; "
    "cd s01.1_p00_ref_MT/ ; "
   
    #1. Annotation 1
    "echo 'testing mitos conda environment' ; "
    "runmitos.py --version ; "
    "mkdir -p annotation_1 ; "
    "runmitos.py -i {input.genome} -c {params.genetic_code} {params.mitos_options} -r {params.refseq_db} -R {params.refseq_dir} -o annotation_1 ; "
    
    #3. Orient based in Annotation set trnF start
    "perl {params.scripts}orient_mitogenome_v1.pl  -f {input.genome} -b annotation_1/result.bed ; "
    "mkdir -p out/annotation_2 ; "
    "mv MT_REF.scaffolds.oriented.fa out/ ; "
    #4. rename and format to ensure multiline fasta
     "cd out ; "
     "{params.scripts}rename_fasta_seq.pl -f  {output.oriented} -n MT_REF_oriented | {params.scripts}FastaToTbl | {params.scripts}TblToFasta > {output.formatted}; "
    # 5. Reannotate formatted mitogenome
      "runmitos.py -i {output.formatted} -c {params.genetic_code} {params.mitos_options} -r {params.refseq_db} -R {params.refseq_dir} -o annotation_2 ; "
    # 6. symbolic link to out file  
    "cd .. ; "
    "ln -s {output.formatted} . ; "
  
rule filter_ont_max_len:
  input:
    reads = config["data"]["ont"],
  output:
    reads = work_dir + "s02.1_p01.1_filtered_long_reads/ont_maxlen.fastq.gz",
  params:
    scripts = scripts_dir,
    max_len = config["parameters"]["max_read_length"]
  conda:
   "../envs/pigz-2.6.yaml"
  threads: 4
  log:
    logs_dir+ str(date) + ".j%j.filter_ont_max_len.out",
    logs_dir+ str(date) + ".j%j.filter_ont_max_len.err",
  shell:
    #example
    "echo \"mitochondrial reference name is {params.max_len}\";"
    "zcat -f {input.reads} | perl {params.scripts}filter_by_max_length_fastq.pl -max_len {params.max_len} | pigz -p {threads} -c > {output.reads} ;"

rule align_ont_to_mitoref:
  input:
    genome = rules.format_mito_ref.output.mitoref,
    reads = rules.filter_ont_max_len.output.reads
  output:
    ont_mito_fq = work_dir + "s03.1_p02.1_map_long_reads/ont_mito.fastq"
  params:
    align_opts = config["parameters"]["minimap2_opts"],
    tmp = work_dir + "tmp/minimap2.sam",
    scripts = scripts_dir,
    min_match = config["parameters"]["ontfilter_minmatch"],
    min_qual = config["parameters"]["ontfilter_minq"],
  conda:
   "../envs/fa5f7adcb4c59485a3a283d809e1402e.yaml"
  threads: 8
  log:
    logs_dir+ str(date) + ".j%j.align_ont_to_mito_ref.out",
    logs_dir+ str(date) + ".j%j.align_ont_to_mito_ref.err",
  shell:
    "mkdir -p tmp; "
    "mkdir -p s03.1_p02.1_map_long_reads/; "
    "cd  s03.1_p02.1_map_long_reads/; "
    #getting REFNAME for fishing reads
    "REFNAME=$( grep \'>\' {input.genome} | sed \'s/>//g\'); "
    "echo \"mitochondrial reference name is $REFNAME\";"

    # pipe to get a sam with the aligned reads to the mitogenome reference calling the md tag
    "minimap2 {params.align_opts} -t {threads} {input.genome} {input.reads} | samtools calmd -S - {input.genome} > {params.tmp};" # specify the output sam is now in the tmpdir

    # CAPTURE READS 1. Using defaults here. my $minmatch = 800 (hard-coded in perl-script) and $minqual = 12 (parameter in perlscript) 
    "perl --version ;"
    "ls -l {params.scripts}sam_id_filter_v02.pl ;"
    "{params.scripts}sam_id_filter_v02.pl -ref $REFNAME -m {params.min_match} -q {params.min_qual} <  {params.tmp} > ont.clean.fastq 2> {output.ont_mito_fq} ;"
    
rule flye_assembly_meta:
  input:
    reads = rules.align_ont_to_mitoref.output.ont_mito_fq
  output:
    assembly = work_dir+"s04.1_p03.1_flye_assembly/flye_meta.fa",
    info = work_dir+"s04.1_p03.1_flye_assembly/flye_meta_info.txt",
    png = work_dir+"s04.1_p03.1_flye_assembly/flye_meta.png",
    circular_contig = work_dir+"s04.1_p03.1_flye_assembly/circular_contig.fa",
    circular_bed = work_dir+"s04.1_p03.1_flye_assembly/circular_contig.bed",
  params:
    scripts = scripts_dir,
    out_dir = work_dir+"s04.1_p03.1_flye_assembly/",
    assembly_opts = config["parameters"]["flye_opts"],
    ref_len = config["parameters"]["reference_length"],
    len_var = config["parameters"]["mito_length_variation"],
    read_type = config["parameters"]["flye_reads"],
  conda:
    "../envs/36b7263cd61f32e4e9f7db777993cd52.yaml"   
  threads: 16
  log:
    logs_dir+ str(date) + ".j%j.flye_assembly.out",
    logs_dir+ str(date) + ".j%j.flye_assembly.err",
  shell:
    "mkdir -p {params.out_dir}out; "
    "cd {params.out_dir}; "
    #added --resume options to avoid redoing flye if computed successfuly (see assembly_opts under rule params)
    "echo 'Running command: flye {params.assembly_opts} -g {params.ref_len} -o {params.out_dir}out -t {threads} --nano-raw {input.reads}'; "
    "flye {params.assembly_opts} -g {params.ref_len} -o {params.out_dir}out -t {threads} {params.read_type} {input.reads}; "
    "ln -s {params.out_dir}out/assembly.fasta {output.assembly}; "
    "ln -s {params.out_dir}out/assembly_info.txt {output.info}; "
    "export XDG_RUNTIME_DIR=$PWD ;" # to avoid Bandage execution error
    "Bandage image {params.out_dir}out/assembly_graph.gfa {output.png}; "    
    "echo 'Selecting the circular contig with length between {params.ref_len}-{params.len_var} and {params.ref_len}+{params.len_var}'; "
    "circular=$(cat {output.info} | gawk \'{{ if ( ($4 == \"Y\") &&  ($2 >= ({params.ref_len}-{params.len_var})) && ($2 <= ({params.ref_len}+{params.len_var})) )   print $1}}\') ; " #double brackets to make it work
    "echo circular contig\(s\)\: $circular ;"
    "{params.scripts}select_contig.pl -f {output.assembly} -t $circular > {output.circular_contig} ; "
    # want to store the coordinates and name of coordinate contig to pass it to samtools in next rule align_illumina_flye_meta
    "{params.scripts}fastalength {output.circular_contig} | gawk \'{{ print $2\"\t\"0\"\t\"$1}}\'  > {output.circular_bed} ; "

# new rule to orient the polsihed mitogenome based on the mitos annotation. uses script orient_circular_contig_v1.pl 
rule orient_circular_contig:
  input:
    genome = rules.flye_assembly_meta.output.circular_contig,
  output:
    oriented = work_dir + "s05.1_p04.1_orient_circular_contig/out/circular_contig.oriented.fa",
    annotation = work_dir + "s05.1_p04.1_orient_circular_contig/annotation/result.gff",
    reannotation = work_dir + "s05.1_p04.1_orient_circular_contig/out/annotation_oriented/result.gff"
#    formatted = work_dir + "s03.1_p02.1_orient_circular_contig/out/formatted.fa",
  params:
    scripts = scripts_dir,
    refseq_dir = utils_dir + "refseq_dir/",
    genetic_code = config["annotation"]["genetic_code"],
    mitos_options = config["annotation"]["mitos_options"],
    refseq_db = config["annotation"]["refseq_database"],
   # tolid = config["parameters"]["tolid"],
  conda:
    "../envs/mitos.yaml"
  threads: 8
  log:
    logs_dir+ str(date) + ".j%j.orient_circular_ctg.out",
    logs_dir+ str(date) + ".j%j.orient_circular_ctg.err",
  shell:
    "cd  s05.1_p04.1_orient_circular_contig ; "
   
    #1. link circular contig polished with flye (default is 1 iteration) 
    "ln -s {input.genome} circular_contig.fa ; "

    #2. Annotation 1
    "echo 'testing mitos conda environment' ; "
    "runmitos.py --version ; "
    "mkdir -p annotation ; "
    "runmitos.py -i circular_contig.fa -c {params.genetic_code} {params.mitos_options} -r {params.refseq_db} -R {params.refseq_dir} -o annotation ; "
    
    #3. Orient based in Annotation set trnF start
    "perl {params.scripts}orient_mitogenome_v1.pl  -f circular_contig.fa -b annotation/result.bed ; "
    "mkdir -p out/ ; "
    "mv circular_contig.oriented.fa out/; "
    #4. rename and format to ensure multiline fasta
     "cd out/; "
    # "{params.scripts}rename_fasta_seq.pl -f  {output.oriented} -n {params.tolid}_MT | {params.scripts}FastaToTbl | {params.scripts}TblToFasta > {output.formatted} ;"
    #5. Reannotate formatted mitogenome
    "mkdir -p out/annotation_oriented ; "
    "runmitos.py -i {output.oriented} -c {params.genetic_code} {params.mitos_options} -r {params.refseq_db} -R {params.refseq_dir} -o annotation_oriented ; "
    # 6. format annotation gff and copy final fasta 
    #  "cd ../ ; "
    #  "cp {output.formatted}  {params.tolid}_MT.fa ; "
     # "{params.scripts}parse_mitos_gff.pl -g {output.annotation} -p MT{params.tolid} > MT{params.tolid}.gff3 ; "

rule flye_polishing:
  input:
    circular = rules.orient_circular_contig.output.oriented,
    reads = rules.align_ont_to_mitoref.output.ont_mito_fq
  output:
    assembly = work_dir+"s06.1_p05.1_flye_polishing/oriented.polished.fa",
  params:
    scripts = scripts_dir,
    out_dir = work_dir+"s06.1_p05.1_flye_polishing/",
    polish_opts = config["parameters"]["polishing_flye"],
    ref_len = config["parameters"]["reference_length"],
    read_type = config["parameters"]["flye_reads"],
  conda:
    "../envs/36b7263cd61f32e4e9f7db777993cd52.yaml"   
  threads: 16
  log:
    logs_dir+ str(date) + ".j%j.flye_polishing.out",
    logs_dir+ str(date) + ".j%j.flye_polishing.err",
  shell:
    "mkdir -p {params.out_dir}out; "
    "cd {params.out_dir}; "
    "ln -s {input.circular} . ; "
    #added --resume options to avoid redoing flye if computed successfuly (see assembly_opts under rule params)
    "echo 'Running command: flye {params.polish_opts} -g {params.ref_len} -o {params.out_dir}out -t {threads} {params.read_type} {input.reads} '; "
    "flye {params.polish_opts} -g {params.ref_len} -o {params.out_dir}out -t {threads} {params.read_type} {input.reads} --polish-target circular_contig.oriented.fa; "
    "ln -s {params.out_dir}out/polished_1.fasta {output.assembly}; "
    
# Finalize mitogenome rename using  ToLID and annotating it

rule final_mitogenome:
  input:
    genome = rules.flye_polishing.output.assembly,
  output:
    annotation = work_dir + "s07.1_p06.1_final_MT/out/annotation/result.gff",
    formatted = work_dir + "s07.1_p06.1_final_MT/out/formatted.fa",
    done=touch("final_mitogenome.done")
  params:
    scripts = scripts_dir,
    refseq_dir = utils_dir + "refseq_dir/",
    genetic_code = config["annotation"]["genetic_code"],
    mitos_options = config["annotation"]["mitos_options"],
    refseq_db = config["annotation"]["refseq_database"],
    tolid = config["parameters"]["tolid"],
  conda:
    "../envs/mitos.yaml"
  threads: 8
  log:
    logs_dir+ str(date) + ".j%j.final_mitogenome.out",
    logs_dir+ str(date) + ".j%j.final_mitogenome.err",
  shell:
    "mkdir -p  s07.1_p06.1_final_MT/out/annotation ; "
    "cd  s07.1_p06.1_final_MT ; "

    # 1. link polished contig 
    "ln -s {input.genome} . ; "

    # 2. rename and format to ensure multiline fasta
     "cd out/; "
     "{params.scripts}rename_fasta_seq.pl -f  {input.genome} -n {params.tolid}_MT | {params.scripts}FastaToTbl | {params.scripts}TblToFasta > {output.formatted} ;"
   
    # 3. Reannotate formatted mitogenome
      "runmitos.py -i {output.formatted} -c {params.genetic_code} {params.mitos_options} -r {params.refseq_db} -R {params.refseq_dir} -o annotation ; "
   
    # 4. format annotation gff and copy final fasta 
      "cd ../ ; "
      "cp {output.formatted}  {params.tolid}_MT.fa ; "
      "{params.scripts}parse_mitos_gff.pl -g {output.annotation} -p MT{params.tolid} > MT{params.tolid}.gff3 ; "

    

# PENDING TO ADD MERQURY SUBWORKFLOW USING A INPUT ASSEMBLY DICTIONARY...