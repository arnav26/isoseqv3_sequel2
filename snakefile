localrules: get_chr13_bam, bed2fasta, tama_orf_seeker
WORKING_DIR = "/cluster/work/pausch/arnav/ecoli/isoseq_2023"

configfile: "config.yaml"

rule all:
    input: expand("{dir}/lima/{smrt}/{sample}/lima.Teloprime_primer_5p--Teloprime_primer_3p.bam", smrt = "SMRT1", sample = ["hom_sus", "hom_res"], dir = WORKING_DIR),
        expand("{dir}/refined/{smrt}/{sample}/refined_flnc.bam", smrt = "SMRT1", sample = ["hom_sus", "hom_res"], dir = WORKING_DIR),
        expand("{dir}/refined/{smrt}/{sample}/refined_flnc.fa", smrt = "SMRT1", sample = ["hom_sus", "hom_res"], dir = WORKING_DIR),
        expand("{dir}/refined_polya_trim/{smrt}/{sample}/polya_flnc.fa", smrt = "SMRT1", sample = ["hom_sus", "hom_res"], dir = WORKING_DIR),
        expand("{dir}/refined_polya_aligned/{smrt}/{sample}/{sample}_{smrt}_ensembl.bam", smrt = "SMRT1", sample = ["hom_sus", "hom_res"], dir = WORKING_DIR),
        expand("{dir}/refined_polya_aligned/{smrt}/{sample}/{sample}_{smrt}_chr13_ensembl.bam", smrt = "SMRT1", sample = ["hom_sus", "hom_res"], dir = WORKING_DIR),
        expand("{dir}/tama_collapse/{smrt}/{sample}/{smrt}_{sample}_collapsed.bed", smrt = "SMRT1", sample = ["hom_sus", "hom_res"], dir = WORKING_DIR),
        expand("{dir}/tama_collapse/{smrt}/{sample}/{smrt}_{sample}_collapsed.fasta", smrt = "SMRT1", sample = ["hom_sus", "hom_res"], dir = WORKING_DIR),
        expand("{dir}/tama_collapse/{smrt}/{sample}/{smrt}_{sample}_ORF.fasta", smrt = "SMRT1", sample = ["hom_sus", "hom_res"], dir = WORKING_DIR),
        expand("{dir}/tama_collapse/{smrt}/{sample}/{sample}_{smrt}_collapsed_ensembl.bam", smrt = "SMRT1", sample = ["hom_sus", "hom_res"], dir = WORKING_DIR),
        expand("{dir}/orf_blast_results/{smrt}/{sample}/{smrt}_{sample}_blast.txt", smrt = "SMRT1", sample = ["hom_sus", "hom_res"], dir = WORKING_DIR),
        expand("{dir}/orf_blast_results/{smrt}/{sample}/{smrt}_{sample}_blast_parsed.txt", smrt = "SMRT1", sample = ["hom_sus", "hom_res"], dir = WORKING_DIR),
        expand("{dir}/tama_cds_regions/{smrt}/{sample}/{smrt}_{sample}_CDS.bed", smrt = "SMRT1", sample = ["hom_sus", "hom_res"], dir = WORKING_DIR)

rule lima:
    input:
        bam_file = lambda wildcards: config["input_files"][wildcards.smrt][wildcards.sample]["reads"][0],
        primer = "{dir}/primer.fasta"
    output:
        processed_file = "{dir}/lima/{smrt}/{sample}/lima.Teloprime_primer_5p--Teloprime_primer_3p.bam"
    conda:  
        "pacbiotools"
    params:
        outputprefix = "{dir}/lima/{smrt}/{sample}",
        flags = "--isoseq --peek-guess --dump-clips -j 24"
    threads: 24
    shell: "lima {input.bam_file} {input.primer} {params.outputprefix}/lima.bam {params.flags}"

rule refine:
    input: 
        lima_output = rules.lima.output,
        primer = rules.lima.input.primer
    output:
        refined_file = "{dir}/refined/{smrt}/{sample}/refined_flnc.bam"
    conda:  
        "pacbiotools"
    params:
        directory = "{dir}/refined/{smrt}/{sample}"
    threads: 24
    shell:'''mkdir -p {params.directory}
        isoseq3 refine {input.lima_output} {input.primer} {output.refined_file} -j 24
        '''

rule bam2fasta:
    input:
        refine_output = rules.refine.output
    output:
        fasta = "{dir}/refined/{smrt}/{sample}/refined_flnc.fa"
    shell:'''
    module load bamtools
    bamtools convert -format fasta -in {input.refine_output}  > {output}
    '''

rule tama_polya_cleanup:
    input: 
        refined_fasta = rules.bam2fasta.output
    output:
        clean_fasta = "{dir}/refined_polya_trim/{smrt}/{sample}/polya_flnc.fa"
    conda: 
        "tama"
    params:
        prefix = "{dir}/refined_polya_trim/{smrt}/{sample}"
    resources:
        mem_mb = 20000
    shell: "tama_flnc_polya_cleanup.py -f {input.refined_fasta} -p {params.prefix}/polya_flnc"

rule minimap2:
    input: 
        fasta = rules.tama_polya_cleanup.output,
        ref = "/cluster/work/pausch/arnav/ecoli/isoseq_rsii/reference/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa"
    output:
        aligned_file = "{dir}/refined_polya_aligned/{smrt}/{sample}/{sample}_{smrt}_ensembl.bam"
    threads: 50
    shell:'''
        module load samtools
        minimap2 -t 50 -ax splice:hq -uf {input.ref} {input.fasta} | samtools view -@ 50 -b - | samtools sort -@ 50 -o {output.aligned_file} - && samtools index {output.aligned_file}
        '''

rule get_chr13_bam:
    input:
        rules.minimap2.output
    output:
        chr13_bam = "{dir}/refined_polya_aligned/{smrt}/{sample}/{sample}_{smrt}_chr13_ensembl.bam"
    shell:'''
        module load samtools
        samtools view -h -b {input} 13 | samtools sort -o {output.chr13_bam} - && samtools index {output.chr13_bam}
        '''

rule tama_collapse:
    input:
        bam = rules.get_chr13_bam.output,
        ref = rules.minimap2.input.ref
    output:
        results = "{dir}/tama_collapse/{smrt}/{sample}/{smrt}_{sample}_collapsed.bed"
    params:
        prefix = "{dir}/tama_collapse/{smrt}/{sample}/{smrt}_{sample}"
    resources:
        mem_mb=25000
    threads: 2
    conda:
        "tama"
    shell: "tama_collapse.py -s {input.bam} -f {input.ref} -p {params.prefix} -x capped -b BAM"

rule bed2fasta:
    input: 
        collapsed_bed = rules.tama_collapse.output,
        ref = rules.minimap2.input.ref
    output: 
        fasta = "{dir}/tama_collapse/{smrt}/{sample}/{smrt}_{sample}_collapsed.fasta"
    conda: "tama"
    shell:"bedtools getfasta -name -split -s -fi {input.ref} -bed {input.collapsed_bed} -fo {output.fasta}"


rule align_collapsed_fasta:
    input: 
        fasta = rules.bed2fasta.output,
        ref = "/cluster/work/pausch/arnav/ecoli/isoseq_rsii/reference/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa"
    output:
        aligned_file = "{dir}/tama_collapse/{smrt}/{sample}/{sample}_{smrt}_collapsed_ensembl.bam"
    threads: 50
    shell:'''
        module load samtools
        minimap2 -t 50 -ax splice:hq -uf {input.ref} {input.fasta} | samtools view -@ 50 -b - | samtools sort -@ 50 -o {output.aligned_file} - && samtools index {output.aligned_file}
        '''
    

rule tama_orf_seeker:
    input: rules.bed2fasta.output
    output: 
        ORF = "{dir}/tama_collapse/{smrt}/{sample}/{smrt}_{sample}_ORF.fasta"
    conda: "tama"
    shell: "tama_orf_seeker.py -f {input} -o {output}"

rule blast:
    input: 
        ORF = rules.tama_orf_seeker.output,
    output: "{dir}/orf_blast_results/{smrt}/{sample}/{smrt}_{sample}_blast.txt"
    threads: 99
    params:
        proteome_prefix = "/cluster/work/pausch/arnav/ecoli/isoseq_rsii/reference/pig_proteome/sus_scrofa_proteome"
    shell: '''
            mkdir -p {wildcards.dir}/orf_blast_results/{wildcards.smrt}/{wildcards.sample}
            module load blast-plus
            blastp -evalue 1e-10 -num_threads 50 -db {params.proteome_prefix} -query {input.ORF}  > {output}
            '''

rule parse_blast:
    input: rules.blast.output
    output: "{dir}/orf_blast_results/{smrt}/{sample}/{smrt}_{sample}_blast_parsed.txt"
    conda: "tama"
    threads: 25
    shell: "tama_orf_blastp_parser.py -b {input} -o {output}"

rule cds_bed_file:
    input:
        orf = rules.parse_blast.output,
        bed = rules.tama_collapse.output,
        fasta = rules.bed2fasta.output
    conda: "tama"
    threads: 25
    output: "{dir}/tama_cds_regions/{smrt}/{sample}/{smrt}_{sample}_CDS.bed"
    shell: "tama_cds_regions_bed_add.py -p {input.orf} -a {input.bed} -f {input.fasta} -o {output}"