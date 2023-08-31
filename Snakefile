# Snakefile
# Workflow: Snakemake-GATKv4.4-Snpeff-Workflow
# Version: 1.0
# Author: Dr. Ravi Prabhakar More (Email: ravipmore7@gmail.com)
# Date: 31-08-2023
# USAGE: snakemake -c config.yaml --cores 16 --dry-run
configfile: "config.yaml"

samples = config["samples"]

#hg38 and gatk bundle resources
genome = config["genome"]
idx = genome

known_sites_dbsnp = config["known_sites_dbsnp"]
idx_dbsnp = known_sites_dbsnp

known_sites_Mills_1000G = config["known_sites_Mills_1000G"]
idx_Mills_1000G = known_sites_Mills_1000G

known_sites_indels = config["known_sites_indels"]
idx_indels = known_sites_indels

snpeff_config = config["snpeff_config"]
idx_snpeff_config = snpeff_config

snpeff_jar = config["snpeff_jar"]
idx_snpeff_jar = snpeff_jar

java_path = config["java_path"]
idx_java_path = java_path

temp_dir = config["temp_dir"]
idx_temp_dir = temp_dir

# variant calling rules 
rule all:
    input: expand("output/{sample}.merged_filtered_snpeff.vcf", sample=samples)
    params: genome=genome

rule bwa_mem:
    input: reads=["input/{sample}.1.fastq.gz", "input/{sample}.2.fastq.gz"],
           idx=idx
    output: "output/{sample}.sam"
    params:
        readgroup=lambda wildcards: wildcards.sample
    shell: 
        """
        bwa mem -R "@RG\\tID:{params.readgroup}\\tLB:library1\\tPL:illumina\\tPU:unit1\\tSM:{params.readgroup}" -t 2 -K 10000000 {idx} {input.reads} > {output}
        """

rule sam_to_bam:
    input:
        sam="output/{sample}.sam"
    output:
        bam="output/{sample}.bam"
    shell:
        """
        samtools view -bS {input.sam} > {output}
        """

rule sort_bam:
    input:
        bam="output/{sample}.bam"
    output:
        sorted_bam="output/{sample}.sorted.bam"
    params:
        threads=4
    shell:
        """
        samtools sort -@ {params.threads} -o {output.sorted_bam} {input.bam}
        """

rule mark_duplicates:
    input:
        bam="output/{sample}.sorted.bam",
        idx_temp_dir = idx_temp_dir
    output:
        dedup_bam="output/{sample}.dedup.bam",
        metrics="output/{sample}.dedup.metrics"
    params:
        java_options="-Xmx4G"
    shell:
        """
        gatk MarkDuplicates --java-options "{params.java_options}" \
            -I {input.bam} \
            -O {output.dedup_bam} \
            -M {output.metrics} \
            --TMP_DIR {idx_temp_dir}
        """

rule base_recalibrator:
    input:
        dedup_bam="output/{sample}.dedup.bam",
        idx=idx,
        idx_dbsnp=idx_dbsnp,
        idx_Mills_1000G=idx_Mills_1000G,
        idx_indels=idx_indels
    output:
        recal_table="output/{sample}.recal.table"
    params:
        java_options="-Xmx4G"
    shell:
        """
        gatk BaseRecalibrator --java-options "{params.java_options}" -I {input.dedup_bam} -R {idx} -O {output.recal_table} \
        --known-sites {idx_dbsnp} --known-sites {idx_Mills_1000G} --known-sites {idx_indels}
        """

rule apply_bqsr:
    input:
        dedup_bam="output/{sample}.dedup.bam",
        recal_table="output/{sample}.recal.table",
        idx=idx
    output:
        recalibrated_bam="output/{sample}.recalibrated.bam"
    shell:
        """
        gatk ApplyBQSR -R {idx} -I {input.dedup_bam} --bqsr-recal-file {input.recal_table} -O {output.recalibrated_bam}
        """

rule haplotype_caller:
    input:
        recalibrated_bam="output/{sample}.recalibrated.bam",
        idx=idx
    output:
        hc_vcf="output/{sample}.hc.vcf"
    shell:
        """
        gatk HaplotypeCaller -R {idx} -I {input.recalibrated_bam} -O {output.hc_vcf}
        """

rule select_variants_snp:
    input:
        hc_vcf="output/{sample}.hc.vcf",
        idx=idx
    output:
        selvar_snp_vcf="output/{sample}.selvarSNP.vcf"
    shell:
        """
        gatk SelectVariants -R {idx} -V {input.hc_vcf} --select-type-to-include SNP -O {output.selvar_snp_vcf}
        """

rule variant_filtration_snp:
    input:
        selvar_snp_vcf="output/{sample}.selvarSNP.vcf",
        idx=idx
    output:
        snp_filtered_vcf="output/{sample}.snp_filtered.vcf"
    shell:
        """
        gatk VariantFiltration -V {input.selvar_snp_vcf} --filter-expression "QD < 2.0 || QUAL < 30.0 || SOR > 3.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "HardFilter" -O {output.snp_filtered_vcf} -R {idx}
        """

rule select_variants_indel:
    input:
        hc_vcf="output/{sample}.hc.vcf",
        idx=idx
    output:
        selvar_indel_vcf="output/{sample}.selvarINDEL.vcf"
    shell:
        """
        gatk SelectVariants -R {idx} -V {input.hc_vcf} --select-type-to-include INDEL -O {output.selvar_indel_vcf}
        """

rule variant_filtration_indel:
    input:
        selvar_indel_vcf="output/{sample}.selvarINDEL.vcf",
        idx=idx
    output:
        indel_filtered_vcf="output/{sample}.indel_filtered.vcf"
    shell:
        """
        gatk VariantFiltration -V {input.selvar_indel_vcf} --filter-expression "QD < 2.0 || QUAL < 30.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filter-name "HardFilter" -O {output.indel_filtered_vcf} -R {idx}
        """

rule merge_vcfs:
    input:
        snp_filtered_vcf="output/{sample}.snp_filtered.vcf",
        indel_filtered_vcf="output/{sample}.indel_filtered.vcf"
    output:
        merged_filtered_vcf="output/{sample}.merged_filtered.vcf"
    shell:
        """
        gatk MergeVcfs -I {input.snp_filtered_vcf} -I {input.indel_filtered_vcf} -O {output.merged_filtered_vcf}
        """

rule snpEff:
    input:
        snpeff_vcf="output/{sample}.merged_filtered.vcf",
        idx_snpeff_jar=idx_snpeff_jar,
        idx_snpeff_config=idx_snpeff_config,
        idx_java_path=idx_java_path
    output:
        snpeff_annot_vcf="output/{sample}.merged_filtered_snpeff.vcf"
    params:
        readgroup=lambda wildcards: wildcards.sample
    shell:
        """
        export PATH={idx_java_path}:$PATH

        java -jar {idx_snpeff_jar} -v -stats output/{params.readgroup}.merged_filtered.stats -c {idx_snpeff_config} -chr chr -classic -o vcf hg38 {input.snpeff_vcf} > {output.snpeff_annot_vcf}
        """

ruleorder:
    bwa_mem > sam_to_bam > sort_bam > mark_duplicates > base_recalibrator > apply_bqsr > haplotype_caller > select_variants_snp > variant_filtration_snp > select_variants_indel > variant_filtration_indel > merge_vcfs > snpEff
