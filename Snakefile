"""
Dans ce qui suit, nous présenterons la syntaxe Snakemake en créant un exemple de workflow. 
Le flux de travail provient du domaine de l’analyse du génome. Il mappe les lectures de 
séquençage sur un génome de référence et appelle des variantes sur les lectures mappées.

SNAKEMAKE
Regle = nom + directives
execution du workflow : ca genere les fichiers cibles
Snakemake ne réexécute les tâches que si l'un des fichiers d'entrée est plus récent que l'un des fichiers de sortie ou si l'un des fichiers d'entrée sera mis à jour par une autre tâche .
"""

# ALL
rule all:
    input:
        "plots/quals.svg"

#----------------------------------------------------------------------------------#
# Mappage des reads
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"

# tri des alignements de reads
rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"

# Indexation des alignements de reads et visualisation du DAG des taches
rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"


# Appel a des variantes gnomiques
SAMPLES = ["A", "B"]
rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
    output:
        "calls/all.vcf"
    shell:
        "bcftools mpileup -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"

# Utilisation de script perso
"""
Nous générerons éventuellement un histogramme des scores de qualité qui ont été attribués aux appels de variantes dans le fichier calls/all.vcf
"""
rule plot_quals:
    input:
        "calls/all.vcf"
    output:
        "plots/quals.svg"
    script:
        "scripts/plot-quals.py"
