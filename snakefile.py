import os
import glob from glob

ID, = glob_wildcards("MO128/raw_reads/{id}.fastq")

rule all:
    input:
        expand("MO128/raw_reads_QC/{id}_fastqc.{extension}", id = ID, extension = ["zip", "html"], allow_missing = True),
        expand("MO128/filtered_reads/{id}_filtered.fastq", id = ID),
        expand("MO128/filtered_reads_QC/{id}_filtered_fastqc.{extension}", id = ID, extension = ["zip", "html"], allow_missing = True),
        directory(expand("MO128/assembly_result/{id}/", id = ID)),
        expand("MO128/mapping_result/{id}.sam", id = ID),
        expand("MO128/consensus_sequence/{id}.fasta", id = ID),
        expand("MO128/mapping_result/{id}.bam", id = ID),
        expand("MO128/sorted_mapping_result/{id}_sorted.bam" , id = ID),
        expand("MO128/sorted_mapping_result/{id}_sorted.bam.bai" , id = ID),
        expand("MO128/sorted_fastq/{id}.fastq", id = ID),
        expand("MO128/spike_reference/reference_mapping/{id}.bam", id = ID),
        expand("MO128/rbd_reference/reference_mapping/{id}.bam", id = ID),
        expand("MO128/spike_reference/sorted_reference_mapping_result/{id}_sorted.bam", id = ID),
        expand("MO128/rbd_reference/sorted_reference_mapping_result/{id}_sorted.bam", id = ID),
        expand("MO128/spike_reference/sorted_reference_mapping_result/{id}_sorted.bam.bai", id = ID),
        expand("MO128/rbd_reference/sorted_reference_mapping_result/{id}_sorted.bam.bai", id = ID),
        expand("MO128/spike_reference/basecounts_result/{id}_Spike", id = ID),
        expand("MO128/rbd_reference/basecounts_result/{id}_RBD", id = ID),
        expand("MO128/rbd_reference/basecounts_ref_result/{id}_RBD.tsv", id = ID),
        expand("MO128/spike_reference/basecounts_ref_result/{id}_Spike.tsv", id = ID)

rule rawFastqc:
    message:
        'Performing Quality control check...'
    input:
        rawread = "MO128/raw_reads/{id}.fastq"
    output:
        zip = "MO128/raw_reads_QC/{id}_fastqc.zip",
        html = "MO128/raw_reads_QC/{id}_fastqc.html"
    params:
        path = "MO128/raw_reads_QC"
    shell:
        """
        fastqc {input.rawread} -o {params.path}
        """
rule nano_trim:
    message:
        'Trimming with NanoFilt...'
    input:
        "MO128/raw_reads/{id}.fastq"
    output:
        "MO128/filtered_reads/{id}_filtered.fastq"      
    shell:
        """
        NanoFilt -l 650 --headcrop 10 -q 10 {input} > {output}
        """
        
rule filtFastqc:
    message:
        'Performing Quality control check...'
    input:
        filtread = "MO128/filtered_reads/{id}_filtered.fastq"
    output:
        zip = "MO128/filtered_reads_QC/{id}_filtered_fastqc.zip",
        html = "MO128/filtered_reads_QC/{id}_filtered_fastqc.html"
    params:
        path = "MO128/filtered_reads_QC"
    shell:
        """
        fastqc {input.filtread} -o {params.path}
        """
rule canu_run:
    message:
        'Running CANU for assembly...'
    input:
        "MO128/raw_reads/{id}.fastq"
    output:
        dir = directory("MO128/assembly_result/{id}/")
    shell:
        """
        canu useGrid=false -p {wildcards.id} -d {output.dir} genomeSize=4k -nanopore {input}
        """
rule mapping_run:
    message:
        'Running minimap2 for mapping...'
    input:
        contigs = "MO128/assembly_result/{id}/{id}.contigs.fasta",
        filtreads = "MO128/filtered_reads/{id}_filtered.fastq"
    output:
        "MO128/mapping_result/{id}.sam"
    shell:
        """
        minimap2 -a {input.contigs} {input.filtreads} > {output}
        """
rule racon_run:
    message:
        'Running racon for creating consensus sequence...'
    input:
        filtreads = "MO128/filtered_reads/{id}_filtered.fastq",
        samfile = "MO128/mapping_result/{id}.sam",
        contigs = "MO128/assembly_result/{id}/{id}.contigs.fasta"
    output:
        "MO128/consensus_sequence/{id}.fasta"
    shell:
        """
        racon -m 10 -t 14 {input.filtreads} {input.samfile} {input.contigs} > {output}
        """				
rule creating_bam:
    message:
        'Running minimap2 for mapping...'
    input:
        "MO128/mapping_result/{id}.sam"
    output:
        "MO128/mapping_result/{id}.bam"
    shell:
        """
        samtools view -bh {input} > {output}
        """
rule samtools_sort_run:
    message:
        'Sort mapping by samtools...'
    input:
        unsorted = "MO128/mapping_result/{id}.bam"
    output:
        sorted = "MO128/sorted_mapping_result/{id}_sorted.bam" 
    shell:
        """
        samtools sort -o {output.sorted} {input.unsorted} 
        """
rule samtools_indexing_run:
    message:
        'Index mapping by samtools...'
    input:
        bamfile = "MO128/sorted_mapping_result/{id}_sorted.bam"
    output:
        indexfile = "MO128/sorted_mapping_result/{id}_sorted.bam.bai" 
    shell:
        """
        samtools index {input.bamfile}
        """
rule sam_to_fasq:
    message:
        'Converting sam to fastq...'
    input:
        bamfile = "MO128/sorted_mapping_result/{id}_sorted.bam"
    output:
        fastqfile = "MO128/sorted_fastq/{id}.fastq" 
    shell:
        """
        java -jar picard/picard.jar SamToFastq I={input.bamfile} FASTQ={output.fastqfile}
        """
rule spike_reference_mapping:
    message:
        'Alignment of fastq reads to the spike reference...'
    input:
        reference = "MO128/reference/Spike_WT.fa",
        fastqfile = "MO128/sorted_fastq/{id}.fastq"
    output:
        "MO128/spike_reference/reference_mapping/{id}.bam"
    shell:
        """
        minimap2 -ax map-ont {input.reference} {input.fastqfile} | samtools view -Sb - > {output}
        """
rule rbd_reference_mapping:
    message:
        'Alignment of fastq reads to the rbd reference......'
    input:
        reference = "MO128/reference/RBD_WT.fa",
        fastqfile = "MO128/sorted_fastq/{id}.fastq"
    output:
        "MO128/rbd_reference/reference_mapping/{id}.bam"
    shell:
        """
        minimap2 -ax map-ont {input.reference} {input.fastqfile} | samtools view -Sb - > {output}
        """
rule spike_samtools_sorting_by_reference:
    message:
        'Sort reference mapping by samtools...'
    input:
        unsorted = "MO128/spike_reference/reference_mapping/{id}.bam"
    output:
        sorted = "MO128/spike_reference/sorted_reference_mapping_result/{id}_sorted.bam" 
    shell:
        """
        samtools sort -o {output.sorted} {input.unsorted} 
        """
rule rbd_samtools_sorting_by_reference:
    message:
        'Sort reference mapping by samtools...'
    input:
        unsorted = "MO128/rbd_reference/reference_mapping/{id}.bam"
    output:
        sorted = "MO128/rbd_reference/sorted_reference_mapping_result/{id}_sorted.bam" 
    shell:
        """
        samtools sort -o {output.sorted} {input.unsorted} 
        """
rule spike_index_sorting_by_reference:
    message:
        'Index mapping by samtools...'
    input:
        bamfile = "MO128/spike_reference/sorted_reference_mapping_result/{id}_sorted.bam"
    output:
        indexfile = "MO128/spike_reference/sorted_reference_mapping_result/{id}_sorted.bam.bai" 
    shell:
        """
        samtools index {input.bamfile}
        """
rule rbd_index_sorting_by_reference:
    message:
        'Index mapping by samtools...'
    input:
        bamfile = "MO128/rbd_reference/sorted_reference_mapping_result/{id}_sorted.bam"
    output:
        indexfile = "MO128/rbd_reference/sorted_reference_mapping_result/{id}_sorted.bam.bai" 
    shell:
        """
        samtools index {input.bamfile}
        """   
rule basecount_against_spike_run:
    message:
        'Basecounting against a reference...'
    input:
        reftsv = "MO128/reference/Spike_WT.tsv",
        bamfile = "MO128/spike_reference/sorted_reference_mapping_result/{id}_sorted.bam",
        reference = "MO128/reference/Spike_WT.fa"
    output:
        "MO128/spike_reference/basecounts_result/{id}_Spike" 
    shell:
        """
        ./basecounts/basecounts -f {input.reference} {input.reftsv} {input.bamfile} -o {output}
         """
rule basecount_against_rbd_run:
    message:
        'Basecounting against a reference...'
    input:
        reftsv = "MO128/reference/RBD_WT.tsv",
        bamfile = "MO128/rbd_reference/sorted_reference_mapping_result/{id}_sorted.bam",
        reference = "MO128/reference/RBD_WT.fa"
    output:
        "MO128/rbd_reference/basecounts_result/{id}_RBD" 
    shell:
        """
        ./basecounts/basecounts -f {input.reference} {input.reftsv} {input.bamfile} -o {output}
         """
rule modifying_rbd_counts:
    input:
        counts = "MO128/rbd_reference/basecounts_result/{id}_RBD.tsv",
        reference = "MO128/reference/RBD_WT_fa_ref.txt"
    output:
        "MO128/rbd_reference/basecounts_ref_result/{id}_RBD.tsv"
    shell:
        """
        paste <(cut -f 1 {input.counts}) <(cut -f 1 {input.reference}) <(cut -f 3- {input.counts}) > {output}
         """
rule modifying_spike_counts:
    input:
        counts = "MO128/spike_reference/basecounts_result/{id}_Spike.tsv",
        reference = "MO128/reference/Spike_WT_fa_ref.txt"
    output:
        "MO128/spike_reference/basecounts_ref_result/{id}_Spike.tsv"
    shell:
        """
        paste <(cut -f 1 {input.counts}) <(cut -f 1 {input.reference}) <(cut -f 3- {input.counts}) > {output}
         """












