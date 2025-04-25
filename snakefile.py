import os
from glob import glob

ID, = glob_wildcards("directory/input/fastq_replicate_1/raw_reads/{id}.fastq")

rule all:
    input:
        expand("directory/input/fastq_replicate_1/raw_reads_QC/{id}_fastqc.{extension}", id = ID, extension = ["zip", "html"], allow_missing = True),
        expand("directory/input/fastq_replicate_1/filtered_reads/{id}_filtered.fastq", id = ID),
        expand("directory/input/fastq_replicate_1/filtered_reads_QC/{id}_filtered_fastqc.{extension}", id = ID, extension = ["zip", "html"], allow_missing = True),
        expand("directory/input/fastq_replicate_1/assembly_result/{id}/{id}.contigs.fasta", id = ID),
        expand("directory/input/fastq_replicate_1/mapping_result/{id}.sam", id = ID),
        expand("directory/input/fastq_replicate_1/consensus_sequence/{id}.fasta", id = ID),
        expand("directory/input/fastq_replicate_1/mapping_result/{id}.bam", id = ID),
        expand("directory/input/fastq_replicate_1/sorted_mapping_result/{id}_sorted.bam" , id = ID),
        expand("directory/input/fastq_replicate_1/sorted_mapping_result/{id}_sorted.bam.bai" , id = ID),
        expand("directory/input/fastq_replicate_1/sorted_fastq/{id}.fastq", id = ID),
        expand("directory/input/fastq_replicate_1/rbd_reference/reference_mapping/{id}.bam", id = ID),
        expand("directory/input/fastq_replicate_1/rbd_reference/sorted_reference_mapping_result/{id}_sorted.bam", id = ID),
        expand("directory/input/fastq_replicate_1/rbd_reference/sorted_reference_mapping_result/{id}_sorted.bam.bai", id = ID),
        expand("directory/input/fastq_replicate_1/rbd_reference/basecounts_ref_result/{id}_RBD", id = ID)
        
        

rule rawFastqc:
    message:
        'Performing Quality control check...'
    input:
        rawread = "directory/input/fastq_replicate_1/raw_reads/{id}.fastq"
    output:
        zip = "directory/input/fastq_replicate_1/raw_reads_QC/{id}_fastqc.zip",
        html = "directory/input/fastq_replicate_1/raw_reads_QC/{id}_fastqc.html"
    params:
        path = "directory/input/fastq_replicate_1/raw_reads_QC"
    shell:
        """
        fastqc {input.rawread} -o {params.path}
        """

rule nano_trim:
    message:
        'Trimming with NanoFilt...'
    input:
        "directory/input/fastq_replicate_1/raw_reads/{id}.fastq"
    output:
        "directory/input/fastq_replicate_1/filtered_reads/{id}_filtered.fastq"      
    shell:
        """
        NanoFilt -l 650 --headcrop 10 -q 10 {input} > {output}
        """

rule filtFastqc:
    message:
        'Performing Quality control check...'
    input:
        filtread = "directory/input/fastq_replicate_1/filtered_reads/{id}_filtered.fastq"
    output:
        zip = "directory/input/fastq_replicate_1/filtered_reads_QC/{id}_filtered_fastqc.zip",
        html = "directory/input/fastq_replicate_1/filtered_reads_QC/{id}_filtered_fastqc.html"
    params:
        path = "directory/input/fastq_replicate_1/filtered_reads_QC"
    shell:
        """
        fastqc {input.filtread} -o {params.path}
        """

rule canu_run:
    message:
        'Running CANU for assembly...'
    input:
        "directory/input/fastq_replicate_1/raw_reads/{id}.fastq"
    output:
        contigs = "directory/input/fastq_replicate_1/assembly_result/{id}/{id}.contigs.fasta",
        outdir = directory("directory/input/fastq_replicate_1/assembly_result/{id}")
    shell:
        """
        canu useGrid=false -p {wildcards.id} -d {output.outdir} genomeSize=4k -nanopore {input}
        """

rule mapping_run:
    message:
        'Running minimap2 for mapping...'
    input:
        contigs = "directory/input/fastq_replicate_1/assembly_result/{id}/{id}.contigs.fasta",
        filtreads = "directory/input/fastq_replicate_1/filtered_reads/{id}_filtered.fastq"
    output:
        "directory/input/fastq_replicate_1/mapping_result/{id}.sam"
    shell:
        """
        minimap2 -a {input.contigs} {input.filtreads} > {output}
        """

rule racon_run:
    message:
        'Running racon for creating consensus sequence...'
    input:
        filtreads = "directory/input/fastq_replicate_1/filtered_reads/{id}_filtered.fastq",
        samfile = "directory/input/fastq_replicate_1/mapping_result/{id}.sam",
        contigs = "directory/input/fastq_replicate_1/assembly_result/{id}/{id}.contigs.fasta"
    output:
        "directory/input/fastq_replicate_1/consensus_sequence/{id}.fasta"
    shell:
        """
        racon -m 10 -t 14 {input.filtreads} {input.samfile} {input.contigs} > {output}
        """
				
rule creating_bam:
    message:
        'Running minimap2 for mapping...'
    input:
        "directory/input/fastq_replicate_1/mapping_result/{id}.sam"
    output:
        "directory/input/fastq_replicate_1/mapping_result/{id}.bam"
    shell:
        """
        samtools view -bh {input} > {output}
        """

rule samtools_sort_run:
    message:
        'Sort mapping by samtools...'
    input:
        unsorted = "directory/input/fastq_replicate_1/mapping_result/{id}.bam"
    output:
        sorted = "directory/input/fastq_replicate_1/sorted_mapping_result/{id}_sorted.bam" 
    shell:
        """
        samtools sort -o {output.sorted} {input.unsorted} 
        """

rule samtools_indexing_run:
    message:
        'Index mapping by samtools...'
    input:
        bamfile = "directory/input/fastq_replicate_1/sorted_mapping_result/{id}_sorted.bam"
    output:
        indexfile = "directory/input/fastq_replicate_1/sorted_mapping_result/{id}_sorted.bam.bai" 
    shell:
        """
        samtools index {input.bamfile}
        """

rule sam_to_fasq:
    message:
        'Converting sam to fastq...'
    input:
        bamfile = "directory/input/fastq_replicate_1/sorted_mapping_result/{id}_sorted.bam"
    output:
        fastqfile = "directory/input/fastq_replicate_1/sorted_fastq/{id}.fastq" 
    shell:
        """
        picard SamToFastq I={input.bamfile} FASTQ={output.fastqfile}
        """

rule rbd_reference_mapping:
    message:
        'Alignment of fastq reads to the rbd reference......'
    input:
        reference = "directory/input/fastq_replicate_1/reference/RBD_WT.fa",
        fastqfile = "directory/input/fastq_replicate_1/sorted_fastq/{id}.fastq"
    output:
        "directory/input/fastq_replicate_1/rbd_reference/reference_mapping/{id}.bam"
    shell:
        """
        minimap2 -ax map-ont {input.reference} {input.fastqfile} | samtools view -Sb - > {output}
        """

rule rbd_samtools_sorting_by_reference:
    message:
        'Sort reference mapping by samtools...'
    input:
        unsorted = "directory/input/fastq_replicate_1/rbd_reference/reference_mapping/{id}.bam"
    output:
        sorted = "directory/input/fastq_replicate_1/rbd_reference/sorted_reference_mapping_result/{id}_sorted.bam" 
    shell:
        """
        samtools sort -o {output.sorted} {input.unsorted} 
        """

rule rbd_index_sorting_by_reference:
    message:
        'Index mapping by samtools...'
    input:
        bamfile = "directory/input/fastq_replicate_1/rbd_reference/sorted_reference_mapping_result/{id}_sorted.bam"
    output:
        indexfile = "directory/input/fastq_replicate_1/rbd_reference/sorted_reference_mapping_result/{id}_sorted.bam.bai" 
    shell:
        """
        samtools index {input.bamfile}
        """   

rule basecount_against_rbd_run:
    message:
        'Basecounting against a reference...'
    input:
        reftsv = "directory/input/fastq_replicate_1/reference/RBD_WT.tsv",
        bamfile = "directory/input/fastq_replicate_1/rbd_reference/sorted_reference_mapping_result/{id}_sorted.bam",
        reference = "directory/input/fastq_replicate_1/reference/RBD_WT.fa" 
    output:
    	"directory/input/fastq_replicate_1/rbd_reference/basecounts_ref_result/{id}_RBD"
    shell:
        """
        ./basecounts/basecounts -f {input.reference} {input.reftsv} {input.bamfile} -o {output}
         """






