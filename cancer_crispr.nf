/*
 * Parameters
 */

params.reads = "$baseDir/data/library.fa"
params.index = "$baseDir/data/index/"
params.outdir = "results"

log.info """\
         Running Cancer-CRISPR
         ===================================
         index        : ${params.index}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()

/*
 * This process takes library.fa converts it to .fq using seqtk
 * maps the library.fq to pre-built hg38 index using bowtie2
 * outputs mapped sgRNAs into read_locs.tsv using samtools
 */

process bowtie2 {

    input:
    path reads from params.reads
    path index from params.index

    output:
    path 'read_locs.tsv' into bowtie2_out_ch

    script:
    """
    seqtk seq -F '#' ${reads} > reads.fq
    idx_base=\$(find ${index}/ -name '*.bt2' | awk -F \".\" '{print \$1 | \"sort -u\"}')
    bowtie2 -x \${idx_base} -U reads.fq -S mapped.sam 
    samtools fixmate -O bam,level=1 mapped.sam fixmate.bam
    samtools sort -l 1 -@8 -o pos.srt.bam -T /tmp/example_prefix fixmate.bam
    samtools markdup -O bam,level=1 pos.srt.bam markdup.bam
    samtools view -@8 markdup.bam -o final.bam
    samtools view final.bam > read_locs.tsv
    """

}

/*
 * This process takes mapped reads from read_locs.tsv
 * calculates end position of reads by the read lenght
 * removes unnecessary columns and checks if
 * provided sgRNA target matches the mapped gene using an Rscript
 */

process get_mapped_reads {

    publishDir params.outdir, mode:'copy'

    input:
    file bowtie2_out from bowtie2_out_ch

    output:
    path 'reads_mapped_table.tsv' into get_mapped_reads_out_ch

    script:
    """
    Rscript $baseDir/bin/process_map_data.R $bowtie2_out
    """
}

/*
 * This process takes mapped reads from read_locs.txt
 * calculates end position of reads by the read lenght
 * removes unnecessary columns and checks if
 * provided sgRNA target matches the mapped gene using an Rscript
 */

process get_tcga_data {

    publishDir params.outdir, mode:'copy'

    input:
    file get_mapped_reads_out from get_mapped_reads_out_ch

    output:
    path 'tcga_data_table.tsv'

    script:
    """
    Rscript $baseDir/bin/get_tcga_data.R $get_mapped_reads_out
    """
}
