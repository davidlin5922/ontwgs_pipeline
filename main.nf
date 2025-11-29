#!/usr/bin/env nextflow

// Include modules
// include { qc } from './modules/qc.nf'
// include { telLength } from './modules/telLength.nf'
// include { alignment } from './modules/alignment.nf'
// include { sv } from './modules/sv.nf'

def reference_genome_link

// get reference genome
if (params.reference_genome in ['hg19', 'GRCh37']) {
    reference_genome_link = "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz"
} else if (params.reference_genome in ['hg38', 'GRCh38']) {
    reference_genome_link = "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz"
} else if (params.reference_genome in ['CHM13', 'T2T', 't2t']) {
    reference_genome_link = "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz"
} else {
    error "Error: unknown reference genome name: ${params.reference_genome} (try: hg19|GRCh37 / hg38|GRCh38 / CHM13|T2T|t2t)"
}

// Read sample names and fastq paths
sample_names_ch = Channel.fromPath(params.sample_names).splitCsv().map{ line -> line[0] }
fastqs_ch = Channel.fromPath(params.fastq_paths).splitCsv().map{ line -> line[0] }
samples_ch = sample_names_ch.merge(fastqs_ch)

process download_ref {

    input:
        val reference_genome_link

    output:
        path "${params.reference_genome}.fa"

    script:
    """
    curl ${reference_genome_link} > ${params.reference_genome}.fa.gz
    gunzip ${params.reference_genome}.fa.gz
    """

}

process qc {

    publishDir "results/qc", mode: 'copy'
    conda 'bioconda::nanoplot==1.46.1'

    input:
        tuple val(sample_id), path(fastq)

    output:
        path "${sample_id}"

    script:
    """
    NanoPlot -o ./${sample_id} --fastq ${fastq} --huge --tsv_stats -c '#1f77b4' --N50 --title ' '
    mv ./${sample_id}/*.html ./
    """

}

process alignment {

    publishDir "results/bams", mode: 'copy'
    conda 'bioconda::minimap2==2.30 bioconda::samtools==1.22.1'

    input:
        tuple val(sample_id), path(fastq)
        path ref_fa

    output:
        tuple val(sample_id), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai")

    script:
    """
    minimap2 -ax map-ont ${ref_fa} ${fastq} | samtools view -bS | samtools sort -o ${sample_id}_sorted.bam
    samtools index ${sample_id}_sorted.bam
    """

}

process bam_qc {
    
    conda 'bioconda::samtools==1.22.1'

    input:
        tuple val(sample_id), path(bam), path(bai)
        path ref_fa

    output:
        tuple val(sample_id), path("summary.tsv"), path("readLength.tsv"), path("mapQuality.tsv"), path("depth.tsv")

    script:
    """
    samtools stats -r ${ref_fa} ${sample_id}_sorted.bam > ${sample_id}_bamqc.txt
    grep ^SN ${sample_id}_bamqc.txt | cut -f 2- > summary.tsv
    grep ^RL ${sample_id}_bamqc.txt | cut -f 2- > readLength.tsv
    grep ^MAPQ ${sample_id}_bamqc.txt | cut -f 2- > mapQuality.tsv
    grep ^COV ${sample_id}_bamqc.txt | cut -f 2- > depth.tsv
    """

}

process plotBam_qc {

    cache false
    publishDir "results/bamqc", mode: 'copy'
    conda "anaconda::python=3.12 anaconda::pandas==2.3.3 anaconda::matplotlib==3.10.6"
    
    input:
        tuple val(sample_id), path(summary), path(readLength), path(mapQuality), path(depth)

    output:
        path("${sample_id}")

    script:
    """
    mkdir ${sample_id}
    cp ${summary} ./${sample_id}
    python3 ${projectDir}/scripts/plot_bamqc.py ${readLength} ${mapQuality} ${depth} ${sample_id}
    """
}

process telLength {

    publishDir "results/telogator2", mode: 'copy'
    conda "${projectDir}/scripts/install_telogator2.yaml"

    input:
        tuple val(sample_id), path(fastq)

    output:
        path("${sample_id}")

    script:
    """
    telogator2 -i ${fastq} -o . --minimap2 \$(which minimap2) -r ont
    mkdir ${sample_id}
    mv ./qc/qc_readlens.png ./${sample_id}/readlens.png
    mv all_final_alleles.png ./${sample_id}/all_final_alleles.png
    mv violin_atl.png ./${sample_id}/violin_atl.png
    mv tlens_by_allele.tsv ./${sample_id}/tlens.tsv
    """

}

process sniffles2 {

    publishDir "results/sniffles2", mode: 'copy'
    conda 'bioconda::sniffles==2.7.1'

    input:
        tuple val(sample_id), path(bam), path(bai)
        path ref_fa

    output:
        path "${sample_id}.vcf"
    
    script:
    """
    sniffles -i ${bam} -v ${sample_id}.vcf --reference ${ref_fa}
    """

}

process sniffles2_snf {

    conda 'bioconda::sniffles==2.7.1'

    input:
        tuple val(sample_id), path(bam), path(bai)
        path ref_fa

    output:
        path "${sample_id}.snf"
    
    script:
    """
    sniffles -i ${bam} --snf ${sample_id}.snf --reference ${ref_fa}
    """

}

process plotSV {

    publishDir "results", mode: 'copy'
    conda 'hcc::sniffles2-plot==0.2.1'

    input:
        path(vcfs)

    output:
        path("sniffles2")

    script:
    """
    mkdir sniffles2
    mv ${vcfs} ./sniffles2
    python3 -m sniffles2_plot -i ./sniffles2
    """
    
}

process combine_snf {

publishDir "results/sniffles2", mode: 'copy'
    conda 'bioconda::sniffles==2.7.1'

    input:
        path(snfs)

    output:
        path("all_sample.vcf")

    script:
    """
    sniffles --input ${snfs} -v all_sample.vcf
    """

}

process plotSV_all {

    publishDir "results/sniffles2", mode: 'copy'
    conda 'hcc::sniffles2-plot==0.2.1'

    input:
        path(snfs)

    output:
        path("all")

    script:
    """
    python3 -m sniffles2_plot -i all_sample.vcf -o ./all
    """

}

process summary {

    cache false
    publishDir "results", mode: 'copy'
    conda "conda-forge::r-base==4.5.2 conda-forge::r-rmarkdown==2.30"

    input:
        val(sample_id)
        val(qc)
        val(bamqc)
        val(svSummary)
        val(allSV)
        val(reference)

    output:
        path "report.html"

    script:
    """
    Rscript -e '
        rmarkdown::render("${projectDir}/scripts/report.rmd",
                          params = list(sample_name = c("${sample_id.join('", "')}"),
                                        qc = c("${qc.join('", "')}"),
                                        bamqc = c("${bamqc.join('", "')}"),
                                        svSummary = c("${svSummary.join('", "')}"),
                                        allSV = "${allSV}",
                                        reference = "${reference}",
                                        processDir = "${workDir}"),
                          output_file= "report.html",
                          output_dir = getwd())
    '
    """

}

workflow {

    // download reference
    ref_ch = download_ref(reference_genome_link)

    // QC
    qc_ch = qc(samples_ch)

    // alignment
    aligned_ch = alignment(samples_ch, ref_ch)

    // bam QC
    bamQC_ch = bam_qc(aligned_ch, ref_ch)
    plotBamQC_ch = plotBam_qc(bamQC_ch)

    // telomere length
    // telLength_ch = telLength(samples_ch)

    // sv calling
    sv_ch = sniffles2(aligned_ch, ref_ch)
    sv_snf_ch = sniffles2_snf(aligned_ch, ref_ch)
    
    plotSV_ch = plotSV(sv_ch.collect())
    combine_snf_ch = combine_snf(sv_snf_ch.collect())
    plotSV_all_ch = plotSV_all(combine_snf_ch)

    id = sample_names_ch.collect()
    qc_results = qc_ch.collect()
    mapqc_results = plotBamQC_ch.collect()
    // telLength_ch.collect()
    sv_results = plotSV_ch.collect()
    allSV_results = plotSV_all_ch
    // summary
    summary(id, qc_results, mapqc_results, sv_results, allSV_results, params.reference_genome)

}
