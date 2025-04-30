#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.vpot_params_abs= file(params.vpot_params).toAbsolutePath().toString()

vcfFile = file(params.vcf)
if( !vcfFile.exists() ) {
    exit 1, "The specified VCF file does not exist: ${params.vcf}"
}
sampleName = vcfFile.baseName

genePanelFile = file(params.gene_panels)
if( !genePanelFile.exists() ) {
    exit 1, "The specified gene panel file does not exist: $params.gene_panels"
}

if (params.column_file != "None") {
    columnFile = file(params.column_file)
} else {
    columnFile = "None"
}

annotatedVcf = file("${params.output_dir}/${sampleName}.${params.ref_name}_multianno.vcf")

process splitVariants { 
    debug true
    input:
    path vcf

    output:
    path "${sampleName}_split.vcf"

    shell:
    """
    bcftools norm -m-both -o ${sampleName}_split.vcf ${vcf} --fasta-ref ${params.ref_fasta} --check-ref w
    """
}

process annotateGene {
    debug true
    input:
    path vcf

    output:
    path "${sampleName}.${params.ref_name}_multianno.vcf"

    script:
    def filterArg = params.enable_pass_filter ? '--convertarg "--filter \'pass\'"' : ''

    """
    ${params.annovar_dir}/table_annovar.pl ${vcf} ${params.humandb_dir} -buildver ${params.ref_name} \
    -out ${sampleName} -remove -protocol refgene -operation g \
    -nastring . ${filterArg} -vcfinput -thread 12
    """
}

process filterByGene {
    input:
    path anno_vcf

    output:
    path "${sampleName}_genefiltered.vcf"

    script:
    """
    cut -d, -f1 ${genePanelFile} | tail -n +2 > genelist.txt
    grep "^#" ${anno_vcf} > ${sampleName}_genefiltered.vcf
    grep -f genelist.txt ${anno_vcf} >> ${sampleName}_genefiltered.vcf
    """
}

process cleanAnnovarAnnotations {
    input:
    path filtered_vcf

    output:
    path "${sampleName}_genefiltered_clean.vcf"

    script:
    """
    sed -e 's/ANNOVAR_DATE=20[0-9][0-9]-[0-9][0-9]-[0-9][0-9];//' \
        -e 's/;ALLELE_END//' ${filtered_vcf} > ${sampleName}_genefiltered_clean.vcf
    """
}

process annotateAll {
    debug true
    input:
    path filtered_vcf

    output:
    path "${sampleName}.${params.ref_name}_multianno.vcf"

    publishDir params.output_dir, mode: 'copy', pattern: '{*_multianno.vcf}'

    shell:
    """
    ${params.annovar_dir}/table_annovar.pl ${filtered_vcf} ${params.humandb_dir} -buildver ${params.ref_name} \
        -out ${sampleName} -remove ${params.annovar_params} -nastring . -vcfinput
    """
}

process vpotPrioritise {
    debug true
    input:
    path vcf

    output:
    path "${sampleName}_final_output_file.txt"

    publishDir params.output_dir, mode: 'copy', pattern: '*_final_output_file.txt'

    shell:
    '''
    echo "!{vcf}	$(grep "#CHROM" !{vcf} | awk '{print $NF}')" > vpot_input.txt
    python !{params.vpot_dir}/VPOT.py priority !{sampleName}_ \
        vpot_input.txt !{params.vpot_params_abs}
    '''
}

process vpotGenePanel {
    debug true
    input:
    path vpol

    output:
    path "${sampleName}_output_genepanels.xlsx"

    publishDir params.output_dir, mode: 'copy', pattern: "${sampleName}_output_genepanels.xlsx"

    shell:
    """
    python ${params.vpot_dir}/VPOT.py genepanelf "${sampleName}_" $vpol $genePanelFile $params.panel_name $columnFile
    """
}

workflow {
    if( annotatedVcf.exists()) {
        println "Already annotated VCF!"
        vpotPrioritise(annotatedVcf)
    } else {
        if ( params.pre_annotated ) {
            println "Gene-filtering pre-annotated vcf"
            filterByGene(vcfFile)
            vpotPrioritise(filterByGene.out)
        } else {
            if ( params.normalise ) {
                splitVariants(vcfFile)
                annotateGene(splitVariants.out)
            } else {
                annotateGene(vcfFile)
            }
            filterByGene(annotateGene.out)
            cleanAnnovarAnnotations(filterByGene.out)
            annotateAll(cleanAnnovarAnnotations.out)
            vpotPrioritise(annotateAll.out)
        }
    }
    vpotGenePanel(vpotPrioritise.out)
}
