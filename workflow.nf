#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/* Script to filter and prioritise variants from WGS/WES data */

params.vpot_params_abs= file(params.vpot_params).toAbsolutePath().toString()

vcfFile = file(params.vcf)
if( !vcfFile.exists() ) {
    exit 1, "The specified VCF file does not exist: ${params.vcf}"
}
sampleName = vcfFile.baseName

genePanelFile = file(params.genepanel)
if( !genePanelFile.exists() ) {
    exit 1, "The specified gene panel file does not exist: $params.genepanel"
}

if (params.column_file != "None") {
    columnFile = file(params.column_file)
} else {
    columnFile = "None"
}

annotatedVcf = file("${params.output_dir}/${sampleName}.${params.buildver}_multianno.vcf")

process splitVariants {
    debug true
    /*
     Split variants with multiple alleles into separate lines.
    */
    input:
    path vcf

    output:
    path "${sampleName}_split.vcf"

    shell:
    """
    bcftools norm -m-both -o ${sampleName}_split.vcf ${vcf}
    """
}

process annotateGene {
    debug true
    /*
     Function to annotate variants using refGene, filtering out variants
     which did not pass quality control.
    */
    input:
    path vcf

    output:
    path "${sampleName}.${params.buildver}_multianno.vcf"

    shell:
    """
    ${params.annovar_dir}/table_annovar.pl ${vcf} ${params.humandb_dir} -buildver ${params.buildver} \
    -out ${sampleName} -remove -protocol refgene -operation g \
    -nastring . --convertarg "--filter 'pass'" -vcfinput
    """
}

process filterByGene {
    /*
     Function to extract the list of genes from the supplied gene panel CSV file
     and use it to select all variants which are associated with those genes.
     The INFO column of the remaining variants can then be cleaned up using sed,
     to remove the ANNOVAR_DATE and ALLELE_END annotations which will be added
     again in the second annotation step.
    */
    input:
    path anno_vcf

    output:
    path "${sampleName}_genefiltered.vcf"

//  publishDir params.output_dir, mode: 'copy', pattern: '{*_genefiltered.vcf}'

    shell:
    """
    cut -d, -f1 ${genePanelFile} | tail -n +2 > genelist.txt
    (grep "^#" ${anno_vcf}; grep -f genelist.txt ${anno_vcf}) | \
    sed -e 's/ANNOVAR_DATE=20[0-9][0-9]-[0-9][0-9]-[0-9][0-9];//' \
        -e 's/;ALLELE_END//' > ${sampleName}_genefiltered.vcf
    """
}

process annotateAll {
    debug true
    /*
     Function to annotate variants using a variety of user-supplied Annovar
     databases.
    */
    input:
    path filtered_vcf

    output:
    path "${sampleName}.${params.buildver}_multianno.vcf"

    publishDir params.output_dir, mode: 'copy', pattern: '{*_multianno.vcf}'

    shell:
    """
    ${params.annovar_dir}/table_annovar.pl ${filtered_vcf} ${params.humandb_dir} -buildver ${params.buildver} \
        -out ${sampleName} -remove ${params.annovar_params} -nastring . -vcfinput
    """
}

process vpotPrioritise {
    debug true
    /*
     Prioritise the annotated VCF using VPOT.
    */
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
    /*
     Output the prioritisation results as a spreadsheet using VPOT-nf's gene
    panel function.
    */
    input:
    path vpol

    output:
    path "${sampleName}_output_genepanels.xlsx"

    publishDir params.output_dir, mode: 'copy', pattern: "${sampleName}_output_genepanels.xlsx"

    shell:
    """
    python ${params.vpot_dir}/VPOT.py genepanelf "${sampleName}_" $vpol $genePanelFile $params.cancer_type $columnFile
    """
}

workflow {
    if( annotatedVcf.exists() ) {
        println "Already annotated VCF!"
        vpotPrioritise(annotatedVcf)
    } else {
        splitVariants(vcfFile)
        annotateGene(splitVariants.out)
        filterByGene(annotateGene.out)
        annotateAll(filterByGene.out)
        vpotPrioritise(annotateAll.out)
    }
    vpotGenePanel(vpotPrioritise.out)
}
