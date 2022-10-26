#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/* Script to filter and prioritise variants from WGS/WES data */

params.output_dir = "results"
params.humandb_dir = "~/bin/annovar/humandb"
params.vpot_dir = "~/bin/VPOT-nf"
params.buildver = "hg19"
params.annovar_params = " -protocol avsnp147,clinvar_20220320,dbnsfp42a,gnomad_exome,1000g2015aug_all,gerp++gt2,fathmm,dann,eigen,caddgt10\
    -operation f,f,f,f,f,f,f,f,f,f"
params.vpot_params = "${params.vpot_dir}/default_params/default_ppf.txt"
params.cancer_type = "BreastCancer"
params.intermediate_files = false
vcfFile = file(params.vcf)
if( !vcfFile.exists() ) {
    exit 1, "The specified VCF file does not exist: ${params.vcf}"
}
sampleName = vcfFile.baseName

genePanelFile = file(params.genepanel)
if( !genePanelFile.exists() ) {
    exit 1, "The specified gene panel file does not exist: $params.genepanel"
}

process annotateGene {
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
    table_annovar.pl ${vcf} ${params.humandb_dir} -buildver ${params.buildver} \
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

    if (params.intermediate_files == true) {
        publishDir params.output_dir, mode: 'copy', pattern: '{*_genefiltered.vcf}'
    }

    shell:
    """
    cut -d, -f1 ${genePanelFile} | tail -n +2 > genelist.txt
    (grep "^#" ${anno_vcf}; grep -f genelist.txt ${anno_vcf}) | \
    sed -e 's/ANNOVAR_DATE=20[0-9][0-9]-[0-9][0-9]-[0-9][0-9];//' \
        -e 's/;ALLELE_END//' > ${sampleName}_genefiltered.vcf
    """
}

process annotateAll {
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
    table_annovar.pl ${filtered_vcf} ${params.humandb_dir} -buildver ${params.buildver} \
        -out ${sampleName} -remove ${params.annovar_params} -nastring . -vcfinput
    """
}

process makeSampleList {
    /*
     VPOT requires the input data to be listed in a text file to support batch samples.
     The first column must be the path of the VCF and the second column is the
    name of the sample header in the VCF file. This function makes the sample
    list file for our annotated VCF.
    */
    input:
    path vcf

    output:
    path "vpot_input.txt"

    if (params.intermediate_files == true) {
        publishDir params.output_dir, mode: 'copy', pattern: 'vpot_input.txt'
    }

    shell:
    '''
    echo "!{launchDir}/!{params.output_dir}/!{vcf}	$(grep "#CHROM" !{vcf} | awk '{print $NF}')" > vpot_input.txt
    '''
}

process vpotPrioritise {
    /*
     Prioritise the annotated VCF using VPOT.
    */
    input:
    path vpot_input

    output:
    path "${sampleName}_final_output_file.txt"

    if (params.intermediate_files == true) {
        publishDir params.output_dir, mode: 'copy', pattern: '*_final_output_file.txt'
    }

    shell:
    """
    python ${params.vpot_dir}/VPOT.py priority ${sampleName}_ \
        $vpot_input $params.vpot_params
    """
}

process vpotGenePanel {
    /*
     Output the prioritisation results as a spreadsheet using VPOT-nf's gene
    panel function.
    */
    input:
    path vpol

    output:
    path "${sampleName}_output_genepanels.xlsx"

    publishDir params.output_dir, mode: 'copy', pattern: '*_output_genepanels.xlsx'

    shell:
    """
    python ${params.vpot_dir}/VPOT.py genepanelf "${sampleName}_" $vpol $genePanelFile $params.cancer_type
    """
}

workflow {
    annotateGene(vcfFile)
    filterByGene(annotateGene.out)
    annotateAll(filterByGene.out)
    makeSampleList(annotateAll.out)
    vpotPrioritise(makeSampleList.out)
    vpotGenePanel(vpotPrioritise.out)
}
