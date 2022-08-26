#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/* Script to filter and prioritise variants from WGS/WES data */

params.output_dir = "results"
params.humandb_dir = "/mnt/d/annovar/humandb"
params.buildver = "hg19"
params.annovar_params = " -protocol avsnp147,dbnsfp42a,gnomad_exome,1000g2015aug_all,gerp++gt2,fathmm,dann,eigen,caddgt10 \
    -operation f,f,f,f,f,f,f,f,f"

vcfFile = file(params.vcf)
if( !vcfFile.exists() ) {
    exit 1, "The specified VCF file does not exist: ${params.vcf}"
}
sampleName = vcfFile.baseName

genePanelFile = file(params.genepanel)
if( !genePanelFile.exists() ) {
    exit 1, "The specified gene panel file does not exist: ${params.genepanel}"
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

    publishDir params.output_dir, mode: 'copy', pattern: '{*_genefiltered.vcf}'

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

workflow {
    annotateGene(vcfFile)
    filterByGene(annotateGene.out)
    annotateAll(filterByGene.out)
}
