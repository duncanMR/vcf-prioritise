#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.output_folder = "results"
params.humandb_folder = "/mnt/d/annovar/humandb"
params.buildver = "hg19"

vcfFile = file(params.vcf)
if( !vcfFile.exists() ) {
    exit 1, "The specified VCF file does not exist: ${params.vcf}"
}
sampleName = vcfFile.baseName
println sampleName

process annotateGenes {
    input:
    path vcf

    output:
    path "${sampleName}.${params.buildver}_multianno.vcf"

    shell:
    """
    table_annovar.pl ${vcf} ${params.humandb_folder} -buildver ${params.buildver} \
    -out ${sampleName} -remove -protocol refgene -operation g \
    -nastring . --convertarg "--filter 'pass'" -vcfinput 
    """
}

workflow {
    annotateGenes(vcfFile)
}
