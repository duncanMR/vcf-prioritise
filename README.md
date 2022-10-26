# vcf-prioritise
Nextflow pipeline for annotating, filtering and prioritising variants from a VCF file, using Annovar (Wang et al., 2010) and VPOT (Ip et al., 2019). Intended for use with Whole Genome Sequencing/Whole Exome Sequencing (WGS/WES) data. The pipeline involves the following steps:

- Annotate a VCF file with RefGene, using Annovar
- Filter out variants which failed QC or are not in a user-defined gene list
- Label filtered VCF files using user-defined annotation databases
- Prioritise variants with VPOT using user-defined parameters
- Output prioritised variants in an Excel spreadsheet for clinical interpretation

## Requirements

- Python 3.6.+ with modules `NumPy` and `xslxwriter`
- Java 11+

## Installation

- [Request a download link for Annovar](https://www.openbioinformatics.org/annovar/annovar_download_form.php) and extract the file into a directory (e.g. `~/bin`)
```
tar -xvf annovar.latest.tar.gz -C ~/bin/
```
- Install Nextflow as per the [official instructions](https://www.nextflow.io/docs/latest/getstarted.html), ensuring to add the nextflow executable to a folder in `$PATH`.
- Clone VPOT-nf and this repository into a directory of your choice (e.g. `~/bin`)
```
git clone https://github.com/duncanMR/VPOT-nf ~/bin/VPOT-nf
git clone https://github.com/duncanMR/vcf-prioritise ~/bin/vcf-prioritise
```
- Add Annovar to `$PATH` (replace .bashrc with .zshrc if you use zsh as your shell)
``` 
$echo "export PATH=$PATH:~/bin/annovar" >> ~/.bashrc
$exec "$SHELL" #restart shell
$which table_annovar.pl  #check that Annovar's executable is in $PATH
/home/duncan/bin/annovar/table_annovar.pl
```
- Install Annovar annotation databases which you would like to use; for a list of all databases for a certain reference genome (e.g. hg19), run

```
annotate_variation.pl -webfrom annovar -downdb avdblist -buildver hg19 humandb/
less ~/bin/annovar/humandb/hg19_avdblist.txt #change path as appropriate
```
For the toy example below, we will need `refGene` and `avsnp147`:
```
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20220320 ~/bin/annovar/humandb/
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene ~/bin/annovar/humandb/
```
## Required inputs
| Argument      | Description                                                                   |
|:--------------|:------------------------------------------------------------------------------|
| --vcf         | Unannotated VCF file to be analysed                                           |
| --panelset   | CSV file with the list of genes and associated panels                         |

## Optional inputs
| Argument             | Description                                                                                                  | Default                                             |
| --output_dir         | Directory in which to output results                                                                         | "results"                                           |
| --humandb_dir        | Directory in which annovar databases are stored                                                              | "~/bin/annovar/humandb"                             |
| --vpot_dir           | Directory in which VPOT-nf is installed                                                                      | "~/bin/VPOT-nf"                                     |
| --buildver           | Reference genome used for alignment                                                                          | "hg19"                                              |
| --annovar_params     | Argument to be passed to ANNOVAR which specifies which annotation databases to use                           | Default annovar parameters (see below)              |
| --vpot_params        | VPOT parameter file location                                                                                 | "${params.vpot_dir}/default_params/default_ppf.txt" |
| --column_file        | Optional CSV file which specifies which headings to use in the excel file and their order. Ignored if "None" | "None"                                              |
| --intermediate_files | True/false option for whether to export all intermediate files or not                                        | false                                               |

The default Annovar parameters are `-protocol avsnp147,1000g2015aug_all,clinvar_20220320,dbnsfp42a,gnomad_exome,gerp++gt2,caddgt10 -operation f,f,f,f,f,f,f`. Note that all the mentioned annotation databases must be available in your humandb directory for the pipeline to work.

## Usage

The panelset CSV should contain a column for genes and at least one other column for which panels each gene belongs to. You may wish to have multiple configurations of gene panels: for example, in my research group we have one gene panel set for colon cancer and another for breast cancer. As you can see in the example table below, BRCA1 is in the BreastHighRisk panel for breast cancer and in the MiscCPG (Cancer-Predisposing Genes) panel for colon cancer.

| Gene   | ColonPanel    | BreastPanel    |
|--------|---------------|----------------|
| TP53   | ColonHighRisk | BreastHighRisk |
| BRCA1  | MiscCPG       | BreastHighRisk |
| RAD51C | MiscCPG       | MiscCPG        |
| CYP2C9 | WarfarinPGX   | WarfarinPGX    |

A trivial example is included in the `test_data` folder. Run the following commands:

```
$ cd ~/bin/vcf-prioritise
$ nextflow workflow.nf --vcf test_data/test.vcf --genepanel test_data/gene_list.csv --annovar_params "-protocol clinvar_20220320 -operation f"
```
If successful, you should find  `test.hg19_multianno.vcf` and `test_output_genepanels.xlsx` in the results folder.

## Acknowledgements

I would like to thank my supervisors, Prof Maritha Kotze, Prof Gerard Tromp and Prof Craig Kinnear for their valuable input and Dr Brigitte Glanzmann for helpful advice.

## References
- Ip, E., Chapman, G., Winlaw, D., Dunwoodie, S.L., Giannoulatou, E., 2019. VPOT: A Customizable Variant Prioritization Ordering Tool for Annotated Variants. Genomics Proteomics Bioinformatics 17, 540–545. https://doi.org/10.1016/j.gpb.2019.11.001
- Wang, K., Li, M., Hakonarson, H., 2010. ANNOVAR: functional annotation of genetic variants from high-throughput sequencing data. Nucleic Acids Research 38, e164–e164. https://doi.org/10.1093/nar/gkq603
