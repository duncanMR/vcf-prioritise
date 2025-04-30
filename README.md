# vcf-prioritise

Nextflow pipeline for annotating, filtering, and prioritising variants from a VCF file, using Annovar (Wang et al., 2010) and VPOT (Ip et al., 2019). Intended for use with Whole Genome Sequencing/Whole Exome Sequencing (WGS/WES) data. The pipeline involves the following main steps:

- Annotate a VCF file with RefGene using Annovar
- Filter out variants which failed QC or are not in a user-defined gene list
- Label filtered VCF files using user-defined annotation databases
- Prioritise variants with VPOT using user-defined parameters
- Output prioritised variants in an Excel spreadsheet for clinical interpretation

## Requirements

- Python 3.6.+ with `numpy`, `xlsxwriter` and `pandas`, as listed in `requirements.txt`.
- Java 11+

## Installation

- [Request a download link for Annovar](https://www.openbioinformatics.org/annovar/annovar_download_form.php) and extract the file into a directory (e.g. `~/bin`):

    ```bash
    tar -xvf annovar.latest.tar.gz -C ~/bin/
    ```

- Install Nextflow as per the [official instructions](https://www.nextflow.io/docs/latest/getstarted.html), ensuring to add the Nextflow executable to a folder in `$PATH`.

- Clone VPOT-nf and this repository into a directory of your choice (e.g. `~/bin`):

    ```bash
    git clone https://github.com/duncanMR/VPOT-nf ~/bin/VPOT-nf
    git clone https://github.com/duncanMR/vcf-prioritise ~/bin/vcf-prioritise
    ```

- Install Annovar annotation databases which you would like to use. For a list of all databases for a certain reference genome (e.g., hg19), run:

    ```bash
    annotate_variation.pl -webfrom annovar -downdb avdblist -buildver hg19 humandb/
    less ~/bin/annovar/humandb/hg19_avdblist.txt # Change path as appropriate
    ```

    For the toy example below, we will need `refGene` and `avsnp147`:

    ```bash
    annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20220320 ~/bin/annovar/humandb/
    annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene ~/bin/annovar/humandb/
    ```

## Configuration options

The parameters of the pipeline must be set in a `.config` file before usage. Examples for both hg19 and hg38 reference sequences are provided in the `example_configs` folder, which have been used in practice. Here are explanations for each option:

### Primary inputs

| Argument         | Description                                                                 |
|------------------|-----------------------------------------------------------------------------|
| `vcf`            | Unannotated VCF file to be analysed                                        |
| `gene_panels`    | CSV file with genes and their assignments to panels                        |
| `panel_name`     | Name of a panel (column) in the provided CSV to be used for prioritisation |
| `vpot_params`    | VPOT parameter file location                                              |
| `annovar_params` | Argument to be passed to ANNOVAR which specifies which annotation databases to use |
| `ref_name`       | Name of the reference genome used for alignment (e.g., hg19, hg38)        |

### Directories

| Argument        | Description                                                                 |
|-----------------|-----------------------------------------------------------------------------|
| `output_dir`    | Directory in which to output results                                       |
| `annovar_dir`   | Directory where Annovar is installed                                       |
| `humandb_dir`   | Directory where Annovar's databases are installed                          |
| `vpot_dir`      | Directory where VPOT-nf is installed                                       |

### Optional inputs

| Argument         | Description                                                                                   | 
|------------------|-----------------------------------------------------------------------------------------------|
| `normalise`      | Enables normalisation of the VCF with BCF tools before annotating.                            |
| `ref_fasta`      | Reference genome FASTA file: required only when `--normalise` is enabled.                     |
| `column_file`    | Optional CSV file specifying which headings to use in the Excel file and their order. Ignored if "None". |

### Niche options

These shouldn't be changed from the recommended defaults below without good reason.

| Argument             | Description                                                                                   | Recommended default |
|----------------------|-----------------------------------------------------------------------------------------------|---------------------|
| `intermediate_files` | True/false option for whether to export all intermediate files or not.                        | `false`             |
| `pre_annotated`      | Indicates whether the VCF has already been annotated before analysis.                         | `false`             |
| `enable_pass_filter` | If true, enables the filter that excludes variants without a PASS on variant call quality.    | `true`              |

## Usage example

A trivial example is included in the `test_data` folder. An extensive example that has been used in practice is also available in the `example_configs` folder. Run the following commands:

```bash
$ cd ~/bin/vcf-prioritise
$ nextflow workflow.nf -c test_data/test.config
```

If successful, you should find `test.hg19_multianno.vcf` and `test_output_genepanels.xlsx` in the results folder.

## Further notes

### Normalisation

The pipeline can optionally normalise and left-align indels in the input VCF if the `normalise` option is enabled. This requires that BCFtools has been installed and is in `PATH`. If so, the reference FASTA file that was used to generate the VCF must be provided.

### Gene panels

The panelset CSV should contain a column for genes and at least one other column for which panels each gene belongs to. You may wish to have multiple configurations of gene panels: for example, in my research group, we have one gene panel set for colon cancer and another for breast cancer. As you can see in the example table below, BRCA1 is in the BreastHighRisk panel for breast cancer and in the MiscCPG (Cancer-Predisposing Genes) panel for colon cancer. If a gene is not assigned a panel, it is automatically moved to the MiscCPG panel.

| Gene   | ColonPanel    | BreastPanel    |
|--------|---------------|----------------|
| TP53   | ColonHighRisk | BreastHighRisk |
| BRCA1  | MiscCPG       | BreastHighRisk |
| RAD51C | MiscCPG       | MiscCPG        |
| CYP2C9 | WarfarinPGX   | WarfarinPGX    |

An extensive gene panel configuration is available in `example_configs/gene_panels.csv` that has been applied to patient data.

### Column ordering

Since Annovar can add many annotations, it can be difficult to find pertinent columns in the output spreadsheet of the pipeline. To remedy this, the `column_file` option can be provided, which is a list of all the columns in the order you want them to appear in the results. Any columns not on the list will be appended in arbitrary order.

## Acknowledgements

I would like to thank my supervisors, Prof Maritha Kotze, Prof Gerard Tromp, and Prof Craig Kinnear for their valuable input and Dr Brigitte Glanzmann for helpful advice.

## References

- Ip, E., Chapman, G., Winlaw, D., Dunwoodie, S.L., Giannoulatou, E., 2019. VPOT: A Customizable Variant Prioritization Ordering Tool for Annotated Variants. Genomics Proteomics Bioinformatics 17, 540–545. https://doi.org/10.1016/j.gpb.2019.11.001
- Wang, K., Li, M., Hakonarson, H., 2010. ANNOVAR: functional annotation of genetic variants from high-throughput sequencing data. Nucleic Acids Research 38, e164–e164. https://doi.org/10.1093/nar/gkq603
