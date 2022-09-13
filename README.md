* vcf-prioritise
Nextflow pipeline for annotating, filtering and prioritising variants from a VCF file, using Annovar (Wang et al., 2010) and VPOT (Ip et al., 2019). Intended for use with Whole Genome Sequencing/Whole Exome Sequencing (WGS/WES) data. The pipeline involves the following steps:

- Annotate a VCF file with RefGene, using Annovar
- Filter out variants which failed QC or are not in a user-defined gene list
- Label filtered VCF files using user-defined annotation databases
- Prioritise variants with VPOT using user-defined parameters
- Output prioritised variants in an Excel spreadsheet for clinical interpretation

** Requirements

- Python 3.6.+ and NumPy
- Java 11+

** Installation

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
** Input files

** Parameters

** Usage

** Output

** Acknowledgements

** References
- Ip, E., Chapman, G., Winlaw, D., Dunwoodie, S.L., Giannoulatou, E., 2019. VPOT: A Customizable Variant Prioritization Ordering Tool for Annotated Variants. Genomics Proteomics Bioinformatics 17, 540–545. https://doi.org/10.1016/j.gpb.2019.11.001
- Wang, K., Li, M., Hakonarson, H., 2010. ANNOVAR: functional annotation of genetic variants from high-throughput sequencing data. Nucleic Acids Research 38, e164–e164. https://doi.org/10.1093/nar/gkq603
