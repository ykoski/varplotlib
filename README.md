# varplotlib
A collection of Python and R scripts used to combine and plot variant call data.

Author: Yrjö Koski, first version August 2018

Contact information:  
email: **ymkoski[at]gmail.com**  
github: **ykoski**

Varplotlib is collection of Python and R scripts used to give quick and elegant visual presentation of variant call data. 
At the moment varplotlib works best with [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/)-annotated variant
call files (.vcf). Two plot types are available: prevalence matrix, which shows mutations in interesting genes per each sample, 
and base substitution frequency plot, which shows the frequencies of six base substitutions and indels. The prevalence matrix
also shows how many functional variants are called per sample.

Running the varplotlib.py script will produce both of these plots and also .csv files used by the R-plotting scripts. 
You can also run R-Scripts separatedly once you have created the .csv-files (this has probably no real use since the Python runs
somewhat fast).

<h2>1. Setting up varplotlib</h2>

Since varplotlib runs on [Python](http://www.python.org/) and [R](http://www.r-project.org/), both of them are obviously needed. This was written on Python 3.6 and it isn't compatible
with Python 2.

R-libraries that are needed:
* ggplot2
* reshape2
* RColorBrewer
* dplyr
* gtable
* gridExtra

Setting up is simple. Just clone into git-repository with:  
`git clone git@github.com:ykoski/varplotlib.git`

After this you're ready to run the script!

<h2>2. Running varplotlib</h2>

Running varplotlib is simple! The basic call for the script is here:

`python3 varplotlib.py -i "my_vcfs.list" -o "output_name" -gl "interesting_genes.txt"`

<h4>Essential arguments:</h4>  

* -i (or --input) is a list of vcf-files with sample names. Each file can belong to a group (for example healthy control group and group with
disease of interest). List is tab separated file where first column is path to file, second column is sample name and third column is sample group.
* -o (or --output) is string that acts as a prefix for all the output files. This can include a path for example 
/path/to/variant_files/output.
* -gl (or --genelist) is a list of genes of interest. The script will plot variants on these genes. This is where 
ANNOVAR-annotations really come in! The script doesn't work with other annotators. Only ANNOVAR's "Gene.RefGene" annotation is supported.
This list can also have gene-groups for grouping genes on the output plot (separated by tab). So make sure every gene
is on a new row and gene groups are separated from the gene name by tab.

<h4>Optional arguments:</h4>

* --remove_zero_variants will remove samples that have no variants in the genes of interest from the prevalence matrix.
* --rename_samples will rename the samples using the sample group given in the input.list. For example if you don't want to
include private patient codes in the final plot you don't have to fix them in the input.list.
* --cw custom width of the prevalence matrix. If the default settings aren't working you can use this parameter to overwrite 
the default width of the plot.
* --ch custom height of the prevalence matrix.

Please notice that the RScript-folder must be located under varplotlib-folder. Otherwise the python script won't be able to call the R-Scripts.
