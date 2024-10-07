suppressPackageStartupMessages(library(rmarkdown))
suppressPackageStartupMessages(library(argparse))


html_automator <- function(template, report_name, output_directory, output_file, cutoff_lfc, cutoff_p, interaction, sample_info_file, data_source) {
  
  rmarkdown::render(
    input = template,
    output_file = output_file,
    output_format = html_document(),
    output_dir = output_directory,
    params = list(
      show_code = FALSE,
      cutoff_p = cutoff_p,
      cutoff_lfc = cutoff_lfc,
      interaction = interaction,
      sample_info_file = sample_info_file,
      data_source = data_source,
      output_dir = output_directory
      
    )
    
  )
  
}


parser <- ArgumentParser(description = "From Gene Level Counts, create a differential expression analysis html report")

parser$add_argument("--template", help= "RMD template to use for creating html report")
parser$add_argument("--report_name", help= "Title of the final html report")
parser$add_argument("--output_directory", help = "Output Directory")
parser$add_argument("--output_file", help = "Name of the final html report file")
parser$add_argument("--lfc_cutoff", help = "Program makes calculations and creates charts with respect to this log fold change value")
parser$add_argument("--p_cutoff", help = "Cutoff p value - i.e 0.05, or 0.01")
parser$add_argument("--interaction", help = "Find out interaction between the factors if design is multifactorial (TRUE or FALSE)")
parser$add_argument("--sample_info_file", help = "File containing sample information (names, treatments and directories)")
parser$add_argument("--data_source", help = "Data source (Can be either ncbi or ensembl)")

args <- parser$parse_args()


html_automator(args$template, args$report_name, args$output_directory, args$output_file, 
               args$lfc_cutoff, args$p_cutoff, args$interaction, args$sample_info_file, 
               args$data_source)



