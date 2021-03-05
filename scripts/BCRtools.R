####################################################################################################
#
## BCRtools
##
## R function to import and collate output from VDJPuzzle
##
##
##
## January 2020
##
##
## EJC
#
#
####################################################################################################

#' VDJPuzzle importer
#' 
#' 
#' This needs to import the sequence results and a separate 'mutations' file.
#' 
#' @param path Path to the vdjpuzzle output folder. Defaults to: ~/nf_bcr_output/VDJPuzzle/
#' @param pattern End of filenames to list. Defaults to ".tsv". Of note the script, selects the tsvs within the 'final_receptor_results' folder. There is no option to change this.
#' 
#' @return list of tibbles, one for each chain
#' 
#' @examples
#' VDJ.chain.results <- VDJPuzzleimporter()
#' 
#' 
#' @export 
#' 

VDJPuzzleimporter <- function(path = "~/nf_bcr_output/VDJPuzzle/", pattern = "*.tsv$"){
  require(readr)
  require(dplyr)
  
  # Get filenames
  vdj.out <- list.files(path = path, pattern = pattern, recursive = T, full.names = T)
  
  ## Remove any NoIndex files:
  vdj.out <- vdj.out[ ! grepl(vdj.out, pattern = "NoCode|NoIndex")]
  
  ## Select only those with /final_receptor_results
  vdj.out.to.read <- vdj.out[ ! grepl(vdj.out, pattern = "final_receptor_results")]
  
  cat("\nFound ",
  length(vdj.out.to.read [ grepl(vdj.out.to.read, pattern = "IGH")]),
  " heavy chain files to import\n")
  
  cat("\nFound ",
      length(vdj.out.to.read [ grepl(vdj.out.to.read, pattern = "IGL|IGK")]),
      " light chain files to import\n")

 ## Read in each file:
  vdj.heavy.results <- lapply(vdj.out.to.read[ grepl(vdj.out.to.read, pattern = "IGH")], function(x) {
    suppressWarnings(
    read_tsv(file = x,col_types = cols(.default = "c"))
    )
  })
  
  vdj.lambda.results <- lapply(vdj.out.to.read[ grepl(vdj.out.to.read, pattern = "IGL")], function(x) {
    suppressWarnings(
    read_tsv(file = x,col_types = cols(.default = "c"))
    )
  })
  
  vdj.kappa.results <- lapply(vdj.out.to.read[ grepl(vdj.out.to.read, pattern = "IGK")], function(x) {
    suppressWarnings(
    read_tsv(file = x,col_types = cols(.default = "c"))
    )
  })
  
  ## Parse the list of tbls into useful objects
  vdj.heavy.results <- vdj.heavy.results %>%
    # Bind rows
    bind_rows() %>%
    group_by(CellID) %>%
    # Take the highest expressed IgH (by kallisto):
    top_n(1, Expression) %>%
    # There are some cells where kallisto value is 0
    # So cannot separate on expression
    # These seem to be duplicate rows.
    dplyr::slice(1) %>%
    # Sort out a 'label' column, where _L001 is taken off the end, and IGx_ off the start of CellID
    mutate(label = gsub(x = gsub(CellID,pattern = "_L001*$", replacement = ""), pattern = "^IGH_", replacement = ""))
  
  
  vdj.kappa.results <- vdj.kappa.results %>% 
    bind_rows() %>%
    group_by(CellID) %>%
    top_n(1, Expression) %>%
    # There are some cells where kallisto value is 0
    # So cannot separate on expression
    # These seem to be duplicate rows.
    dplyr::slice(1) %>%
    # Sort out a 'label' column, where _L001 is taken off the end, and IGx_ off the start of CellID
    mutate(label = gsub(x = gsub(CellID,pattern = "_L001*$", replacement = ""), pattern = "^IGK_", replacement = ""))  
  
  
  
  vdj.lambda.results <- vdj.lambda.results %>% 
    bind_rows() %>%
    group_by(CellID) %>%
    top_n(1, Expression) %>%
    # There are some cells where kallisto value is 0
    # So cannot separate on expression
    # These seem to be duplicate rows.
    dplyr::slice(1) %>%
    # Sort out a 'label' column, where _L001 is taken off the end, and IGx_ off the start of CellID
    mutate(label = gsub(x = gsub(CellID,pattern = "_L001*$", replacement = ""), pattern = "^IGL_", replacement = ""))
  
  # return these objects as a list: 1=heavy chain, 2=kappa, 3=lambda
  vdj.results <- list(
    as.data.frame(vdj.heavy.results),
    as.data.frame(vdj.kappa.results),
    as.data.frame(vdj.lambda.results))

 return(vdj.results)
}
