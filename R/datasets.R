#' Human-mouse ortholog conversion table
#'
#' A conversion table of stable orthologs between Hs and Mm.
#' 
#' @format A dataframe containing gene ortholog mapping.
#' @source \url{https://www.ensembl.org/Mus_musculus/Info/Index}
'Hs2Mm.convert.table'

#' Cell cycling signatures
#'
#' A list of cell cycling signatures (G1.S and G2.M phases),
#' for mouse and human.
#' 
#' @format A list of cycling signatures.
#' @source \doi{10.1126/science.aad0501}
'cell.cycle.obj'

#' Test dataset for ProjecTILs
#'
#' A small dataset of CD8 T cells, to test the ProjecTILs installation
#' 
#' @format A Seurat object
#' @source \url{https://pmc.ncbi.nlm.nih.gov/articles/PMC6673650/}
'query_example_seurat'