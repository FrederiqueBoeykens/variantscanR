#' @title chromosomenameR
#' @description Changes 'NC_xxx' chromosome nomenclature to 'Chr' nomenclature.
#'
#' @param a VCF file that needs a different chromosome nomenclature
#' @param b Conversion file that contains the original nomenclature in the first column and the desired nomenclature in the second.
#'
#' @return Returns a VCF file with different chromosome nomenclature
#' @export
#'
#' @examples
#' \donttest{}

chromosomenameR <- function(a,b){
  new <- a
  new[] <- lapply(a, function(x) b[,2][match(x, b[,1])])
  a[,1] <- new[,1]
  a[] <- lapply(a, as.character)
  return(a)
}
