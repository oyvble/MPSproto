##### Function to obtain the reverse complement of the original sequence

#' @title getComplement
#'
#' @author Oyvind Bleka <oyvble.at.hotmail.com>
#'
#' @description Outputs the reverse complement sequence
#'
#' @details Helpfunction 1. If marker is positioned in antisense strand, the function
#' obtains the complement function and subsequently reverts into 5' to 3'
#'
#' @param seqs Original sequence
#'
#' @return seqs Reverse complement of the original seqs
#'
#' @export
#'
#' @importFrom stringi stri_reverse
#'
#' @examples
#' \dontrun{
#' getComplement("ATCGATCG")
#' }
getComplement <- function(seqs) {
 if(!all(is.character(seqs))) stop("Please only input character vectors")

 seqs = toupper(seqs) #make sure that letters are upper case
 seqs = gsub("A","t",seqs)
 seqs = gsub("T","a",seqs)
 seqs = gsub("C","g",seqs)
 seqs = gsub("G","c",seqs)
 return( stringi::stri_reverse(toupper(seqs)) )
}


