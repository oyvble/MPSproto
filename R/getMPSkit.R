#' @title getMPSkit
#' @description Function to get kit information
#' @details Returns kit information. Modified by Oyvind Bleka
#' @param kit Short name of kit: See supported kits with getKit()
#' @param fileName Name of kit file used to extract kitinformation
#' @param folderName Name of folder containing the kit file 'fileName'. Using package folder if not specified.
#' @return res A data frame with kit information
#' @export

getMPSkit <- function(kit=NULL, fileName = "kit.txt", folderName=NULL) {  
  .separator <- .Platform$file.sep # Platform dependent path separator. 
  
  if(is.null(folderName)) {
    packagePath <- path.package("MPSproto", quiet = FALSE) # Get package path.
    folderName <- paste(packagePath,"extdata",sep=.separator) #get folder containing the filename
  }
  filePath <- paste(folderName, fileName, sep=.separator) #get full pathname of kit file
  .kitInfo <- read.delim(file=filePath, header = TRUE, sep = "\t", quote = "\"",dec = ".", fill = TRUE, stringsAsFactors=FALSE)
 
  # Available kits. Must match else if construct.
  kits<-unique(.kitInfo$Short.Name)
  
	if (is.null(kit)) {	# Check if NULL
		res <- kits
	} else {	# String provided.
		# Check if number or string.
		if (is.numeric(kit)) {
			index<-kit # Set index to number.
		} else {
			index<-match(toupper(kit),toupper(kits)) # Find matching kit index (case insensitive)
		}
		
		if (any(is.na(index))) { 		# No matching kit.
			return(NA)
		# Assign matching kit information.
		} else {
		  currentKit <- .kitInfo[.kitInfo$Short.Name==kits[index], ]
              res <- data.frame(
                        Marker = currentKit$Marker,
                        Allele = currentKit$Allele,
                        Size = currentKit$Size,
                        stringsAsFactors = FALSE)
		  res$Marker <- factor(res$Marker, levels=unique(res$Marker))# Create useful factors. 
		} 
	}
          
  return(res)  
} #end function

