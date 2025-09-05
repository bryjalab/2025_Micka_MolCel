library("idpr")
library("dplyr")

### Define functions
## custom hendersonHasselbach equation - we added X (for phosphoserine) and Z (for phosphothreonine)
phospho_hasselbalch <- function(
    pKa,
    pH = 7.0,
    residue) {
  
  validAcidicResidues <- c("D", "E", "C", "Y",
                           "acid", "negative",
                           "COOH", "COO")
  validBasicResidues <-  c("H", "K", "R",
                           "base", "positive",
                           "NH2", "NH3")
  validPhosphoResidues <- c("X","Z")
  #---- Validating Input
  if (!residue %in% c(validAcidicResidues, validBasicResidues, validPhosphoResidues)) {
    stop("Please set residue equal to an accepted value.")
  }
  if (is.numeric(pKa) == FALSE || is.numeric(pH) == FALSE) {
    stop("pKa and pH each require one numeric value")
  }
  #------ Calculations depending on residue type
  if (residue %in% validAcidicResidues) {
    charge <- -1 / (1 + 10 ^ (pKa - pH))
    return(charge)
  }
  if (residue %in% validBasicResidues) {
    charge <- 1 / (1 + 10 ^ (pH - pKa))
    return(charge)
  }
  if (residue %in% validPhosphoResidues) {
    charge <- -2 / (1 + 10 ^ (pKa - pH))
    return(charge)
  }
}
environment(phospho_hasselbalch) <- asNamespace('idpr')
assignInNamespace("hendersonHasselbalch", phospho_hasselbalch, ns = "idpr")

## custom protein sequence - we added X (for phosphoserine) and Z (for phosphothreonine)
mycheck <- function(
    sequence,
    method = "stop",
    outputType = "string",
    nonstandardResidues = NA,
    suppressAAWarning = FALSE,
    suppressOutputMessage = FALSE) {
  if (!all(is.character(outputType), is.character(method))) {
    stop("Error: method and outputType must be character vectors,")
  }
  if (!any(is.character(sequence), 
           (is(sequence)[1]  %in% c("AAString", "BString", 
                                    "AAStringSet", "BStringSet"))
  )){
    stop("Error: sequence must be a character vector or an AAString Object")
  }
  if (!(method %in% c("stop", "warn"))) {
    stop('Error: method is not equal to a valid term.
            Set method equal to "stop" or "warn"')
  }
  #-----
  #This section will confirm what to do with the amino acid sequence
  if(is(sequence)[1] %in% c("AAString", "BString", 
                            "AAStringSet", "BStringSet")){
    sequence <- as.character(sequence)
  }
  if (length(sequence) == 1) {
    #this is to see if the string is a .fasta / .fa file
    if (grepl("\\.fa", sequence, ignore.case = TRUE)) {
      sequence <- Biostrings::readAAStringSet(sequence, format="fasta")
      sequence <- as.character(sequence)
    }
    separatedSequence <- strsplit(sequence, "")
    names(separatedSequence) <- NULL
    separatedSequence <- unlist(separatedSequence)
  } else {
    separatedSequence <- sequence
  }
  #----- Test for valid residues
  aa <- "ACDEFGHIKLMNPQRSTVWYXZ"
  aa <- strsplit(aa, "")
  aa <- unlist(aa)
  if (!is.na(nonstandardResidues)) {
    aa <- c(aa, nonstandardResidues)
    if (!suppressAAWarning) {
      warningMessage <- paste(
        "This validation allows the following non-standard amino acids: ",
        nonstandardResidues,
        ". If this is an error, please set nonstandardResidues = NA . ",
        "If this is not an error, please set suppressAAWarning = T. ",
        sep = "")
      warning(warningMessage)
    }
  }
  #-----
  #This section checks if the amino acid sequence contains invalid residues
  aaError <- FALSE
  #Used for returning messages later. Set to False unless there is an error
  if (all(separatedSequence %in% aa) == FALSE) {
    aaError <- TRUE #Used for returning messages later reporting an error
    invalidResidues <- separatedSequence[!(separatedSequence %in% aa)]
    invalidResidues <- unique(invalidResidues)
    warningMessage <- paste(
      "Protein contains the following invalid residues: ",
      invalidResidues, ". ", sep = "", collapse = "")
    #--- below reports the error
    if (method == "stop") {
      stop(warningMessage)
    } else {
      warning(warningMessage)
    }
    #--- below reports the error
  }
  #------
  #this section creates the output
  if (outputType == "string") {
    outputSequence <- paste(separatedSequence, sep = "", collapse = "")
  }
  if (outputType == "vector") {
    outputSequence <- separatedSequence
  }
  if (suppressOutputMessage == FALSE) {
    if (aaError == FALSE) {
      validMessage <- paste("The sequence contains no invalid residues.")
    }
    if (aaError == TRUE) {
      validMessage <-
        paste("INVALID SEQUENCE! There are invalid residues.")
    }
    message(validMessage)
  }
  if (!outputType == "none") {
    return(outputSequence)
  }
}
environment(mycheck) <- asNamespace('idpr')
assignInNamespace("sequenceCheck", mycheck, ns = "idpr")

## custom pKa - added pKa 5.6 for phosphoserine, 5.9 for phosphothreonine
pka <-pKaData %>%   
  select("AA", "IPC_peptide") %>% 
  filter(AA != "citation")
custom <- data.frame(AA = c("X", "Z"),
                     IPC_peptide = as.character(c(5.6, 5.9)))

custom_pka <- pka %>% 
  rows_insert(custom)
colnames(custom_pka) <- c("AA", "pKa")
custom_pka$pKa <- as.numeric(custom_pka$pKa)
custom_pka$AA <- as.character(custom_pka$AA)

