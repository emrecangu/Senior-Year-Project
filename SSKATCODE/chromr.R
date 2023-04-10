library(vcfR)

vcf <- read.vcfR("e:/bitirme_tezi/maincohort/ob.vcf")

CHROM <- create.chromR <- function(vcf, name="CHROM", seq=NULL, ann=NULL, verbose=TRUE){

  stopifnot( inherits(vcf, "vcfR") )
  
  if( length( unique( getCHROM(vcf) ) ) > 1 ){
    myChroms <- unique( getCHROM(vcf) )
    message('vcfR object includes more than one chromosome (CHROM).')
    message( paste(myChroms, collapse = ", ") )
    message("Subsetting to the first chromosome")
    vcf <- vcf[ getCHROM(vcf) == myChroms[1],]
  }
  
  if( length( names(seq) ) > 1 ){
    mySeqs <- names(seq)
    message('DNAbin object includes more than one chromosome.')
    message( paste(mySeqs, collapse = ", ") )
    message("Subsetting to the first chromosome")
    seq <- seq[ mySeqs[1] ]
  }
  
  if( length( unique( ann[,1] ) ) > 1 ){
    myChroms <- unique( ann[,1] )
    message('Annotations include more than one chromosome.')
    message( paste(myChroms, collapse = ", ") )
    message("Subsetting to the first chromosome")
    ann <- ann[ann[,1] == myChroms[1], , drop = FALSE]
  }

  }
x <- new(Class="chromR")
names(x) <- x
#  setName(x) <- name

# Insert vcf into Chom.
if(length(vcf)>0){
  #    x <- vcf2chromR(x, vcf)
  x@vcf <- vcf
}

# Insert seq into chromR
# Needs to handle lists and matrices of DNAbin
# Matrices are better behaved.
#
if(is.null(seq)){
  POS <- getPOS(x)
  x@len <- POS[length(POS)]
  #    x@len <- x@vcf.fix$POS[length(x@vcf.fix$POS)]
  #} else if (class(seq)=="DNAbin"){
} else if ( inherits(seq, "DNAbin") ){
  x <- seq2chromR(x, seq)
} else {
  #stopifnot( class(seq)=="DNAbin" )
  stopifnot( inherits(seq, "DNAbin") )
}

# Annotations.
if( !is.null(ann) ){
  if( nrow(ann) > 0 ){
    #  if(nrow(ann) > 0){
    #stopifnot(class(ann) == "data.frame")
    stopifnot( inherits(ann, "data.frame") )
    #if(class(ann[,4]) == "factor"){ann[,4] <- as.character(ann[,4])}
    if( inherits(ann[,4], "factor") ){ann[,4] <- as.character(ann[,4])}
    #if(class(ann[,5]) == "factor"){ann[,5] <- as.character(ann[,5])}
    if( inherits(ann[,5], "factor") ){ann[,5] <- as.character(ann[,5])}
    #if(class(ann[,4]) == "character"){ann[,4] <- as.numeric(ann[,4])}
    if( inherits(ann[,4], "character") ){ann[,4] <- as.numeric(ann[,4])}
    #if(class(ann[,5]) == "character"){ann[,5] <- as.numeric(ann[,5])}
    if( inherits(ann[,5], "character") ){ann[,5] <- as.numeric(ann[,5])}
    x@ann <- ann
    
    # Manage length
    if( max(as.integer(as.character(ann[,4]))) > x@len ){
      x@len <- max(as.integer(as.character(ann[,4])))
    }
    if( max(as.integer(as.character(ann[,5]))) > x@len ){
      x@len <- max(as.integer(as.character(ann[,5])))
    }
  }
}

# Report names of objects to user.
if(verbose == TRUE){
  # Print names of elements to see if they match.
  message("Names in vcf:")
  chr_names <- unique(getCHROM(x))
  message(paste('  ', chr_names, sep=""))
  #    message(paste('  ', unique(as.character(x@vcf.fix$CHROM)), sep=""))
  
  #if(class(x@seq) == "DNAbin"){
  if( inherits(x@seq, "DNAbin") ){
    message("Names of sequences:")
    message(paste('  ', unique(labels(x@seq)), sep=""))
    
    #      if(unique(as.character(x@vcf.fix$CHROM)) != unique(labels(x@seq))){
    if(chr_names != unique(labels(x@seq))){
      warning("
        Names in variant data and sequence data do not match perfectly.
        If you choose to proceed, we'll do our best to match the data.
        But prepare yourself for unexpected results.")
      #        message("Names in variant file and sequence file do not match perfectly.")
      #        message("If you choose to proceed, we'll do our best to match data.")
      #        message("But prepare yourself for unexpected results.")
    }
  }
  
  if(nrow(x@ann) > 0){
    message("Names in annotation:")
    message(paste('  ', unique(as.character(x@ann[,1])), sep=""))
    #      if(unique(as.character(x@vcf.fix$CHROM)) != unique(as.character(x@ann[,1]))){
    if( length( unique(as.character(x@ann[,1])) ) > 1 ){
      warning("The annotation data appear to include more than one chromosome.\nUsing only the first.\n")
      firstChrom <- unique(as.character(x@ann[,1]))[1]
      x@ann <- x@ann[ x@ann[,1] == firstChrom, , drop = FALSE]
      myChrom <- unique( x@ann[,1] )
      warning( paste('Using annotation chromosome:', myChrom, '\n') )
    }
    
    if(chr_names != unique(as.character(x@ann[,1]))){
      warning("
        Names in variant data and annotation data do not match perfectly.
        If you choose to proceed, we'll do our best to match the data.
        But prepare yourself for unexpected results.")
      #        message("Names in variant file and annotation file do not match perfectly.\n")
      #        message("If you choose to proceed, we'll do our best to match data.\n")
      #        message("But prepare yourself for unexpected results.\n")
    }
  }
}

# Check to see if annotation positions exceed seq position.
if( nrow(x@ann) > 0 ){
  if( max(as.integer(as.character(x@ann[,4]))) > x@len | max(as.integer(as.character(x@ann[,5]))) > x@len ){
    stop("Annotation positions exceed chromosome positions.  Is this the correct set of annotations?")
  }
}

if( verbose == TRUE ){
  message("Initializing var.info slot.")
}
x@var.info <- data.frame( CHROM = x@vcf@fix[,"CHROM"] , POS = as.integer(x@vcf@fix[,"POS"]) )
#  mq <- getINFO(x, element="MQ")
mq <- extract.info(x, element = 'MQ', as.numeric = TRUE)
if( length(mq) > 0 ){ x@var.info$MQ <- mq }
#  dp <- getDP(x)
dp <- extract.info(x, element = 'DP', as.numeric = TRUE)
if( length(dp) > 0 ){ x@var.info$DP <- dp }
if( nrow(x@var.info) > 0 ){
  x@var.info$mask <- TRUE
}
if( verbose == TRUE ){
  message("var.info slot initialized.")
}
return(x)




##### ##### ##### ##### #####
#
# chromR data loading functions
#
##### ##### ##### ##### #####

#' @rdname create_chromR
#' @export
#' @aliases vcf2chromR
# @aliases chromR-methods vcf2chromR
#'
#'
# @description
# Methods to work with objects of the chromR class
# Reads in a vcf file and stores it in a vcf class.
#'
# @param x an object of class chromR
#'
#'
vcfR2chromR <- function(x, vcf){
  x@vcf.fix <- as.data.frame(vcf@fix)
  #  colnames(x@vcf.fix) <- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO')
  #  x@vcf.fix[,2] <- as.numeric(x@vcf.fix[,2])
  #  x@vcf.fix[,6] <- as.numeric(x@vcf.fix[,6])
  #
  for(i in 1:ncol(vcf@gt)){
    vcf@gt[,i] <- as.character(vcf@gt[,i])
  }
  x@vcf.gt <- vcf@gt
  #
  x@vcf.meta <- vcf@meta
  #
  # Initialize var.info slot
  x@var.info <- data.frame(matrix(ncol=5, nrow=nrow(vcf@fix)))
  names(x@var.info) <- c('CHROM', 'POS', 'mask', 'DP','MQ')
  #  names(x@var.info) <- c('DP','MQ', 'mask')
  #
  x@var.info$CHROM <- x@vcf.fix$CHROM
  x@var.info$POS <- x@vcf.fix$POS
  x@var.info$mask <- rep(TRUE, times=nrow(x@vcf.fix))
  #
  if(length(grep("DP=", vcf@fix[,8])) > 0){
    x@var.info$DP <- unlist(lapply(strsplit(unlist(lapply(strsplit(as.character(vcf@fix[,8]), ";"), function(x){grep("^DP=", x, value=TRUE)})),"="),function(x){as.numeric(x[2])}))
  }
  if(length(grep("MQ=", vcf@fix[,8])) > 0){
    x@var.info$MQ <- unlist(lapply(strsplit(unlist(lapply(strsplit(as.character(vcf@fix[,8]), ";"), function(x){grep("^MQ=", x, value=TRUE)})),"="),function(x){as.numeric(x[2])}))
  }
  #
  # assign may be more efficient.
  return(x)
}


# Needs to handle lists and matrices of DNAbin.
# Matrices appear better behaved.
#
#' @rdname create_chromR
#' @export
#' @aliases seq2chromR
#'
seq2chromR <- function(x, seq=NULL){
  # A DNAbin will store in a list when the fasta contains
  # multiple sequences, but as a matrix when the fasta
  # only contains one sequence.
  if(is.list(seq)){
    #stopifnot(length(seq)==1)
    if( length(seq) != 1 ){
      stop("seq2chromR expects a DNAbin object with only one sequence in it.")
    }
    x@seq <- as.matrix(seq)
    x@len <- length(x@seq)
  } else if (is.matrix(seq)){
    stopifnot(nrow(seq)==1)
    #    x@seq <- ape::as.DNAbin(as.character(seq)[1,])
    #    dimnames(pinf_dna)[[1]][1]
    x@seq <- seq
    x@len <- length(x@seq)
  } else {
    stop("DNAbin is neither a list or matrix")
  }
  return(x)
}


#' @rdname create_chromR
#' @export
#' @aliases ann2chromR
#'
ann2chromR <- function(x, gff){
  x@ann <- as.data.frame(gff)
  colnames(x@ann) <- c('seqid','source','type','start','end','score','strand','phase','attributes')
  x@ann$start <- as.numeric(as.character(x@ann$start))
  x@ann$end   <- as.numeric(as.character(x@ann$end))
  return(x)
}