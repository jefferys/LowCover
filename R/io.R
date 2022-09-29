
#' Load a bed-like file as a GRanges object
#'
#' Loads a bed-like file of space or tab-delimited genomic location data into a
#' genomic ranges object. Columns are identified based on provided column names.
#' Start positions are assumed to be 0 based in the file. Genomic ranges objects
#' use 1-based positions.
#' 
#' @param file The bed file to load
#' @param colnames Character vector that specifies how to map loaded data to a
#'   gr object. Only the specified columns (in order) will be loaded; any
#'   following unspecified columns will be dropped, as will any columns whose
#'   names are NA. Three column are required "seqnames", "start", and
#'   "end". The default is for these to be the first three columns. If "strand"
#'   is named, this column will be used as the strand, otherwise the strand for
#'   all loci will be "*". If other columns are named, those will be used as
#'   metadata columns, w
#' @param oneBased By default, input start values are assumed to be zero-based,
#' as per the bed specification. Set `TRUE` if input start values are one-based.
#' @param ... Extra parameters for read.table when reading the file, e.g.
#' setting `skip= 1` or `header= TRUE` to skip loading headings (file headings
#' are never used; `colnames` must always be specified.)
#'
#' @return GenomicRanges::GRanges object
#'
#' @examples
#' # Sample bed detail file
#' bed6File <- system.file( "extdata", "bed6.bed", package = "LowCover" )
#' readLines(bed6File)
#' 
#' # Load location information only, use a (default) strand of '*'
#' gr <- loadBed( bed6File )
#' gr
#' 
#' # Load all data, including strand and two meta-data columns, 'score' and 'id'
#' mapCols <- c("seqname", "start", "end", "strand", "score", "id")
#' gr <- loadBed( bed6File, colnames= mapCols)
#' gr
#' 
#' # Load required columns + score only. Use (default) strand of '*'
#' mapCols <- c("seqname", "start", "end", NA, "score")
#' gr <- loadBed( bed6File, colnames= mapCols)
#' gr
#' 
#' # Load bed-like file with header and 1 based location column info
#' bed6LikeFile <- system.file( "extdata", "bed6like.txt", package = "LowCover" )
#' readLines(bed6LikeFile)
#' 
#' mapCols <- c("seqname", "start", "end", "strand", "score", "id")
#' gr <- loadBed( bed6File, colnames= mapCols, oneBased= TRUE, skip= 1)
#' gr
#' 
#' @export
loadBed <- function( file, colnames= c("seqname", "start", "end"), oneBased= FALSE, ... ) {
	if (! methods::is(file, "connection") && ! file.exists(file)) {
		stop("Error: Can't see the file '", file, "'.", call.= FALSE)
	}
	bedDat <- utils::read.table( file, row.names= NULL, ... )
	colnames(bedDat) <- colnames
	bedDat <- bedDat[ , ! is.na(colnames(bedDat)), drop=FALSE ]
	GenomicRanges::makeGRangesFromDataFrame( bedDat, keep.extra.columns= TRUE, starts.in.df.are.0based= ! oneBased )
}

#' Save a Granges object as a bed file.
#' 
#' Saves a Granges object to a tab-delimited file with "seqnames", "start",
#' "end", and "strand" as the first four columns, followed by all metadata
#' columns. Will convert to 0 based starts and will not include column names
#' unless parameters indicate otherwise. Strand can be excluded if `strand` is
#' set false. To save compressed output, use a connection like gzcon instead of
#' a file. You will have to close it yourself after use.
#'
#' @param x GenomicRanges::Granges object to write to a bed file
#' @param file Filename or connection to write to.
#' @param strand By default, writes a strand column. Set this TRUE to skip.
#' @param oneBased By default, writes 0-based start locations. Set this TRUE to
#'   use one-based starts
#' @param col.names By default, writes no header. Pass a vector of names or set
#'   this TRUE to use the default "as.data.frame" GRanges column names.
#' @param ... Additional parameters to write.table
#'
#' @return The result of calling write.table(), which is undocumented.
#'
#' @examples
#' gr <- GenomicRanges::GRanges(
#'   seqnames= c( "chr1", "chr1", "chr1", "chr2", "chr10" ),
#'   ranges <- IRanges::IRanges(
#'     start=    c( 1L, 11L, 11L, 1L, 16L ),
#'     end=      c( 10L, 20L, 20L, 249000000L, 15L )),
#'     strand=   c( "*", "+", "+", "*", "-" ),
#'     score=    c( 1.1, 2.2, 3.3, 4.4, 5.5 ),
#'     id=       c( "bob", "sue", "tom", "april", "alice" ))
#' 
#' bedFile <- tempfile( fileext= ".bed" )
#' saveBed( gr, bedFile )
#' readLines( bedFile )
#' 
#' bed3File <- tempfile( fileext= ".bed" )
#' grMin <- gr
#' GenomicRanges::mcols( gr ) <- NULL
#' saveBed( grMin, bed3File, strand= FALSE )
#' readLines( bed3File )
#' 
#' tsvFile <- tempfile( fileext= ".tsv" )
#' saveBed( gr, tsvFile, oneBased=TRUE, col.names = TRUE )
#' 
#' @export
saveBed <- function( x, file, strand= TRUE, oneBased= FALSE, col.names= FALSE, ... ) {
	df <- as.data.frame( x )
	df$width <- NULL
	if ( ! oneBased ) {
		df$start <- df$start - 1
	}
	if ( ! strand ) {
		df$strand <- NULL
	}
	utils::write.table( df, file, quote=FALSE, sep= "\t", row.names=FALSE, col.names= col.names, ... )
}
