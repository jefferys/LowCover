stats  <- function( targets, regions, chrY= NA ) {
	if (! is.na(chrY)) {
		targets <- targets[ GenomicRanges::seqnames(targets) != chrY ]
		regions <- regions[ GenomicRanges::seqnames(regions) != chrY ]
	}
	reducedTargets <- GenomicRanges::reduce( targets )
	reducedRegions <- GenomicRanges::reduce( regions )
	bases= sum(GenomicRanges::width(reducedTargets))
	lowBases= sum(GenomicRanges::width(reducedRegions))
	segments = length(reducedTargets)
	lowSegments = length(reducedRegions)
	genes= length(unique( targets$geneName ))
	lowGenes= length(unique( regions$geneName ))
	
	list(
		bases= bases,
		lowBases= lowBases,
		coveredBaseRat=  (bases - lowBases)/ bases,
		segments= segments,
		lowSegments=  lowSegments,
		coveredSegmentRat = (segments - lowSegments)  / segments,
		genes= genes,
		lowGenes= lowGenes,
		coveredGeneRat = (genes - lowGenes) / genes
	)
}

#' Run LowCover as an application
#'
#' @param args Character vector of command line parameters, by default will
#' get them from the command line used to start the R session in which this
#' runs in the usual way.
#'
#' @return N/A - run as an app that saves everything to files
#'
#' @export
LowCoverApp <- function( args= commandArgs( trailingOnly= TRUE )) {
	opts <- parseCLI( args )
	if ( length(opts$help) > 0) {
		cat( opts$help )
		return()
	}
	if ( length(opts$version) > 0) {
		cat( paste0( "version: ", opts$version, "\n" ))
		return()
	}
	
	# Set parallel plan,
	if (opts$parallel == "sequential") {
		oldPlan <- future::plan( "sequential" )
		on.exit( future::plan( oldPlan ), add= TRUE )
	}
	else if (opts$parallel != "default") {
		oldPlan <- future::plan( opts$parallel, workers= opts$workers )
		on.exit( future::plan( oldPlan ), add= TRUE )
	}
	
	targetGR <- loadBed( opts$targetsFile, colnames = c( "seqname", "start", "end", "geneName" ))
	coverageGR <- loadBed( opts$coverageFile, colnames = c( "seqname", "start", "end", "coverTag" ))
	select <- coverageGR$coverTag %in% opts$keep
	LowCoverageGR <- coverageGR[ select ]
	
	regionsGR <- groupedIntersect(targetGR, LowCoverageGR, "geneName" )
	
	saveBed( regionsGR, opts$regionsFile, strand= FALSE )
	badGenes <- unique( regionsGR$geneName )
	if (is.null(badGenes) || (length(badGenes) == 1 && is.na(badGenes)) || length(badGenes) < 1) {
		file.create( opts$badGenesFile )
	}
	else {
		writeLines( badGenes, opts$badGenesFile )
	}
	goodGenes <- setdiff( unique( targetGR$geneName), badGenes )
	if (is.null(goodGenes) || (length(goodGenes) == 1 && is.na(goodGenes)) || length(goodGenes) < 1) {
		file.create( opts$goodGenesFile )
	}
	else {
		writeLines( goodGenes, opts$goodGenesFile )
	}
	
	stats <- stats( targetGR, regionsGR )
	writeLines( paste( names(stats), stats, sep="\t" ),  opts$summaryFile )
	if ( ! is.na(opts$chrY)) {
		stats <- stats( targetGR, regionsGR, chrY= opts$chrY )
		writeLines( paste( names(stats), stats, sep="\t" ),  opts$summaryFileNoY )
	}
	invisible(NULL)
}
