lowcoverStats  <- function( targets, regions, chrY= NA ) {
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

lowcoverApp <- function( args= commandArgs( trailingOnly= TRUE )) {
	opts <- parseCLI( args )
	if ( length(opts$help) > 0) {
		cat( opts$help )
		return()
	}
	if ( length(opts$version) > 0) {
		cat( paste0( "version: ", opts$version, "\n" ))
		return()
	}
	
	
	targetGR <- loadBed( opts$targetsFile, colnames = c( "seqname", "start", "end", "geneName" ))
	coverageGR <- loadBed( opts$coverageFile, colnames = c( "seqname", "start", "end", "coverTag" ))
	select <- coverageGR$coverTag %in% opts$keep
	lowCoverageGR <- coverageGR[ select ]
	
	regionsGR <- groupedIntersect(targetGR, lowCoverageGR, "geneName" )
	
	saveBed( regionsGR, opts$regionsFile, strand= FALSE )
	badGenes <- unique( regionsGR$geneName )
	if (is.null(badGenes) || is.na(badGenes) || length(badGenes) < 1) {
		file.create( opts$badGenesFile )
	}
	else {
		writeLines( badGenes, opts$badGenesFile )
	}
	goodGenes <- setdiff( unique( targetGR$geneName), badGenes )
	if (is.null(goodGenes) || is.na(goodGenes) || length(goodGenes) < 1) {
		file.create( opts$goodGenesFile )
	}
	else {
		writeLines( goodGenes, opts$goodGenesFile )
	}
	
	stats <- lowcoverStats( targetGR, regionsGR )
	writeLines( paste( names(stats), stats, sep="\t" ),  opts$summaryFile )
	if ( ! is.na(opts$chrY)) {
		stats <- lowcoverStats( targetGR, regionsGR, chrY= opts$chrY )
		writeLines( paste( names(stats), stats, sep="\t" ),  opts$summaryFileNoY )
	}
}
