makeDat <- function( name ) {
	dir <- tempfile( name )
	dir.create(dir)
	list(
		runDir= dir,
		targetFile= file.path( dir, "targets.bed" ),
		coverageFile= file.path( dir, "coverage.bed" ),
		regionsFile= file.path( dir, "lowcover.regions.bed" ),
		badGenesFile= file.path( dir, "lowcover.badgenes.txt" ),
		goodGenesFile= file.path( dir, "lowcover.goodgenes.txt" ),
		summaryFile= file.path( dir, "lowcover.summary.tsv"),
		summaryFileNoY= file.path( dir, "lowcover.summaryNoY.tsv")
	)
}

coverageGR <- loadBed( system.file("extdata", "coverage.bed", package= packageName() ),
					   c("seqnames", "start", "end", "tag"))
targetsGR <- loadBed( system.file("extdata", "targets.bed", package= packageName() ),
					   c("seqnames", "start", "end", "genes"))

describe( "lowcoverApp() test runs", {
	describe( "Complex data, keep Y", {
		dat <- makeDat( "test_lowcover_complex_withY" )
		withr::with_dir( dat$runDir, {
			saveBed( coverageGR, dat$coverageFile, strand= FALSE )
			saveBed( targetsGR, dat$targetFile, strand= FALSE )
			
			cli <- c("-t", dat$targetFile, "--coverageFile", dat$coverageFile, "--chrY", "NA")
			expect_silent(lowcoverApp( cli ))
			
			it( "Creates expected files", {
				expect_true( file.exists( dat$regionsFile ))
				expect_true( file.exists( dat$badGenesFile ))
				expect_true( file.exists( dat$goodGenesFile ))
				expect_true( file.exists( dat$summaryFile ))
				expect_false( file.exists( dat$summaryFileNoY ))
			})
			it( "The regions file contents are correct", {
				got <- readLines( dat$regionsFile )
				want <- c(
					"chr1\t10\t25\tgene_A",
					"chr1\t200\t250\tgene_A",
					"chr1\t200\t250\tgene_B",
					"chr1\t290\t300\tgene_B",
					"chr1\t290\t310\tgene_C",
					"chr2\t10\t45\tgene_D",
					"chrY\t0\t100\tgene_YE"
				)
				expect_equal(got, want)
			})
			it( "The good genes file contents are correct", {
				got <- readLines( dat$goodGenesFile )
				want <- c( "gene_YF" )
				expect_equal(got, want)
			})
			it( "The bad genes file contents are correct", {
				got <- readLines( dat$badGenesFile )
				want <- c( "gene_A", "gene_B", "gene_C", "gene_D", "gene_YE" )
				expect_equal(got, want)
			})
			it( "The summary file contents are correct", {
				got <- read.delim( dat$summaryFile, header = FALSE,)
				wantNames <- c( "bases",    "lowBases",    "coveredBaseRat",
								"segments", "lowSegments", "coveredSegmentRat",
								"genes",    "lowGenes",    "coveredGeneRat" )
				expect_equal( got[,1], wantNames )
				wantValues <- c( 630,220,(630-220)/630,5,5,0,6,5,(6-5)/6 )
				expect_equal( got[,2], wantValues )
			})
		})
	})
	describe( "Complex data, drop Y", {
		dat <- makeDat( "test_lowcover_complex_dropY" )
		withr::with_dir( dat$runDir, {
			saveBed( coverageGR, dat$coverageFile, strand= FALSE )
			saveBed( targetsGR, dat$targetFile, strand= FALSE )
			
			cli <- c("-t", dat$targetFile, "--coverageFile", dat$coverageFile, "--chrY", "chrY")
			expect_silent(lowcoverApp( cli ))
			
			it( "Creates expected files", {
				expect_true( file.exists( dat$regionsFile ))
				expect_true( file.exists( dat$badGenesFile ))
				expect_true( file.exists( dat$goodGenesFile ))
				expect_true( file.exists( dat$summaryFile ))
				expect_true( file.exists( dat$summaryFileNoY ))
			})
			it( "The regions file contents are correct", {
				got <- readLines( dat$regionsFile )
				want <- c(
					"chr1\t10\t25\tgene_A",
					"chr1\t200\t250\tgene_A",
					"chr1\t200\t250\tgene_B",
					"chr1\t290\t300\tgene_B",
					"chr1\t290\t310\tgene_C",
					"chr2\t10\t45\tgene_D",
					"chrY\t0\t100\tgene_YE"
				)
				expect_equal(got, want)
			})
			it( "The good genes file contents are correct", {
				got <- readLines( dat$goodGenesFile )
				want <- c( "gene_YF" )
				expect_equal(got, want)
			})
			it( "The bad genes file contents are correct", {
				got <- readLines( dat$badGenesFile )
				want <- c( "gene_A", "gene_B", "gene_C", "gene_D", "gene_YE" )
				expect_equal(got, want)
			})
			it( "The summary file contents are correct", {
				got <- read.delim( dat$summaryFile, header = FALSE,)
				wantNames <- c( "bases",    "lowBases",    "coveredBaseRat",
								"segments", "lowSegments", "coveredSegmentRat",
								"genes",    "lowGenes",    "coveredGeneRat" )
				expect_equal( got[,1], wantNames )
				wantValues <- c( 630,220,(630-220)/630,5,5,0,6,5,(6-5)/6 )
				expect_equal( got[,2], wantValues )
			})
			it( "The summary file excluding Y contents are correct", {
				got <- read.delim( dat$summaryFileNoY, header = FALSE,)
				wantNames <- c( "bases",    "lowBases",    "coveredBaseRat",
								"segments", "lowSegments", "coveredSegmentRat",
								"genes",    "lowGenes",    "coveredGeneRat" )
				expect_equal( got[,1], wantNames )
				wantValues <- c( 430,120,(430-120)/430,3,4,(3-4)/3,4,4,0 )
				expect_equal( got[,2], wantValues )
			})
		})
	})
	describe( "No low genes", {
		dat <- makeDat( "test_lowcover_noLow" )
		withr::with_dir( dat$runDir, {
			mcols(coverageGR)["tag"] <- rep.int("HI_COVERAGE", length(coverageGR))
			saveBed( coverageGR, dat$coverageFile, strand= FALSE )
			saveBed( targetsGR, dat$targetFile, strand= FALSE )

			cli <- c("-t", dat$targetFile, "--coverageFile", dat$coverageFile, "--chrY", "chrY")
			expect_silent(lowcoverApp( cli ))
			
			it( "Creates expected files", {
				expect_true( file.exists( dat$regionsFile ))
				expect_true( file.exists( dat$badGenesFile ))
				expect_true( file.exists( dat$goodGenesFile ))
				expect_true( file.exists( dat$summaryFile ))
			})
			it( "The regions file contents are correct", {
				got <- readLines( dat$regionsFile )
				expect_true( length(got) == 0 )
			})
			it( "The good genes file contents are correct", {
				got <- readLines( dat$goodGenesFile )
				want <- c( "gene_A", "gene_B", "gene_C", "gene_D", "gene_YE", "gene_YF" )
				expect_equal(got, want)
			})
			it( "The bad genes file contents are correct", {
				got <- readLines( dat$badGenesFile )
				expect_true( length(got) == 0 )
			})
			it( "The summary file contents are correct", {
				got <- read.delim( dat$summaryFile, header = FALSE,)
				wantNames <- c( "bases",    "lowBases",    "coveredBaseRat",
								"segments", "lowSegments", "coveredSegmentRat",
								"genes",    "lowGenes",    "coveredGeneRat" )
				expect_equal( got[,1], wantNames )
				wantValues <- c( 630,0,1,5,0,1,6,0,1 )
				expect_equal( got[,2], wantValues )
			})
			it( "The summary file excluding Y contents are correct", {
				got <- read.delim( dat$summaryFileNoY, header = FALSE,)
				wantNames <- c( "bases",    "lowBases",    "coveredBaseRat",
								"segments", "lowSegments", "coveredSegmentRat",
								"genes",    "lowGenes",    "coveredGeneRat" )
				expect_equal( got[,1], wantNames )
				wantValues <- c( 430,0,1,3,0,1,4,0,1 )
				expect_equal( got[,2], wantValues )
			})
		})
	})
	describe( "No good genes", {
		dat <- makeDat( "test_lowcover_all_low" )
		withr::with_dir( dat$runDir, {
			saveBed( targetsGR, dat$targetFile, strand= FALSE )
			coverageGR <- targetsGR
			mcols(coverageGR)["genes"] <- rep.int( "NO_COVERAGE", length( coverageGR ))
			saveBed( coverageGR, dat$coverageFile, strand= FALSE )
			
			cli <- c("-t", dat$targetFile, "--coverageFile", dat$coverageFile, "--chrY", "chrY")
			expect_silent(lowcoverApp( cli ))
			
			it( "Creates expected files", {
				expect_true( file.exists( dat$regionsFile ))
				expect_true( file.exists( dat$badGenesFile ))
				expect_true( file.exists( dat$goodGenesFile ))
				expect_true( file.exists( dat$summaryFile ))
				expect_true( file.exists( dat$summaryFileNoY ))
			})
			it( "The regions file contents are correct", {
				got <- readLines( dat$regionsFile )
				want <- c(
					"chr1\t10\t100\tgene_A",
					"chr1\t150\t250\tgene_A",
					"chr1\t200\t300\tgene_B",
					"chr1\t290\t400\tgene_C",
					"chr2\t10\t100\tgene_D",
					"chrY\t0\t100\tgene_YE",
					"chrY\t600\t700\tgene_YF"
				)
				expect_equal(got, want)
			})
			it( "The good genes file contents are correct", {
				got <- readLines( dat$goodGenesFile )
				expect_true( length(got) == 0 )
			})
			it( "The bad genes file contents are correct", {
				got <- readLines( dat$badGenesFile )
				want <- c( "gene_A", "gene_B", "gene_C", "gene_D", "gene_YE", "gene_YF" )
				expect_equal(got, want)
			})
			it( "The summary file contents are correct", {
				got <- read.delim( dat$summaryFile, header = FALSE,)
				wantNames <- c( "bases",    "lowBases",    "coveredBaseRat",
								"segments", "lowSegments", "coveredSegmentRat",
								"genes",    "lowGenes",    "coveredGeneRat" )
				expect_equal( got[,1], wantNames )
				wantValues <- c( 630,630,0,5,5,0,6,6,0 )
				expect_equal( got[,2], wantValues )
			})
			it( "The summary file excluding Y contents are correct", {
				got <- read.delim( dat$summaryFileNoY, header = FALSE,)
				wantNames <- c( "bases",    "lowBases",    "coveredBaseRat",
								"segments", "lowSegments", "coveredSegmentRat",
								"genes",    "lowGenes",    "coveredGeneRat" )
				expect_equal( got[,1], wantNames )
				wantValues <- c( 430,430,0,3,3,0,4,4,0 )
				expect_equal( got[,2], wantValues )
			})
		})
	})
})