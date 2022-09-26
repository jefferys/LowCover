
runDir <- tempfile("test-cli")
dir.create(runDir)
targetFile <- file.path(runDir, "targets.bed")
file.create( targetFile )
coverageFile <- file.path(runDir, "coverage.bed")
file.create( coverageFile )
withr::with_dir( runDir, {
	describe( "Command line parsing", {
		describe( "default behavior", {
			it( "will parse without error", {
				expect_silent( parseCLI() )
			})
			it( "Has correct defaults for all options, including rawOpts", {
				want <- list(
					targetsFile="targets.bed",
					coverageFile="coverage.bed",
					force=FALSE,
					regionsFile="lowcover.regions.bed",
					goodGenesFile="lowcover.goodgenes.txt",
					badGenesFile="lowcover.badgenes.txt",
					summaryFile="lowcover.summary.tsv",
					keep=c("NO_COVERAGE", "LOW_COVERAGE")
				)
				wantNames <- sort(names(want))
				want <- want[wantNames]
				got <- parseCLI()
				gotRaw <- got$rawOpts
				got$rawOpts <- NULL
				gotNames <- sort(names(got))
				got <- got[gotNames]
				expect_equal(got,want)
				gotRaw <- gotRaw[sort(names(gotRaw))]
				wantRaw <- want
				wantRaw$keep <- paste0( wantRaw$keep, collapse= ",")
				expect_equal(gotRaw,wantRaw)
			})
		})
		describe( "--targetsFile (-t)", {
			it( "Long opt changes option value", {
				aFile <- file.path( runDir, "aFile.bed" )
				cli <- c('--targetsFile', aFile )
				file.create( aFile )
				got <- parseCLI( cli )
				expect_equal( got$targetsFile, aFile )
				unlink( aFile )
			})
			it( "Short opt changes option value", {
				aFile <- file.path( runDir, "aFile.bed" )
				cli <- c('-t', aFile )
				file.create( aFile )
				got <- parseCLI( cli )
				expect_equal( got$targetsFile, aFile )
				unlink( aFile )
			})
			it( "File must exist", {
				cli <- c('--targetsFile', "noSuchFile.bed" )
				wantErrRe  <- "--targetsFile noSuchFile\\.bed\n\tCan't find file\\."
				expect_error( parseCLI( cli ), wantErrRe )
			})
		})
		describe( "--coverageFile", {
			it( "Long opt changes option value", {
				aFile <- file.path( runDir, "aFile.bed" )
				cli <- c('--coverageFile', aFile )
				file.create( aFile )
				got <- parseCLI( cli )
				expect_equal( got$coverageFile, aFile )
				unlink( aFile )
			})
			it( "Short opt changes option value", {
				aFile <- file.path( runDir, "aFile.bed" )
				cli <- c('-c', aFile )
				file.create( aFile )
				got <- parseCLI( cli )
				expect_equal( got$coverageFile, aFile )
				unlink( aFile )
			})
			it( "File must exist", {
				cli <- c('-c', "noSuchFile.bed" )
				wantErrRe  <- "--coverageFile noSuchFile\\.bed\n\tCan't find file\\."
				expect_error( parseCLI( cli ), wantErrRe )
			})
		})
		describe( "--regionsFile", {
			it( "Long opt changes option value", {
				noSuchFile <- file.path( runDir, "noSuchFile.bed" )
				cli <- c('--regionsFile', noSuchFile )
				got <- parseCLI( cli )
				expect_equal( got$regionsFile, noSuchFile )
			})
			it( "Short opt changes option value", {
				noSuchFile <- file.path( runDir, "noSuchFile.bed" )
				cli <- c('-r', noSuchFile )
				got <- parseCLI( cli )
				expect_equal( got$regionsFile, noSuchFile )
			})
			it( "File can't exist without force", {
				aFile <- file.path( runDir, "aFile.bed" )
				cli <- c('--regionsFile', aFile )
				file.create( aFile )
				wantErrRe  <- paste0("--regionsFile ", aFile, "\n\tFile exists \\(use --force to overwrite\\.\\)" )
				expect_error( parseCLI( cli ), wantErrRe )
				unlink( aFile )
			})
			it( "File deleted with warning if use --force", {
				aFile <- file.path( runDir, "aFile.bed" )
				cli <- c('--regionsFile', aFile, "--force" )
				file.create( aFile )
				wantWarnRe  <- paste0("--regionsFile ", aFile, "\n\tExisting file deleted; will be replaced \\(as --force used\\.\\)" )
				expect_warning( parseCLI( cli ), wantWarnRe )
				unlink( aFile )
			})
			it( "Missing diectories created if --force", {
				aFile <- file.path( runDir, "aFile.bed" )
				cli <- c('--regionsFile', aFile, "--force" )
				file.create( aFile )
				wantWarnRe  <- paste0("--regionsFile ", aFile, "\n\tExisting file deleted; will be replaced \\(as --force used\\.\\)" )
				expect_warning( parseCLI( cli ), wantWarnRe )
				unlink( aFile )
			})
		})
		describe( "--goodGenesFile", {
			it( "Long opt changes option value", {
				noSuchFile <- file.path( runDir, "noSuchFile.txt" )
				cli <- c('--goodGenesFile', noSuchFile )
				got <- parseCLI( cli )
				expect_equal( got$goodGenesFile, noSuchFile )
			})
			it( "Short opt changes option value", {
				noSuchFile <- file.path( runDir, "noSuchFile.txt" )
				cli <- c('-g', noSuchFile )
				got <- parseCLI( cli )
				expect_equal( got$goodGenesFile, noSuchFile )
			})
			it( "File can't exist without force", {
				aFile <- file.path( runDir, "aFile.bed" )
				cli <- c('-g', aFile )
				file.create( aFile )
				wantErrRe  <- paste0("--goodGenesFile ", aFile, "\n\tFile exists \\(use --force to overwrite\\.\\)" )
				expect_error( parseCLI( cli ), wantErrRe )
				unlink( aFile )
			})
			it( "File deleted with warning if use --force", {
				aFile <- file.path( runDir, "aFile.bed" )
				cli <- c('--goodGenesFile', aFile, "--force" )
				file.create( aFile )
				wantWarnRe  <- paste0("--goodGenesFile ", aFile, "\n\tExisting file deleted; will be replaced \\(as --force used\\.\\)" )
				expect_warning( parseCLI( cli ), wantWarnRe )
				unlink( aFile )
			})
			
		})
		describe( "--badGenesFile", {
			it( "Long opt changes option value", {
				noSuchFile <- file.path( runDir, "noSuchFile.txt" )
				cli <- c('--badGenesFile', noSuchFile )
				got <- parseCLI( cli )
				expect_equal( got$badGenesFile, noSuchFile )
			})
			it( "Short opt changes option value", {
				noSuchFile <- file.path( runDir, "noSuchFile.txt" )
				cli <- c('-b', noSuchFile )
				got <- parseCLI( cli )
				expect_equal( got$badGenesFile, noSuchFile )
			})
			it( "File can't exist without force", {
				aFile <- file.path( runDir, "aFile.bed" )
				cli <- c('--badGenesFile', aFile )
				file.create( aFile )
				wantErrRe  <- paste0("--badGenesFile ", aFile, "\n\tFile exists \\(use --force to overwrite\\.\\)" )
				expect_error( parseCLI( cli ), wantErrRe )
				unlink( aFile )
			})
			it( "File deleted with warning if use --force", {
				aFile <- file.path( runDir, "aFile.bed" )
				cli <- c('-b', aFile, "--force" )
				file.create( aFile )
				wantWarnRe  <- paste0("--badGenesFile ", aFile, "\n\tExisting file deleted; will be replaced \\(as --force used\\.\\)" )
				expect_warning( parseCLI( cli ), wantWarnRe )
				unlink( aFile )
			})
		})
		describe( "--summaryFile", {
			it( "Long opt changes option value", {
				noSuchFile <- file.path( runDir, "noSuchFile.tsv" )
				cli <- c('--summaryFile', noSuchFile )
				got <- parseCLI( cli )
				expect_equal( got$summaryFile, noSuchFile )
			})
			it( "Short opt changes option value", {
				noSuchFile <- file.path( runDir, "noSuchFile.tsv" )
				cli <- c('-s', noSuchFile )
				got <- parseCLI( cli )
				expect_equal( got$summaryFile, noSuchFile )
			})
			it( "File can't exist without force", {
				aFile <- file.path( runDir, "aFile.bed" )
				cli <- c('-s', aFile )
				file.create( aFile )
				wantErrRe  <- paste0("--summaryFile ", aFile, "\n\tFile exists \\(use --force to overwrite\\.\\)" )
				expect_error( parseCLI( cli ), wantErrRe )
				unlink( aFile )
			})
			it( "File deleted with warning if use --force", {
				aFile <- file.path( runDir, "aFile.bed" )
				cli <- c('--summaryFile', aFile, "--force" )
				file.create( aFile )
				wantWarnRe  <- paste0("--summaryFile ", aFile, "\n\tExisting file deleted; will be replaced \\(as --force used\\.\\)" )
				expect_warning( parseCLI( cli ), wantWarnRe )
				unlink( aFile )
			})
		})
		describe( "--force", {
			it( "Long opt changes option value", {
				cli <- c('--force' )
				got <- parseCLI( cli )
				expect_equal( got$force, TRUE )
			})
			it( "Short opt changes option value", {
				cli <- c('-f' )
				got <- parseCLI( cli )
				expect_equal( got$force, TRUE )
			})
			
		})
		describe( "--keep", {
			it( "Long opt changes option value", {
				cli <- c('--keep', "this,or that,or the other")
				got <- parseCLI( cli )
				expect_equal( got$keep, c("this", "or that", "or the other") )
			})
			it( "Short opt changes option value", {
				cli <- c('-k', "this,or that,or the other")
				got <- parseCLI( cli )
				expect_equal( got$keep, c("this", "or that", "or the other") )
			})
		})
	})
})
