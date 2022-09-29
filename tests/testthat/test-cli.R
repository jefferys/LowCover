
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
					regionsFile="LowCover.regions.bed",
					goodGenesFile="LowCover.goodgenes.txt",
					badGenesFile="LowCover.badgenes.txt",
					summaryFile="LowCover.summary.tsv",
					keep=c("NO_COVERAGE", "LOW_COVERAGE")
				)
				wantNames <- sort(names(want))
				want <- want[wantNames]
				got <- parseCLI()
				gotRaw <- got$rawOpts
				
				# Can compare with NA values, so check and remove
				# (both opt and opt$rawOpts have to be checked)
				expect_true( is.na( got$chrY ))
				expect_equal( gotRaw$chrY, "NA" )
				expect_true( is.na( got$summaryFileNoY ))
				expect_equal( gotRaw$summaryFileNoY, "LowCover.summaryNoY.tsv" )
				
				# clean and sort for direct  comparison
				got$chrY <- NULL
				gotRaw$chrY <- NULL
				got$summaryFileNoY <- NULL
				gotRaw$summaryFileNoY <- NULL
				got$rawOpts <- NULL
				gotNames <- sort(names(got))
				got <- got[gotNames]
				gotRaw <- gotRaw[sort(names(gotRaw))]
				
				expect_equal(got,want)

				# Fix raw content				
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
				aFile <- file.path( runDir, "deleteMe", "nestedDir", "aFile.bed" )
				expect_false( dir.exists( "deleteMe" ))
				cli <- c('--regionsFile', aFile, "--force" )
				expect_silent( parseCLI( cli ))
				expect_true( dir.exists( "deleteMe" ))
				expect_true( dir.exists( file.path( "deleteMe", "nestedDir" )))
				unlink( "deleteMe", recursive = TRUE )
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
			it( "Missing diectories created if --force", {
				aFile <- file.path( runDir, "deleteMe", "nestedDir", "aFile.bed" )
				expect_false( dir.exists( "deleteMe" ))
				cli <- c('-g', aFile, "--force" )
				expect_silent( parseCLI( cli ))
				expect_true( dir.exists( "deleteMe" ))
				expect_true( dir.exists( file.path( "deleteMe", "nestedDir" )))
				unlink( "deleteMe", recursive = TRUE )
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
			it( "Missing diectories created if --force", {
				aFile <- file.path( runDir, "deleteMe", "nestedDir", "aFile.bed" )
				expect_false( dir.exists( "deleteMe" ))
				cli <- c('--badGenesFile', aFile, "--force" )
				expect_silent( parseCLI( cli ))
				expect_true( dir.exists( "deleteMe" ))
				expect_true( dir.exists( file.path( "deleteMe", "nestedDir" )))
				unlink( "deleteMe", recursive = TRUE )
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
			it( "Missing diectories created if --force", {
				aFile <- file.path( runDir, "deleteMe", "nestedDir", "aFile.bed" )
				expect_false( dir.exists( "deleteMe" ))
				cli <- c('-s', aFile, "--force" )
				expect_silent( parseCLI( cli ))
				expect_true( dir.exists( "deleteMe" ))
				expect_true( dir.exists( file.path( "deleteMe", "nestedDir" )))
				unlink( "deleteMe", recursive = TRUE )
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
			it( "Handles single value", {
				cli <- c('--keep', "this or that or the other")
				got <- parseCLI( cli )
				expect_equal( got$keep, c("this or that or the other") )
			})
			it( "Strips spaces around commas", {
				cli <- c('-k', "this, or that ,or , the other")
				got <- parseCLI( cli )
				expect_equal( got$keep, c("this", "or that", "or", "the other") )
			})
		})
		describe( "--chrY", {
			it( "Long opt changes option value", {
				cli <- c('--chrY', "chrY")
				got <- parseCLI( cli )
				expect_equal( got$chrY, "chrY" )
			})
			it( "Short opt changes option value", {
				cli <- c('-y', "Y")
				got <- parseCLI( cli )
				expect_equal( got$chrY, "Y" )
			})
			it( "Handles setting to NA", {
				cli <- c('--chrY', "NA")
				got <- parseCLI( cli )
				expect_true( is.na( got$chrY ) )
			})
			it( "Handles empty string as NA", {
				cli <- c('-y', "")
				got <- parseCLI( cli )
				expect_true( is.na( got$chrY ) )
			})
		})
		describe( "--summaryFileNoY", {
			it( "Ignored if set without setting --chrY", {
				noSuchFile <- file.path( runDir, "noSuchFile.tsv" )
				cli <- c('--summaryFileNoY', noSuchFile )
				got <- parseCLI( cli )
				expect_true( is.na( got$summaryFileNoY ))
			})
			it( "Long opt changes option value", {
				noSuchFile <- file.path( runDir, "noSuchFile.tsv" )
				cli <- c('--summaryFileNoY', noSuchFile, '--chrY', 'chrY' )
				got <- parseCLI( cli )
				expect_equal( got$summaryFileNoY, noSuchFile )
			})
			it( "Short opt changes option value", {
				noSuchFile <- file.path( runDir, "noSuchFile.tsv" )
				cli <- c('-S', noSuchFile, '-y', 'chrY' )
				got <- parseCLI( cli )
				expect_equal( got$summaryFileNoY, noSuchFile )
			})
			it( "File can't exist without force", {
				aFile <- file.path( runDir, "aFile.bed" )
				cli <- c('--summaryFileNoY', aFile, '--chrY', 'chrY' )
				file.create( aFile )
				wantErrRe  <- paste0("--summaryFileNoY ", aFile, "\n\tFile exists \\(use --force to overwrite\\.\\)" )
				expect_error( parseCLI( cli ), wantErrRe )
				unlink( aFile )
			})
			it( "File deleted with warning if use --force", {
				aFile <- file.path( runDir, "aFile.bed" )
				cli <- c('-S', aFile, '-y', 'chrY', "--force" )
				file.create( aFile )
				wantWarnRe  <- paste0("--summaryFileNoY ", aFile, "\n\tExisting file deleted; will be replaced \\(as --force used\\.\\)" )
				expect_warning( parseCLI( cli ), wantWarnRe )
				unlink( aFile )
			})
			it( "Missing diectories created if --force", {
				aFile <- file.path( runDir, "deleteMe", "nestedDir", "aFile.bed" )
				expect_false( dir.exists( "deleteMe" ))
				cli <- c('--summaryFileNoY', aFile, '--chrY', 'chrY',  "--force" )
				expect_silent( parseCLI( cli ))
				expect_true( dir.exists( "deleteMe" ))
				expect_true( dir.exists( file.path( "deleteMe", "nestedDir" )))
				unlink( "deleteMe", recursive = TRUE )
			})
		})
		describe( "--version", {
			desc <- packageDescription( packageName() )
			it( "Long opt changes option value", {
				cli <- c('--version' )
				got <- parseCLI( cli )
				expect_equal( got$version, paste( packageName(), desc$Version ))
			})
			it( "Short opt changes option value", {
				cli <- c('-V' )
				got <- parseCLI( cli )
				expect_equal( got$version, paste( packageName(), desc$Version ))
			})
		})
		describe( "--help", {
			desc <- packageDescription( packageName() )
			it( "Long opt changes option value", {
				cli <- c( '--help' )
				got <- parseCLI( cli )
				expect_match( got$help, desc$Title)
			})
			it( "Short opt changes option value", {
				cli <- c( '-h' )
				got <- parseCLI( cli )
				expect_match( got$help, desc$Description )
			})
		})
	})
})
