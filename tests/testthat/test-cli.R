
runDir <- tempfile("test-cli")
dir.create(runDir)
targetFile <- file.path(runDir, "targets.bed")
file.create( targetFile )
coverageFile <- file.path(runDir, "coverage.bed")
file.create( coverageFile )
withr::with_dir( runDir, {
	describe( "Command line parsing", {
		it( "identifies unknown paramaters when specified", {
			cli <- c( '--whatIsThis', "some value" )
			wantErrRE <- "Undefined argument labels supplied"
			wantMsgRE <- "Argument '--whatIsThis' is not a defined optional argument or flag"
			expect_message(
				expect_error( parseCLI( cli ), wantErrRE),
				wantMsgRE
			)
		})
		describe( "default behavior", {
			it( "command line with no options parses without error", {
				expect_silent( parseCLI() )
			})
			it( "Has the expected default 'raw' options", {
				wantRaw <- list(
					targetsFile=    "targets.bed",
					coverageFile=   "coverage.bed",
					force=          FALSE,
					regionsFile=    "LowCover.regions.bed",
					goodGenesFile=  "LowCover.goodgenes.txt",
					badGenesFile=   "LowCover.badgenes.txt",
					summaryFile=    "LowCover.summary.tsv",
					summaryFileNoY= "LowCover.summaryNoY.tsv",
					keep=           "NO_COVERAGE,LOW_COVERAGE",
					chrY=           "NA",
					parallel=       "default",
					workers=        NA_integer_
				)
				got <- parseCLI()
				gotRaw <- got$rawOpts
				expect_equal( sort(names(gotRaw)), sort(names(wantRaw)))
				for (opt in names(gotRaw)) {
					if (opt %in% c("workers")) {
						expect_true( is.na( gotRaw[[opt]] ), label= opt)
					}
					else {
						expect_equal( gotRaw[[opt]], wantRaw[[opt]], label= opt )
					}
				}
			})
			it( "Has the expected default 'parsed' options", {
				want <- list(
					targetsFile=    "targets.bed",
					coverageFile=   "coverage.bed",
					force=          FALSE,
					regionsFile=    "LowCover.regions.bed",
					goodGenesFile=  "LowCover.goodgenes.txt",
					badGenesFile=   "LowCover.badgenes.txt",
					summaryFile=    "LowCover.summary.tsv",
					summaryFileNoY= NA, # as chrY is NA.
					keep=           c("NO_COVERAGE", "LOW_COVERAGE"),
					chrY=           NA,
					parallel=       "default",
					workers=        NA_integer_
				)
				got <- parseCLI()
				got$rawOpts <- NULL
				expect_equal( sort(names(got)), sort(names(want)))
				for (opt in names(got)) {
					if (anyNA(got[[opt]])) {
						expect_true( is.na( want[[opt]] ), label = opt)
					}
					else {
						expect_equal( got[[opt]], want[[opt]], label = opt)
					}
				}
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
		describe(  "--parallel", {
			it( "Long opt changes option value", {
				cli <- c('--parallel', "sequential")
				got <- parseCLI( cli )
				expect_equal( got$parallel, "sequential" )
			})
			it( "Short opt changes option value", {
				cli <- c('-p', "multisession")
				got <- parseCLI( cli )
				expect_equal( got$parallel, "multisession" )
			})
			it( "Converts upper-case to lower case", {
				cli <- c('-p', "MultiSession")
				got <- parseCLI( cli )
				expect_equal( got$parallel, "multisession" )
			})
			it( "Error (with original case) if invalid strategies are suggested.", {
				cli <- c('-p', "NotTheStrategy")
				wantErrRE <- " --parallel NotTheStrategy.+Not a valid value for this option; see --help\\."
				expect_error( parseCLI( cli ), wantErrRE )
			})
			it( "It replaces 'guess' with a strategy based on 'future::supportsMulticore", {
				cli <- c('-p', "guess")
				mockthat::with_mock( `future::supportsMulticore` = function(...) TRUE, {
					got <- parseCLI( cli )
					expect_equal( got$parallel, "multicore" )
				})
				mockthat::with_mock( `future::supportsMulticore` = function(...) FALSE, {
					got <- parseCLI( cli )
					expect_equal( got$parallel, "multisession" )
				})
			})
		})
		describe(  "--workers", {
			mockthat::with_mock( `future::availableCores` = function(...) 2, {
				it( "Long opt changes option value", {
					cli <- c('--parallel', 'multisession', '--workers', '1')
					got <- parseCLI( cli )
					expect_equal( got$workers, 1 )
				})
				it( "Short opt changes option value", {
					cli <- c('-p', 'multisession', '-w', '1')
					got <- parseCLI( cli )
					expect_equal( got$workers, 1 )
				})
				it( "warning that it is ignored when --parallel is unset, 'default' or 'sequential'", {
					cli <- c('-w', '10')
					wantWarnRE <- " --workers 10.+Ignored with --parallel 'default'\\."
					expect_warning( parseCLI( cli ), wantWarnRE )

					cli <- c('-p', 'default', '--workers', '1')
					wantWarnRE <- " --workers 1.+Ignored with --parallel 'default'\\."
					expect_warning( parseCLI( cli ), wantWarnRE )
					
					cli <- c('--parallel', 'sequential', '-w', '-1')
					wantWarnRE <- " --workers -1.+Ignored with --parallel 'sequential'\\."
					expect_warning( parseCLI( cli ), wantWarnRE )
				})
				it( "replaces values <= 0 or > max with max = 'future::availableCores()'", {
					cli <- c('--parallel', 'multisession', '--workers', '0')
					got <- parseCLI( cli )
					expect_equal( got$workers, 2 )  # availableCores() is mocked
					
					cli <- c('-w', '-102', '-p', 'multicore')
					got <- parseCLI( cli )
					expect_equal( got$workers, 2 )  # availableCores() is mocked
					
					cli <- c('-w', '102', '-p', 'multisession')
					got <- parseCLI( cli )
					expect_equal( got$workers, 2 )  # availableCores() is mocked
				})
			})
		})
		
	})
})
