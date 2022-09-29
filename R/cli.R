initParser <- function( desc ) {
	argparser::arg_parser( paste(desc$Title, desc$Description, sep="\n\n"), desc$Package, hide.opts= TRUE )
}

defineInterface <- function( parser ) {
	parser <- argparser::add_argument(
		parser, "--help", short= "-h", flag= TRUE,
		help= "Show this help message and exit."
	)
	parser <- argparser::add_argument(
		parser, "--version", short= "-V", flag= TRUE,
		help= "Print version info to stderr and exit."
	)
	parser <- argparser::add_argument(
		parser, "--targetsFile", short="-t", default="targets.bed",
		help= "Bed4 file of target regions to limit coverage reporting to (colum #4 = gene name), gzipped ok."
	)
	parser <- argparser::add_argument(
		parser, "--coverageFile", short="-c", default="coverage.bed",
		help= "Bed4 file reporting regions and their coverage class (column #4 = coverage tag), gzipped ok."
	)
	parser <- argparser::add_argument(
		parser, "--keep", short="-k", default="NO_COVERAGE,LOW_COVERAGE",
		help= "Comma separated string of tags identifying low coverage regions."
	)
	parser <- argparser::add_argument(
		parser, "--chrY", short="-y", default="NA",
		help= "Name of 'Y' chromosome, if want to filter by sex."
	)
	parser <- argparser::add_argument(
		parser, "--force", short= "-f", flag= TRUE,
		help= "Create any non-existing directories and overwrite any existing files."
	)
	parser <- argparser::add_argument(
		parser, "--regionsFile", short="-r", default="LowCover.regions.bed",
		help= "Bed4 file of regions with low coverage (column #4 = gene name)"
	)
	parser <- argparser::add_argument(
		parser, "--badGenesFile", short= "-b", default="LowCover.badgenes.txt",
		help= "Text file listing  genes with at least 1 base with low coverage (one per line)"
	)
	parser <- argparser::add_argument(
		parser, "--goodGenesFile", short= "-g", default="LowCover.goodgenes.txt",
		help= "Text file listing genes with no low covage (one per line)"
	)
	parser <- argparser::add_argument(
		parser, "--summaryFile", short= "-s", default="LowCover.summary.tsv",
		help= "Stats table"
	)
	parser <- argparser::add_argument(
		parser, "--summaryFileNoY", short= "-S", default="LowCover.summaryNoY.tsv",
		help= "Stats table ignoring chrY. Ignored if --chrY not set"
	)
}

getRawOpt <- function( opts, optName ) {
	if ( ! "rawOpts" %in% names(opts)) {
		stop( "Can't find raw options list `rawOpts`.",  .Call= FALSE )
	}
	opts <- opts$rawOpts
	if ( ! optName %in% names(opts)) {
		stop( "Can't find option '", optName, "' in raw options.",  .Call= FALSE )
	}
	opts[[optName]]
}

getOpt <- function( opts, optName ) {
	if ( ! optName %in% names(opts)) {
		stop( "Can't find option '", optName, "' in parsed options.",  .Call= FALSE )
	}
	opts[[optName]]
}

setOpt <- function( opts, optName, optValue ) {
	opts[[optName]] <- optValue
	opts
}

optError <- function( optName, optValue, message ) {
	stop( "--", optName, " ", optValue, "\n\t", message, call. = FALSE )
}

optWarning <- function( optName, optValue, message ) {
	warning( "--", optName, " ", optValue, "\n\t", message, call. = FALSE )
}

assertOptInFile <- function( optName,  optValue ) {
	if (! checkmate::testFileExists( optValue )) {
		optError( optName, optValue, "Can't find file." )
	}
	if (! checkmate::testFileExists( optValue, access = "r" )) {
		optError( optName, optValue, "File is not readable." )
	}
	TRUE
}

assertOptOutDir <- function( optName, dir, force= NA ) {
	if (! dir.exists( dir )) {
		if (is.na(force)) {
			optError( optName, dir,
					  "Output directory must exist." )
		}
		else if (! force) {
			optError( optName, dir,
					  "Output directory must exist (unless --force creates it.)" )
		}
		else {
			dir.create( dir, showWarnings = TRUE, recursive = TRUE )
			if (! dir.exists( dir )) {
				optError( optName, dir,
						  "Failed to create output directory (with --force).")
			}
		}
	}
	if (! checkmate::testDirectory( dir, access = "w" )) {
		optError( optName, dir,
				  "Output directory must be writeable." )
	}
	TRUE
}


#' Assert option is valid output file name
#' 
#' Checks that a value provided as an option is suitable for use as an output
#' file path. Can be relative or absolute. It is an error if the file exists
#' or its directory does not,  unless "force" is set TRUE, in which case the
#' directory will be (silently) created or the pre-existing file will be
#' deleted  (with a warning). Failure to create a needed directory or delete
#' a pre-existing file also results in an error.
#'
#' @param optName The name of the option (used in error messages).
#' @param file The file name (path) to check for useability.
#' @param force If set TRUE, will delete `file` if it exists (with a warning)
#'   and create any needed directories (silently). If FALSE or NA, these will be
#'   errors, with only the text of the error differing. NA, the default, means
#'   no 'force' option is mentioned, while FALSE implies `--force` is available
#'   as an option, but was not set.
#'
#' @return Exits with error or returns TRUE
#' @export
#'
assertOptOutFile <- function( optName, file, force= NA ) {
	if (file.exists( file )) {
		if (is.na(force)) {
			optError( optName, file,
					  "File exists." )
		}
		else if (! force) {
			optError( optName, file,
					  "File exists (use --force to overwrite.)" )
		}
		else { # force == TRUE
			if (! checkmate::testPathForOutput( file, overwrite = TRUE )) {
				optError( optName, file,
						  "File exists and appears impossible to overwrite." )
			}
			optWarning( optName, file,
						"Existing file deleted; will be replaced (as --force used.)" )
			file.remove(file)
		}
	}
	pathTo <- dirname(file)
	if ( ! dir.exists( pathTo ) || ! isSamePath( getwd(), pathTo)) {
		assertOptOutDir( optName, pathTo, force )
	}
	TRUE
}

ensureOpt_targetsFile <- function( opts ) {
	optName <- "targetsFile"
	optValue  <- getRawOpt( opts, optName )
	assertOptInFile( optName, optValue)

	setOpt( opts, optName, optValue)
}

ensureOpt_coverageFile <- function( opts ) {
	optName <- "coverageFile"
	optValue  <- getRawOpt( opts, optName )
	assertOptInFile( optName, optValue)
	
	setOpt( opts, optName, optValue)
}

ensureOpt_force <- function( opts ) {
	optName <- "force"
	optValue <- getRawOpt( opts, optName )
	# Intrface ensures flags are True or False values
	setOpt( opts, optName, optValue)
}

ensureOpt_regionsFile <- function( opts ) {
	optName <- "regionsFile"
	optValue  <- getRawOpt( opts, optName )
	force <- getOpt( opts, "force" )
	assertOptOutFile( optName, optValue, force)
	
	setOpt( opts, optName, optValue)
}

ensureOpt_badGenesFile <- function( opts ) {
	optName <- "badGenesFile"
	optValue  <- getRawOpt( opts, optName )
	force <- getOpt( opts, "force" )
	assertOptOutFile( optName, optValue, force)
	
	setOpt( opts, optName, optValue)
}

ensureOpt_goodGenesFile <- function( opts ) {
	optName <- "goodGenesFile"
	optValue  <- getRawOpt( opts, optName )
	force <- getOpt( opts, "force" )
	assertOptOutFile( optName, optValue, force)
	
	setOpt( opts, optName, optValue)
}

ensureOpt_summaryFile <- function( opts ) {
	optName <- "summaryFile"
	optValue  <- getRawOpt( opts, optName )
	force <- getOpt( opts, "force" )
	assertOptOutFile( optName, optValue, force)
	
	setOpt( opts, optName, optValue)
}

ensureOpt_keep <- function( opts ) {
	optName <- "keep"
	optValue  <- getRawOpt( opts, optName )
	
	tags <- strsplit( optValue, "[\\s]*[,][\\s]*", perl= TRUE )[[1]]
	setOpt( opts, optName, tags)
}

ensureOpt_chrY <- function( opts ) {
	optName <- "chrY"
	optValue <- getRawOpt( opts, optName )
	
	if (optValue == "NA" || nchar(optValue) < 1) {
		optValue <- NA
	}
	setOpt( opts, optName, optValue )
}

ensureOpt_summaryFileNoY <- function( opts ) {
	optName <- "summaryFileNoY"
	optValue  <- getRawOpt( opts, optName )
	chrY <- getOpt( opts, "chrY" )
	if (is.na(chrY)) {
		optValue <- NA
	}
	else {
		force <- getOpt( opts, "force" )
		assertOptOutFile( optName, optValue, force)
	}
	setOpt( opts, optName, optValue)
}

#' Validate CLI options
#'
#' Validates (and possibly transforms) command line options. Halts with error on
#' the first invalid option found, or returns a list of valid options and their
#' post-validation values. The list of raw options, including any invalid opts,
#' is included in the returned list as `rawOpts` and returns a list of valid
#' options.How unknown options are handled depends on `unvalidated`.
#' @param rawOpts The unvalidated options list
#' @param unvalidated What to do with unvalidated options:
#'
#' - `warning` - By default if an option is not validated, it is used as is, but
#' a warning will be generated.
#' - `error`  - Unvalidated options will trigger an error.
#' - `ok` - Unvalidated options will silently be used as is.
#' - `ignore` - Unvalidated options will silently be ignored.
#'
#' @return The list of validated options, including any unvalidated options if
#'   `unvalidated` is "warning" or "ok".
#' 
validateOptions <- function(
	rawOpts, unvalidated= c("warning", "error", "ok", "ignore" )
) {
	unvalidated <- match.arg(unvalidated)
	# Ensure not overwriting existing option "rawOpts"
	if ( "rawOpts" %in% names( rawOpts )) {
		stop( "Can't have a CLI element named 'rawOpts'.", call. = FALSE )
	}
	
	# If unvalidated options are ok, just add them all as "validated".
	# Otherwise, options are added one by one as validated.
	if (unvalidated == "ok") {
		validOpts <- rawOpts
	}
	else {
		validOpts <- list()
	}
	
	# Want to track the original options as presented on the command line for
	# use in error messages and to identify any unvalidated options if needed.
	validOpts$rawOpts <- rawOpts
	
	# Validate and format options
	validOpts <- ensureOpt_force(          validOpts ) 
	validOpts <- ensureOpt_targetsFile(    validOpts )
	validOpts <- ensureOpt_coverageFile(   validOpts )
	validOpts <- ensureOpt_regionsFile(    validOpts )
	validOpts <- ensureOpt_badGenesFile(   validOpts )
	validOpts <- ensureOpt_goodGenesFile(  validOpts )
	validOpts <- ensureOpt_summaryFile(    validOpts )
	validOpts <- ensureOpt_keep(           validOpts )
	validOpts <- ensureOpt_chrY(           validOpts )
	validOpts <- ensureOpt_summaryFileNoY( validOpts )
	
	# Report on unvalidated options, if required.
	# (Extra "validOpts" added during validation are fine)
	if ( unvalidated == "warning" || unvalidated  == "error" ) {
		select <- names(rawOpts) %in% names(validOpts)
		extraOpts <- names(rawOpts)[ ! select ]
		if (length(extraOpts) > 0)  {
			optsString <- paste0("'", paste( extraOpts, collapse = "', '"), "'" )
			if (unvalidated == "error") {
				stop( "The following unexpected command line options were present: ",
					  optsString, ".", call. = FALSE )
			}
			if (unvalidated == "warning") {
				warning( "Ignoring unused command line options (",
								 optsString, ").", call. = FALSE )
			}
		}
	}
	
	validOpts
}

#' Parse the command line options and arguments.
#'
#' Returns an unvalidated list of options by long name. (-h,--help) and
#' (-V,--version) are parsed manually, so only help=TRUE will be present if help
#' specified, or if no help asked for but versions is specified, only
#' version=TRUE will be present.
#'
#' @param args The vector or command line tokens to parse, defaults to
#'   retrieving them itself.
#'
#' @return A list of the parsed command line options and arguments. Contains an
#'   unnamed option as the first element due to undocumented behavior in
#'   argparser.
#'   
parseCLI <- function( args= commandArgs( trailingOnly= TRUE )) {
	myPackage <- utils::packageName()
	desc <- utils::packageDescription( myPackage )
	
	parser <- initParser( desc )
	parser <- defineInterface( parser )
	
	# Hack: Don't want to call "quit()" after printing help message as argparse
	# will if allowed to handle help itself. So have to preparse and handle "-h"
	# or "--help". If specified, the help message that argparse would normally
	# generate is saved as the value of the `help` option, and this is returned.
	wantHelp <- any(args %in% c("-h", "--help"))
	if (wantHelp) {
		helpMessage <- paste( utils::capture.output( print(parser), type="message" ), collapse= "\n" ) 
		return( list( "help"= helpMessage ))
	}

	# Since have to handle help manually here, also handle version manually here,
	# printing a custom version message to STDOUT.
	wantVersion <- any(args %in% c("-V", "--version"))
	if (wantVersion) {
		versionMessage <- paste0( utils::packageName(),  " ", desc$Version ) 
		return( list( "version"= versionMessage ))
	}
	
	rawOpts <- argparser::parse_args(parser, args)
	
	# Remove weird initial unnamed option and
	# already handled "--help" and "--version"
	rawOpts[nchar(names(rawOpts)) < 1] <- NULL
	rawOpts[c("help", "version")] <- NULL
	validateOptions( rawOpts )
}
