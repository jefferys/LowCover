### Data for testing

#  Minimal 3-column genomic location data
bed3df <- data.frame(
	seqnames= c( "chr1", "chr1", "chr1", "chr2",     "chr10" ),
	start=    c( 0L,     10L,    10L,    0L,         15L     ),
	end=      c( 10L,     20L,    20L,    249000000L, 15L     )
)

# 6-column genomic location data with strand and 2 extra columns, score and id
bed6df <- data.frame(
	seqnames= c( "chr1", "chr1", "chr1", "chr2",     "chr10" ),
	start=    c( 0L,     10L,    10L,    0L,         15L     ),
	end=      c( 10L,     20L,    20L,    249000000L, 15L     ),
	strand=   c( "*",    "+",    "+",    "*",        "-"     ),
	score=    c( 1.1,    2.2,    3.3,    4.4,        5.5     ),
	id=       c( "bob",  "sue",  "tom",  "april",    "alice" )
)

# GRanges object equivalent to bed6df
grPlus2 <- GRanges(
	seqnames= bed6df$seqnames,
	ranges= IRanges::IRanges( start= bed6df$start + 1, end= bed6df$end ),
	strand= bed6df$strand,
	score= bed6df$score,
	id= bed6df$id
)

describe( "loadBed()", {
	describe("Load a bed file with no extra columns", {
		
		bed3File <- tempfile( fileext = ".bed" )
		write.table( bed3df, bed3File, quote=FALSE, row.names = FALSE, col.names = FALSE,  sep="\t")
		
		describe("Using default parameters", {
			it( "Runs silently if no error", {
				expect_silent( loadBed( bed3File ))
			})
			describe("It returns the correct result", {
				gr <- loadBed( bed3File )
				it( "is a Genomic Ranges object", {
					expect_true( is( gr, "GenomicRanges" ))
				})
				it( "has location columns matched to the bed file", {
					expect_equal( as.character(seqnames(gr)), bed3df$seqnames )
					expect_equal( as.integer(start(gr)),    bed3df$start + 1 )
					expect_equal( as.integer(end(gr)),      bed3df$end )
				})
				it( "has default '*' values for strandedness",  {
					expect_equal( as.character(strand(gr)), rep.int ("*", nrow( bed3df )))
				})
				it( "has no metadata columns", {	
					expect_equal( length(mcols(gr)), 0 )
				})
			})
		})
	})
	describe("Load a bed file with extra columns", {
		
		bed6File <- tempfile( fileext = ".bed" )
		write.table( bed6df, bed6File, quote=FALSE, row.names = FALSE, col.names = FALSE,  sep="\t")
		
		describe("Using default paramters", {
			it( "Runs silently if no error", {
				expect_silent( loadBed( bed6File ))
			})
			describe( "It returns the correct result", {
				gr <- loadBed( bed6File )
				it( "Retuns a Genomic Ranges object", {
					expect_true( is( gr, "GenomicRanges" ))
				})
				it( "has location columns matched to the bed file", {
					expect_equal( as.character(seqnames(gr)), bed6df$seqnames )
					expect_equal( as.integer(start(gr)),    bed6df$start + 1 )
					expect_equal( as.integer(end(gr)),      bed6df$end )
				})
				it( "has default '*' values for strandedness",  {
					expect_equal( as.character(strand(gr)),   rep.int ("*", nrow( bed6df )))
				})
				it( "has no metadata columns", {
					expect_equal( length(mcols(gr)), 0 )
				})
			})
		})
		describe("With metadata columns specified (including strandedness)", {
			colnames = names(bed6df)
			it( "Runs silently if no error", {
				expect_silent( loadBed( bed6File, colnames = colnames ))
			})
			describe( "It returns the correct result", {
				gr <- loadBed( bed6File, colnames = colnames )
				it( "Retuns a Genomic Ranges object", {
					expect_true( is( gr, "GenomicRanges" ))
				})
				it( "has location columns matched to the bed file", {
					expect_equal( as.character(seqnames(gr)), bed6df$seqnames )
					expect_equal( as.integer(start(gr)),    bed6df$start + 1 )
					expect_equal( as.integer(end(gr)),      bed6df$end )
				})
				it( "has specified values for strandedness",  {
					expect_equal( as.character(strand(gr)), bed6df$strand )
				})
				it( "has the expected metadata columns", {
					meta <- mcols(gr)
					expect_equal( names(meta), c("score", "id") )
					expect_equal( mcols(gr)[,"score"], bed6df$score )
					expect_equal( mcols(gr)[,"id"], bed6df$id )
				})
			})
		})
		describe( "With some columns not named or named as NA", {
			colnames <- c( "seqnames", "start", "end", NA, "score" ) 
			it( "Runs silently if no error", {
				expect_silent( loadBed( bed6File, colnames = colnames ))
			})
			describe( "It returns the correct result", {
				gr <- loadBed( bed6File, colnames = colnames )
				
				it( "Retuns a Genomic Ranges object", {
					expect_true( is( gr, "GenomicRanges" ))
				})
				it( "has location columns matched to the bed file", {
					expect_equal( as.character(seqnames(gr)), bed6df$seqnames )
					expect_equal( as.integer(start(gr)),    bed6df$start + 1 )
					expect_equal( as.integer(end(gr)),      bed6df$end )
				})
				it( "has default '*' values for NA (skipped) strand column",  {
					expect_equal( as.character(strand(gr)), rep.int ("*", nrow( bed6df )))
				})
				it( "has the expected metadata columns", {	
					meta <- mcols(gr)
					expect_equal( names(meta), c("score" ) )
					expect_equal( mcols(gr)[,"score"], bed6df$score )
				})
			})
		})
		describe( "Loading from a connection", {
			colnames = names(bed6df)
			
			it( "Runs silently if no error", {
				con <- file(bed6File)
				expect_silent( loadBed( con, colnames = colnames ))
			})
			describe( "It returns the correct result", {
				con <- file(bed6File)
				gr <- loadBed( con, colnames = colnames )
				it( "Retuns a Genomic Ranges object", {
					expect_true( is( gr, "GenomicRanges" ))
				})
				it( "has location columns matched to the bed file", {
					expect_equal( as.character(seqnames(gr)), bed6df$seqnames )
					expect_equal( as.integer(start(gr)),    bed6df$start + 1 )
					expect_equal( as.integer(end(gr)),      bed6df$end )
				})
				it( "has specified values for strandedness",  {
					expect_equal( as.character(strand(gr)), bed6df$strand )
				})
				it( "has the expected metadata columns", {
					meta <- mcols(gr)
					expect_equal( names(meta), c("score", "id") )
					expect_equal( mcols(gr)[,"score"], bed6df$score )
					expect_equal( mcols(gr)[,"id"], bed6df$id )
				})
			})
			
		})
	})
	describe("Can load a near-bed file with a header and with one based start values", {
		
		# Alt format bed-like file - 1 based starts, header, and space delimited
		bedLikeDf <- bed6df
		bedLikeDf$start <- bedLikeDf$start + 1
		bedLikeFile <- tempfile( fileext = ".bed" )
		write.table( bedLikeDf, bedLikeFile, quote=FALSE, row.names = FALSE, col.names = TRUE,  sep=" ")
		
		describe("With metadata columns specified (including strandedness)", {
			colnames = names(bedLikeDf)
			it( "Runs silently if no error", {
				expect_silent( loadBed( bedLikeFile, colnames = colnames, oneBased= TRUE, skip= 1 ))
			})
			describe( "It returns the correct result", {
				gr <- loadBed( bedLikeFile, colnames = colnames, oneBased= TRUE, skip= 1 )
				it( "Retuns a Genomic Ranges object", {
					expect_true( is( gr, "GenomicRanges" ))
				})
				it( "has location columns matched to the bed file", {
					expect_equal( as.character(seqnames(gr)), bedLikeDf$seqnames )
					expect_equal( as.integer(start(gr)),    bedLikeDf$start )
					expect_equal( as.integer(end(gr)),      bedLikeDf$end )
				})
				it( "has specified values for strandedness",  {
					expect_equal( as.character(strand(gr)), bedLikeDf$strand )
				})
				it( "has the expected metadata columns", {
					meta <- mcols(gr)
					expect_equal( names(meta), c("score", "id") )
					expect_equal( mcols(gr)[,"score"], bedLikeDf$score )
					expect_equal( mcols(gr)[,"id"], bedLikeDf$id )
				})
			})
		})
		
	})
	describe("Loads a G-zipped bed transparently", {
		
		# gzipped file (6-col)
		bed6GzFile <- tempfile( fileext = ".bed.gz" )
		gzCon <- gzfile(bed6GzFile)
		write.table( bed6df, gzCon, quote=FALSE, row.names = FALSE, col.names = FALSE,  sep="\t")

		describe("With metadata columns specified (including strandedness)", {
			colnames = names(bed6df)
			it( "Runs silently if no error", {
				expect_silent( loadBed( bed6GzFile, colnames = colnames ))
			})
			describe( "It returns the correct result", {
				gr <- loadBed(  bed6GzFile, colnames = colnames )
				it( "Retuns a Genomic Ranges object", {
					expect_true( is( gr, "GenomicRanges" ))
				})
				it( "has location columns matched to the bed file", {
					expect_equal( as.character(seqnames(gr)), bed6df$seqnames )
					expect_equal( as.integer(start(gr)),    bed6df$start + 1 )
					expect_equal( as.integer(end(gr)),      bed6df$end )
				})
				it( "has specified values for strandedness",  {
					expect_equal( as.character(strand(gr)), bed6df$strand )
				})
				it( "has the expected metadata columns", {
					meta <- mcols(gr)
					expect_equal( names(meta), c("score", "id") )
					expect_equal( mcols(gr)[,"score"], bed6df$score )
					expect_equal( mcols(gr)[,"id"], bed6df$id )
				})
			})
		})
	})
	describe( "errors", {
		it( "Throws an error if file does not exist", {
			expect_false( file.exists( "noSuchFile" ))
			wantErrorRE <- "Error: Can't see the file 'noSuchFile'\\."
			expect_error( loadBed( "noSuchFile" ), wantErrorRE )
		})
	})
})

describe( "saveBed()", {
	describe( "Save a granges object with no meta-data to a bed-like file", {
		# Make minimal gr object
		grMin <- grPlus2
		mcols(grMin) <- NULL

		it( "Using default parameters, creates the correct file without error", {
			file <- tempfile(fileext = ".bed")
			expect_silent( saveBed( grMin, file ))
			got <- loadBed( file, colnames= names(bed6df)[1:4])
			want <- grMin
			expect_equal( got, want )
			unlink(file)
		})
		it( "Can exclude strand when creating a minimal bed file", {
			file <- tempfile(fileext = ".bed")
			saveBed( grMin, file, strand= FALSE )
			got <- loadBed( file, colnames= names(bed6df)[1:3])
			want <- grMin
			strand( want ) <- "*"
			expect_equal( got, want )
			unlink(file)
		})
	})
	describe( "Save a granges object with meta-data to a bed-like file", {
		it( "Using default parameters, creates the correct file without error", {
			file <- tempfile(fileext = ".bed")
			expect_silent( saveBed( grPlus2, file ))
			got <- loadBed( file, colnames= names(bed6df))
			want <- grPlus2
			expect_equal( got, want )
			unlink(file)
		})
		describe( "Can exclude strand", {
			file <- tempfile(fileext = ".bed")
			saveBed( grPlus2, file, strand= FALSE )
			got <- loadBed( file, colnames= names(bed6df)[-4])
			want <- grPlus2
			strand( want ) <- "*"
			expect_equal( got, want )
			unlink(file)
		})
	})
	describe( "Can save to gzipped connection transparently", {
		file <- tempfile(fileext = ".bed.is_gz_file")
		con <- gzfile( file )
		expect_silent( saveBed( grPlus2, con ))

		con <- file( file )
		expect_equal( summary(con)$class, "gzfile" )
		
		got <- loadBed( con, colnames= names(bed6df))
		want <- grPlus2
		expect_equal( got, want )
		unlink(file)
	})
	describe( "Can save to a tsv file with header and 1-s based starts", {
		it( "Using default parameters, creates the correct file without error", {
			file <- tempfile(fileext = ".tsv")
			expect_silent( saveBed( grPlus2, file, oneBased = TRUE, col.names = TRUE ))
			got <- loadBed( file, colnames= names(bed6df), oneBased = TRUE, skip= 1 )
			want <- grPlus2
			expect_equal( got, want )
			headers <- unlist(strsplit( readLines(file, n= 1), "\t", fixed = TRUE))
			expect_equal( headers, names( bed6df ))
			unlink(file)
		})
		it( "Can skip strand", {
			colNames <- names( bed6df )[-4]
			file <- tempfile(fileext = ".tsv")
			expect_silent( saveBed( grPlus2, file, strand= FALSE, oneBased = TRUE, col.names = TRUE ))
			got <- loadBed( file, colnames= colNames, oneBased = TRUE, skip= 1 )
			want <- grPlus2
			strand(want) <- "*" # Default if not set
			expect_equal( got, want )
			headers <- unlist(strsplit( readLines(file, n= 1), "\t", fixed = TRUE))
			expect_equal( headers, colNames )
			unlink(file)
		})
	})
})