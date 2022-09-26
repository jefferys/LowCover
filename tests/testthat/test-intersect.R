
singleQ <- GenomicRanges::GRanges(
	seqnames= rep.int( "chr1", 5), 	strand="*",
	kind= rep.int( "S", 5),
	ranges= IRanges::IRanges(
		start= c( 11, 23, 34, 42, 51 ),
		end=   c( 15, 28, 37, 49, 53 )),
	tests= c( "l-over", "r-over", "inside", "outside", "distant" )
)

singleT <- GenomicRanges::GRanges(
	seqnames= rep.int( "chr1", 5), 	strand="*",
	ranges= IRanges::IRanges(
		start= c( 13, 21, 32, 44, 56 ),
		end=   c( 18, 25, 39, 47, 59 ))
)

singleI <- GenomicRanges::GRanges(
	seqnames= rep.int( "chr1", 4), 	strand="*",
	ranges= IRanges::IRanges(
		start= c( 13, 23, 34, 44 ),
		end=   c( 15, 25, 37, 47 )),
	kind= rep.int( "S", 4)
)

multiQ <- GenomicRanges::GRanges(
	seqnames= c( rep.int( "chr1", 3 ), rep.int( "chr2", 3 ),
				 rep.int( "chr3", 3 ), rep.int( "chr4", 3 ),
				 rep.int( "chr5", 3 ), rep.int( "chr6", 3 ), rep.int( "chr7", 4 )),
	strand= "*",
	ranges= IRanges::IRanges(
		start= c( c( 1, 2,  3 ), c( 1, 2, 3 ), c( 5, 4, 3 ), c( 3, 2,  1 ), c( 1, 5,  8 ), c( 3, 5, 7), c(3, 5, 7, 10 )),
		end=   c( c( 8, 9, 10 ), c( 5, 6, 7 ), c( 6, 7, 8 ), c( 8, 9, 10 ), c( 3, 6, 10 ), c( 4, 6, 8), c(3, 5, 7, 10 ))),
	kind= c( rep.int( "M", 3 ), rep.int( "M", 3 ),
			 rep.int( "M", 3 ), rep.int( "M", 3 ),
			 rep.int( "M", 3 ), rep.int( "M", 3 ), rep.int( "M", 4 ))
)
multiT <- GenomicRanges::GRanges(
	seqnames= c( rep.int( "chr1", 3 ), rep.int( "chr2", 3 ),
				 rep.int( "chr3", 3 ), rep.int( "chr4", 3 ),
				 rep.int( "chr5", 3 ), rep.int( "chr6", 3 ), rep.int( "chr7", 4 )),
	strand= "*",
	ranges= IRanges::IRanges(
		start= c( c( 5, 6,  7 ), c( 1, 2,  3 ), c(  1, 2, 3 ), c( 3, 4, 5 ), c( 3, 4, 5 ), c(  8, 5, 1 ), c( 1, 3, 6, 8 )),
		end=   c( c( 8, 9, 10 ), c( 8, 9, 10 ), c( 10, 9, 8 ), c( 8, 7, 6 ), c( 6, 7, 8 ), c( 10, 6, 3 ), c( 2, 4, 7, 9 )))
)

multiI <- GenomicRanges::GRanges(
	seqnames= c( "chr1", "chr2", "chr3", "chr4", rep.int( "chr5", 3 ), rep.int( "chr6", 3 ), rep.int( "chr7", 2 )),
	strand="*",
	ranges= IRanges::IRanges(
		start= c(  5, 1, 3, 3, c( 3, 5, 8 ), c( 3, 5, 8 ), c( 3, 7 )),
		end=   c( 10, 7, 8, 8, c( 3, 6, 8 ), c( 3, 6, 8 ), c( 3, 7 ))),
	kind= rep.int( "M", 12)
)

grNoneQ <- GenomicRanges::GRanges(
	seqnames= c( "chr1", "chr1", "chr2" ),
	kind= c("one", "one", "two"),
	ranges= IRanges::IRanges(
		start= c( 10, 16, 10 ),
		end=   c( 15, 20, 15 )
	)
)
grNoneT <- GenomicRanges::GRanges(
	seqnames= c( "chr1", "chr1", "chr2" ),
	ranges= IRanges::IRanges(
		start= c( 21, 1, 16 ),
		end=   c( 99, 9, 20 )
	)
)

describe( "groupedIntersect()", {
	describe( "With only one group", {
		describe( "and simplify= FALSE", {
			describe( "and group provided by name", {
				groupBy= "kind"
				it( "Runs quietly", {
					expect_silent( groupedIntersect( singleQ,  singleT, groupBy, simplify= FALSE ))
					expect_silent( groupedIntersect(  multiQ,   multiT, groupBy, simplify= FALSE ))
				})
				describe( "Returns correct list of GRanges objects", {
					gotSingGRL  <- groupedIntersect( singleQ,  singleT, groupBy, simplify= FALSE )
					gotMultiGRL <- groupedIntersect(  multiQ,   multiT, groupBy, simplify= FALSE )
					it( "Returns a GRangesList object", {
						expect_s4_class( gotSingGRL, "GRangesList" )
						expect_s4_class( gotMultiGRL, "GRangesList" )
					})
					describe( "Intersections are correctly identified", {
						it( "single overlaps are correct", {
							wantSingGR <- singleI
							GenomicRanges::mcols(wantSingGR)$kind <- NULL
							wantSingGRL <- GenomicRanges::GRangesList( S=wantSingGR )
							GenomicRanges::mcols(wantSingGRL)$kind <- NULL
							expect_equal( wantSingGRL, gotSingGRL )
						})
						it( "multiple overlaps are correct", {
							wantMultiGR <- multiI
							GenomicRanges::mcols(wantMultiGR)$kind <- NULL
							wantMultiGR <- GenomicRanges::GRangesList( M=wantMultiGR )
							GenomicRanges::mcols(wantMultiGR)$kind <- NULL
							expect_equal( wantMultiGR, gotMultiGRL )
						})
					})
				})
			})
			describe( "With group provided by list", {
				groupBySing= list( kind= GenomicRanges::mcols(singleQ)[, "kind"])
				groupByMulti= list( kind= GenomicRanges::mcols(multiQ)[, "kind"])
				it( "Runs quietly", {
					expect_silent( groupedIntersect( singleQ,  singleT,  groupBySing, simplify= FALSE ))
					expect_silent( groupedIntersect(  multiQ,   multiT, groupByMulti, simplify= FALSE ))
				})
				describe( "Returns correct list of GRanges objects", {
					gotSingGRL  <- groupedIntersect( singleQ,  singleT, groupBySing,  simplify= FALSE )
					gotMultiGRL <- groupedIntersect(  multiQ,   multiT, groupByMulti, simplify= FALSE )
					it( "Returns a GRangesList object", {
						expect_s4_class( gotSingGRL, "GRangesList" )
						expect_s4_class( gotMultiGRL, "GRangesList" )
					})
					describe( "Intersections are correctly identified", {
						it( "single overlaps are correct", {
							wantGR <- singleI
							GenomicRanges::mcols(wantGR)$kind <- NULL
							wantGRL <- GenomicRanges::GRangesList( S=wantGR )
							expect_equal( wantGRL, gotSingGRL )
						})
						it( "multiple overlaps are correct", {
							wantGR <- multiI
							GenomicRanges::mcols(wantGR)$kind <- NULL
							wantGRL <- GenomicRanges::GRangesList( M=wantGR )
							expect_equal( wantGRL, gotMultiGRL )
						})
					})
				})
			})
			describe( "With group provided by vector", {
				groupBySing=  GenomicRanges::mcols(singleQ)[, "kind"]
				groupByMulti= GenomicRanges::mcols(multiQ)[, "kind"]
				it( "Runs quietly", {
					expect_silent( groupedIntersect( singleQ,  singleT,  groupBySing, simplify= FALSE ))
					expect_silent( groupedIntersect(  multiQ,   multiT, groupByMulti, simplify= FALSE ))
				})
				describe( "Returns correct list of GRanges objects", {
					gotSingGRL <-  groupedIntersect( singleQ,  singleT, groupBySing, simplify= FALSE )
					gotMultiGRL <- groupedIntersect( multiQ,  multiT, groupByMulti, simplify= FALSE )
					it( "Returns a GRangesList object", {
						expect_s4_class( gotSingGRL, "GRangesList" )
						expect_s4_class( gotMultiGRL, "GRangesList" )
					})
					describe( "Intersections are correctly identified", {
						it( "single overlaps are correct", {
							wantGR <- singleI
							GenomicRanges::mcols(wantGR)$kind <- NULL
							wantGRL <- GenomicRanges::GRangesList( S=wantGR )
							expect_equal( wantGRL, gotSingGRL )
						})
						it( "multiple overlaps are correct", {
							wantGR <- multiI
							GenomicRanges::mcols(wantGR)$kind <- NULL
							wantGRL <- GenomicRanges::GRangesList( M=wantGR )
							expect_equal( wantGRL, gotMultiGRL )
						})
					})
				})
			})
		})
		describe( "With simplify= TRUE", {
			describe( "With group provided by name", {
				groupBy <- "kind"
				it( "Runs quietly", {
					expect_silent( groupedIntersect( singleQ,  singleT, groupBy, simplify= TRUE ))
					expect_silent( groupedIntersect(  multiQ,   multiT, groupBy, simplify= TRUE ))
				})
				describe( "Returns correct merged GRanges object", {
					gotSingGR <- groupedIntersect( singleQ,  singleT, groupBy, simplify= TRUE )
					gotMultiGR <- groupedIntersect( multiQ,  multiT, groupBy, simplify= TRUE )
					it( "Returns a GRanges object", {
						expect_s4_class( gotSingGR, "GRanges" )
						expect_s4_class( gotMultiGR, "GRanges" )
					})
					describe( "Intersections are correctly identified", {
						it( "single overlaps are correct", {
							wantGR <- singleI
							expect_equal( wantGR, gotSingGR )
						})
						it( "multiple overlaps are correct", {
							wantGR <- multiI
							expect_equal( wantGR, gotMultiGR )
						})
					})
				})
			})
			describe( "group provided by list", {
				groupBySing= list( kind= GenomicRanges::mcols(singleQ)[, "kind"])
				groupByMulti= list( kind= GenomicRanges::mcols(multiQ)[, "kind"])
				it( "Runs quietly", {
					expect_silent( groupedIntersect( singleQ,  singleT, groupBySing,  simplify= TRUE ))
					expect_silent( groupedIntersect(  multiQ,   multiT, groupByMulti, simplify= TRUE ))
				})
				describe( "Returns correct merged GRanges object", {
					gotSingGR  <- groupedIntersect( singleQ,  singleT,  groupBySing, simplify= TRUE )
					gotMultiGR <- groupedIntersect(  multiQ,   multiT, groupByMulti, simplify= TRUE )
					it( "Returns a GRanges object", {
						expect_s4_class( gotSingGR, "GRanges" )
						expect_s4_class( gotMultiGR, "GRanges" )
					})
					describe( "Intersections are correctly identified", {
						it( "single overlaps are correct", {
							wantGR <- singleI
							expect_equal( wantGR, gotSingGR )
						})
						it( "multiple overlaps are correct", {
							wantGR <- multiI
							expect_equal( wantGR, gotMultiGR )
						})
					})
				})
			})
			describe( "group provided by vector", {
				groupSingBy= GenomicRanges::mcols(singleQ)[, "kind"]
				groupMultiBy= GenomicRanges::mcols(multiQ)[, "kind"]
				it( "Runs quietly", {
					expect_silent( groupedIntersect( singleQ,  singleT,  groupSingBy, simplify= TRUE ))
					expect_silent( groupedIntersect(  multiQ,   multiT, groupMultiBy, simplify= TRUE ))
				})
				describe( "Returns correct merged GRanges object", {
					gotSingGR  <- groupedIntersect( singleQ,  singleT, groupSingBy, simplify= TRUE )
					gotMultiGR <- groupedIntersect(  multiQ,   multiT, groupMultiBy, simplify= TRUE )
					it( "Returns a GRanges object", {
						expect_s4_class( gotSingGR, "GRanges" )
						expect_s4_class( gotMultiGR, "GRanges" )
					})
					describe( "Intersections are correctly identified", {
						it( "single overlaps are correct", {
							wantGR <- singleI
							names(GenomicRanges::mcols(wantGR))[ names(GenomicRanges::mcols(wantGR)) == "kind" ] <- "group"
							expect_equal( wantGR, gotSingGR )
						})
						it( "multiple overlaps are correct", {
							wantGR <- multiI
							names(GenomicRanges::mcols(wantGR))[ names(GenomicRanges::mcols(wantGR)) == "kind" ] <- "group"
							expect_equal( wantGR, gotMultiGR )
						})
					})
				})
			})
		})
	})
	describe( "With multiple groups", {
		groupedQ <- c( singleQ, multiQ )
		groupedT <- c( singleT, multiT )
		describe( "and simplify= FALSE", {
			describe( "and group provided by name", {
				groupBy= "kind"
				it( "Runs quietly", {
					expect_silent( groupedIntersect( groupedQ,  groupedT, groupBy, simplify= FALSE ))
				})
				describe( "The list of GRanges objects returned", {
					got <- groupedIntersect( groupedQ,  groupedT, groupBy, simplify= FALSE )
					it( "Is a GRangesList object", {
						expect_s4_class( got, "GRangesList" )
					})
					it( "Reports the correct overlaps grouped correctly", {
						wantMI <- multiI
						GenomicRanges::mcols(wantMI)$kind <- NULL
						wantSI <- singleI
						GenomicRanges::mcols(wantSI)$kind <- NULL
						
						want <- GenomicRanges::GRangesList( M=wantMI, S=wantSI )
						expect_equal( want, got )
					})
				})
			})
			describe( "and group provided by list", {
				groupByList= list( kind= GenomicRanges::mcols(groupedQ)[, "kind"])
				it( "Runs quietly", {
					expect_silent( groupedIntersect(  groupedQ,  groupedT,  groupByList, simplify= FALSE ))
				})
				describe( "Returns correct list of GRanges objects", {
					got <- groupedIntersect( groupedQ,  groupedT, groupByList, simplify= FALSE )
					it( "Returns a GRangesList object", {
						expect_s4_class( got, "GRangesList" )
					})
					it( "Reports the correct overlaps grouped correctly", {
						wantMI <- multiI
						GenomicRanges::mcols(wantMI)$kind <- NULL
						wantSI <- singleI
						GenomicRanges::mcols(wantSI)$kind <- NULL
						
						want <- GenomicRanges::GRangesList( M=wantMI, S=wantSI )
						expect_equal( want, got )
					})
				})
			})
			describe( "group provided by vector", {
				groupByVec=  GenomicRanges::mcols(groupedQ)[, "kind"]
				it( "Runs quietly", {
					expect_silent( groupedIntersect(  groupedQ,  groupedT,  groupByVec, simplify= FALSE ))
				})
				describe( "Returns correct list of GRanges objects", {
					got <- groupedIntersect( groupedQ,  groupedT, groupByVec, simplify= FALSE )
					it( "Returns a GRangesList object", {
						expect_s4_class( got, "GRangesList" )
					})
					it( "Reports the correct overlaps grouped correctly", {
						wantMI <- multiI
						GenomicRanges::mcols(wantMI)$kind <- NULL
						wantSI <- singleI
						GenomicRanges::mcols(wantSI)$kind <- NULL
						
						want <- GenomicRanges::GRangesList( M=wantMI, S=wantSI )
						expect_equal( want, got )
					})
				})
			})
		})
		describe( "and simplify= TRUE (default)", {
			describe( "group provided by name", {
				groupBy= "kind"
				it( "Runs quietly", {
					expect_silent( groupedIntersect( groupedQ,  groupedT, groupBy ))
				})
				describe( "The GRanges object returned", {
					got <- groupedIntersect( groupedQ,  groupedT, groupBy, simplify= TRUE )
					it( "Is a GRanges object", {
						expect_s4_class( got, "GRanges" )
					})
					it( "Reports the correct overlaps grouped correctly", {
						want <- c( multiI, singleI )
						expect_equal( want, got )
					})
				})
			})
			describe( "group provided by list", {
				groupByList= list( kind= GenomicRanges::mcols(groupedQ)[, "kind"])
				it( "Runs quietly", {
					expect_silent( groupedIntersect(  groupedQ,  groupedT,  groupByList, simplify= TRUE ))
				})
				describe( "The list of GRanges objects returned", {
					got <- groupedIntersect( groupedQ,  groupedT, groupByList, simplify= TRUE )
					it( "Is a GRanges object", {
						expect_s4_class( got, "GRanges" )
					})
					it( "Reports the correct overlaps grouped correctly", {
						want <- c( multiI, singleI )
						expect_equal( want, got )
					})
				})
			})
			describe( "group provided by vector", {
				groupByVec=  GenomicRanges::mcols(groupedQ)[, "kind"]
				it( "Runs quietly", {
					expect_silent( groupedIntersect(  groupedQ,  groupedT,  groupByVec ))
				})
				describe( "The list of GRanges objects returned", {
					got <- groupedIntersect( groupedQ,  groupedT, groupByVec, simplify= TRUE )
					it( "Is a GRanges object", {
						expect_s4_class( got, "GRanges" )
					})
					it( "Reports the correct overlaps grouped correctly", {
						want <- c( multiI, singleI )
						names(GenomicRanges::mcols(want))[ names(GenomicRanges::mcols(want)) == "kind" ] <- "group"
						expect_equal( want, got )
					})
				})
			})
		})
	})
	describe( "Corner cases and error handling", {
		describe( "small GRanges objects", {
			emptyGR <- GenomicRanges::GRanges()
			emptyGRL <- GenomicRanges::GRangesList()
			it( "returns empty ganges if either or both inputs empty", {
				got <- groupedIntersect( singleQ, emptyGR, "kind",      simplify=  TRUE )
				expect_equal( got, emptyGR )
				got <- groupedIntersect( emptyGR, singleQ, list("A"=1), simplify=  TRUE )
				expect_equal( got, emptyGR )
				got <- groupedIntersect( emptyGR, emptyGR, c("A", "B"), simplify= FALSE )
				expect_equal( got, emptyGRL )
			})
			it( "returns NA if either or both inputs are NA", {
				expect_true(is.na( groupedIntersect( NA,      emptyGR, "kind",      simplify= TRUE  )))
				expect_true(is.na( groupedIntersect( emptyGR, NA,      list("A"=1), simplify= FALSE )))
				expect_true(is.na( groupedIntersect( NA,      NA,      c("A", "B"), simplify= FALSE )))
			})
			it( "returns empty ganges if no overlapping sequence", {
				grouping <- c( "one", "one", "two")
				got <- groupedIntersect( grNoneQ, grNoneT, "kind", simplify= TRUE )
				expect_true( length(got) == 0 )
				expect_s4_class( got, "GRanges" )
				got <- groupedIntersect( grNoneQ, grNoneT, "kind", simplify= FALSE )
				expect_true( length(got) == length(unique(grouping)))
				expect_true( all( sort( names( got )) == sort( unique( grouping ))))
				for ( x in names( got )) {
					
					expect_equal( length( got[[x]] ), 0, label= x )
					expect_s4_class( got[[x]], "GRanges" )
				}
			})
		})	
		describe( "Bad group information", {
			it( "returns an error if named group not found", {
				wantErrRE <- "Length of 'group' must be same as length of 'x'\\. Perhaps you specified an unknown data column\\?"
				expect_error( groupedIntersect( grNoneQ, grNoneT, group="noSuchGroup", simplify= FALSE ), wantErrRE )
				expect_error( groupedIntersect( grNoneQ, grNoneT, group="noSuchGroup", simplify= TRUE  ), wantErrRE )
			})
			it( "errors if group values length not same as GRanges query items", {
				wantErrRE <- "Length of 'group' must be same as length of 'x'\\."
				expect_error( groupedIntersect( grNoneQ, grNoneT, group= 1:5, simplify= FALSE ), wantErrRE )
				expect_error( groupedIntersect( grNoneQ, grNoneT, group= c("A", "B"), simplify= TRUE ), wantErrRE )
				expect_error( groupedIntersect( grNoneQ, grNoneT, group= list( "kind"= 1:5 ), simplify= FALSE ), wantErrRE )
				expect_error( groupedIntersect( grNoneQ, grNoneT, group= list( "bob"= c( "A", "B" )), simplify= TRUE ), wantErrRE )
			})
			it( "errors if group is provided as an unnamed list or a list with other than one element", {
				wantErrRE <-  "If 'group' is specified as a list, it must have exactly one element which must be named\\."
				expect_error( groupedIntersect( grNoneQ, grNoneT, list(), simplify= FALSE ), wantErrRE )
				expect_error( groupedIntersect( grNoneQ, grNoneT, list( x="a", y="b", z="c"), simplify= TRUE ), wantErrRE )
				expect_error( groupedIntersect( grNoneQ, grNoneT, list( c( "a", "b", "c" )), simplify= FALSE ), wantErrRE )
			})
		})
	})
})
