#' Intersect groups of ranges with common target range set
#'
#' @param x The grouped (query) ranges
#' @param y The common target (subject) ranges
#' @param group The groupings, must be a vector of values the same length
#'   as x (query ranges). Can be specified as a vector, as a list with one
#'   element, or as the name of a meta-data column in x (query ranges). See
#'   `simplify` for additional info.
#' @param simplify Set FALSE to return a GRangesList, by group value. By default
#'   simplify is TRUE and a concatenated GRanges object will be returned, with
#'   a metadata column giving the group value named based on the format of 'group'
#'   * vector: column will be named "group".
#'   * list: column will be named `names(group)[1]`.
#'   * column name: column will have the same name.
#' @return A GRanges object or GRangesList object giving the intersection of
#'   all 'y' (subject) ranges with each group of 'x' (query ranges) as specified
#'   by "group". If either is NA, returns NA. If either is empty, returns an
#'   empty GRangesList object, or an empty GRanges object if `simplify` is `TRUE`.
#'
#' @export
groupedIntersect <- function( x, y, group, simplify= TRUE ) {
	if (( length( x ) == 1 && is.na( x )) || ( length( y ) == 1 && is.na( y ))) {
		return( NA )
	}
	if (! methods::is(x, "GRanges") || ! methods::is(y, "GRanges")) {
		stop( "Don't know how to intersect anything except 'GRanges' objects.")
	}
	if (length(x) == 0 || length(y) == 0 ) {
		if (simplify) {
			return( GenomicRanges::GRanges() )
		}
		else {
			return( GenomicRanges::GRangesList() )
		}
	}
	if ( is.character( group ) && length(group) == 1 && group %in% names(GenomicRanges::mcols(x))) {
		groupName <- group
		group <- GenomicRanges::mcols(x)[, groupName]
	}
	else if (is.list(group)) {
		if ( length(group) != 1 || is.null(names(group))) {
			stop( "If 'group' is specified as a list, it must have exactly one element which must be named.")
		}
		groupName <- names(group)
		group <- group[[1]]
	}
	else {
		groupName <- "group"
	}
	if (length(group) != length(x)) {
		message <- "Length of 'group' must be same as length of 'x'."
		if (length(group) == 1) {
			message <- paste0(message, " Perhaps you specified an unknown data column?" )
		}
		stop( message )
	}
	result <- GenomicRanges::GRangesList(unlist(tapply( x, group, intersect, y)))
	if (simplify) {
		result <- unlist( result, recursive= TRUE, use.names= TRUE )
		GenomicRanges::mcols(result)[, groupName] <- names(result)
		names(result) <- NULL
	}
	result
}
