
# Convienience so don't have to remember this implementation detail
isSamePath <- function( testPath, path ) {
	return( normalizePath( testPath ) == normalizePath( path ))
}

# # Basic idea from  xfun::is_abs_path (Yihui Xie, MIT licensed, cran 2022-09-23)
# isAbsPath <- function( path ) {
# 	if (.Platform$OS.type == "unix") {
# 		return( grepl( "^[~/]", path ))
# 	}
# 	else {
# 		testPath <- file.path(".", path)
# 		return( ! isSamePath( testPath, path ))
# 	}
# }
# 
# isRelPath <- function( path ) {
# 	! isAbsPath( path )
# }
# 
