#' check_vec
#'
#' @description
#'  Run checks on input vector
#' @param nm_x name of the object to check
#' @param nm_y name of the reference object
#' @param nm_fcn optional name of the function the check is being run in
#' @param test_1 logical test function 1
#' @param test_2 logical test function 2
#'
#' @details
#' Used inside bamExtras functions to check if input vectors are numeric and are the same length as a reference object nm_y. Often vectors need to be numeric and the same length as a vector of ages, and the function shouldn't run unless these criteria are met.
#'
#' @keywords bam stock assessment fisheries
#' @author Nikolai ?warn
#' @export
#' @examples
#' aser <- rdat_BlackSeaBass$a.series
#' M <- aser$M
#' age <- aser$age
#' check_vec("M","age") # Should invisibly return TRUE
#'
#' # Generate an error
#' M <- M[-1]
#' check_vec("M","age",nm_fcn="pop_demo")
#'

# Check if a vector is numeric and of a certain length
# check_vec <- function(nm_x,
#                       nm_y,
#                       nm_fcn=NULL,
#                       test_1=function(x){is.numeric(x)},
#                       test_2=function(x,y){length(x)==length(y)}
#                       )
#                       {
#   x <- get(nm_x )
#   y <- get(nm_y)
#
#   res <- (test_1(x)&test_2(x,y))
#
#   if(!res){
#     a <- ifelse(!is.null(nm_fcn),paste(nm_fcn,"\n"),"")
#     stop(paste(a,nm_x,"must be a numeric vector the same length as",nm_y),call.=FALSE)
#   }
#   invisible(res)
# }

check_vec <- function(x,
                      y,
                      nm_fcn=NULL,
                      test_1=function(x){is.numeric(x)},
                      test_2=function(x,y){length(x)==length(y)}
)
{
  nm_x <- deparse(substitute(x))
  nm_y <- deparse(substitute(y))

  res <- (test_1(x)&test_2(x,y))

  if(!res){
    a <- ifelse(!is.null(nm_fcn),paste(nm_fcn,"\n"),"")
    stop(paste(a,nm_x,"must be a numeric vector the same length as",nm_y),call.=FALSE)
  }
  invisible(res)
}
