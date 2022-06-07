#' @importFrom Rcpp evalCpp
#' @useDynLib cdt, .registration=TRUE
#' @noRd
NULL

#' Title
#'
#' Description
#' 
#' By default holes are removed, as are triangles that are outside the boundary
#' but inside the convex hull. Internally this is done by CGAL 'marking domains' by
#' nesting level. This corresponds to the even-odd rule for holes. If 'mark_domains' is set to `FALSE`
#' all triangles are returned (this can be handy for meshes). 
#' 
#' @param points xx
#' @param edges xx
#' @param mark_domains logical, holes are removed (see Details)
#' @return xx
#' @export
#'
#' @importFrom Rvcg vcgGetEdge
#' 
#' @examples
#' nsides <- 12L
#' angles <- seq(0, 2*pi, length.out = nsides+1L)[-1L]
#' points <- cbind(cos(angles), sin(angles))
#' points <- rbind(points, points/1.5)
#' # constraint edges
#' indices <- 1L:nsides
#' edges_outer <- cbind(
#'   indices, c(indices[-1L], indices[1L])
#' )
#' edges_inner <- edges_outer + nsides
#' edges <- rbind(edges_outer, edges_inner)
#' cdel(points, edges)
cdel <- function(points, edges, mark_domains = TRUE){
   del <- t(tridel(points, edges, mark_domains))
  rglfake <- list(vb = rbind(tpoints, 1), it = del)
  class(rglfake) <- "mesh3d"
  edges <- `colnames<-`(
    as.matrix(vcgGetEdge(rglfake))[, -3L], c("v1", "v2", "border")
  )
  edges
}


tridel <- function(points, edges, mark_domains = TRUE) {
    if(!is.matrix(points) || !is.numeric(points)){
    stop(
      "The `points` argument must be a numeric matrix with two or three columns.", 
      call. = TRUE
    )
  }
  dimension <- ncol(points)
  if(!is.element(dimension, c(2L, 3L))){
    stop(
      "The `points` argument must be a numeric matrix with two or three columns.", 
      call. = TRUE
    )
  }
  if(any(is.na(points))){
    stop("Missing values are not allowed.", call. = TRUE)
  }
  if(anyDuplicated(points)){
    stop("There are some duplicated points.", call. = TRUE)
  }
  storage.mode(points) <- "double"
  tpoints <- t(points)
  if(!is.matrix(edges) || !is.numeric(edges) || ncol(edges) != 2L){
    stop(
      "The `edges` argument must be an integer matrix with two columns.", 
      call. = TRUE
    )
  }
  if(any(is.na(edges))){
    stop("Missing values in `edges` are not allowed.", call. = TRUE)
  }
  storage.mode(edges) <- "integer"
  stopifnot(all(edges >= 1L))
  stopifnot(all(edges <= nrow(points)))

  ## bit faster than apply()
  edges <- cbind(pmin(edges[,1L], edges[,2L]), pmax(edges[,1L], edges[,2L]))
  if(anyDuplicated(edges)){
    stop("There are some duplicated constraint edges.", call. = TRUE)
  }
  if(any(edges[, 1L] == edges[, 2L])){
    stop("There are some invalid constraint edges.", call. = TRUE)
  }
  t(del2d_constrained_cpp(tpoints, t(edges), mark_domains))
}