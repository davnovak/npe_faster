
# require: FNN, distr, distrEx, emdist

.get_likeness_distributions <- function(
  coords, annot, knn, k, normalise
) {
  ## Normalise data ----
  if (normalise) coords <- scale(coords)
  
  ## Find nearest neighbours ----
  if (is.null(knn)) knn <- FNN::knn.index(coords, k) else knn <- knn[, 1:k]
  
  ## Fill in NAs ----
  coords[is.na(coords)] <- 0
  
  ## Get counts of like neighbours for each point ----
  reference <- as.integer(annot)
  neighbours <- matrix(as.integer(annot[knn]), ncol = k)
  neighbours <- split(neighbours, row(neighbours))
  like <- mapply(function(x, y) sum(x == y), neighbours, reference)
  like <- as.vector(like)
  
  ## Split counts by population ----
  like_by_pop <- split(like, reference)
  
  ## Compute distribution over these counts for each population ----
  dist <- do.call(cbind, lapply(like_by_pop, function(x) as.vector(table(factor(x[x != 0], levels = 1:k)))))
  dist[1, colSums(dist) == 0] <- 1
  dist <- apply(dist, 2, function(x) x / sum(x))
  
  ## Convert to list of vectors ----
  dist <- split(dist, col(dist))
  names(dist) <- levels(annot)
  dist
}

npe <- function(
  hd, ld, annot, knn = NULL, k = NULL, normalise = FALSE, exclude_pops = c(), metric = 'total_variation_distance', reduce = 'sum', plot_path = NULL
) {
  
  ## Prepare inputs ----
  annot <- as.factor(annot)
  stopifnot(nrow(hd) == nrow(ld))
  stopifnot(!(is.null(knn) && is.null(k)))
  n <- nrow(hd)
  if (is.null(k))
    k <- ncol(knn)
  
  ## Resolve distance function between likeness distributions in HD and LD ----
  metric <- match.arg(arg = metric, choices = c(
    'total_variation_distance', 'tvd',
    'earth_movers_distance', 'emd')
  )
  
  ## Resolve reduction function (whether to take sum or mean of dissimilarities between likeness distributions) ----
  reduce <- match.arg(arg = reduce, choices = c(
    'sum', 'mean', 'average', 'avg', 'none'
  ))
  redf <- if (reduce == 'sum') sum else if (reduce %in% c('mean', 'average', 'avg')) mean else if (reduce == 'none') function(x) x
  
  ## Compute distribution over the counts of like neighbours for each population ----
  dist_hd <- .get_likeness_distributions(hd, annot, knn, k, normalise) # in high dimension
  dist_ld <- .get_likeness_distributions(ld, annot, NULL, k, normalise) # in low dimension
  
  ## Exclude select populations ----
  exclude <- levels(annot) %in% exclude_pops
  dist_hd <- dist_hd[!exclude]
  dist_ld <- dist_ld[!exclude]
  
  ## Plot distributions per population ----
  if (!is.null(plot_path)) {
    if (tolower(substr(plot_path, nchar(plot_path) - 3, nchar(plot_path))) != '.pdf')
      plot_path <- paste0(plot_path, '.pdf')
    p <- par(no.readonly = TRUE)
    pdf(file = plot_path, width = 3*nlevels(annot), height = 6)
    par(mfrow = c(2, nlevels(annot)), mar = c(2, 2, 2, 2))
    bounds <- sapply(levels(annot), function(pop) max(dist_hd[[pop]], dist_ld[[pop]]))
    for (pop in levels(annot))
      barplot(dist_hd[[pop]], ylim = c(0, bounds[pop]), col = 'darkblue', border = 'darkblue', main = pop)
    for (pop in levels(annot))
      barplot(dist_ld[[pop]], ylim = c(0, bounds[pop]), col = 'darkred', border = 'darkred', main = pop)
    dev.off()
    do.call(par, p)
  }
  
  ## Change format of distributions to suit the selected distance function ----
  if (metric %in% c('total_variation_distance', 'tvd')) {
    format <- function(d) distr::DiscreteDistribution(1:k, d)
    distf <- distrEx::TotalVarDist
  } else if (metric %in% c('earth_movers_distance', 'emd')) {
    format <- function(d) matrix(c(d, 1:k), ncol = 2)
    distf <- emdist::emd
  }
  dist_hd <- lapply(dist_hd, format)
  dist_ld <- lapply(dist_ld, format)
  
  ## Compute dissimilarities between likeness distributions and reduce ----
  vardist_by_pop <- mapply(distf, dist_hd, dist_ld)
  redf(vardist_by_pop)
}



