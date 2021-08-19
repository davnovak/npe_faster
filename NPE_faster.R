
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
  
  ## Convert to list of DiscreteDistribution objects ----
  dist <- apply(dist, 2, function(x) distr::DiscreteDistribution(1:k, x))
  names(dist) <- levels(annot)
  dist
}

npe <- function(
  hd, ld, annot, knn = NULL, k = NULL, normalise = FALSE, exclude_pops = c()
) {
  require(FNN)
  require(distr)
  require(distrEx)
  
  annot <- as.factor(annot)
  stopifnot(nrow(hd) == nrow(ld))
  stopifnot(!(is.null(knn) && is.null(k)))
  n <- nrow(hd)
  if (is.null(k))
    k <- ncol(knn)
  
  ## Compute distribution over the counts of like neighbours for each population ----
  dist_hd <- .get_likeness_distributions(hd, annot, knn, k, normalise) # in high dimension
  dist_ld <- .get_likeness_distributions(ld, annot, NULL, k, normalise) # in low dimension
  
  ## Exclude select populations ----
  exclude <- annot %in% exclude_pops
  dist_hd <- dist_hd[!exclude]
  dist_ld <- dist_ld[!exclude]
  
  ## Compute the sum of total variation distances between distributions ----
  vardist_by_pop <- mapply(distrEx::TotalVarDist, dist_hd, dist_ld)
  sum(vardist_by_pop)
}



