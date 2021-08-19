#### Improved Neighbourhood Proportion Error implementation
#### David Novak, 2021 (davidnovak@ugent.be)

# This script showcases an improved version of the Neighbourhood Proportion Error (NPE).
# It is faster, allows for re-using a previously generated k-NN graph and can exclude some
# populations (eg. unassigned cells) from the final result. Also, normalisation of data
# is not hard-coded.

# Original reference: Konstorum et al. (2018) https://doi.org/10.1101/273862
# Original authors: Anna Konstorum, Nathan Jekel (2019)

source('./NPE.R')
source('./NPE_faster.R')

library(HDCytoData)
library(uwot)

## Get 39-dimensional CyTOF dataset and preprocess it ----
data <- Samusik_01_SE()
exprs <- data@assays$data$exprs
exprs <- exprs[, data@colData$marker_class == 'type']
exprs <- asinh(exprs/5)
annot <- data@elementMetadata$population_id
rm(data)

## Get a sub-sample to speed things up ----
set.seed(1)
idcs_samp <- sample(x = 1:nrow(exprs), size = 5000)
exprs <- exprs[idcs_samp, ]
annot <- as.factor(as.character(annot[idcs_samp])) # re-level if populations vanished

## Create the 2-dimensional projection ----
proj <- umap(exprs, n_components = 2)

## Plot the projection, if interested ----
palette <- c('#1B9E77', '#D95F02', '#EBEBED', '#E7298A', '#66A61E', '#E6AB02', '#A6761D',
             '#EDD8D8', '#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99', '#E31A1C',
             '#FDBF6F', '#FF7F00', '#CAB2D6', '#6A3D9A', '#FFFF99', '#B15928', '#8DD3C7',
             '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5',
             '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F')
col <- palette[as.integer(annot)]
p <- par(no.readonly = TRUE) # save original graphical parameters
par(bg = 'black', mar = c(0.1, 0.1, 2, 0.1))
plot(proj, col = col, pch = 20, cex = 0.5, axes = FALSE, main = 'with unassigneds', col.main = 'white')
do.call(par, p)

## Load all required libraries beforehand to make the timing comparison fair ----
library(entropy)
library(FNN)
library(som)
library(topicmodels)
library(distrEx)
library(distr)
library(emdist)

## Compute NPE using the original function and time it ----

t0 <- system.time(
  res0 <- npe_error(
    sample_data = data.frame(exprs, as.integer(annot), annot),
    reduced_data = data.frame(proj, as.integer(annot), annot),
    population_key = data.frame(1:nlevels(annot), levels(annot)),
    k = 100
  )
)

t1 <- system.time(
  res1 <- npe(
    hd = exprs,
    ld = proj,
    annot = annot,
    k = 100,
    normalise = TRUE
  )
)

cat('The original function takes', t0['elapsed'], 'seconds\n')
cat('The improved function takes', t1['elapsed'], 'seconds\n')
if (res0 == res1) cat('Their results are the same\n')
