
### 29 July 2026 Michael Greenacre, Marina Martínez-Álvaro
### Requires package easyCODA and vegan
packages <- c("easyCODA", "vegan")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

### ====================================================================================
### function to compute concordance measure (closeness to a 45 degree line through zero) 
### i.e., closeness to exact numerical agreement, or scaled concordance (after a slope)

CONCORD <- function(x1,x2, scaled=FALSE) {
  
  if (length(x1) != length(x2)) {
    stop("x1 and x2 must have the same length.")
  }
  
  # Scaled x2 when scaled = TRUE
  if (scaled) {
    beta <- sum(x1 * x2) / sum(x2^2) # regression coeff. for x1 = b* x2
    x2 <- beta * x2
  }
  
  # Lin's coefficient of concordance
  s1 <- sd(x1)
  s2 <- sd(x2)
  corr <- cor(x1,x2)
  concord <- 2*corr*s1*s2 / ((mean(x1) - mean(x2))^2 + s1^2 + s2^2)
  concord
}



### ====================================================================================
unsuper_subco <- function(data, nsub=100, propsub=1/3, weight=TRUE, chi=FALSE, seed=123)  {             
  
### Input checks
  data <- as.matrix(data)
  if (!is.numeric(data)) {
    stop("data must be numeric.")
  }
  
  if (any(!is.finite(data))) {
    stop("data contains missing or non-finite values.")
  }
  
  if (any(data < 0)) {
    stop("data must contain non-negative values.")
  }
  
  if (any(rowSums(data) == 0)) {
    stop("All rows must have a positive sum.")
  }
  
  if (nsub < 1 || nsub != as.integer(nsub)) {
    stop("nsub must be a positive integer.")
  }
  
  if (propsub <= 0 || propsub > 1) {
    stop("propsub must be greater than 0 and no greater than 1.")
  }
  
  if (is.null(colnames(data))) {
    colnames(data) <- paste0("Feature_", seq_len(ncol(data)))
  }
  

### Create output matrices
concord.self <- matrix(NA, nrow = nsub, ncol = ncol(data), dimnames = list(NULL, colnames(data))) # Feature coherence as values concordance
concord.self.chi <- matrix(NA, nrow = nsub, ncol = ncol(data), dimnames = list(NULL, colnames(data))) # Feature coherence as values concordance (with chi-normalization)
scaledconcord.self <- matrix(NA, nrow = nsub, ncol = ncol(data), dimnames = list(NULL, colnames(data))) #  Feature coherence as values scaled concordance
scaledconcord.self.chi <- matrix(NA, nrow = nsub, ncol = ncol(data), dimnames = list(NULL, colnames(data))) #  Feature coherence as values scaled concordance (with chi-normalization)
pearson.self <- matrix(NA, nrow = nsub, ncol = ncol(data), dimnames = list(NULL, colnames(data))) #  Feature coherence as values correlation 
pearson.self.chi <- matrix(NA, nrow = nsub, ncol = ncol(data), dimnames = list(NULL, colnames(data))) #  Feature coherence as values correlation(with chi-normalization)
procr.sub <- rep(NA,nsub)                                  # Procrustes correlations 
procr.chi.sub <- rep(NA,nsub)                              # Procrustes correlations (with chi-normalization)
concord.dist.sub <- rep(NA,nsub)                           # concordance of distances
concord.dist.chi.sub <- rep(NA,nsub)                        # concordance of distances (with chi-normalization)
concord.cor.sub <- rep(NA,nsub)                            # concordance of correlations (chi-normalization makes no difference)
concord.cov.sub <- rep(NA,nsub)                            # concordance of covariances
concord.cov.chi.sub <- rep(NA,nsub)                        # concordance of covariances (with chi-normalization)


### Check and close the data
data.pro <- easyCODA::CLOSE(data)
data.cm <- colMeans(data.pro) 

### Transform data
data.trans <- data.pro
if(chi) data.trans <- sweep(data.pro, 2, sqrt(data.cm), FUN="/")

### The following are needed for comparing with full compositions
### column principal coordinates for Procrustes
data.pro.PCA   <- easyCODA::PCA(data.pro, weight=FALSE)
if(chi) data.trans.PCA <- easyCODA::PCA(data.trans, weight=FALSE)
data.pro.cpc   <- data.pro.PCA$colpcoord
if(chi) data.trans.cpc <- data.trans.PCA$colpcoord

### distances between columns of data.pro and data.trans
data.pro.dist   <- as.matrix(dist(t(data.pro)))
if(chi) data.trans.dist <- as.matrix(dist(t(data.trans)))

### correlations between data, i.e. normalized by standard deviations (so chi-normalization makes no difference)
data.pro.corr <- cor(data.pro)

### covariances between data, using chi-square normalization on relative abundances
data.pro.cov <- cov(data.pro)
if(chi) data.trans.cov <- cov(data.trans)


### sampling will use weights proportional to average data relative abundances (in data.cm)
### -------------------------------------------------------------------------------------
### start of loop on nsub subcompositions and storing results 
### we use subcompositions with 1/3 the number of data
### sub.rand contains column numbers of selected subcompositional data 

set.seed(seed)
for(isub in 1:nsub) {	
  if(weight)  sub.rand <- sample(1:ncol(data.pro), round(propsub*ncol(data.pro)), prob=data.cm)   
  if(!weight) sub.rand <- sample(1:ncol(data.pro), round(propsub*ncol(data.pro)))

  # create the closed subcomposition
  data.sub.pro <- easyCODA::CLOSE(data.pro[,sub.rand, drop = FALSE])
  data.sub.cm <- colMeans(data.sub.pro)
  if(chi) data.sub.trans <- sweep(data.sub.pro, 2, sqrt(data.sub.cm), FUN="/")
 
  # Feature coherence
  common_cols <- intersect(colnames(data.sub.pro), colnames(data.pro))
  for (feature in common_cols) {
    x1 <- data.sub.pro[,feature]
    if(chi) x1_chi<- data.sub.trans[,feature]
    x2 <- data.pro[,feature]
    if(chi) x2_chi <- data.trans[,feature]
    
    concord.self[isub, feature] <- CONCORD(x1, x2) 
    if(chi) concord.self.chi[isub, feature] <- CONCORD(x1_chi, x2_chi) 
    scaledconcord.self[isub, feature] <- CONCORD(x1, x2, scaled=TRUE)
    if(chi) scaledconcord.self.chi[isub, feature] <- CONCORD(x1_chi, x2_chi, scaled=TRUE)
    pearson.self[isub, feature] <- cor(x1, x2)
    if(chi) pearson.self.chi[isub, feature] <- cor(x1_chi, x2_chi)
  }
 
  
  # Procrustes, omit last dimension corresponding to 0 variance  
  data.sub.pro.PCA <- easyCODA::PCA(data.sub.pro, weight=FALSE)
  if(chi) data.sub.trans.PCA <- easyCODA::PCA(data.sub.trans, weight=FALSE)
  data.sub.pro.cpc <- data.sub.pro.PCA$colpcoord
  if(chi) data.sub.trans.cpc <- data.sub.trans.PCA$colpcoord
  procr.sub[isub] <- vegan::protest(data.pro.cpc[sub.rand, -ncol(data.pro.cpc)], data.sub.pro.cpc[,-ncol(data.sub.pro.cpc)], permutations=0)$t0                              
  if(chi) procr.chi.sub[isub] <- vegan::protest(data.trans.cpc[sub.rand, -ncol(data.trans.cpc)], data.sub.trans.cpc[,-ncol(data.sub.trans.cpc)], permutations=0)$t0
  
  
  # Distances between features, and then is scaled concordance between compositional and subcompositional distances. As_dist considers only lower.tri
  data.sub.pro.dist <- as.matrix(dist(t(data.sub.pro)))
  if(chi) data.sub.trans.dist <- as.matrix(dist(t(data.sub.trans)))
  concord.dist.sub[isub] <- CONCORD(as.dist(data.pro.dist[sub.rand,sub.rand]), as.dist(data.sub.pro.dist), scaled=TRUE)
  if(chi) concord.dist.chi.sub[isub] <- CONCORD(as.dist(data.trans.dist[sub.rand,sub.rand]), as.dist(data.sub.trans.dist), scaled=TRUE)
  
  # Correlations, and then concordance between compositional and subcompositional correlations. As_dist considers only lower.tri. ( chi-normalization makes no difference)
  data.sub.pro.corr <- cor(data.sub.pro)
  concord.cor.sub[isub] <- CONCORD(as.dist(data.pro.corr[sub.rand,sub.rand]),as.dist(data.sub.pro.corr), scaled=FALSE)
  
  # Var and Covariances (with chi-square normalization), and then scaled concordance between var-covariances
  data.sub.pro.cov <- cov(data.sub.pro)
  if(chi) data.sub.trans.cov <- cov(data.sub.trans)
  
  data.pro.cov.vec <- as.numeric(data.pro.cov[sub.rand,sub.rand][lower.tri(data.pro.cov[sub.rand,sub.rand], diag = TRUE)]) # Vector var-cov in the composition
  if(chi) data.trans.cov.vec <- as.numeric(data.trans.cov[sub.rand,sub.rand][lower.tri(data.trans.cov[sub.rand,sub.rand], diag = TRUE)])
  data.sub.pro.cov.vec <- as.numeric(data.sub.pro.cov[lower.tri(data.sub.pro.cov, diag = TRUE)]) # Vector var-cov in the subcomposition
  if(chi) data.sub.trans.cov.vec <- as.numeric(data.sub.trans.cov[lower.tri(data.sub.trans.cov, diag = TRUE)])
  
  concord.cov.sub[isub] <- CONCORD(data.pro.cov.vec, data.sub.pro.cov.vec, scaled=TRUE)
  if(chi) concord.cov.chi.sub[isub] <- CONCORD(data.trans.cov.vec, data.sub.trans.cov.vec, scaled=TRUE)
  
}

### ---------end of loop-----------------------------------------------------------------
  
### Feature coherence summaries across the subcompositions in which
### each feature was selected

feature.coherence.summary <- data.frame(
  feature = colnames(data.pro),
  n_selected = colSums(!is.na(concord.self)),
  concordance = colMeans(
    concord.self,
    na.rm = TRUE
  ),
  scaled_concordance = colMeans(
    scaledconcord.self,
    na.rm = TRUE
  ),
  pearson_correlation = colMeans(
    pearson.self,
    na.rm = TRUE
  ),
  row.names = NULL
)

if (chi) {
  feature.coherence.summary.chi <- data.frame(
    feature = colnames(data.pro),
    n_selected = colSums(!is.na(concord.self.chi)),
    concordance = colMeans(
      concord.self.chi,
      na.rm = TRUE
    ),
    scaled_concordance = colMeans(
      scaledconcord.self.chi,
      na.rm = TRUE
    ),
    pearson_correlation = colMeans(
      pearson.self.chi,
      na.rm = TRUE
    ),
    row.names = NULL
  )
}

## Prepare the output:
out <- list(call=match.call(),
  feature.coherence = list(
    concordance = concord.self,
    scaled_concordance = scaledconcord.self,
    pearson_correlation = pearson.self,
    summary = feature.coherence.summary
  ),
  procr = procr.sub,
  concord.dist = concord.dist.sub,
  concord.cor = concord.cor.sub,
  concord.cov = concord.cov.sub,
  settings = list(
    nsub = nsub,
    propsub = propsub,
    weight = weight,
    transformation = if (chi) "chi-squared relative abundances" else "relative abundances",
    seed = seed
  )
)

if (chi) {
  out$feature.coherence$chi <- list(
    summary = feature.coherence.summary.chi,
    concordance = concord.self.chi,
    scaled_concordance = scaledconcord.self.chi,
    pearson_correlation = pearson.self.chi
  )
  out$procr.chi <- procr.chi.sub
  out$concord.dist.chi <- concord.dist.chi.sub
  out$concord.cov.chi <- concord.cov.chi.sub
}

class(out) <- "unsuper_subco"
return(out)


}


print.unsuper_subco <- function(object, digits = 5) {
  
  # Function to summarize one metric
  metric_summary <- function(x) {
    x <- x[is.finite(x)]
    c(
      Mean   = mean(x),
      Median = median(x),
      Min    = min(x),
      Max    = max(x)
    )
  }
  
  print_metrics <- function(metrics, title) {
    
    cat(title, "\n")
    cat("--------------------------------------------------------------------------------------------\n")
    cat(sprintf(
      "%-48s | %7s | %7s | %-21s\n",
      "Metric Descriptor", "Mean", "Median", "Range"
    ))
    cat("--------------------------------------------------------------------------------------------\n")
    
    
    for(i in seq_len(nrow(metrics))) {
      
      cat(sprintf(
        "%-48s | %7.*f | %7.*f | [%.*f, %.*f]\n",
        rownames(metrics)[i],
        digits, metrics[i,"Mean"],
        digits, metrics[i,"Median"],
        digits, metrics[i,"Min"],
        digits, metrics[i,"Max"]
      ))
    }
    
    cat("\n")
  }
  
  cat("\n")
  cat("=========================================\n")
  cat("SUBCOMPOSITION COHERENCE ANALYSIS SUMMARY\n")
  cat("=========================================\n\n")
  
  cat("Call:\n")
  print(object$call)
  
  cat("\n")
  cat("Number of tested subcompositions:", object$settings$nsub, "\n")
  cat("Proportion of parts sampled:     ", object$settings$propsub, "\n\n")
  
  ## Relative abundances
  metrics <- rbind(
    "Scaled Concordance of feature values"   = metric_summary(object$feature.coherence$summary$scaled_concordance),
    "Procrustes Correlation"                 = metric_summary(object$procr),
    "Scaled Concordance of pair-wise distances" = metric_summary(object$concord.dist),
    "Concordance of pair-wise correlations"  = metric_summary(object$concord.cor),
    "Scaled Concordance of pair-wise var-cov"= metric_summary(object$concord.cov)
  )
  
  print_metrics(metrics, "Relative abundances:")
  
  ## Chi-normalized abundances
  if (!is.null(object$feature.coherence$chi)) {
    
    metrics.chi <- rbind(
      "Scaled Concordance of feature values"   = metric_summary(object$feature.coherence$chi$summary$scaled_concordance),
      "Procrustes Correlation"                 = metric_summary(object$procr.chi),
      "Scaled Concordance of pair-wise distances" = metric_summary(object$concord.dist.chi),
      "Concordance of pair-wise correlations"  = metric_summary(object$concord.cor),
      "Scaled Concordance of pair-wise var-cov"= metric_summary(object$concord.cov.chi)
    )
    
    print_metrics(metrics.chi, "Chi-normalized relative abundances:")
  }
  
  invisible(object)
}