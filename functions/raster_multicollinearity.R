raster_multicollinearity <- function (raster.stack, multicollinearity.cutoff = 0.7,
                                            select.variables = FALSE, sample.points = FALSE, nb.points = 10000,
                                            plot = FALSE, method = "pearson") 
{
  if (sample.points) {
    if (!is.numeric(nb.points)) {
      stop("nb.points must be a numeric value corresponding to the number of pixels to sample from raster.stack")
    }
    env.df <- sampleRandom(raster.stack, size = nb.points, 
                           na.rm = TRUE)
  }
  else {
    env.df <- getValues(raster.stack)
    if (any(is.na(env.df))) {
      env.df <- env.df[-unique(which(is.na(env.df), arr.ind = T)[, 
                                                                 1]), ]
    }
  }
  if (!is.numeric(multicollinearity.cutoff)) {
    stop("You must provide a numeric cutoff between 0 and 1 in multicollinearity.cutoff")
  }
  else if (multicollinearity.cutoff > 1 | multicollinearity.cutoff < 
           0) {
    stop("You must provide a numeric cutoff between 0 and 1 in multicollinearity.cutoff")
  }
  cor.matrix <- matrix(data = 0, nrow = nlayers(raster.stack), 
                       ncol = nlayers(raster.stack), dimnames = list(names(raster.stack), 
                                                                     names(raster.stack)))
  cor.matrix <- 1 - abs(stats::cor(env.df, method = method))
  dist.matrix <- stats::as.dist(cor.matrix)
  ahc <- stats::hclust(dist.matrix, method = "complete")
  groups <- stats::cutree(ahc, h = 1 - multicollinearity.cutoff)
  if (length(groups) == max(groups)) {
    message(paste("  - No multicollinearity detected in your data at threshold ", 
                  multicollinearity.cutoff, "\n", sep = ""))
    mc <- FALSE
  }
  else {
    mc <- TRUE
  }
  if (plot) {
    op <- par(no.readonly = TRUE)
    method_list <- c("pearson", "spearman", "kendall")
    ylabel <- c("Pearson's r", "Spearman's rho", "Kendall's tau")[which(method_list %in% method)]
    graphics::par(mar = c(5.1, 5.1, 4.1, 3.1))
    plot(ahc, hang = -1, xlab = "", ylab = paste0("Distance (1 - ", ylabel, ")"), 
         main = "", las = 1, sub = "", axes = F)
    graphics::axis(2, at = seq(0, 1, length = 6), las = 1)
    if (mc) {
      graphics::title(paste("Groups of intercorrelated variables at cutoff", 
                            multicollinearity.cutoff))
      par(xpd = T)
      rect.hclust(ahc, h = 1 - multicollinearity.cutoff)
    }
    else {
      graphics::title(paste("No intercorrelation among variables at cutoff", 
                            multicollinearity.cutoff))
    }
    par(op)
  }
  if (select.variables) {
    sel.vars <- NULL
    for (i in 1:max(groups)) {
      sel.vars <- c(sel.vars, sample(names(groups[groups == 
                                                    i]), 1))
    }
  }
  else {
    if (mc) {
      sel.vars <- list()
      for (i in groups) {
        sel.vars[[i]] <- names(groups)[groups == i]
      }
    }
    else {
      sel.vars <- names(raster.stack)
    }
  }
  return(sel.vars)
}
