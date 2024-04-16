#' Detect channel combinations that potentially contain aggregates
#'
#' Channels will be subset to positive events only, since aggregates can only be detected in positive events. Pearson correlation is performed for every channel combination and a threshold is applied to determine which channel combinations may contain aggregates.
#'
#' @param seu Seurat object.
#' @param threshold Threshold for Pearson R. Channel combinations with a Pearson R higher than this threshold will be tagged potentially containing aggregates.
#' @param pos_threshold Treshold that distinguishes positive and negative events. Default: 1.
#' @param percent_cells Only tag channels as potentially containing aggregates if at least a certain percentage of all cells are double positive. This is too avoid high Pearson R in presence of very limited number of cells. Default: 0.5% (= 0.005)
#' @return List of correlation matrix with pearson R, correlation matrix with 1 for above threshold and 0 for below and data.frame with channel combinations that may contain aggregates.
#' @export
#' @examples
#' # identify channel combinations that potentially contain aggregates
#' problem_channels <- detect_aggregate_channels(seu)
#' # pheatmap can be used to quickly plot correlation matrices
#' library(pheatmap)
#' pheatmap(problem_channels[[1]])
#' pheatmap(problem_channels[[2]])
detect_aggregate_channels <- function(seu,
                                      threshold = 0.7,
                                      pos_threshold = 1,
                                      percent_cells = 0.005) {
    # function for correlation computation
    do_correlate <- function(column, row) {
        # get data.frame of only 2 markers
        marker <- data[, c(row, column)]
        # subset events to only double positive (> 1)
        marker <- as.matrix(marker[which(marker[, row] > pos_threshold & marker[, column] > pos_threshold), ])
        # only proceed if enough cells
        if(nrow(marker) > (nrow(data) * percent_cells)) {
             # write pearson coefficient into correlation matrix
             return(cor(marker[,1], marker[,2], method = "pearson"))
         }else{
             return(0)
         }
    }
    # get data from Seurat Object
    data <- t(GetAssayData(seu))
    # generate correlation matrix
    corr_mtx <- sapply(colnames(data), FUN = function(row) 
        sapply(colnames(data), FUN = function(column) do_correlate(column, row)))
    # name rows and columns
    row.names(corr_mtx) <- colnames(data)
    colnames(corr_mtx) <- colnames(data)    
    # make correlation matrix with 0 and 1 if threshold was passed
    corr_num <- ifelse(corr_mtx > threshold, 1, 0)
    # make list of channel combinations that are above threshold
    corr_logical <- corr_mtx > threshold
    # indices of correlations above threshold
    indices <- which(corr_logical, arr.ind = TRUE)
    # exclude combinations with the same variable on both sides
    indices <- indices[indices[, 1] != indices[, 2], ]
    if(nrow(indices) < 1) stop("No channels containing aggregates found. Either your dataset does not contain aggregates, or you need to adjust the parameters.", call. = FALSE)
    # sort the indices to ensure non-redundant pairs
    sorted_indices <- t(apply(indices, 1, sort))
    # remove duplicated indices (to exclude vice versa combinations)
    unique_indices <- unique(sorted_indices)
    # create a data frame with the variable pairs and their correlations
    result_df <- data.frame(
      Channel_1 = rownames(corr_mtx)[unique_indices[, 1]],
      Channel_2 = rownames(corr_mtx)[unique_indices[, 2]],
      Pearson_R = corr_mtx[unique_indices]
    )
    # make list of correlation matrix, correlation binary matrix and data.frame of channel combinations that may contain aggregates
    corr_all <- list(corr_mtx, corr_num, result_df)
    return(corr_all)
}


#' Modified Random Sample Consensus (RANSAC) to model aggregate events in channels
#'
#' Channels will be subset to positive events only, since aggregates can only be detected in positive events. A modified RANSAC is run that samples on events in the upper top right quadrant (events > 50% of fluorescence in both channels) and only allows linear models with a slope between 0.75 - 1.25 since aggregates follow an almost x = y model. Function returns the best fitting model.
#'
#' @param formula Formula for linear model containing dependent and independent variable.
#' @param data Data to run model on (events that pass the pos_threshold).
#' @param min_to_fit Minimum number of events to sample and fit the linear models. Default: 10.
#' @param max_iterations Number of iterations RANSAC should run.
#' @param fit_threshold Threshold for residuals (y - predicted) to be counted as inliers.
#' @param data_close Minimum data points counted as inliers to consider the model.
#' @param seed Set seed for reproducible results. Default: 42. 
#' @return Best fitting model for formula and data.
#' @export
run_modified_ransac <- function(formula,
                                data,
                                min_to_fit = 10,
                                max_iteration = 10,
                                fit_threshold = 0.5,
                                data_close = 10,
                                seed = 42) {
  # set seet
  set.seed(seed)
  # initial fit and error
  best_fit <- NULL
  best_error <- 1e6
  # variables from formula
  dependent_variable <- as.character(formula)[[2]]
  independent_variable <- as.character(formula)[[3]]
  # subset data to upper right quadrant; use this for sampling
  x50 <- (max(data[independent_variable]) - min(data[independent_variable])) / 2
  y50 <- (max(data[dependent_variable]) - min(data[dependent_variable])) / 2
  sample_dat <- data[which(data[independent_variable] > x50 & data[dependent_variable] > y50),]
  # iterate to generate models
  for (iteration in 1:max_iteration) {
    # take sample from sample data and build model
    possible_inliers <- sample(1:nrow(sample_dat), size = min_to_fit, replace = FALSE)
    initial_model <- lm(formula, data = data[possible_inliers, ])
    # calculate residuals and identify other inliers
    residuals <- abs(data[, dependent_variable] - predict(initial_model, newdata = data))
    also_inliers <- which(residuals < fit_threshold)
    # evaluate model and replace old model with new one if the sigma is better
    if (length(also_inliers) > data_close) {
      inliers <- c(possible_inliers, also_inliers)
      better_model <- lm(formula, data = data[inliers, ])
      slope <- summary(better_model)$coef[2,1]
      # only use model if slope is between 0.75 and 1.25
      if (slope > 0.75 & slope < 1.25) {
          this_error <- summary(better_model)$sigma
          if (this_error < best_error) {
            best_fit <- better_model
            best_error <- this_error
          }
      }
    }
  }
  if(is.null(best_fit)) warning("Could not find a suitable model. Try reduced min_to_fit or higher max_iteration.")
  return(best_fit)
}


#' Labels potential aggregates based on modified RANSAC model in Seurat object (seu$aggregate)
#'
#' Problematic channels that may contain aggregates can be identified with [detect_aggregate_channels()] function. This function runs a modified RANSAC algorithm (see [run_modified_ransac()]) on data of each channel combinations that may contain aggregates and scores the probability of events being aggregates. If cells are labelled as aggregates in several channels (score_threshold), cells are regarded as aggregates, which is stored in Seurat Object (seu$aggregate). The number of channels in which a cells is labelled as aggregate is stored in seu$aggregate_score.
#'
#' @param seu Seurat object.
#' @param aggregate_channels A data.frame with columns giving the channel combinations that potentially contain aggregates. This can be determine with [detect_aggregate_channels()] which returns a list. Element 3 of the list is this data.frame (aggregate_channels[[3]]).
#' @param score_threshold Number of channel combinations in which a cells should at be labeled as an aggregate (aggregate_score) to regard it as an aggregate. The higher the number, the less likely false positives will occur, but the more likely aggregates could be missed. Default: >=2.
#' @param fit_threshold Threshold for residuals (y - predicted) to be counted as aggregates in modified RANSAC model. Also see [run_modified_ransac()].
#' @param pos_threshold Treshold that distinguishes positive and negative events. Should be the value as provided for [detect_aggregate_channels()]. Default: 1.
#' @param ... Arguments that will be passed on to [run_modified_ransac()].
#' @return Seurat object with new cell-level meta.data "aggregate_score" and "aggregate".
#' @export
#' @examples
#' # detect aggregates in problematic channels
#' seu <- detect_aggregates(seu, problem_channels[[3]])
#' # remove aggregate events from Seurat object
#' seu <- subset(seu, subset = aggregate == FALSE)
detect_aggregates <- function(seu,
                              aggregate_channels,
                              score_threshold = 2,
                              fit_threshold = 0.25,
                              pos_threshold = 1,
                              ...) {
    # function for correlation computation
    do_ransac <- function(x, y) {
        # get data.frame of only 2 markers
        data_sub <- na.omit(data[, c(x, y)])
        # subset events to only double positive (> 1)
        data_sub <- as.data.frame(data_sub[which(data_sub[, x] > pos_threshold & data_sub[, y] > pos_threshold), ])
        # rename columns to avoid issues with LM in ransac function
        colnames(data_sub) <- c("x", "y")
        # only proceed if enough cells
        best_model <- suppressWarnings({run_modified_ransac(formula = formula("y ~ x"), data = data_sub, ...)})
        # get cells above fit_threshold of predicted value
        if(!is.null(best_model)) {
            # predict y values and calculate residuals
            data_sub$predicted <- predict(best_model, data_sub)
            data_sub$residuals <- abs(data_sub[, "y"] - data_sub$predicted)
            # get cells above threshold and return cell names
            potential_aggregates <- row.names(data_sub[which(data_sub$residuals < fit_threshold), ])
            return(potential_aggregates)
        }else{
            warning(paste0("Could not find suitable model for channel combination ", x, " vs ", y, ". ",
                           "Try reduced min_to_fit or higher max_iteration if this warning occurs in too many channel combinations."), call. = FALSE)
            return(NULL)
        }
    }
    # get data from Seurat Object
    data <- t(GetAssayData(seu))
    # run ransac for all channels and get a list of cell IDs
    potential_aggregates_list <- lapply(1:nrow(aggregate_channels), FUN = function(x) do_ransac(aggregate_channels[x, 1], aggregate_channels[x, 2]))
    potential_aggregates_df <- data.frame(table(unlist(potential_aggregates_list)))     
    # make new Seurat meta.data with aggregate scores
    seu$aggregate_score <- 0
    seu@meta.data[as.character(potential_aggregates_df[, 1]), "aggregate_score"] <- potential_aggregates_df[, 2]
    seu$aggregate <- seu$aggregate_score >= score_threshold
    # feedback on how many aggregate events found
    message(paste("Found",
                  sum(seu$aggregate),
                  "aggregate events, which is",
                  round((sum(seu$aggregate) / ncol(seu)) * 100, 2),
                  "percent of all cells."))
    return(seu)
}
