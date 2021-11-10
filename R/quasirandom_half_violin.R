#' @noRd
#' @keywords internal
#' @importFrom graphics points
#' @export

# Taken from https://github.com/teunbrand/ggh4x/tree/master/R/utils.R
# Function for grabbing internal functions of ggplot2 that are also used here
.grab_ggplot_internals <- function() {
  objects <- c(
    "collide",
    "expand_range4",
    "pos_dodge"
  )

  objects <- stats::setNames(objects, objects)
  out <- lapply(objects, function(i) utils::getFromNamespace(i, "ggplot2"))
}

# Store the needed ggplot internals here
.beeint <- .grab_ggplot_internals()

# This is the primary modification made
crop_half <- function(df, right, nudge) {
  x <- NULL # to appease R CMD check
  if (right) {
    df <- dplyr::filter(df, x >= mean(df$x)) # right by default (as this --> top)
    if (!is.null(nudge)) df$x <- df$x + nudge
  }
  else {
    df <- dplyr::filter(df, x <= mean(df$x))
    if (!is.null(nudge)) df$x <- df$x - nudge
  }
  df
}

#### Main plotting functions ---------------------------------------------------
# Adapted from https://github.com/eclarke/ggbeeswarm (change is to combine a quasirandom
# and halfviolin plot)

position_quasirandom <- function(method = "quasirandom",
                                 width = NULL, varwidth = FALSE,
                                 bandwidth = 0.5, nbins = NULL,
                                 dodge.width = NULL, half = FALSE,
                                 right = NULL, nudge = NULL) {
  ggplot2::ggproto(NULL, PositionQuasirandom,
                   method = method,
                   width = width,
                   varwidth = varwidth,
                   bandwidth = bandwidth,
                   nbins = nbins,
                   dodge.width = dodge.width,
                   half = half,
                   right = right,
                   nudge = nudge
  )
}

PositionQuasirandom <- ggplot2::ggproto(
 "PositionQuasirandom",
 ggplot2::Position,
 required_aes = c("x", "y"),
 setup_params = function(self, data) {
   flipped_aes <- ggplot2::has_flipped_aes(data)
   data <- ggplot2::flip_data(data, flipped_aes)

   # get number of points in each x axis group and
   # find the largest group
   max.length <- max(data.frame(table(data$x))$Freq)

   list(
     method = self$method,
     width = self$width,
     varwidth = self$varwidth,
     bandwidth = self$bandwidth,
     nbins = self$nbins,
     max.length = max.length,
     dodge.width = self$dodge.width,
     half = self$half,
     right = self$right,
     nudge = self$nudge,
     flipped_aes = flipped_aes
   )
 },
 compute_panel = function(data, params, scales) {
   data <- ggplot2::flip_data(data, params$flipped_aes)

   # set width if not specified
   if (is.null(params$width)) {
     params$width <- ggplot2::resolution(
       data$x, zero = FALSE) * 0.05
   }

   data <- .beeint$collide(
     data,
     params$dodge.width,
     name = "position_quasirandom",
     strategy = .beeint$pos_dodge,
     check.width = FALSE
   )

   # split data.frame into list of data.frames
   if(!is.null(params$dodge.width)) {
     data <- split(data, data$group)
   } else {
     data <- split(data, data$x)
   }

   # perform swarming separately for each data.frame
   data <- lapply(
     data,
     pos_quasirandom,
     method = params$method,
     width = params$width,
     vary.width = params$varwidth,
     adjust = params$bandwidth,
     nbins = params$nbins,
     max.length = params$max.length
   )

   if (params$half) {
     data <-
       lapply(
         seq_along(data),
         function(i) crop_half(data[[i]],
                               params$right,
                               params$nudge)
       )
   }

   # recombine list of data.frames into one
   data <- Reduce(rbind, data)

   ggplot2::flip_data(data, params$flipped_aes)
 }
)

pos_quasirandom <- function(df, width = 0.4, vary.width = FALSE,
                            max.length = NULL,...) {
  x.offset <- vipor::aveWithArgs(
    df$y, df$x,
    FUN = vipor::offsetSingleGroup,
    maxLength = if (vary.width) {max.length} else {NULL},
    ...
  )

  x.offset <- x.offset * width
  df$x <- df$x + x.offset
  df
}

geom_quasirandom <- function(
  mapping = NULL,
  data = NULL,
  stat = "identity", ...,
  method = "quasirandom",
  width = NULL,
  varwidth = FALSE,
  bandwidth = 0.5, #
  nbins = NULL,
  dodge.width = NULL,
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  half = FALSE,
  right = TRUE,
  nudge = 0) {
    position <- position_quasirandom(
      method = method,
      width = width,
      varwidth = varwidth,
      bandwidth = bandwidth,
      nbins = nbins,
      dodge.width = dodge.width,
      half = half,
      right = right,
      nudge = nudge
    )

    ggplot2::layer(
      data = data,
      mapping = mapping,
      stat = stat,
      geom = ggplot2::GeomPoint,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(
        na.rm = na.rm,
        ...
      )
    )
  }
