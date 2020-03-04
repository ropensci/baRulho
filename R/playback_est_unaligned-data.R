#' Extended selection table with re-recorded playbacks before alignment
#' 
#' The data contains a subset of the selections in the example data 'playback_est' but in this subset the re-recorded signals are not aligned in time with the corresponding reference signals (see \code{\link{spcc_align}} for more details on aligning signals). This data set is intended mostly for using as an example in \code{\link{spcc_align}}. The data contains recordings of \emph{Phaethornis longirostris} (Long-billed Hermit) songs from different song types (column 'signal.type') that were broadcast and re-recorded at 4 distances (1m, 5m, 10m, 15m, column 'distance'). The data was created by the function \code{\link[warbleR]{selection_table}} from the \code{\link{warbleR}} package.
#' @format Extended selection table object in the \code{\link{warbleR}} format, which contains annotations and acoustic data
#' 
#' @usage data(playback_est_unaligned)
#' 
#' @source Marcelo Araya-Salas 
"playback_est_unaligned"
