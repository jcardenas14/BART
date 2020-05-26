aheatmap_circle <- function (x, color = "-RdYlBu2:100", breaks = NA, border_color = NA, 
                             circleSize = 1, cellwidth = NA, cellheight = NA, scale = "none", 
                             Rowv = TRUE, Colv = TRUE, revC = identical(Colv, "Rowv") || is_NA(Rowv) || 
                               (is.integer(Rowv) && length(Rowv) > 1) || is(Rowv, "silhouette"), 
                             distfun = "euclidean", hclustfun = "complete", 
                             reorderfun = function(d, w) reorder(d, w), treeheight = 50, 
                             legend = TRUE, annCol = NA, annRow = NA, annColors = NA, 
                             annLegend = TRUE, labRow = NULL, labCol = NULL, subsetRow = NULL, 
                             subsetCol = NULL, txt = NULL, fontsize = 10, cexRow = min(0.2 + 
                                                                                         1/log10(nr), 1.2), cexCol = min(0.2 + 1/log10(nc), 1.2), 
                             filename = NA, width = NA, height = NA, main = NULL, sub = NULL, 
                             info = NULL, verbose = getOption("verbose"), gp = gpar()) 
{
  ol <- lverbose(verbose)
  on.exit(lverbose(ol))
  if (is(x, "ExpressionSet")) {
    requireNamespace("Biobase")
    if (isTRUE(annCol)) 
      annCol <- atrack(x)
    x <- Biobase::exprs(x)
  }
  mat <- x
  if (!is.null(txt)) {
    if (!all(dim(mat), dim(x))) {
      stop("Incompatible data and text dimensions: arguments x and txt must have the same size.")
    }
  }
  res <- list()
  if (length(treeheight) == 1) 
    treeheight <- c(treeheight, treeheight)
  treeheight_row <- treeheight[1]
  treeheight_col <- treeheight[2]
  if (!is.null(subsetRow)) {
    if (verbose) 
      message("Compute row subset indexes")
    subsetRow <- subset_index(mat, 1L, subsetRow)
  }
  if (!is.null(subsetCol)) {
    if (verbose) 
      message("Compute column subset indexes")
    subsetCol <- subset_index(mat, 2L, subsetCol)
  }
  if (is.null(labRow) && is.null(rownames(mat))) 
    labRow <- 1L
  if (!is.null(labRow)) {
    if (verbose) 
      message("Process labRow")
    rownames(mat) <- generate_dimnames(labRow, nrow(mat), 
                                       rownames(mat))
  }
  if (is.null(labCol) && is.null(colnames(mat))) 
    labCol <- 1L
  if (!is.null(labCol)) {
    if (verbose) 
      message("Process labCol")
    colnames(mat) <- generate_dimnames(labCol, ncol(mat), 
                                       colnames(mat))
  }
  if (!is.null(subsetRow)) {
    mat <- mat[subsetRow, ]
  }
  if (!is.null(subsetCol)) {
    mat <- mat[, subsetCol]
  }
  tree_row <- if (!is_NA(Rowv)) {
    if (verbose) 
      message("Cluster rows")
    if (isReal(Rowv)) {
      treeheight_row <- Rowv
      Rowv <- TRUE
    }
    cluster_mat(mat, Rowv, distfun = distfun, hclustfun = hclustfun, 
                reorderfun = reorderfun, subset = subsetRow, verbose = verbose)
  }
  else NA
  if (identical(Rowv, FALSE) || !is_treedef(tree_row)) 
    treeheight_row <- 0
  tree_col <- if (!is_NA(Colv)) {
    if (identical(Colv, "Rowv")) {
      if (ncol(mat) != nrow(mat)) 
        stop("aheatmap - Colv='Rowv' but cannot treat columns and rows symmetrically: input matrix is not square.")
      treeheight_col <- treeheight_row
      tree_row
    }
    else {
      if (isReal(Colv)) {
        treeheight_col <- Colv
        Colv <- TRUE
      }
      if (verbose) 
        message("Cluster columns")
      cluster_mat(t(mat), Colv, distfun = distfun, hclustfun = hclustfun, 
                  reorderfun = reorderfun, subset = subsetCol, 
                  verbose = verbose)
    }
  }
  else NA
  if (identical(Colv, FALSE) || !is_treedef(tree_col)) 
    treeheight_col <- 0
  if (!is_NA(tree_row)) {
    if (revC) {
      if (verbose) 
        message("Reverse row clustering")
      tree_row <- rev(tree_row)
    }
    if (is_treedef(tree_row)) {
      res$Rowv <- tree_row$dendrogram
      res$rowInd <- order.dendrogram(tree_row$dendrogram)
      if (length(res$rowInd) != nrow(mat)) 
        stop("aheatmap - row dendrogram ordering gave index of wrong length (", 
             length(res$rowInd), ")")
    }
    else {
      res$rowInd <- tree_row
      tree_row <- NA
    }
  }
  else if (revC) {
    res$rowInd <- nrow(mat):1L
  }
  res$rowInd <- subset2orginal_idx(res$rowInd, subsetRow)
  if (!is.null(res$rowInd)) {
    if (!is.integer(res$rowInd) || length(res$rowInd) != 
        nrow(mat)) 
      stop("aheatmap - Invalid row ordering: should be an integer vector of length nrow(mat)=", 
           nrow(mat))
    if (verbose) 
      message("Order rows")
    subInd <- attr(res$rowInd, "subset")
    ri <- if (is.null(subInd)) 
      res$rowInd
    else subInd
    mat <- mat[ri, , drop = FALSE]
    if (!is.null(txt)) 
      txt <- txt[ri, , drop = FALSE]
  }
  if (!is_NA(tree_col)) {
    if (is_treedef(tree_col)) {
      res$Colv <- tree_col$dendrogram
      res$colInd <- order.dendrogram(tree_col$dendrogram)
      if (length(res$colInd) != ncol(mat)) 
        stop("aheatmap - column dendrogram ordering gave index of wrong length (", 
             length(res$colInd), ")")
    }
    else {
      res$colInd <- tree_col
      tree_col <- NA
    }
  }
  res$colInd <- subset2orginal_idx(res$colInd, subsetCol)
  if (!is.null(res$colInd)) {
    if (!is.integer(res$colInd) || length(res$colInd) != 
        ncol(mat)) 
      stop("aheatmap - Invalid column ordering: should be an integer vector of length ncol(mat)=", 
           ncol(mat))
    if (verbose) 
      message("Order columns")
    subInd <- attr(res$colInd, "subset")
    ci <- if (is.null(subInd)) 
      res$colInd
    else subInd
    mat <- mat[, ci, drop = FALSE]
    if (!is.null(txt)) 
      txt <- txt[, ci, drop = FALSE]
  }
  if (isTRUE(info) || is.character(info)) {
    if (verbose) 
      message("Compute info")
    if (!is.character(info)) 
      info <- NULL
    linfo <- NULL
    if (is_treedef(tree_row) && !is.null(tree_row$dist.method)) 
      linfo <- paste("rows:", tree_row$dist.method, 
                     "/", tree_row$method)
    if (is_treedef(tree_col) && !is.null(tree_col$dist.method)) 
      linfo <- paste(linfo, paste(" - cols:", tree_col$dist.method, 
                                  "/", tree_col$method))
    info <- c(info, linfo)
  }
  if (is_treedef(tree_col)) 
    tree_col <- tree_col$dendrogram
  if (is_treedef(tree_row)) 
    tree_row <- tree_row$dendrogram
  if (verbose) 
    message("Scale matrix")
  mat = as.matrix(mat)
  mat = scale_mat(mat, scale)
  color <- ccRamp(color)
  if (is_NA(breaks) || isNumber(breaks)) {
    if (verbose) 
      message("Generate breaks")
    cbreaks <- if (isNumber(breaks)) 
      breaks
    else NA
    breaks = generate_breaks(as.vector(mat), length(color), 
                             center = cbreaks)
  }
  if (isTRUE(legend)) {
    if (verbose) 
      message("Generate data legend breaks")
    legend = grid.pretty(range(as.vector(breaks)))
  }
  else {
    legend = NA
  }
  mat = scale_colours(mat, col = color, breaks = breaks)
  annotation_legend <- annLegend
  annotation_colors <- annColors
  annCol_processed <- atrack(annCol, order = res$colInd, .SPECIAL = specialAnnotation(2L), 
                             .DATA = amargin(x, 2L), .CACHE = annRow)
  annRow_processed <- atrack(annRow, order = res$rowInd, .SPECIAL = specialAnnotation(1L), 
                             .DATA = amargin(x, 1L), .CACHE = annCol)
  specialAnnotation(clear = TRUE)
  annTracks <- renderAnnotations(annCol_processed, annRow_processed, 
                                 annotation_colors = annotation_colors, verbose = verbose)
  nr <- nrow(mat)
  nc <- ncol(mat)
  res$vp <- heatmap_motor(mat, border_color = border_color, circleSize = circleSize,
                          cellwidth = cellwidth, cellheight = cellheight, treeheight_col = treeheight_col, 
                          treeheight_row = treeheight_row, tree_col = tree_col, 
                          tree_row = tree_row, filename = filename, width = width, 
                          height = height, breaks = breaks, color = color, legend = legend, 
                          annTracks = annTracks, annotation_legend = annotation_legend, 
                          txt = txt, fontsize = fontsize, fontsize_row = cexRow * 
                            fontsize, fontsize_col = cexCol * fontsize, main = main, 
                          sub = sub, info = info, verbose = verbose, gp = gp)
  invisible(res)
}

verbose <- NMF:::verbose
lverbose <- NMF:::lverbose
cluster_mat <- NMF:::cluster_mat
is_NA <- pkgmaker:::is_NA
isReal <- pkgmaker:::isReal
isNumber <- pkgmaker:::isNumber
grid.pretty <- grid:::grid.pretty
atrack <- NMF:::atrack
gpar <- grid:::gpar
convertUnit <- grid:::convertUnit
unit <- grid:::unit
grid.rect <- grid:::grid.rect
pushViewport <- grid:::pushViewport
viewport <- grid:::viewport
grid.text <- grid:::grid.text

renderAnnotations <- NMF:::renderAnnotations
generate_annotation_colours <- NMF:::generate_annotation_colours

ccPalette <- NMF:::ccPalette
ccSpec <-  NMF:::ccSpec

revPalette <- NMF:::revPalette


sequential_hcl <- colorspace:::sequential_hcl

hex <- colorspace:::hex

polarLUV <- colorspace:::polarLUV


round.pretty <- NMF:::round.pretty

convert_annotations <-  NMF:::convert_annotations


ccRamp <- NMF:::ccRamp


as_treedef <- NMF:::as_treedef


isLogical <- NMF:::isLogical

cutdendro <- NMF:::cutdendro


cutheight <- NMF:::cutheight


digest <- digest:::digest



subset_index <- NMF:::subset_index

generate_dimnames <- NMF:::generate_dimnames



is_treedef <- NMF:::is_treedef

subset2orginal_idx <-  NMF:::subset2orginal_idx

scale_mat <- NMF:::scale_mat
generate_breaks <-  NMF:::generate_breaks

scale_colours <- NMF:::scale_colours


scale_vec_colours <- NMF:::scale_vec_colours

amargin <- NMF:::amargin

adata <- NMF:::adata


specialAnnotation <-  NMF:::specialAnnotation

is.grob <- grid:::is.grob

textGrob <- grid:::textGrob

grid.draw <- grid:::grid.draw




heatmap_motor <-  function (matrix, border_color, circleSize, cellwidth, cellheight, tree_col, 
                            tree_row, treeheight_col, treeheight_row, filename = NA, 
                            width = NA, height = NA, breaks, color, legend, txt = NULL, 
                            annTracks, annotation_legend = TRUE, new = TRUE, fontsize, 
                            fontsize_row, fontsize_col, main = NULL, sub = NULL, info = NULL, 
                            verbose = getOption("verbose"), gp = gpar()) 
{
  annotation_colors <- annTracks$colors
  row_annotation <- annTracks$annRow
  annotation <- annTracks$annCol
  writeToFile <- !is.na(filename)
  if (writeToFile) {
    gfile(filename)
    on.exit(dev.off())
  }
  vpp <- current.vpPath_patched()
  if (is.null(vpp)) {
    if (verbose) 
      message("Detected path: [ROOT]")
    mf <- par("mfrow")
    new <- if (!identical(mf, c(1L, 1L))) {
      if (verbose) 
        message("Detected mfrow: ", mf[1], " - ", 
                mf[2], " ... MIXED")
      opar <- grid.base.mix(trace = verbose > 1)
      on.exit(grid.base.mix(opar))
      FALSE
    }
    else {
      if (verbose) {
        message("Detected mfrow: ", mf[1], " - ", 
                mf[2])
        message("Honouring ", if (missing(new)) 
          "default ", "argument `new=", new, 
          "` ... ", if (new) 
            "NEW"
          else "OVERLAY")
      }
      new
    }
  }
  else {
    if (verbose) 
      message("Detected path: ", vpp)
    if (missing(new)) {
      if (verbose) 
        message("Missing argument `new` ... OVERLAY")
      new <- FALSE
    }
    else if (verbose) 
      message("Honouring argument `new=", new, "` ... ", 
              if (new) 
                "NEW"
              else "OVERLAY")
  }
  if (new) {
    if (verbose) 
      message("Call: plot.new")
    plot.new()
  }
  mainGrob <- if (!is.null(main) && !is.grob(main)) 
    textGrob(main, gp = c_gpar(gp, fontsize = 1.2 * fontsize, 
                               fontface = "bold"))
  subGrob <- if (!is.null(sub) && !is.grob(sub)) 
    textGrob(sub, gp = c_gpar(gp, fontsize = 0.8 * fontsize))
  infoGrob <- if (!is.null(info) && !is.grob(info)) {
    grobTree(gList(rectGrob(gp = gpar(fill = "grey80")), 
                   textGrob(paste(info, collapse = " | "), x = unit(5, 
                                                                    "bigpts"), y = 0.5, just = "left", 
                            gp = c_gpar(gp, fontsize = 0.8 * fontsize))))
  }
  glo = lo(coln = colnames(matrix), rown = rownames(matrix), 
           nrow = nrow(matrix), ncol = ncol(matrix), cellwidth = cellwidth, 
           cellheight = cellheight, treeheight_col = treeheight_col, 
           treeheight_row = treeheight_row, legend = legend, annTracks = annTracks, 
           annotation_legend = annotation_legend, fontsize = fontsize, 
           fontsize_row = fontsize_row, fontsize_col = fontsize_col, 
           main = mainGrob, sub = subGrob, info = infoGrob, gp = gp)
  if (writeToFile) {
    if (verbose) 
      message("Compute size for file graphic device")
    m <- par("mar")
    if (is.na(height)) 
      height <- glo$height
    if (is.na(width)) 
      width <- glo$width
    dev.off()
    if (verbose) 
      message("Resize file graphic device to: ", 
              width, " - ", height)
    gfile(filename, width = width, height = height)
    if (new) {
      if (verbose) 
        message("Call again plot.new")
      op <- par(mar = c(0, 0, 0, 0))
      plot.new()
      par(op)
    }
    if (verbose) 
      message("Push again top viewport")
    pushViewport(glo$vp)
    if (verbose) 
      grid.rect(width = unit(glo$width, "inches"), 
                height = unit(glo$height, "inches"), gp = gpar(col = "blue"))
  }
  mindim <- glo$mindim
  if (mindim < 3) 
    border_color = NA
  if (!is_NA(tree_col) && treeheight_col != 0) {
    vplayout("ctree")
    draw_dendrogram(tree_col, horizontal = T)
    upViewport()
  }
  if (!is_NA(tree_row) && treeheight_row != 0) {
    vplayout("rtree")
    draw_dendrogram(tree_row, horizontal = F)
    upViewport()
  }
  fontsize_row <- convertUnit(min(unit(fontsize_row, "points"), 
                                  unit(0.6 * glo$cellheight, "bigpts")), "points")
  fontsize_col <- convertUnit(min(unit(fontsize_col, "points"), 
                                  unit(0.6 * glo$cellwidth, "bigpts")), "points")
  vplayout("mat")
  draw_matrix(matrix, circleSize = circleSize, border_color, txt = txt, gp = gpar(fontsize = fontsize_row))
  upViewport()
  if (length(colnames(matrix)) != 0) {
    vplayout("cnam")
    draw_colnames(colnames(matrix), gp = c_gpar(gp, fontsize = fontsize_col))
    upViewport()
  }
  if (length(rownames(matrix)) != 0) {
    vplayout("rnam")
    draw_rownames(rownames(matrix), gp = c_gpar(gp, fontsize = fontsize_row))
    upViewport()
  }
  if (!is_NA(annotation)) {
    vplayout("cann")
    draw_annotations(annotation, border_color)
    upViewport()
  }
  if (!is_NA(row_annotation)) {
    vplayout("rann")
    draw_annotations(row_annotation, border_color, horizontal = FALSE)
    upViewport()
  }
  if (annotation_legend && !is_NA(annotation_colors)) {
    vplayout("aleg")
    draw_annotation_legend(annotation_colors, border_color, 
                           gp = c_gpar(gp, fontsize = fontsize))
    upViewport()
  }
  if (!is_NA(legend)) {
    vplayout("leg")
    draw_legend(color, breaks, legend, gp = c_gpar(gp, fontsize = fontsize))
    upViewport()
  }
  if (!is.null(mainGrob)) {
    vplayout("main")
    grid.draw(mainGrob)
    upViewport()
  }
  if (!is.null(subGrob)) {
    vplayout("sub")
    grid.draw(subGrob)
    upViewport()
  }
  if (!is.null(infoGrob)) {
    vplayout("info")
    grid.draw(infoGrob)
    upViewport()
  }
  upViewport()
  NULL
}

current.vpPath_patched <-  NMF:::current.vpPath_patched
.use.grid.patch <- NMF:::.use.grid.patch


grid.base.mix <- NMF:::grid.base.mix

upViewport <- grid:::upViewport

current.vpPath <- grid:::current.vpPath 


grid.Call <- grid:::grid.Call

rootVP <- grid:::rootVP

vpPathFromVector <- grid:::vpPathFromVector

grid.Call.graphics <- grid:::grid.Call.graphics

record <- grid:::record
baseViewports <- gridBase:::baseViewports 

lo <- NMF:::lo

c_gpar <- NMF:::c_gpar

vplayout <- NMF:::vplayout

gfile <- NMF:::gfile


draw_dendrogram <-  NMF:::draw_dendrogram

draw_matrix <- function(matrix, circleSize = circleSize, border_color,txt = NULL, gp = gpar()) 
{ 
  n = nrow(matrix)
  m = ncol(matrix)
  x = (1:m)/m - 1/2/m
  y = (1:n)/n - 1/2/n
  if (!is.null(txt)) 
    txt[is.na(txt)] <- ""
  for (i in 1:m) {
    grid::grid.circle(x = x[i], y = y, r = .005*circleSize, default.units = "native", 
                      gp = gpar(fill = matrix[, i], col = border_color))
    if (!is.null(txt)) {
      grid.text(label = txt[, i], x = x[i], y = y, rot = 0, 
                check.overlap = FALSE, default.units = "npc", 
                gp = gp, )
    }
  }
}

draw_colnames <- NMF:::draw_colnames
draw_rownames <- NMF:::draw_rownames

draw_annotations <- NMF:::draw_annotations
draw_annotation_legend <- NMF:::draw_annotation_legend 
draw_legend <- NMF:::draw_legend  

gridPLT <- gridBase:::gridPLT
current.transform <- grid:::current.transform
