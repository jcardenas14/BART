library(shiny)
library(shinydashboard)
library(shinyjs)
library(ggplot2)
library(RColorBrewer)
library(fastcluster)
library(NMF)
library(grid)
library(clValid)
library(VennDiagram)
library(gtools)
library(scales)
library(reshape2)
library(data.table)
library(tidyverse)
library(janitor)
library(DT)

source("R/heatmap.R")
source("R/helpers.R")
source("R/modules.R")

helpPopup <- function(title, content,
                      placement=c('right', 'top', 'left', 'bottom'),
                      trigger=c('click', 'hover', 'focus', 'manual')) {
  tagList(
    singleton(
      tags$head(
        tags$script("$(function() { $(\"[data-toggle='popover']\").popover(); })")
      )
    ),
    tags$a(
      href = "#", class = "btn btn-mini", `data-toggle` = "popover",`data-html`="true",
      title = title, `data-content` = content, `data-animation` = TRUE,
      `data-placement` = match.arg(placement, several.ok=TRUE)[1],
      `data-trigger` = match.arg(trigger, several.ok=TRUE)[1],
      
      tags$i(class="fa fa-question-circle")
    )
  )
}

infoPopup <- function(title, content,
                      placement=c('right', 'top', 'left', 'bottom'),
                      trigger=c('click', 'hover', 'focus', 'manual')) {
  tagList(
    singleton(
      tags$head(
        tags$script("$(function() { $(\"[data-toggle='popover']\").popover(); })")
      )
    ),
    tags$a(
      href = "#", class = "btn btn-mini", `data-toggle` = "popover",`data-html`="true",
      title = title, `data-content` = content, `data-animation` = TRUE,
      `data-placement` = match.arg(placement, several.ok=TRUE)[1],
      `data-trigger` = match.arg(trigger, several.ok=TRUE)[1],
      
      tags$i(class="fa fa-info-circle")
    )
  )
}

old <- setwd(tempdir())

mycorrection<-function(dat,alpha,correction){
  k<-1:length(dat)
  index<-alpha*k/length(dat)
  sort.dat<-dat[order(dat)]
  
  if(correction=="FDR") maxk <- ifelse(length(which(sort.dat<=index))==0,0,max(which(sort.dat<=index)))
  if(correction=="BONF") maxk <- ifelse(length(which(sort.dat<=(alpha/length(sort.dat))))==0,0,max(which(sort.dat<=(alpha/length(sort.dat)))))
  if(correction=="RAW") maxk <- max(which(sort.dat<=alpha))
  
  percentage <- maxk/length(k)
  return(percentage)
}

getNorm <- function(x, y, colname, id, mynames, index_sid, index_refvar,
                    ref_level, keep = TRUE) {
  dat <- x[which(x[, index_sid] == id), mynames]
  lev_name <- dat[dat[, 3] == ref_level, colname]
  if (keep == TRUE) {
    lev_norm <- y[, as.character(dat[, colname])] - y[, as.character(lev_name)]
  }
  if (keep == FALSE) {
    index <- setdiff(as.character(dat[, colname]), as.character(lev_name))
    lev_norm <- y[, index] - y[, as.character(lev_name)]
    if (length(index) == 1) {
      lev_norm <- data.frame(lev_norm)
      names(lev_norm) <- as.character(dat[, colname][dat[, colname] %in% index])
    }
  }
  return(lev_norm)
}

getNorm2 <- function(y, colnames_lev, colnames_lev_rm, keep = TRUE) {
  index_lev <- which(names(y) %in% colnames_lev)
  index_lev_rm <- which(names(y) %in% colnames_lev_rm)
  if (keep == FALSE) {
    normdat <- y[, index_lev_rm] - apply(y[, index_lev], 1, mean)
  }
  if (keep == TRUE) {
    normdat <- y[, c(index_lev, index_lev_rm)] - apply(y[, index_lev], 1, mean)
  }
  return(normdat)
}

manipulateData <- function(y, x, colname, norm.method = "mean", ref.var = NULL, 
                           ref.val = NULL, long = FALSE, subject.id = NULL, 
                           keep.ref = TRUE) {
  y <- y[, match(x[, colname], colnames(y), nomatch = 0)]
  x <- x[match(colnames(y), x[, colname], nomatch = 0), ]
  if (!is.null(ref.var) & !is.null(ref.val)) {
    if (!long) {
      y.norm <- y - apply(y[, x[, ref.var] == ref.val], 1, norm.method, 
                          na.rm = TRUE)
      x.norm <- x
      if (!keep.ref) {
        y.norm <- y.norm[, x[, ref.var] != ref.val]
        x.norm <- x[match(colnames(y.norm), x[, colname], nomatch = 0), ]
      }
    }
    if (long) {
      if (!is.null(subject.id)) {
        subjects.ref <- as.character(
          unique(x[, subject.id][x[, ref.var] == ref.val])
        )
        subjects.not.ref <- as.character(
          unique(x[, subject.id][x[, ref.var] != ref.val])
        )
        subjects <- intersect(subjects.ref, subjects.not.ref)
        y.norm <- list()
        for (i in 1:length(subjects)) {
          index <- x[, subject.id] %in% subjects[i]
          index.ref <- index & x[, ref.var] == ref.val
          y.norm[[i]] <- y[, index] - y[, index.ref]
        }
        y.norm <- data.frame(do.call("cbind", y.norm))
        x.norm <- x[match(colnames(y.norm), x[, colname], nomatch = 0), ]
        if (!keep.ref) {
          y.norm <- y.norm[, x.norm[, ref.var] != ref.val]
          x.norm <- x[match(colnames(y.norm), x[, colname], nomatch = 0), ]
        }
      }
      if (is.null(subject.id)) {
        return(
          warning("Must specify subject.id when long = TRUE")
        )
      }
    }
  } else {
    y.norm <- y - apply(y, 1, norm.method, na.rm = TRUE)
    x.norm <- x
  }
  normData <- list(exprs.norm = y.norm, design.norm = x.norm)
  return(normData)
}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid::grid.newpage()
    pushViewport(viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
