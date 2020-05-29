## Shiny Modules

# Filtering Modules
filterOptsUI <- function(id){
  # Create a namespace function using the provided id
  ns <- NS(id)
  condCall1 <- paste0("input['",ns("TSoptions"),"']")
  condCall2 <- paste0("input['",ns("showfc"),"'] == true")
  
  tagList(
    div(style = "display:inline-block", checkboxInput(ns("TSoptions"), strong("P-value/log2 FC filtering", style = "color:#456dae"), FALSE)),
    div(style = "display:inline-block", helpPopup("P-value/log2 FC filtering", "These options allow the user to filter on signifance threshold and fold change.",
                                                  trigger = "hover")),
    conditionalPanel(condition = condCall1,
                     selectInput(ns("correction_method"),"Multiple testing correction:",
                                 c("FDR","Bonferroni","Raw"),"Raw"),
                     numericInput(ns('alphalevel'),"Significance threshold:", min = 0, max = 1, value = 0.05, step = 0.025),
                     checkboxInput(inputId = ns("showfc"),label = strong("Filter on fold change"),value = FALSE),
                     conditionalPanel(condition = condCall2,
                                      numericInput(inputId = ns("fcval"), label = "Log2 FC cut off:",
                                                   min = 0, max = 10, value = 0, step = 0.25),
                                      selectInput(ns("sign"),"Log2 FC sign:", c("+", "-", "Both"), c("Both")))
    )
  )
}

filterOpts <- function(input, output, session, data, comparison, data.type = NULL, paloOrFirst = NULL, geneList = NULL){
  results <- data()
  sum.custom <- function(x){
    result <- sum(x, na.rm = TRUE)
    if(all(is.na(x))){
      result <- 0
    }
    return(result)
  }
  if(data.type %in% c("genes", "flow", "metab")){
    results <- results[,c(1,2,which(colnames(results) %in% paste0("Estimate.", comparison()) | colnames(results) %in% paste0("Test.statistic.", comparison()) | 
                                      colnames(results) %in% paste0("P.Value.", comparison()) | colnames(results) %in% paste0("FDR.P.Value.", comparison())))]
    bonf <- p.adjust(results[,5], method = "bonferroni")
    results <- do.call("cbind", list(results, Bonf = bonf))
    if(data.type == "genes"){
      colnames(results) <- c("Transcript.ID", "Gene.Symbol", "Log2FC", "Test.Statistic", "P.Value", "FDR", "Bonferroni")
    } else if(data.type == "flow"){
      results[,2] <- NULL
      colnames(results) <- c("Flow.Variables","Log2FC", "Test.Statistic", "P.Value", "FDR", "Bonferroni")
    } else if(data.type == "metab"){
      results[,2] <- NULL
      colnames(results) <- c("Metab.Variables","Log2FC", "Test.Statistic", "P.Value", "FDR", "Bonferroni")
    }
    results <- data.frame(results)
    results <- results %>%
      filter(
        case_when(input$correction_method == "Raw" ~ P.Value <= input$alphalevel,
                  input$correction_method == "FDR" ~ FDR <= input$alphalevel,
                  TRUE ~ Bonferroni <= input$alphalevel)
      ) %>%
      arrange(P.Value)
    if(input$showfc){
      results <- results %>%
        filter(
          case_when(input$sign == "Both" ~ abs(Log2FC) >= input$fcval,
                    input$sign == "+" ~ Log2FC >= input$fcval,
                    TRUE ~ Log2FC <= -input$fcval)
        )
    }
    if(nrow(results) == 0){
      if(data.type == "genes"){
        results <- data.frame("No Genes Present")
        names(results) <- "Transcript.ID"
      } else if(data.type == "flow"){
        results <- data.frame("No Flow Variables Present")
        names(results) <- "Flow.Variables"
      } else if(data.type == "metab"){
        results <- data.frame("No Metabolites Present")
        names(results) <- "Metab.Variables"
      }
    }
  } else if(data.type == "geneSets") {
    if(is.null(results)){return(NULL)}
    if(length(comparison()) < 1){return(NULL)}
    results <- results %>% filter(Comparison %in% comparison())
    first_comp <- results %>% filter(Comparison %in% comparison()[1])
    if(input$showfc){
      first_cut <- unique(first_comp %>% 
                            filter(
                              case_when(input$sign == "Both" ~ abs(Log2FC) >= input$fcval,
                                        input$sign == "+" ~ Log2FC >= input$fcval,
                                        TRUE ~ Log2FC <= -input$fcval)
                            ) %>% pull(pathway.name))
      palo_cut <- unique(results %>% 
                           filter(
                             case_when(input$sign == "Both" ~ abs(Log2FC) >= input$fcval,
                                       input$sign == "+" ~ Log2FC >= input$fcval,
                                       TRUE ~ Log2FC <= -input$fcval)
                           ) %>% pull(pathway.name))
      results <- results %>% 
        filter(
          case_when(length(comparison()) == 1 || (length(paloOrFirst()) > 0 & paloOrFirst() == 1) ~ pathway.name %in% palo_cut,
                    TRUE ~ pathway.name %in% first_cut)
        )
    }
    
    pval <- case_when(input$correction_method == "Raw" ~ "p.Value", input$correction_method == "FDR" ~ "FDR", TRUE ~ "Bonf")
    keep <- unique(results %>%
                     filter(
                       case_when(length(comparison()) == 1 || (length(paloOrFirst()) > 0 & paloOrFirst() == 1) ~ !!sym(pval) <= input$alphalevel,
                                 TRUE ~ (Comparison == comparison()[1] & !!sym(pval) <= input$alphalevel))
                     ) %>% pull(pathway.name))
    if(length(keep) == 0){
      return(NULL)
    }
    results <- results %>% filter(pathway.name %in% keep)
    results.temp <- vector("list", length(comparison()))
    for(i in 1:length(comparison())){
      results.temp[[i]] <- results %>% filter(Comparison %in% comparison()[i])
      if(i == 1){
        results.temp[[i]] <- results.temp[[i]] %>% arrange(desc(Log2FC))
      } else {
        results.temp[[i]] <- results.temp[[i]][match(results.temp[[1]]$pathway.name, results.temp[[i]]$pathway.name, nomatch = 0),]
      }
      Index <- 1:nrow(results.temp[[i]])
      results.temp[[i]] <- cbind(Index, results.temp[[i]])
    }
    results <- do.call("rbind", results.temp)
    results[,2] <- factor(results[,2], levels = unique(results[,2]))
    results$SIG = "Not Significant"
    results$SIG[which(results[[pval]] <= input$alphalevel)] = "Significant"
  } else if(data.type == "percents"){
    dataset <- data()
    index <- grep("^P.Value",colnames(dataset))
    p.names <- gsub("P.Value.","",colnames(dataset)[index])
    y <- as.matrix(dataset[,index, drop = FALSE])
    if(input$correction_method == "Raw"){
      dat <- y
    } else if(input$correction_method == "FDR"){
      dat <- as.matrix(dataset[,grep("FDR.P.Value", colnames(dataset)), drop = FALSE])
    } else{
      dat <- apply(y, 2, p.adjust, method="bonferroni")
    }
    dat[dat <= input$alphalevel] <- 2
    dat[dat != 2] <- 0
    dat[dat == 2] <- 1
    
    fcindex <- grep("Estimate",colnames(dataset),fixed=T)
    dataset.fcindex <- dataset[,fcindex, drop = FALSE]
    if(input$showfc == TRUE){
      dataset.fcindex[dataset.fcindex < input$fcval & dataset.fcindex > -input$fcval] <- 0
      dat <- sign(dataset.fcindex)*dat
      if(input$sign == "+"){
        dat[dat == -1] <- 0
      } else{
        dat[dat == 1] <- 0
      }
    } else {
      dat <- sign(dataset.fcindex)*dat
    }
    dat <- data.frame(cbind(Transcript.ID = as.character(dataset$Transcript.ID),dat))
    dat <- inner_join(dat, geneList(), by = "Transcript.ID")
    dat$Transcript.ID <- NULL
    count_matrix <- data.frame(dat %>% group_by(Module) %>% summarise_all(sum.custom))
    rownames(count_matrix) <- count_matrix$Module
    count_matrix$Module <- NULL
    freq <- data.frame(table(geneList()$Module))$Freq
    names(freq) <- data.frame(table(geneList()$Module))$Var1
    freq.matched <- freq[match(rownames(count_matrix), names(freq), nomatch = 0)]
    prop_matrix <- count_matrix/freq.matched
    colnames(prop_matrix) <- p.names
    prop_matrix2 <- prop_matrix
    prop_matrix2[prop_matrix < .1 & prop_matrix > -.1] <- 0
    results <- list(prop_matrix = prop_matrix, prop_matrix2 = prop_matrix2)
  }
  return(results)
}




## Create data to be plotted in heatmap

# Sample Subset and Ordering Modules
subsetAndOrderUI <- function(id){
  ns <- NS(id)
  condCall1 <- paste0("input['",ns("subsetCols"),"']")
  condCall2 <- paste0("input['",ns("orderCols"),"'] == true")
  
  tagList(
    div(style = "display:inline-block", checkboxInput(ns("subsetCols"),strong("Subset samples", style = "color:#456dae"),FALSE)),
    div(style = "display:inline-block", helpPopup("Subset samples", "These options allow the user to plot a subset of the samples within the expression data set.",
                                                  trigger = "hover")),
    conditionalPanel(condition = condCall1,
                     uiOutput(ns('subsetVars')),
                     uiOutput(ns('subsetVals'))
    ),
    br(),
    div(style = "display:inline-block", checkboxInput(ns('orderCols'), strong("Order samples", style = "color:#456dae"), FALSE)),
    div(style = "display:inline-block", helpPopup("Order samples", "These options allow the user to specify column ordering by column variable.",
                                                  trigger = "hover")),
    conditionalPanel(condition = condCall2,
                     uiOutput(ns('orderingVars'))
    )
  )
}

subsetAndOrderRenderUI <- function(input, output, session, des){
  return(
    tagList(
      output$subsetVars <- renderUI({
        ns <- session$ns
        selectizeInput(ns("subsetVars"),"Variable(s) to subset heatmap:",c(colnames(des())), selected = NULL,multiple = TRUE)
      }),
      output$subsetVals <- renderUI({
        ns <- session$ns
        if(is.null(input$subsetVars)){
          return(NULL)
        }
        subsetVals <- lapply(1:length(input$subsetVars), function(i) {
          selectizeInput(ns(paste0("subsetVals",i)),paste0("Select ",input$subsetVars[i], " values"), 
                         choices = unique(des()[,input$subsetVars[i]]), selected = unique(des()[,input$subsetVars[i]]), multiple = TRUE)
        })
        subsetVals <- do.call(tagList, subsetVals)
        return(subsetVals)
      }),
      output$orderingVars <- renderUI({
        ns <- session$ns
        selected <- colnames(des())[1]
        selectizeInput(ns("orderingVars"), "Select variables to order heatmap:", choices = colnames(des()), selected = selected, multiple = TRUE)
      })
    )
  )
}

subsetAndOrder <- function(input, output, session, des, data, sampleAnnot){
  x <- data()
  if(!is.null(sampleAnnot())){
    colnames(x) <- as.character(des()[,sampleAnnot()])
  }
  if(input$subsetCols){
    if(!is.null(input$subsetVars)){
      subsetVals <- keep <- list()
      design <- des()
      for(i in 1:length(input$subsetVars)){
        subsetVals[[i]] <- eval(parse(text = paste0("input$subsetVals",i)))
        keep[[i]] <- which(design[,input$subsetVars[i]] %in% subsetVals[[i]])
        design[,input$subsetVars[i]] <- as.character(design[,input$subsetVars[i]])
        design[,input$subsetVars[i]] <- factor(design[,input$subsetVars[i]], levels = unique(subsetVals[[i]]))
      }
      keep <- Reduce(intersect, keep)
      x <- x[,keep]
      design <- design[keep,]
      order.vars <- do.call(order, design[,input$subsetVars, drop = FALSE])
      x <- x[,order.vars]
      design <- design[order.vars,]
    }
    if(input$orderCols){
      order.vars <- do.call(order, design[,input$orderingVars, drop = FALSE])
      x.ord <- x[,order.vars]
      colAnnot <- design[order.vars,input$orderingVars,drop = FALSE]
      design <- design[order.vars,]
    } else {
      orderingVars <- input$subsetVars
      order.vars <- do.call(order, design[,orderingVars, drop = FALSE])
      x.ord <- x[,order.vars]
      colAnnot <- design[order.vars,orderingVars, drop = FALSE]
      if(is.null(input$subsetVars)){
        orderingVars <- colnames(des())[1]
        x.ord <- x[,order(design[,orderingVars])]
        colAnnot <- design[order(design[,orderingVars]),orderingVars,drop = FALSE]
      }
    }
  } else {
    if(input$orderCols){
      order.vars <- do.call(order, des()[,input$orderingVars, drop = FALSE])
      x.ord <- x[,order.vars]
      colAnnot <- des()[order.vars,input$orderingVars,drop=FALSE]
      design <- des()[order.vars,]
    } else {
      orderingVars <- colnames(des())[1]
      x.ord <- x[,order(des()[,orderingVars])]
      colAnnot <- des()[order(des()[,orderingVars]),orderingVars,drop = FALSE]
      design <- des()
    }
  }
  for(i in 1:ncol(colAnnot)){
    if(length(unique(colAnnot[,i])) > 10){
      if(is.numeric(colAnnot[,i])){
        colAnnot[,i] <- colAnnot[,i]
      } else{
        colAnnot[,i] <- as.numeric(as.factor(as.character(colAnnot[,i])))
      }
    }
  }
  return(list(dat = x.ord, colAnnot = colAnnot, design = design))
}

annColors <- function(input, output, session, data){
  colors <- colorRampPalette(c("navy","yellow","firebrick3"))(length(unique(data()[,1])))
  eval(parse(text=paste(colnames(data())[1],"=","colors",sep="")))
  eval(parse(text=paste("first_color=list(",colnames(data())[1],"=",colnames(data())[1],")")))
  return(first_color)
}


# Sample Clustering Modules
colClusterUI <- function(id){
  ns <- NS(id)
  condCall1 <- paste0("input['",ns("clusterOptions"),"']")
  condCall2 <- paste0("input['",ns("clusterGroups"),"'] == true")
  
  tagList(
    div(style = "display:inline-block", checkboxInput(ns('clusterOptions'), strong("Cluster samples", style = "color:#456dae"), FALSE)),
    div(style = "display:inline-block", helpPopup("Cluster column options", "These options allow the user to cluster the samples and label the groups they are clustered in.",
                                                  trigger = "hover")),
    conditionalPanel(condition = condCall1,
                     checkboxInput(ns("colCluster"),"Cluster samples (columns)",FALSE),
                     checkboxInput(ns('clusterGroups'), "Show cluster groups", FALSE),
                     conditionalPanel(condition = condCall2,
                                      div(style = "display:inline-block", uiOutput(ns('clusterCuts')))),
                     div(style = "display:inline-block", infoPopup("Cluster Cuts", "Optimal number calculated using Dunn's Index.",
                                                                   placement = "right", trigger = "hover")))
  )
}

colCluster <- function(input, output, session, des = NULL, data){
  x <- data()
  colAnnot <- NULL
  if(!is.null(des)){
    colAnnot <- des()
  }
  colddm <- NA
  dist <- dist(t(x))
  hcl <- fastcluster::hclust(dist)
  if(input$colCluster == TRUE){
    colddm <- as.dendrogram(hcl)
    Rowv <- rowMeans(t(x), na.rm = T)
    colddm <- reorder(colddm, Rowv)
  }
  if(input$clusterGroups == TRUE){
    clusters = cutree(hcl, input$clusterCuts)
    colAnnot$Clusters = as.character(clusters)
    if(length(unique(colAnnot$Clusters)) > 10){
      colAnnot$Clusters <- as.numeric(colAnnot$Clusters)
    }
  }
  d <- sapply(2:round(nrow(colAnnot)/2), function(y) clValid::dunn(dist, cutree(hcl,y)))
  opt_num <- which(d == max(d)) + 1
  output$clusterCuts <- renderUI({
    ns <- session$ns
    numericInput(ns('clusterCuts'), "Number of clusters", min = 2, value = opt_num, step = 1)
  })
  return(list(opt_num = opt_num, colddm = colddm, hcl = hcl, d = d, colAnnot = colAnnot))
}


#Range of Values Plotted Module

maxValuesUI <- function(id){
  ns <- NS(id)
  numericInput(ns("setCutOff"),"Max value on color key:",2)
}

maxValues <- function(input, output, session, data){
  x <- data()
  if(input$setCutOff !=0 ){
    cut1 <- as.numeric(input$setCutOff)
    x[x < -cut1] <- -cut1
    x[x > cut1] <- cut1
  }
  return(x)
}

# Cluster Association Modules

clusterAssociationUI <- function(id){
  ns <- NS(id)
  tagList(
    uiOutput(ns("selectGroup")),
    uiOutput(ns("clusterNumber"))
  )
}

clusterAssociationRenderUI <- function(input, output, session, data){
  tagList(
    output$selectGroup <- renderUI({
      ns <- session$ns
      selected <- data()[,1]
      selectizeInput(ns("selectGroup"), "Choose a Group for Association Test", choices = colnames(data()), selected = NULL, multiple = TRUE, options = list(maxItems = 1))
    }),
    output$clusterNumber <- renderUI({
      ns <- session$ns
      numericInput(ns("clusterNumber"), "Number of clusters:", min = 2, value = 2, step = 1)
    })
  )
}

clusterAssociation <- function(input, output, session, des, hclObj){
  if(is.null(input$selectGroup)){return(NULL)}
  clusters <- cutree(hclObj(), input$clusterNumber)
  tab <- aggregate(des()[[input$selectGroup]] ~ clusters, data = des(), FUN = table)
  tab <- as.data.frame(tab[,2])
  table1 <- cbind("Cluster" = rownames(tab), tab)
  low_exp_count <- sapply(1:length(tab), function(y) if(sum(tab[,y]) == 0) y)
  low_exp_count <- unlist(low_exp_count)
  table2 <- tab
  if(!is.null(low_exp_count)){
    table2 <- as.data.frame(tab[, -low_exp_count, drop = FALSE])
  }
  if(ncol(table2) > 1){
    chi_sqr <- chisq.test(table2)
    fish_exact <- fisher.test(table2, simulate.p.value = T, B = 10000)
    percents <- list()
    for(i in 1:nrow(table2)){
      percents[[i]] <- paste(table2[i,], " ", "(", round(100*(table2[i,]/sum(table2[i,])), 2), "%",")", sep = "")
    }
    percents <- do.call("rbind", percents)
    colnames(percents) <- colnames(table2)
    table2 <- cbind("Cluster" = rownames(table2), percents)
  }
  else{
    chi_sqr <- ""
    fish_exact <- ""
  }
  table1[,1] <- as.character(table1[,1])
  table1_total <- rbind(data.frame(table1[,1,drop = F], stringsAsFactors = FALSE), "Total")
  table1_colsums <- rbind(table1[,-1], colSums(table1[,-1]))
  table1 <- cbind(table1_total, table1_colsums)
  table1$Total <- rowSums(table1[,-1])
  return(list(table1 = table1, table2 = table2, chi_sqr = chi_sqr, fish_exact = fish_exact))
}


# Upload Transcripts Modules

uploadVarsUI <- function(id, varType){
  ns <- NS(id)
  condCall1 <- paste0("input['",ns("uploadVars"),"']")
  uploadType <- paste0("Upload ", varType, ":")
  uploadDescription <- paste0("Allows the user to provide their own list of ", varType, " (CSV) to plot. The CSV file should contain a single column named 'Transcript.ID' or 'Gene.Symbol', 
                               depending on whether the list provided are transcript ids or gene symbols") 
  
  tagList(
    div(style = "display:inline-block", checkboxInput(ns("uploadVars"), strong(uploadType, style = "color:#456dae"), FALSE)),
    div(style = "display:inline-block", infoPopup(uploadType, uploadDescription, placement = "right", trigger = "hover")),
    conditionalPanel(condition = condCall1,
                     fileInput(ns('varSelect'), '', multiple = FALSE, accept=c(".csv")),
                     div(style = "display:inline-block",checkboxInput(ns("rowCluster"), strong("Row cluster"), FALSE)),
                     div(style = "display:inline-block", infoPopup("Row cluster", 'If unchecked, the order of the probes will be exactly the same as the order
                                                                   of the probes you upload. If checked, the rows will cluster')))
  )
}

uploadVarsRowCluster <- function(input, output, session, data, dendro){
  if(input$uploadVars){
    varnames <- read.csv(input$varSelect$datapath, header = TRUE)
    rowAnnot <- NA
    labelRows <- NULL
    x <- data()[match(varnames[,1], rownames(data()), nomatch = 0), ]
    if(ncol(varnames) > 1){
      rowAnnot <- varnames[,-1,drop = FALSE]
      rowAnnot <- rowAnnot[match(rownames(x), varnames[,1], nomatch = 0),,drop = FALSE]
    }
    if(input$rowCluster){
      dist <- dist(x)
      hcl <- fastcluster::hclust(dist)
      ddm <- as.dendrogram(hcl)
      Rowv <- rowMeans(x, na.rm = TRUE)
      ddm <- reorder(ddm, Rowv)
    }
    if(!input$rowCluster){
      ddm <- NA
    }
  } else {
    x <- data()
    x <- x[order.dendrogram(dendro()),]
    rowAnnot <- NA
    ddm <- NA
    labelRows = NA
  }
  return(list(ddm = ddm, x = x, labelRows = labelRows, rowAnnot = rowAnnot))
}

# Heatmap Graphing Options

graphOptionsUI <- function(id, varType = "transcripts"){
  ns <- NS(id)
  conCall1 <- condCall1 <- paste0("input['",ns("graphOptions"),"']")
  
  if(varType == "modules"){
    return(tagList(
      div(style = "display:inline-block", checkboxInput(ns("graphOptions"), strong("Graphing options", style = "color:#456dae"), FALSE)),
      div(style = "display:inline-block", helpPopup(ns("Graphing options"), "These options allow the user to make adjustments to the plot (i.e. width and height).",
                                                    trigger = "hover")),
      conditionalPanel(condition = condCall1,
                       sliderInput(ns('graphWidth'), "Plot width", min = 400, max = 2000, value = 750, step = 25),
                       sliderInput(ns('graphHeight'), "Plot height", min = 400, max = 2000, value = 800, step = 25),
                       numericInput(ns('circleSize'),"Circle size (.25 to 4):",min=.25,max=4,value=1,step=.25),
                       numericInput(ns('fontSize'), "Font size:", min = 0, value = 8, step = 2),
                       numericInput(ns('legendSize'), "Legend size:", min = 1, max = 5, value = 1, step = .50),
                       numericInput(ns('treeHeight'), "Tree height:", min = 0, value = 50, step = 5),
                       numericInput(ns('plotResolution'), "Plot resolution (for downloaded plot):", min = 72, max = 400, value = 72, step = 1))
    ))
  }
  return(tagList(
    div(style = "display:inline-block", checkboxInput(ns("graphOptions"), strong("Graphing options", style = "color:#456dae"), FALSE)),
    div(style = "display:inline-block", helpPopup(ns("Graphing options"), "These options allow the user to make adjustments to the plot (i.e. width and height).",
                                                  trigger = "hover")),
    conditionalPanel(condition = condCall1,
                     sliderInput(ns('graphWidth'), "Plot width", min = 400, max = 2000, value = 750, step = 25),
                     sliderInput(ns('graphHeight'), "Plot height", min = 400, max = 2000, value = 650, step = 25),
                     numericInput(ns('fontSize'), "Font size:", min = 10, value = 10, step = 2),
                     numericInput(ns('legendSize'), "Legend size:", min = 1, max = 5, value = 1, step = .50),
                     numericInput(ns('treeHeight'), "Tree height:", min = 0, value = 50, step = 5),
                     numericInput(ns('plotResolution'), "Plot resolution (for downloaded plot):", min = 72, max = 400, value = 72, step = 1))
  ))
}

graphOptions <- function(input, output, session, varType = "transcripts"){
  width <- input$graphWidth
  height <- input$graphHeight
  fontSize <- input$fontSize
  legendSize <- input$legendSize
  treeHeight <- input$treeHeight
  resolution <- input$plotResolution
  if(varType == "modules"){
    circleSize <- input$circleSize
  } else {
    circleSize <- NA
  }
  return(list(width = width, height = height, fontSize = fontSize, legendSize = legendSize, treeHeight = treeHeight, resolution = resolution, circleSize = circleSize))
}
