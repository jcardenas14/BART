options(shiny.maxRequestSize=1000*1024^2)
shinyServer(function(input,output, session){
  output$test<-renderUI({
    if(is.null(values$rowdend3)){
      if(ctrl()$id ==TRUE){
        try<-list(3,4)
        names(try) <- c(paste0("All Samples ",str_to_sentence(values$norm.method), " Normalized"),"All Samples Healthy Normalized")
      }
      if(ctrl()$id == FALSE){
        try<-list(3)
        names(try) <- paste0("All Samples ",str_to_sentence(values$norm.method), " Normalized")
      }
    }

    if(!is.null(values$rowdend3)){
      if(ctrl()$id == TRUE){
        try<-list(1,2,3,4,5)
        names(try) <- c(paste0("Baseline ",str_to_sentence(values$norm.method), " Normalized"),"Baseline Healthy Normalized",paste0("All Samples ",str_to_sentence(values$norm.method), " Normalized"),
                        "All Samples Healthy Normalized","All Samples Baseline Normalized")
      }

      if(ctrl()$id == FALSE){
        try<-list(1,3,5)
        names(try) <- c(paste0("Baseline ",str_to_sentence(values$norm.method), " Normalized"), paste0("All Samples ",str_to_sentence(values$norm.method), " Normalized"),"All Samples Baseline Normalized")
      }
    }
    selectizeInput("set", "Select heatmap:", choices = as.list(try), multiple = FALSE)
  })

  list.projects <- reactive({
    setwd(old)
    list.files("data")
  })
  
  observe({
    updateSelectizeInput(session,"exampleData", "Select example BART results file for Demo:", choices = list.projects(), selected = NULL)
  })
  
  output$uploadExampleData <- renderUI({
    req(input$exampleData)
    withBusyIndicatorUI(
      actionButton("uploadExampleData", "Upload Demo BART file", class = "btn-primary", style = "color: white; background-color: #363f47")
    )
  })

  fileindex <- reactive({
    req(input$file1)
    index <- which(input$file1$name == "bartResults.rda")
    index
  })

  values <- reactiveValues()

  updateData <- function(path, route){
    vars <- load(file = path, envir = .GlobalEnv)
    setwd(route)
    for (var in vars){
      values[[var]] <- get(var, .GlobalEnv)
    }
    return(values)
  }
  
  baylorMod <- reactive({
   if(is.null(values$dge.gsets)){
    baylorMod <- NULL
   } else{
    baylorMod <- input$baylorMods
   }
   baylorMod
  })
  
  unsupervisedBaylorMod <- reactive({
   baylorMod <- input$baylorModules
   baylorMod
  })

  observe({
    req(input$file1)
    mypath <- input$file1[[fileindex(), 'datapath']]
    updateData(mypath, tempdir())
  })
  
  observeEvent(input$uploadExampleData,{
   withBusyIndicatorServer("uploadExampleData", {
    setwd(old)
    mypath <- "data/Longitudinal Macaque TB Data/bartResults.rda"
    updateData(mypath, tempdir())
   })
  })
  
  design <- reactive({
    req(values$design)
    des <- values$design[[which(names(values$design) %in% c("microarray", "rnaseq"))]]
    flow <- values$design[[which(names(values$design) %in% c("flow"))]]
    metab <- values$design[[which(names(values$design) %in% c("metab"))]]
    z <- list(des = des, flow = flow, metab = metab)
    return(z)
  })
  
  sample.id <- reactive({
    req(values$sample.id)
    id <- values$sample.id[[which(names(values$sample.id) %in% c("microarray", "rnaseq"))]]
    flow <- values$sample.id[[which(names(values$sample.id) %in% c("flow"))]]
    metab <-  values$sample.id[[which(names(values$sample.id) %in% c("metab"))]]
    z <- list(id = id, flow = flow, metab = metab)
    return(z)
  })
  
  subject.id <- reactive({
    req(values$subject.id)
    id <- values$subject.id[[which(names(values$subject.id) %in% c("microarray", "rnaseq"))]]
    flow <- values$subject.id[[which(names(values$subject.id) %in% c("flow"))]]
    metab <- values$subject.id[[which(names(values$subject.id) %in% c("metab"))]]
    z <- list(id = id, flow = flow, metab = metab)
    return(z)
  })
  
  control.var <- reactive({
    req(values$control.var)
    id <- values$control.var[[which(names(values$control.var) %in% c("microarray", "rnaseq"))]]
    flow <- values$control.var[[which(names(values$control.var) %in% c("flow"))]]
    metab <- values$control.var[[which(names(values$control.var) %in% c("metab"))]]
    z <- list(id = id, flow = flow, metab = metab)
    return(z)
  })
  
  control.val <- reactive({
    req(values$control.val)
    id <- values$control.val[[which(names(values$control.val) %in% c("microarray", "rnaseq"))]]
    flow <- values$control.val[[which(names(values$control.val) %in% c("flow"))]]
    metab <- values$control.val[[which(names(values$control.val) %in% c("metab"))]]
    z <- list(id = id, flow = flow, metab = metab)
    return(z)
  })
  
  baseline.var <- reactive({
    req(values$baseline.var)
    id <- values$baseline.var[[which(names(values$baseline.var) %in% c("microarray", "rnaseq"))]]
    flow <- values$baseline.var[[which(names(values$baseline.var) %in% c("flow"))]]
    metab <- values$baseline.var[[which(names(values$baseline.var) %in% c("metab"))]]
    z <- list(id = id, flow = flow, metab = metab)
    return(z)
  })
  
  baseline.val <- reactive({
    req(values$baseline.val)
    id <- values$baseline.val[[which(names(values$baseline.val) %in% c("microarray", "rnaseq"))]]
    flow <- values$baseline.val[[which(names(values$baseline.val) %in% c("flow"))]]
    metab <- values$baseline.val[[which(names(values$baseline.val) %in% c("metab"))]]
    z <- list(id = id, flow = flow, metab = metab)
    return(z)
  })
  
  time.var <- reactive({
    if(is.null(values$time.var)){
     id <- flow <- metab <- NULL
    } else{
     id <- values$time.var[[which(names(values$time.var) %in% c("microarray", "rnaseq"))]]
     flow <- values$time.var[[which(names(values$time.var) %in% c("flow"))]]
     metab <- values$time.var[[which(names(values$time.var) %in% c("metab"))]]
    }
    z <- list(id = id, flow = flow, metab = metab)
    return(z)
  })
  
  ctrl <- reactive({
    if(is.null(values$control.var) || is.null(values$control.val)){return(NULL)}
    id <- flow <- metab <- FALSE
    if(!is.null(control.var()$id) & !is.null(control.val()$id)){
      id <- TRUE
    }
    if(!is.null(control.var()$flow) & !is.null(control.val()$flow)){
      flow <- TRUE
    }
    if(!is.null(control.var()$metab) & !is.null(control.val()$metab)){
      metab <- TRUE
    }
    z <- list(id = id, flow = flow, metab = metab)
    return(z)
  })

  output$versionbox <- renderValueBox({
    valueBox(
      paste0("Software Version: 1.0.0 \nRelease Date: TBD"), version$version.string, icon = icon("th-list"),
      color = "purple"
    )
  })

  output$version <- renderPrint({
    writeLines("Software Version: 1.0.0 \nRelease Date: TBD")
    writeLines(version$version.string)
  })

  output$projectName <- renderText({
    values$project.name
  })


  ############ Design File and Summary Statistics ##############################

  output$SumStat <- renderMenu({
    if(is.null(values$design)){
      return(strong(""))
    }
    if(is.null(design()$des) == FALSE){
      if(length(which(input$file1$name %in% grep(".html", input$file1$name, value = TRUE, fixed = TRUE))) == 0 &
         length(which(input$file1$name %in% grep(".csv", input$file1$name, value = TRUE, fixed = TRUE))) == 0 &
         length(which(input$file1$name %in% grep(".png", input$file1$name, value = TRUE, fixed = TRUE))) == 0){
        return(menuItem("Summary Stats and QC", icon = icon("table"), tabName = "summarystats",
                        menuSubItem("Design file", tabName = "design"),
                        menuSubItem("Summary Statistics", tabName = "summary")))
      }
      if(length(which(input$file1$name %in% grep(".png", input$file1$name, value = TRUE, fixed = TRUE))) > 0 &
         length(which(input$file1$name %in% grep(".html", input$file1$name, value = TRUE, fixed = TRUE))) == 0 &
         length(which(input$file1$name %in% grep(".csv", input$file1$name, value = TRUE, fixed = TRUE))) == 0){
        return(menuItem("Summary Stats and QC", icon = icon("table"), tabName = "summarystats",
                        menuSubItem("Design file", tabName = "design"),
                        menuSubItem("Summary Statistics", tabName = "summary")))
      }
      if(length(which(input$file1$name %in% grep(".png", input$file1$name, value = TRUE, fixed = TRUE))) > 0 &
         length(which(input$file1$name %in% grep(".html", input$file1$name, value = TRUE, fixed = TRUE))) > 0 &
         length(which(input$file1$name %in% grep(".csv", input$file1$name, value = TRUE, fixed = TRUE))) == 0){
        return(menuItem("Summary Stats and QC", icon = icon("table"), tabName = "summarystats",
                        menuSubItem("Design file", tabName = "design"),
                        menuSubItem("Summary Statistics", tabName = "summary"),
                        menuSubItem("Fast QC Report", tabName = "fastqc")))
      }
      if(length(which(input$file1$name %in% grep(".png", input$file1$name, value = TRUE, fixed = TRUE))) > 0 &
         length(which(input$file1$name %in% grep(".html", input$file1$name, value = TRUE, fixed = TRUE))) == 0 &
         length(which(input$file1$name %in% grep(".csv", input$file1$name, value = TRUE, fixed = TRUE))) > 0){
        return(menuItem("Summary Stats and QC", icon = icon("table"), tabName = "summarystats",
                        menuSubItem("Design file", tabName = "design"),
                        menuSubItem("Summary Statistics", tabName = "summary"),
                        menuSubItem("QC files upload", tabName = "qcupload")))
      }
      if(length(which(input$file1$name %in% grep(".png", input$file1$name, value = TRUE, fixed = TRUE))) == 0 &
         length(which(input$file1$name %in% grep(".html", input$file1$name, value = TRUE, fixed = TRUE))) > 0 &
         length(which(input$file1$name %in% grep(".csv", input$file1$name, value = TRUE, fixed = TRUE))) == 0){
        return(menuItem("Summary Stats and QC", icon = icon("table"), tabName = "summarystats",
                        menuSubItem("Design file", tabName = "design"),
                        menuSubItem("Summary Statistics", tabName = "summary"),
                        menuSubItem("Fast QC Report", tabName = "fastqc")))
      }
      if(length(which(input$file1$name %in% grep(".png", input$file1$name, value = TRUE, fixed = TRUE))) == 0 &
         length(which(input$file1$name %in% grep(".html", input$file1$name, value = TRUE, fixed = TRUE))) == 0 &
         length(which(input$file1$name %in% grep(".csv", input$file1$name, value = TRUE, fixed = TRUE))) > 0){
        return(menuItem("Summary Stats and QC", icon = icon("table"), tabName = "summarystats",
                        menuSubItem("Design file", tabName = "design"),
                        menuSubItem("Summary Statistics", tabName = "summary"),
                        menuSubItem("QC files upload", tabName = "qcupload")))
      }
      if(length(which(input$file1$name %in% grep(".png", input$file1$name, value = TRUE, fixed = TRUE))) == 0 &
         length(which(input$file1$name %in% grep(".html", input$file1$name, value = TRUE, fixed = TRUE))) > 0 &
         length(which(input$file1$name %in% grep(".csv", input$file1$name, value = TRUE, fixed = TRUE))) > 0){
        return(menuItem("Summary Stats and QC", icon = icon("table"), tabName = "summarystats",
                        menuSubItem("Design file", tabName = "design"),
                        menuSubItem("Summary Statistics", tabName = "summary"),
                        menuSubItem("Fast QC Report", tabName = "fastqc"),
                        menuSubItem("QC files upload", tabName = "qcupload")))
      }
      else{
        return(menuItem("Summary Stats and QC", icon = icon("table"), tabName = "summarystats",
                        menuSubItem("Design file", tabName = "design"),
                        menuSubItem("Summary Statistics", tabName = "summary"),
                        menuSubItem("Fast QC Report", tabName = "fastqc"),
                        menuSubItem("QC files upload", tabName = "qcupload")))
      }
    }
  })


  output$designTable <- renderDT({
    datatable(design()$des,
              class = "table-condensed",
              style = "bootstrap4",
              rownames = FALSE,
              options = list(autowidth = TRUE, scrollX = TRUE)
    )
  })

  output$downloadDesign <- downloadHandler(

    filename = function() {paste(values$project.name,'_Design','.csv', sep='')  },
    content = function(file) {
      write.csv(design()$des, file, row.names = FALSE)
    }
  )
  
  observe({
    req(designForSummaryStats())
    updateSelectizeInput(session,"summaryVar", "Main Variable for Summarization:", colnames(designForSummaryStats()), selected = NULL)
    updateSelectizeInput(session,"groupingVars", "Grouping variable(s)", colnames(designForSummaryStats()),selected = NULL)
  })
  
  summaryTable <- reactive({
    req(input$summaryVar)
    design <- designForSummaryStats()
    if(is.numeric(design[,input$summaryVar])){
      summaryTable <- design %>% group_by(!!!syms(input$groupingVars)) %>%
        summarise(Min = min(!!! syms(input$summaryVar), na.rm = TRUE), Q25 = quantile(!!! syms(input$summaryVar), na.rm = TRUE)[2], Median = median(!!! syms(input$summaryVar), na.rm = TRUE), 
                  Mean = mean(!!! syms(input$summaryVar), na.rm = TRUE), Q75 = quantile(!!! syms(input$summaryVar), na.rm = TRUE)[4], max = max(!!! syms(input$summaryVar), na.rm = TRUE),
                  SD = sd(!!! syms(input$summaryVar), na.rm = TRUE))
    } else {
      summaryTable <- design %>% 
        mutate_all(as.character) %>% 
        group_by(!!! syms(c(input$groupingVars,input$summaryVar))) %>% 
        summarise(N = n()) %>%
        pivot_wider(names_from = !! sym(input$summaryVar), values_from = N) %>%
        replace(is.na(.), 0) %>%
        adorn_totals("col","-",TRUE,"Total",everything())
      if(length(input$groupingVars) > 0){
        summaryTable <- summaryTable %>%
          adorn_totals("row")
      }
    }
    return(summaryTable)
  })

  output$downloadSummary0 <- downloadHandler(

    filename = function() {'SummaryStats_Table1.csv'},
    content = function(file) {
      write.csv(summaryTable(), file, row.names = FALSE)
    }
  )

  output$summaryText <- renderUI({
    req(input$summaryVar)
    table_description <- paste0(strong("Table: "), "Summary Statistics for ", input$summaryVar)
    if(length(input$groupingVars) > 0){
      for(i in 1:length(input$groupingVars)){
        if(i == 1){
          if(i+1 > length(input$groupingVars)){
            table_description <- paste0(table_description, " by ", input$groupingVars[i]) 
          } else{
            table_description <- paste0(table_description, " by ", input$groupingVars[i],", ")
          }
        } else{
          if(i+1 > length(input$groupingVars)){
            table_description <- paste0(table_description, input$groupingVars[i]) 
          } else{
            table_description <- paste0(table_description, input$groupingVars[i], ", ")
          }
        }
      }
    }
    HTML(table_description)
  })

  designForSummaryStats <- reactive({
    design <- design()$des
    time.var <- time.var()$id
    for(i in 1:ncol(design)){
     if(is.null(time.var)){
      if(is.factor(design[,i])){
       design[,i] <- as.character(design[,i])
      }
     } else{
      if(is.factor(design[,i]) || colnames(design)[i] == time.var){
       design[,i] <- as.character(design[,i])
      }
     }
      if(length(which(design[,i] == "")) > 0){
        design[,i][which(design[,i] == "")] <- NA
      }
      if(is.character(design[,i])){
        design[,i] <- as.factor(design[,i])
      }
    }
    return(design)
  })
  
  dig <- function(){
    as.numeric(input$digits)
  }
  
  output$summary0 <- renderUI({
    if(input$sortableTable){
      return(DTOutput("summaryTable"))
    } else{
      return(tableOutput("summaryTable2"))
    }
  })
  
  output$summaryTable <- renderDT({
    datatable(summaryTable(),
              class = "table-condensed",
              style = "bootstrap4",
              rownames = FALSE,
              options = list(autowidth = TRUE, scrollX = TRUE)
    )
  })
  
  output$summaryTable2 <- renderTable({
    summaryTable()
  },digits = dig, include.rownames = FALSE)

  ######################### Unsupervised Side Menu ###############################

output$Unsupervised <- renderMenu({
  if(is.null(values$exprs)){
      return(strong(""))
  }
  
  if(is.null(values$scores.base) & is.null(values$scores.ctrl)){
    return(menuItem("Unupervised Analysis", icon = icon("line-chart"), tabName = "unsupervised",
                    menuSubItem("Gene Level Heat Maps", tabName = "probeheatmap"))
    )
  }
  
  if(!is.null(values$scores.base) || !is.null(values$scores.ctrl)){
    return(menuItem("Unupervised Analysis", icon = icon("line-chart"), tabName = "unsupervised",
                    menuSubItem("Gene Level Heat Maps", tabName = "probeheatmap"),
                    menuItem("Module Maps", icon = icon("angle-double-right"), tabName = "moduleMap"))
    )
  }
})


  #################################### Module Maps ############################
  
  output$baseOrCtrl <- renderUI({
    if(is.null(values$scores.ctrl) & is.null(values$scores.base)){return(NULL)}
    if(!is.null(values$scores.ctrl) & !is.null(values$scores.base)){
      return(
        selectizeInput("baseOrCtrl", "Proportions With Respect to Baseline or Controls:", choices = c("With Respect to Baseline", "With Respect to Controls"))
      )
    }
    if(!is.null(values$scores.ctrl) & is.null(values$scores.base)){
      return(
        selectizeInput("baseOrCtrl", "Proportions With Respect to Baseline or Controls:", choices = c("With Respect to Controls"))
      )
    } 
    if(is.null(values$scores.ctrl) & !is.null(values$scores.base)){
      return(
        selectizeInput("baseOrCtrl", "Proportions With Respect to Baseline or Controls:", choices = c("With Respect to Baseline"))
      )
    } 
  })
  
  output$modSelection <- renderUI({
    if(!unsupervisedBaylorMod()){return(NULL)}
    selectizeInput("moduleSelection", "Module to include:", c("All", "Only Annotated", "First Round", "First Two Rounds", "First Three Rounds", "First Four Rounds", "First Five Rounds",
                                                              "First Six Rounds", "First Seven Rounds", "First Eight Rounds"), "First Six Rounds")
  })
  
  modData <- reactive({
    req(input$baseOrCtrl)
    if(input$baseOrCtrl == "With Respect to Baseline"){
      dat <- 100*(values$scores.base-1)
    } 
    if(input$baseOrCtrl == "With Respect to Controls"){
      dat <- 100*(values$scores.ctrl-1)
    } 
    if(unsupervisedBaylorMod()){
      if(input$moduleSelection == "Only Annotated"){
        dat <- dat[grep(" ", rownames(dat)),]
      } else if(input$moduleSelection == "First Round"){
        dat <- dat[grep("^M1", rownames(dat)),]
      } else if(input$moduleSelection == "First Two Rounds"){
        dat <- dat[grep("^M1|^M2", rownames(dat)),]
      } else if(input$moduleSelection == "First Three Rounds"){
        dat <- dat[grep("^M1|^M2|^M3", rownames(dat)),]
      } else if(input$moduleSelection == "First Four Rounds"){
        dat <- dat[grep("^M1|^M2|^M3|^M4", rownames(dat)),]
      } else if(input$moduleSelection == "First Five Rounds"){
        dat <- dat[grep("^M1|^M2|^M3|^M4|^M5", rownames(dat)),]
      } else if(input$moduleSelection == "First Six Rounds"){
        dat <- dat[grep("^M1|^M2|^M3|^M4|^M5|^M6", rownames(dat)),]
      } else if(input$moduleSelection == "First Seven Rounds"){
        dat <- dat[grep("^M1|^M2|^M3|^M4|^M5|^M6|^M7", rownames(dat)),]
      } else if(input$moduleSelection == "First Eight Rounds"){
        dat <- dat[grep("^M1|^M2|^M3|^M4|^M5|^M6|^M7|^M8", rownames(dat)),]
      } 
    }
    des <- design()$des
    dist <- dist(dat)
    hcl <- fastcluster::hclust(dist)
    ddm <- as.dendrogram(hcl)
    Rowv <- rowMeans(dat, na.rm = TRUE)
    ddm <- reorder(ddm, Rowv)
    des <- des[match(colnames(dat), des[,sample.id()$id], nomatch = 0),]
    x <- list(dat = dat, ddm = ddm, des = des)
    return(x)
  })
  
  modGraphParams <- eventReactive(input$goMod,{
    params <- callModule(graphOptions, "modGraphOptions", varType = "modules")
    params <- list(width = params$width, height = params$height, fontSize = params$fontSize, legendSize = params$legendSize, 
                   treeHeight = params$treeHeight, resolution = params$resolution, circleSize = params$circleSize)
    return(params)
  })
  
  callModule(subsetAndOrderRenderUI, "modSubOrder", des = reactive(design()$des))
  callModule(clusterAssociationRenderUI, "modAssociation", data = reactive(design()$des))
  
  modRowCluster <- eventReactive(input$goMod,{
    dat <- callModule(uploadVarsRowCluster, "mod", data = reactive(modData()$dat), dendro = reactive(modData()$ddm))
    x <- dat$x
    ddm <- dat$ddm
    labelRows <- dat$labelRows
    y <- list(x = x, ddm = ddm, labelRows = labelRows)
    return(y)
  })
  
  modClusterData <- eventReactive(input$goMod,{
    x <- callModule(colCluster, "modCluster", des = reactive(modOrderedData()$colAnnot), data = reactive(modOrderedData()$x))
    return(x)
  })
  
  modGen_clustTab <- reactive({
    x <- callModule(clusterAssociation, "modAssociation", des = reactive(modOrderedData()$design), hclObj = reactive(modClusterData()$hcl))
    return(x)
  })
  
  modOrderedData <- eventReactive(input$goMod,{
    dat <- callModule(subsetAndOrder, "modSubOrder", des = reactive(modData()$des), data = reactive(modRowCluster()$x), 
                      sampleAnnot = reactive(sample.id()$id))
    x <- dat$dat
    colAnnot <- dat$colAnnot
    design <- dat$design
    z <- list(x = x, colAnnot = colAnnot, design = design)
    return(z)
  })
  
  modColors <- eventReactive(input$goMod,{
    x <- callModule(annColors, "modSubOrder", reactive(modOrderedData()$colAnnot))
    return(x)
  })
  
  output$modHeatmap <- renderPlot({
    withProgress(message = 'Making plot',
                 detail = 'This may take a while...', value = 1,{
                   aheatmap_circle(modOrderedData()$x,Rowv = modRowCluster()$ddm, Colv = modClusterData()$colddm, circleSize = modGraphParams()$circleSize, treeheight = modGraphParams()$treeHeight, fontsize = modGraphParams()$fontSize, cexRow = 1.2, 
                             color = color.heatmap(),annCol = modClusterData()$colAnnot,annColors = modColors(),
                             breaks=seq(-100,100,by=4))
                 }
    )
  }, width = function(){modGraphParams()$width}, height = function(){modGraphParams()$height})
  
  output$downloadModPlot3Other <- downloadHandler(
    filename = function() {paste(values$project.name, '_Longitudinal_ModuleMap','.png', sep = '')},
    content = function(file){
      res <- modGraphParams()$resolution
      height <- modGraphParams()$height
      width <- modGraphParams()$width
      png(file, width = (res/72)*width, height = (res/72)*height, res = res)
      print(aheatmap_circle(modOrderedData()$x,Rowv = modRowCluster()$ddm, Colv = modClusterData()$colddm, circleSize = modGraphParams()$circleSize, treeheight = modGraphParams()$treeHeight, fontsize = modGraphParams()$fontSize, cexRow = 1.2, 
                      color = color.heatmap(),annCol = modClusterData()$colAnnot,annColors = modColors(),
                      breaks=seq(-100,100,by=4)))
      dev.off()
    }
  )
  
  output$downloadModMap3Other <- downloadHandler(
    filename = function() {paste(values$project.name,'_Longitudinal_ModuleMapData','.csv', sep='')  },
    content = function(file) {
      x <- modOrderedData()$x
      if(!is.na(modRowCluster()$ddm)){
       x <- x[order.dendrogram(modRowCluster()$ddm),]
      }
      if(!is.na(modClusterData()$colddm)){
       x <- x[,order.dendrogram(modClusterData()$colddm)]
      }
      write.csv(x, file, row.names = TRUE)
    }
  )


  output$modOptimalNumber <- renderText({
    paste("Optimal number of clusters =", modClusterData()$opt_num )
  })
  
  output$modDunnIndexPlot <- renderPlot({
    barplot(modClusterData()$d, names.arg = 2:round(nrow(modOrderedData()$colAnnot)/2), xlab = "Number of Clusters", ylab = "Dunn's Index",cex.main = 1.5, col = "#4ba9d6")
  }, height = 400)

  output$downloadClusterPlot2Other <- downloadHandler(
    filename = function() {paste(values$project.name,'_','Cluster_Plot_Longitudinal','.png', sep = '')},
    content = function(file){
      png(file, width = 800)
      print(barplot(modClusterData()$d, names.arg = 2:round(nrow(modOrderedData()$colAnnot)/2), xlab = "Number of Clusters", ylab = "Dunn's Index"))
      dev.off()
    }
  )
  
  output$modClusterTableFull <- renderTable({
    req(modGen_clustTab())
    modGen_clustTab()$table1
  }, include.rownames = FALSE, digits = 0)
  
  output$modClusterTableDescription <- renderText({
    req(modGen_clustTab())
    if(ncol(modGen_clustTab()$table2) > 1){
      print("The table below is the table the tests are run on. It is the same as the table above, except the columns with zero counts have been deleted.")
    }
  })
  
  output$modClusterTable <- renderTable({
    req(modGen_clustTab())
    if(ncol(modGen_clustTab()$table2) > 1){
      modGen_clustTab()$table2
    }
  }, include.rownames = FALSE)
  
  
  output$modAssociationTestResults <- renderPrint({
    req(modGen_clustTab())
    if(ncol(modGen_clustTab()$table2) > 1){
      if(modGen_clustTab()[[3]]$p.value < .001){
        if(modGen_clustTab()[[4]]$p.value < .001){
          cat(paste0("Chi_square statistic = ", round(modGen_clustTab()[[3]]$statistic, 2), ", ", "p-value < .001","\nFishers Exact Test: p-value < .001"))
        } else{
          cat(paste0("Chi_square statistic = ", round(modGen_clustTab()[[3]]$statistic, 2), ", ", "p-value < .001",
                     "\nFishers Exact Test: p-value = ", round(modGen_clustTab()[[4]]$p.value, 3)))
        }
      } else{
        if(modGen_clustTab()[[4]]$p.value < .001){
          cat(paste0("Chi-Square Test Statistic = ", round(modGen_clustTab()[[3]]$statistic, 2), ", ", "p-value = ", round(modGen_clustTab()[[3]]$p.value, 3),
                     "\nFishers Exact Test: p-value < .001"))
        } else{
          cat(paste0("Chi-Square Test Statistic = ", round(modGen_clustTab()[[3]]$statistic, 2), ", ", "p-value = ", round(modGen_clustTab()[[3]]$p.value, 3),
                     "\nFishers Exact Test: p-value = ", round(modGen_clustTab()[[4]]$p.value, 3)))
        }
      }
    }
  })

  ##################################### Gene Level Heat Map #####################################

  heatmapData <- reactive({
    if (input$set==1){#heatmapbase1
      if(ctrl()$id == TRUE){
        base_sample_name <- design()$des$columnname[which(design()$des[,baseline.var()$id] == baseline.val()$id & design()$des[,control.var()$id] != control.val()$id)]
        index <- which(colnames(values$exprs) %in% base_sample_name)
      }
      if(ctrl()$id == FALSE){
        base_sample_name <- design()$des$columnname[which(design()$des[,baseline.var()$id] == baseline.val()$id)]
        index <- which(colnames(values$exprs) %in% base_sample_name)
      }
        exp_base_sam <- values$exprs[, index]
        des_base_sam <- design()$des[which(design()$des$columnname %in% colnames(exp_base_sam)),]
        y<-manipulateData(y = exp_base_sam, x = des_base_sam,colname = "columnname")
        ddm<-values$rowdend1b
    }
    if (input$set==2){#heatmapbase2
      base_sample_name <- design()$des$columnname[which(design()$des[,baseline.var()$id]==baseline.val()$id | design()$des[,control.var()$id]==control.val()$id)]
      index <- which(colnames(values$exprs) %in% base_sample_name)
      exp_base_sam <- values$exprs[, index]
      des_base_sam <- design()$des[which(design()$des$columnname %in% colnames(exp_base_sam)), ]
      y<-manipulateData(y=exp_base_sam,x=des_base_sam,colname ="columnname",ref.var=control.var()$id,ref.val=control.val()$id,long=FALSE,keep.ref=TRUE)
      ddm<-values$rowdend2b
    }
    if(input$set==3){#heatmap1
      y<-manipulateData(y=values$exprs,x=design()$des,colname="columnname")
      ddm<-values$rowdend1
    }
    if(input$set==4){#heatmap2
      y <- manipulateData(y = values$exprs, x = design()$des, colname = "columnname", ref.var = control.var()$id, ref.val = control.val()$id, long = FALSE, keep.ref = TRUE)
      ddm<-values$rowdend2
    }
    if(input$set==5){#heatmap3
      if(ctrl()$id==TRUE){
        des_w_controls<-design()$des[which(design()$des$columnname %in% colnames(values$exprs)),]
        des_wo_controls<-design()$des[-which(design()$des[,control.var()$id]==control.val()$id),]
        h5index<-which(colnames(values$exprs) %in% des_wo_controls$columnname)
        y<-manipulateData(y=values$exprs[,h5index],x=des_wo_controls,colname="columnname",ref.var=baseline.var()$id,ref.val=baseline.val()$id,long=TRUE,subject.id=subject.id()$id,keep.ref=FALSE)
      }
      if(ctrl()$id==FALSE){
        y<-manipulateData(y=values$exprs,x=design()$des,colname="columnname",ref.var=baseline.var()$id,ref.val=baseline.val()$id,long=TRUE,subject.id=subject.id()$id,keep.ref=FALSE)
      }
      ddm<-values$rowdend3
    }
    z<-list(y=y,ddm=ddm)
    return(z)
  })
  
  output$plotGeneSymbols <- renderUI({
    req(values$results.file)
    if(all.equal(as.character(values$results.file$Transcript.ID), as.character(values$results.file$Gene.Symbol)) == TRUE){
      return(NULL)
    }
    checkboxInput("plotGeneSymbols", "Plot Gene Symbols?", FALSE)
  })
  
  output$dgePlotGeneSymbols <- renderUI({
    req(values$results.file)
    if(all.equal(as.character(values$results.file$Transcript.ID), as.character(values$results.file$Gene.Symbol)) == TRUE){
      return(NULL)
    }
    checkboxInput("dgePlotGeneSymbols", "Plot Gene Symbols?", FALSE)
  })

  heatNormType <- reactive({
    if(input$set==1) heattxt<-paste0("Baseline ",str_to_sentence(values$norm.method)," Normalized")
    if(input$set==2) heattxt<-"Baseline Healthy Normalized"
    if(input$set==3) heattxt<-paste0("All Samples ",str_to_sentence(values$norm.method), " Normalized")
    if(input$set==4) heattxt<-"All Samples Healthy Normalized"
    if(input$set==5) heattxt<-"All Samples Normalized to each Subjects Baseline"
    return(heattxt)
  })
  
  graphParams <- reactive({
    params <- callModule(graphOptions, "unsupervisedGraphOptions")
    params <- list(width = params$width, height = params$height, fontSize = params$fontSize, legendSize = params$legendSize, 
                   treeHeight = params$treeHeight, resolution = params$resolution)
    return(params)
  })
  
  callModule(subsetAndOrderRenderUI, "unsupervisedSubOrder", des = reactive(heatmapData()$y$design.norm))
  callModule(clusterAssociationRenderUI, "unsupervisedAssociation", data = reactive(heatmapData()$y$design.norm))
  
  rowCluster <- eventReactive(input$go, {
    dat <- callModule(uploadVarsRowCluster, "transcripts", data = reactive(heatmapData()$y$exprs.norm), dendro = reactive(heatmapData()$ddm))
    x <- dat$x
    ddm <- dat$ddm
    labelRows <- dat$labelRows
    rowAnnot <- dat$rowAnnot
    y <- list(x = x, ddm = ddm, labelRows = labelRows, rowAnnot = rowAnnot)
    return(y)
  })
  
  clusterData <- eventReactive(input$go,{
    x <- callModule(colCluster, "unsupervisedCluster", des = reactive(orderedData()$colAnnot), data = reactive(orderedData()$x))
    return(x)
  })
  
  gen_clustTab <- reactive({
    x <- callModule(clusterAssociation, "unsupervisedAssociation", des = reactive(orderedData()$design), hclObj = reactive(clusterData()$hcl))
    return(x)
  })

  orderedData <- reactive({
    dat <- callModule(subsetAndOrder, "unsupervisedSubOrder", des = reactive(heatmapData()$y$design.norm), data = reactive(rowCluster()$x), 
                    sampleAnnot = reactive(sample.id()$id))
    x <- dat$dat
    colAnnot <- dat$colAnnot
    design <- dat$design
    z <- list(x = x, colAnnot = colAnnot, design = design)
    return(z)
  })
  
  maxRangeData <- eventReactive(input$go, {
    x <- callModule(maxValues, "unsupervisedMaxValues", reactive(orderedData()$x))
    if(!is.null(values$results.file) & all.equal(values$results.file$Transcript.ID, values$results.file$Gene.Symbol) != TRUE){
     if(input$plotGeneSymbols){
      rownames(x) <- make.unique(as.character(values$results.file$Gene.Symbol[match(rownames(x), as.character(values$results.file$Transcript.ID))]))
     }
    }
    return(x)
  })
  
  heatColors <- eventReactive(input$go, {
    x <- callModule(annColors, "unsupervisedSubOrder", reactive(orderedData()$colAnnot))
    return(x)
  })
  
  output$heatmap <- renderPlot({
    withProgress(message = 'Making plot',
                 detail = 'This may take a while...', value = 1,{
                   aheatmap(maxRangeData(),Rowv = rowCluster()$ddm,Colv = clusterData()$colddm, treeheight = graphParams()$treeHeight, fontsize = graphParams()$fontSize, cexRow = 1.2, 
                             color = colorRampPalette(c("navy", "yellow", "firebrick3"))(100),annCol = clusterData()$colAnnot,annColors = heatColors(),labRow=rowCluster()$labelRows,
                             annRow = rowCluster()$rowAnnot,breaks=0)
                 }
    )
  }, width = function(){graphParams()$width}, height = function(){graphParams()$height})
  
  output$downloadHeatmap <- downloadHandler(
    filename = function() {paste('Gene_Level_Heatmap','.png', sep = '')},
    content = function(file){
      png(file, width = (graphParams()$resolution/72)*graphParams()$width, height = (graphParams()$resolution/72)*graphParams()$height, res = graphParams()$resolution)
        print(aheatmap(maxRangeData(),Rowv = rowCluster()$ddm,Colv = clusterData()$colddm, treeheight = graphParams()$treeHeight, fontsize = graphParams()$fontSize, cexRow = 1.2, 
                        color = colorRampPalette(c("navy", "yellow", "firebrick3"))(100),annCol = clusterData()$colAnnot,annColors = heatColors(),labRow=rowCluster()$labelRows, 
                        breaks=0))
      dev.off()
    }
  )
  
  output$downloadHeatmapData <- downloadHandler(
    filename = function() {paste(values$project.name,'_',heatNormType(),'.csv', sep='')  },
    content = function(file) {
      write.csv(heatDataDownload(), file, row.names = TRUE)
    }
  )
  
  heatDataDownload <- reactive({
    x <- orderedData()$x
    if(is.na(rowCluster()$ddm)){
      if(is.na(clusterData()$colddm)){
        x <- x
      } else {
        x <- x[,order.dendrogram(clusterData()$colddm)]
      }
    } else {
      if(is.na(clusterData()$colddm)){
        x <- x[order.dendrogram(rowCluster()$ddm),]
      } else {
        x <- x[order.dendrogram(rowCluster()$ddm), order.dendrogram(clusterData()$colddm)]
        x <- x[nrow(x):1,]
      }
    }
    return(x)
  })


  output$dunnIndexPlot <- renderPlot({
    barplot(clusterData()$d, names.arg = 2:round(nrow(orderedData()$colAnnot)/2), xlab = "Number of Clusters", ylab = "Dunn Index",main = "Dunn Index for Cluster Number Selection", cex.main = 1.5, col = "#4ba9d6")
  }, height = 400)

  output$downloadClusterPlot3 <- downloadHandler(
    filename = function() {paste('DunnIndexPlot','.png', sep = '')},
    content = function(file){
      png(file, width = 800)
      print(barplot(clusterData()$d, names.arg = 2:round(nrow(orderedData()$colAnnot)/2), xlab = "Number of Clusters", ylab = "Dunn's Index", col = "#4ba9d6"))
      dev.off()
    }
  )

  output$OptimalNumber <- renderText({
    paste("Optimal number of clusters =", clusterData()$opt_num )
  })

  output$clusterTableFull <- renderTable({
    req(gen_clustTab())
    gen_clustTab()$table1
  }, include.rownames = FALSE, digits = 0)

  output$clusterTableDescription <- renderText({
    req(gen_clustTab())
    if(ncol(gen_clustTab()$table2) > 1){
      print("The table below is the table the tests are run on. It is the same as the table above, except the columns with zero counts have been deleted.")
    }
  })

  output$clusterTable <- renderTable({
    req(gen_clustTab())
    if(ncol(gen_clustTab()$table2) > 1){
      gen_clustTab()$table2
    }
  }, include.rownames = FALSE)
  
  output$associationTestResults <- renderPrint({
    req(gen_clustTab())
    if(ncol(gen_clustTab()$table2) > 1){
      if(gen_clustTab()[[3]]$p.value < .001){
        if(gen_clustTab()[[4]]$p.value < .001){
          cat(paste0("Chi_square statistic = ", round(gen_clustTab()[[3]]$statistic, 2), ", ", "p-value < .001","\nFishers Exact Test: p-value < .001"))
        } else{
          cat(paste0("Chi_square statistic = ", round(gen_clustTab()[[3]]$statistic, 2), ", ", "p-value < .001",
                       "\nFishers Exact Test: p-value = ", round(gen_clustTab()[[4]]$p.value, 3)))
        }
      } else{
        if(gen_clustTab()[[4]]$p.value < .001){
          cat(paste0("Chi-Square Test Statistic = ", round(gen_clustTab()[[3]]$statistic, 2), ", ", "p-value = ", round(gen_clustTab()[[3]]$p.value, 3),
                     "\nFishers Exact Test: p-value < .001"))
        } else{
          cat(paste0("Chi-Square Test Statistic = ", round(gen_clustTab()[[3]]$statistic, 2), ", ", "p-value = ", round(gen_clustTab()[[3]]$p.value, 3),
                     "\nFishers Exact Test: p-value = ", round(gen_clustTab()[[4]]$p.value, 3)))
        }
      }
    }
  })

  ################################# DGE ###########################################

  output$dge <- renderMenu({
    if(is.null(values$results.file)){
      return(strong(""))
    }
    else{
      return(menuItem("DGE Analysis", icon = icon("area-chart"), tabName = "genelistmaker"))
    }
  })

  dgeOverview <- reactive({
    results.file <- values$results.file
    estimates <- results.file[,grep("Estimate", names(results.file)), drop = FALSE]
    pvals <- results.file[,grep("^P.Value", names(results.file)), drop = FALSE]
    fdr.pvals <- results.file[,grep("FDR.P.Value", names(results.file)), drop = FALSE]
    bonf.pvals <- data.frame(apply(pvals, 2, p.adjust, method = "bonferroni"))
    colnames(bonf.pvals) <- gsub("P.Value.","Bonf.P.Value.",colnames(bonf.pvals))
    comparisons <- gsub("Estimate.", "", colnames(estimates))
    raw <- fdr <- bonf <- c()
    if(input$overviewFc){
      if(input$selectFcSign == "+"){
        for(i in 1:ncol(estimates)){
          raw[i] <- sum(pvals[,i] <= input$alphalevel2 & estimates[,i] >= input$selectFcValue, na.rm = TRUE)
          fdr[i] <- sum(fdr.pvals[,i] <= input$alphalevel2 & estimates[,i] >= input$selectFcValue, na.rm = TRUE)
          bonf[i] <- sum(bonf.pvals[,i] <= input$alphalevel2 & estimates[,i] >= input$selectFcValue, na.rm = TRUE)
        }
      } else if(input$selectFcSign == "-"){
        for(i in 1:ncol(estimates)){
          raw[i] <- sum(pvals[,i] <= input$alphalevel2 & estimates[,i] <= -input$selectFcValue, na.rm = TRUE)
          fdr[i] <- sum(fdr.pvals[,i] <= input$alphalevel2 & estimates[,i] <= -input$selectFcValue, na.rm = TRUE)
          bonf[i] <- sum(bonf.pvals[,i] <= input$alphalevel2 & estimates[,i] <= -input$selectFcValue, na.rm = TRUE)
        }
      } else{
        for(i in 1:ncol(estimates)){
          raw[i] <- sum(pvals[,i] <= input$alphalevel2 & (estimates[,i] >= input$selectFcValue | estimates[,i] <= -input$selectFcValue), na.rm = TRUE)
          fdr[i] <- sum(fdr.pvals[,i] <= input$alphalevel2 & (estimates[,i] >= input$selectFcValue | estimates[,i] <= -input$selectFcValue), na.rm = TRUE)
          bonf[i] <- sum(bonf.pvals[,i] <= input$alphalevel2 & (estimates[,i] >= input$selectFcValue | estimates[,i] <= -input$selectFcValue), na.rm = TRUE)
        }
      }
    } else{
      raw <- apply(pvals, 2, function(x) sum(x <= input$alphalevel2, na.rm = TRUE))
      fdr <- apply(fdr.pvals, 2, function(x) sum(x <= input$alphalevel2, na.rm = TRUE))
      bonf <- apply(bonf.pvals, 2, function(x) sum(x <= input$alphalevel2, na.rm = TRUE))
    }
    overview <- data.frame(Comparison = comparisons, Raw = raw, FDR = fdr, Bonf = bonf)
    y <- list(z = overview)
    return(y)
  })

  output$dgeOverviewTable <- renderDT({
    withProgress(message = 'Making the table',
                 detail = 'This may take a while...', value = 1,{
                   datatable(dgeOverview()$z[order(dgeOverview()$z$Raw,decreasing=TRUE),],
                             class = "table-condensed",
                             style = "bootstrap4",
                             rownames = FALSE,
                             options = list(autowidth = TRUE, scrollX = TRUE)
                   )
                 })
  })
  
  index <- reactive({
    which(names(values$results.file) == paste0("P.Value.",input$diagnosticsComparison))
  })
  
  pcomp <- reactive({
    x <- sub("^", "P.Value.", input$comparison)
    x
  })

  plottitle1 <- reactive({
    paste("Distribution of Raw p-values for", input$diagnosticsComparison)
  })

  output$pvalDistPlot<-renderPlot({
    y<-max(hist(values$results.file[, index()])$density)
    hist(values$results.file[,index()], freq=FALSE, xlim=c(0,1), ylim=c(0,y),main=plottitle1(),
         xlab="Raw p-value's", ylab="Density")
    lines(c(0,1),c(1,1),lwd=2,lty=2)
  }, height = 400)
  
  pvalThreshData <-reactive({
    dat <- c()
    thresh <- c(0.001,0.01,1:19/20)
    for(i in thresh){
      dat <- rbind(dat, c(mycorrection(values$results.file[,index()],i,"RAW"),
                          mycorrection(values$results.file[,index()],i,"FDR"),
                          mycorrection(values$results.file[,index()],i,"BONF")))
    }
    dat <- data.frame(cbind(thresh,dat))
    colnames(dat) <- c("Alpha","Raw","FDR","Bonf")
    return(dat)
  })

  output$pvalThreshTable <-renderDT({
    dat <- cbind(pvalThreshData()[,1,drop=FALSE],pvalThreshData()[,2:4]*dim(values$results.file)[1])
    datatable(dat,
              class = "table-condensed",
              style = "bootstrap4",
              rownames = FALSE,
              options = list(autowidth = TRUE, scrollX = TRUE)
    )
  })
  
  multGeneLists <- reactive({
    if(length(grep("All", input$comparisonsDownload)) > 0){
      nam <- names(values$results.file)
      index <- grep("^P.Value",nam)
      p.names <- nam[index]
      comps <- gsub("P.Value.", "", p.names)
    }
    else{
      comps <- input$comparisonsDownload
    }
    x.all <- list()
    for(i in 1:length(comps)){
      x <- callModule(filterOpts, "dgeResultsTable", data = reactive(values$results.file), comparison = reactive(comps[i]),"genes")
      x.all[[i]] <- x
    }
    names(x.all) <- comps
    return(x.all)
  })

  dgeModuleDat <- reactive({
    if(!is.null(values$dge.gsets)){
      genes <- unlist(values$dge.gsets)
      gsetNames <- rep(names(values$dge.gsets), times = lapply(values$dge.gsets, length))
      modinfo <- data.frame(Module = gsetNames, Transcript.ID = genes)
      if(input$dgeResults == "Modular DGE Analysis"){
        dat <- callModule(filterOpts, "modDgeMap", data = reactive(values$results.file), data.type = "percents", geneList = reactive(modinfo))
      }
      if(input$dgeResults == "Gene Lists"){
        dat <- callModule(filterOpts, "dgeResultsTable", data = reactive(values$results.file), data.type = "percents", geneList = reactive(modinfo))
      }
      if(baylorMod() == TRUE){
        modnames <- unique(modinfo$Module)
        modnum <- gsub("M","",modnames)
        modvec <- as.numeric(unlist(strsplit(modnum,".",fixed=TRUE)))
        modmat <- matrix(modvec,ncol=2,byrow=T)
        modordnum <- table(factor(modinfo$Module,levels=modnames[order(modmat[,1],modmat[,2])],ordered=TRUE))
        mod_for_merge <- data.frame(Module=names(modordnum),Size=as.vector(modordnum))
        prop_matrix <- dat$prop_matrix[match(as.character(mod_for_merge$Module), rownames(dat$prop_matrix), nomatch = 0),,drop = FALSE]
        prop_matrix2 <- dat$prop_matrix2[match(as.character(mod_for_merge$Module), rownames(dat$prop_matrix2), nomatch = 0),,drop = FALSE]
        z <- list(prop_matrix = prop_matrix, prop_matrix2 = prop_matrix2)
      } else{
        z <- list(prop_matrix = dat$prop_matrix, prop_matrix2 = dat$prop_matrix2)
      }
      return(z)
    } else{
      return(NULL)
    }
  })
  
  dgeModuleMap <- reactive({
    numgen <- function(x){
      d <- c()
      odd <- 2*(1:20)-1
      if(length(x) == 1){
        return(odd[1:x])
      }
      if(length(x) > 1){
        for(i in 1:length(x)){
          d <- c(d,odd[1:x[i]])
        }
        return(d)
      }
    }
    numgen2 <- function(x){
      d <- c()
      odd <- 2*(1:8)-1
      for(i in 1:length(x)){
        d <- c(d,rep(odd[9-i],x[i]))
      }
      return(d)
    }
    colorvar <- dgeModuleDat()$prop_matrix2[1:97,input$comparison]
    colgrp <- findInterval(colorvar,seq(0,2,length.out=10))
    colfunc <- colorRampPalette(c("blue","white", "red"))
    collist <- colfunc(length(unique(colgrp)))
    mycolors <- collist[colgrp]
    mycolors = colorvar
    df2 <- data.frame(
      x = numgen(c(2,3,6,16,15,20,20,15)),
      y = c(numgen2(c(2,3,6,16,15,20,20,15))),
      proportion = mycolors
    )
    test <- ggplot(df2, aes(xmin=x-1,xmax=x+1,ymin=y-1,ymax=y+1),environment=environment())+
      theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),panel.grid.major = element_blank()
            ,panel.grid.minor = element_blank()
            ,panel.border = element_blank())+geom_rect(fill="white", colour="black")+
      scale_y_continuous(breaks=c(15,13,11,9,7,5,3,1), labels=c("M1","M2","M3","M4","M5","M6","M7","M7 (21-35)"),limits=c(0,40))+
      scale_x_continuous(breaks=(2*(1:20)-1),labels=1:20,limits=c(0,40))
    
    test2<-test + coord_cartesian(xlim = c(0, 40), ylim=c(0, 16))
    test3<-test2 + geom_point(aes(x=x,y=y,colour=proportion,size=500))+
      scale_colour_gradient2(low="blue", high="red", guide="colorbar",midpoint=0, limits = c(-1,1))+ 
      theme(axis.title.x=element_blank(),axis.title.y=element_blank()) + scale_size_continuous(limits=c(1,500)) + 
      guides(size = FALSE) 
    return(test3)
  })
  
  output$moduleMap <- renderUI({
    if(!input$baylorMods){return(NULL)}
    plotOutput('dgeModuleMap')
  })
  
  output$dgeModuleMapResolution <- renderUI({
    if(!input$baylorMods){return(NULL)}
    numericInput('dgeModuleMapResolution', "Plot resolution (for downloaded plot):", min = 72, max = 400, value = 72, step = 1)
  })
  
  output$dgeModuleMap <- renderPlot({
    dgeModuleMap()
  }, width = 700, height = 190)
  
  output$baylorMods <- renderUI({
    req(values$dge.gsets)
    checkboxInput("baylorMods", strong("Are the provided gene sets baylor modules?"), FALSE)
  })
  
  output$downloadDgeModMap <- renderUI({
    if(!input$baylorMods){return(NULL)}
    downloadButton('downloadDgeModuleMap', "Download Figure")
  })
  
  output$downloadDgeModuleMap <- downloadHandler(
    filename = function() {paste('dgeModuleMap','.png', sep = '')},
    content = function(file){
      png(file, width = (input$dgeModuleMapResolution/72)*700, height = (input$dgeModuleMapResolution/72)*190, res = input$dgeModuleMapResolution)
      print(dgeModuleMap())
      dev.off()
    }
  )
  
  output$mergeGeneSets <- renderUI({
    if(!is.null(values$dge.gsets)){
      return(checkboxInput("mergeGeneSets", "Merge gene set information", FALSE))
    } else{
      return(NULL)
    }
  })
  
  mergeDgeGeneSet <- reactive({
    x <- reshape2::melt(values$dge.gsets)
    x <- x[,c(2,1)]
    colnames(x) <- c("Gene.Set", "Gene")
    if(!is.null(values$dge.annots)){
      y <- values$dge.annots
      colnames(y) <- c("Gene.Set", "Annotation")
      x <- full_join(x, y, by = "Gene.Set")
    }
    return(x)
  })

  output$dgeResultTable <- renderDT({
    input$comparison
    y <- callModule(filterOpts, "dgeResultsTable", data = reactive(values$results.file), comparison = reactive(input$comparison), "genes")
    if(!is.null(values$dge.gsets)){
      if(input$mergeGeneSets){
        x <- mergeDgeGeneSet()
        colnames(x)[2] <- "Transcript.ID"
        y <- left_join(y, x, by = "Transcript.ID") %>%
          select(Transcript.ID, Gene.Symbol, Gene.Set, Annotation, everything())
      }
    }
    for(i in 1:ncol(y)){
      if(is.numeric(y[,i])){
        y[,i] <- as.numeric(formatC(y[,i], digits = 4))
      }
    }
    y$Gene.Symbol <- paste("<a href=http://www.genecards.org/cgi-bin/carddisp.pl?gene=",y$Gene.Symbol," target = '_blank'",'>',y$Gene.Symbol,"</a>",sep='')
    datatable(y,
              class = "table-condensed",
              style = "bootstrap4",
              rownames = FALSE,
              escape = FALSE,
              options = list(autowidth = TRUE, scrollX = TRUE)
    )
  })

  dgeData <- reactive({
    input$comparisonDgeHeatmap
    genes <- callModule(filterOpts, "dgeHeatmap", data = reactive(values$results.file), comparison = reactive(input$comparisonDgeHeatmap),"genes")
    dat <- values$exprs[match(genes$Transcript.ID, rownames(values$exprs), nomatch = 0),]
    genes <- genes[match(rownames(dat), genes$Transcript.ID, nomatch = 0),]
    genes$Gene.Symbol <- as.character(genes$Gene.Symbol)
    genes$Transcript.ID <- as.character(genes$Transcript.ID)
    if(length(which(genes$Gene.Symbol == "")) > 0){
     rows <- genes$Gene.Symbol
     rows[which(rows == "")] <- genes$Transcript.ID[which(rows == "")]
    } else{
     rows <- genes$Gene.Symbol
    }
    if(input$dgePlotGeneSymbols){
     rownames(dat) <- make.unique(as.character(rows))
    }
    return(dat)
  })

  output$test1<-renderUI({
   if(is.null(values$rowdend3)){
    if(ctrl()$id ==TRUE){
     try<-list(3,4)
     names(try) <- c(paste0("All Samples ",str_to_sentence(values$norm.method), " Normalized"),"All Samples Healthy Normalized")
    }
    
    if(ctrl()$id == FALSE){
     try<-list(3)
     names(try) <- paste0("All Samples ",str_to_sentence(values$norm.method), " Normalized")
    }
   }
   
   if(!is.null(values$rowdend3)){
    if(ctrl()$id == TRUE){
     try<-list(1,2,3,4,5)
     names(try) <- c(paste0("Baseline ",str_to_sentence(values$norm.method), " Normalized"),"Baseline Healthy Normalized",paste0("All Samples ",str_to_sentence(values$norm.method), " Normalized"),
                     "All Samples Healthy Normalized","All Samples Baseline Normalized")
    }
    if(ctrl()$id == FALSE){
     try<-list(1,3,5)
     names(try) <- c(paste0("Baseline ",str_to_sentence(values$norm.method), " Normalized"),paste0("All Samples ",str_to_sentence(values$norm.method), " Normalized"),"All Samples Baseline Normalized")
    }
   }
   selectizeInput("set1", "Select heatmap:",as.list(try))
  })


  dgeHeatmapData <-reactive({
    dat <- dgeData()
    if (input$set1==1){#heatmapbase1
      if(ctrl()$id == TRUE){
        base_sample_name <- design()$des$columnname[which(design()$des[,baseline.var()$id] == baseline.val()$id & design()$des[,control.var()$id] != control.val()$id)]
        index <- which(colnames(dat) %in% base_sample_name)
      }
      if(ctrl()$id == FALSE){
        base_sample_name <- design()$des$columnname[which(design()$des[,baseline.var()$id] == baseline.val()$id)]
        index <- which(colnames(dat) %in% base_sample_name)
      }
      exp_base_sam <- dat[, index]
      des_base_sam <- design()$des[which(design()$des$columnname %in% colnames(exp_base_sam)),]
      y<-manipulateData(y = exp_base_sam, x = des_base_sam,colname = "columnname")
    }
    if (input$set1==2){#heatmapbase2
      base_sample_name <- design()$des$columnname[which(design()$des[,baseline.var()$id]==baseline.val()$id | design()$des[,control.var()$id]==control.val()$id)]
      index <- which(colnames(dat) %in% base_sample_name)
      exp_base_sam <- dat[, index]
      des_base_sam <- design()$des[which(design()$des$columnname %in% colnames(exp_base_sam)), ]
      y<-manipulateData(y=exp_base_sam,x=des_base_sam,colname ="columnname",ref.var=control.var()$id,ref.val=control.val()$id,long=FALSE,keep.ref=TRUE)
    }
    if(input$set1==3){#heatmap1
      y<-manipulateData(y=dat,x=design()$des,colname="columnname")
    }
    if(input$set1==4){#heatmap2
      y <- manipulateData(y = dat, x = design()$des, colname = "columnname", ref.var = control.var()$id, ref.val = control.val()$id, long = FALSE, keep.ref = TRUE)
    }
    if(input$set1==5){#heatmap3
      if(ctrl()$id==TRUE){
        des_w_controls<-design()$des[which(design()$des$columnname %in% colnames(dat)),]
        des_wo_controls<-design()$des[-which(design()$des[,control.var()$id]==control.val()$id),]
        h5index<-which(colnames(dat) %in% des_wo_controls$columnname)
        y<-manipulateData(y=dat[,h5index],x=des_wo_controls,colname="columnname",ref.var=baseline.var()$id,ref.val=baseline.val()$id,long=TRUE,subject.id=subject.id()$id,keep.ref=FALSE)
      }
      if(ctrl()$id==FALSE){
        y<-manipulateData(y=dat,x=design()$des,colname="columnname",ref.var=baseline.var()$id,ref.val=baseline.val()$id,long=TRUE,subject.id=subject.id()$id,keep.ref=FALSE)
      }
    }
    z<-list(y=y)
    return(z)
  })

  heatNormType1<-reactive({
    if(input$set1==1) heattxt<-paste0("Baseline ",str_to_sentence(values$norm.method), " Normalized")
    if(input$set1==2) heattxt<-"Baseline Healthy Normalized"
    if(input$set1==3) heattxt<-paste0("All Samples ",str_to_sentence(values$norm.method), " Normalized")
    if(input$set1==4) heattxt<-"All Samples Healthy Normalized"
    if(input$set1==5) heattxt<-"All Samples Normalized to each Subjects Baseline"
    heattxt
  })
  
  dgeGraphParams <- reactive({
    params <- callModule(graphOptions, "dgeGraphOptions")
    params <- list(width = params$width, height = params$height, fontSize = params$fontSize, legendSize = params$legendSize, 
                   treeHeight = params$treeHeight, resolution = params$resolution)
    return(params)
  })
  
  callModule(subsetAndOrderRenderUI, "dgeSubOrder", des = reactive(dgeHeatmapData()$y$design.norm))
  
  dgeRowCluster <- eventReactive(input$go2, {
    x <- dgeHeatmapData()$y$exprs.norm 
    ddm <- NA
    if(input$rowCluster){
      dist <- dist(x)
      hcl <- fastcluster::hclust(dist)
      ddm <- as.dendrogram(hcl)
      Rowv <- rowMeans(x, na.rm = TRUE)
      ddm <- reorder(ddm, Rowv)
    }
    labelRows <- NULL
    y <- list(x = x, ddm = ddm, labelRows = labelRows)
    return(y)
  })
  
  dgeClusterData <- eventReactive(input$go2,{
    x <- callModule(colCluster, "dgeCluster", des = reactive(dgeOrderedData()$colAnnot), data = reactive(dgeOrderedData()$x))
    return(x)
  })
  
  dgeOrderedData <- eventReactive(input$go2,{
    dat <- callModule(subsetAndOrder, "dgeSubOrder", des = reactive(dgeHeatmapData()$y$design.norm), data = reactive(dgeRowCluster()$x), 
                      sampleAnnot = reactive(sample.id()$id))
    x <- dat$dat
    colAnnot <- dat$colAnnot
    design <- dat$design
    z <- list(x = x, colAnnot = colAnnot, design = design)
    return(z)
  })
  
  dgeMaxRangeData <- eventReactive(input$go2, {
    x <- callModule(maxValues, "dgeMaxValues", reactive(dgeOrderedData()$x))
    return(x)
  })
  
  dgeHeatColors <- eventReactive(input$go2, {
    x <- callModule(annColors, "dgeSubOrder", reactive(dgeOrderedData()$colAnnot))
    return(x)
  })
  
  output$dgeHeatmap <- renderPlot({
    withProgress(message = 'Making plot',
                 detail = 'This may take a while...', value = 1,{
                   aheatmap(dgeMaxRangeData(),Rowv = dgeRowCluster()$ddm,Colv = dgeClusterData()$colddm, treeheight = dgeGraphParams()$treeHeight, fontsize = dgeGraphParams()$fontSize, cexRow = 1.2, 
                             color = colorRampPalette(c("navy", "yellow", "firebrick3"))(100),annCol = dgeClusterData()$colAnnot,annColors = dgeHeatColors(),labRow=dgeRowCluster()$labelRows,
                             breaks=0)
                 }
    )
  }, width = function(){dgeGraphParams()$width}, height = function(){dgeGraphParams()$height})
  
  
  output$downloadDgeHeatmapFigure <- downloadHandler(
    filename = function() {paste('SigVarsHeatmap','.png', sep = '')},
    content = function(file){
      png(file, width = (dgeGraphParams()$resolution()/72)*dgeGraphParams()$width, height = (dgeGraphParams()$resolution()/72)*dgeGraphParams()$height, res = dgeGraphParams()$resolution())
      print(aheatmap(dgeMaxRangeData(),Rowv = dgeRowCluster()$ddm,Colv = dgeClusterData()$colddm, treeheight = dgeGraphParams()$treeHeight, fontsize = dgeGraphParams()$fontSize, cexRow = 1.2, 
                      color = colorRampPalette(c("navy", "yellow", "firebrick3"))(100),annCol = dgeClusterData()$colAnnot,annColors = dgeHeatColors(),labRow=dgeRowCluster()$labelRows,
                      breaks=0))
      dev.off()
    }
  )
  
  output$downloadDgeHeatmapData <- downloadHandler(
    filename = function() {paste(values$project.name,'_',heatNormType1(),'.csv', sep='')  },
    content = function(file) {
      write.csv(heatDataDownload2(), file, row.names = TRUE)
    }
  )
  
  heatDataDownload2 <- reactive({
    x <- dgeOrderedData()$x
    if(is.na(dgeRowCluster()$ddm)){
      if(is.na(dgeClusterData()$colddm)){
        x <- x
      } else {
        x <- x[,order.dendrogram(dgeClusterData()$colddm)]
      }
    } else {
      if(is.na(dgeClusterData()$colddm)){
        x <- x[order.dendrogram(dgeRowCluster()$ddm),]
      } else {
        x <- x[order.dendrogram(dgeRowCluster()$ddm), order.dendrogram(dgeClusterData()$colddm)]
        x <- x[nrow(x):1,]
      }
    }
    return(x)
  })

  observe({
    req(values$results.file)
    nam <- names(values$results.file)
    index <- grep("^P.Value",nam)
    p.names <- nam[index]
    p.names <- gsub("P.Value.", "", p.names)
    updateSelectizeInput(session,"comparison", "Comparison:", choices = p.names, selected = p.names[1])
    updateSelectizeInput(session,"comparisonModuleMat", "Comparison:", choices = p.names, selected = p.names[1])
    updateSelectizeInput(session,"comparisonsDownload", "Comparison:", choices = c("All",p.names), selected = "All")
    updateSelectizeInput(session,"diagnosticsComparison", "Comparison:", choices = p.names, selected = p.names[1])
    updateSelectizeInput(session,"Vcomparison", "Comparison:", choices = p.names, selected = p.names[1])
  })
  
  observe({
    req(values$results.file)
    nam <- names(values$results.file)
    index <- grep("^P.Value",nam)
    p.names <- nam[index]
    p.names <- gsub("P.Value.", "", p.names)
    updateSelectizeInput(session, "comparisonDgeHeatmap", "Comparison:", choices = p.names, selected = input$comparison)
  })
  
  output$pvalThreshCurve <- renderPlot({
    plot(pvalThreshData()$alpha,pvalThreshData()$raw,type="l",main="Multiple Testing Comparison Plot",col="black",xlim=c(0,1),ylim=c(0,1),xlab=expression(alpha),ylab="% of Total Probes in Gene List")
    lines(pvalThreshData()$alpha,pvalThreshData()$fdr,col="red",lwd=2)
    lines(pvalThreshData()$alpha,pvalThreshData()$bonf,col="green")
    legend("bottomright",legend=c("Raw P-value","FDR","Bonf."), lty=c(1,1,1),col=c("black","red","green") )
    axis(4,at=1:10/10,labels=round(1:10/10*dim(values$results.file)[1],0))
    lines(c(0,1),c(0,1),lty=2)
  }, height = 400)

  output$downloadSC <- downloadHandler(
    filename = function() {paste("Significance_Comparison_Overview", '_', input$alphalevel2, '.csv', sep = '') },
    content = function(file){
      write.csv(dgeOverview()$z, file,row.names = FALSE)
    }
  )

  output$downloadData <- downloadHandler(
    filename = function() {paste(substring(pcomp(),12),'_',input$correction_method1,input$alphalevel1,'.csv', sep='')},
    content = function(file) {
      input$comparison
      y <- callModule(filterOpts, "dgeResultsTable", data = reactive(values$results.file), comparison = reactive(input$comparison),"genes")
      write.csv(y, file,row.names = FALSE)
    }
  )
  
  output$downloadSelComp <- downloadHandler(
    filename = function() {paste(input$ziptext,".zip", sep = "")},
    content = function(file) {
      file.names <- c()
      for(i in 1:length(multGeneLists())){
        file.names[i] <- paste(names(multGeneLists())[i],".csv", sep = "")
        write.csv(multGeneLists()[[i]], file = paste(names(multGeneLists())[i], ".csv", sep = ""), row.names = FALSE)
      }
      zip(zipfile = file, files = file.names)
    },
    contentType = "application/zip"
  )
  
  output$intro <-renderText({paste("Number of probes in Gene list by method and alpha parameter.")})

  modDgeGraphParams <- reactive({
    params <- callModule(graphOptions, "modDgeGraphOptions", varType = "modules")
    params <- list(width = params$width, height = params$height, fontSize = params$fontSize, legendSize = params$legendSize, 
                   treeHeight = params$treeHeight, resolution = params$resolution,circleSize = params$circleSize)
    return(params)
  })
  
  output$dgeModSelection <- renderUI({
    if(!input$baylorMods){return(NULL)}
    if(!is.null(values$dge.annots)){
     return(
      selectizeInput("dgeModuleSelection", "Module to include:", c("All", "Only Annotated", "First Round", "First Two Rounds", "First Three Rounds", "First Four Rounds", "First Five Rounds",
                                                                   "First Six Rounds", "First Seven Rounds", "First Eight Rounds"), "First Six Rounds")
     )
    } else{
     return(
      selectizeInput("dgeModuleSelection", "Module to include:", c("All", "First Round", "First Two Rounds", "First Three Rounds", "First Four Rounds", "First Five Rounds",
                                                                   "First Six Rounds", "First Seven Rounds", "First Eight Rounds"), "First Six Rounds")
     )
    }
  })
  
  modDgeData <- eventReactive(input$goModDge,{
    dat <- dgeModuleDat()$prop_matrix2*100
    if(!is.null(values$dge.annots)){
     annots <- values$dge.annots[match(rownames(dat)[which(rownames(dat) %in% as.character(values$dge.annots[,1]))], 
                                       as.character(values$dge.annots[,1]), nomatch = 0),]
     rownames(dat)[which(rownames(dat) %in% as.character(values$dge.annots[,1]))] <- 
      paste0(rownames(dat)[which(rownames(dat) %in% as.character(values$dge.annots[,1]))], " ", as.character(annots[,2]))
    }
    ddm <- NA
    colddm <- NA
    if(input$compSelect){
      dat <- dat[,match(input$comparisonModuleMat, colnames(dat), nomatch = 0), drop = FALSE]
    }
    if(baylorMod()){
     if(input$dgeModuleSelection == "Only Annotated"){
      dat <- dat[grep(" ", rownames(dat)),, drop = FALSE]
     } else if(input$dgeModuleSelection == "First Round"){
      dat <- dat[grep("^M1", rownames(dat)),, drop = FALSE]
     } else if(input$dgeModuleSelection == "First Two Rounds"){
      dat <- dat[grep("^M1|^M2", rownames(dat)),, drop = FALSE]
     } else if(input$dgeModuleSelection == "First Three Rounds"){
      dat <- dat[grep("^M1|^M2|^M3", rownames(dat)),, drop = FALSE]
     } else if(input$dgeModuleSelection == "First Four Rounds"){
      dat <- dat[grep("^M1|^M2|^M3|^M4", rownames(dat)),, drop = FALSE]
     } else if(input$dgeModuleSelection == "First Five Rounds"){
      dat <- dat[grep("^M1|^M2|^M3|^M4|^M5", rownames(dat)),, drop = FALSE]
     } else if(input$dgeModuleSelection == "First Six Rounds"){
      dat <- dat[grep("^M1|^M2|^M3|^M4|^M5|^M6", rownames(dat)),, drop = FALSE]
     } else if(input$dgeModuleSelection == "First Seven Rounds"){
      dat <- dat[grep("^M1|^M2|^M3|^M4|^M5|^M6|^M7", rownames(dat)),, drop = FALSE]
     } else if(input$dgeModuleSelection == "First Eight Rounds"){
      dat <- dat[grep("^M1|^M2|^M3|^M4|^M5|^M6|^M7|^M8", rownames(dat)),, drop = FALSE]
     } 
    }
    if(input$uploadModules){
      modNames <- read.csv(input$modSelect$datapath, header = TRUE)
      dat <- dat[which(rownames(dat) %in% modNames[,1]),]
    }
    if(input$LMMdeleterows){
      dat <- dat[-which(rowSums(dat) == 0),]
    }
    if(input$modDgeCluster){
      ddm <- FALSE
    }
    if(input$colCluster){
      colddm <- TRUE
    }
    z <- list(dat = dat, ddm = ddm, colddm = colddm)
    return(z)
  })
  
  output$modDgeMap <- renderPlot({
    withProgress(message = 'Making plot',
                 detail = 'This may take a while...', value = 1,{
                   aheatmap_circle(modDgeData()$dat,Rowv = modDgeData()$ddm,Colv = modDgeData()$colddm,circleSize = modDgeGraphParams()$circleSize, treeheight = modDgeGraphParams()$treeHeight, fontsize = modDgeGraphParams()$fontSize, cexRow = 1.2, 
                             color = color.heatmap(),breaks=seq(-100,100,by=4))
                 }
    )
  }, width = function(){modDgeGraphParams()$width}, height = function(){modDgeGraphParams()$height})
  
  output$downloadModDgeMap <- downloadHandler(
    filename = function() {paste(values$project.name,"_LMM_ModuleMap",'_',input$correction_method1,input$alphalevel1,'.png', sep='')},
    content = function(file) {
      png(file, width = (modDgeGraphParams()$resolution/72)*modDgeGraphParams()$width, height = (modDgeGraphParams()$resolution/72)*modDgeGraphParams()$height, 
          res = modDgeGraphParams()$resolution)
      print(aheatmap_circle(modDgeData()$dat,Rowv = modDgeData()$ddm,Colv = modDgeData()$colddm,circleSize = modDgeGraphParams()$circleSize, treeheight = modDgeGraphParams()$treeHeight, fontsize = modDgeGraphParams()$fontSize, cexRow = 1.2, 
                      color = color.heatmap(),breaks=seq(-100,100,by=4)))
      dev.off()
    }
  )

  output$downloadModDgeData <- downloadHandler(

    filename = function() {paste(values$project.name,"_LMM_ModuleData",'_',input$correction_method1,input$alphalevel1,'.csv', sep='')},
    content = function(file) {
      write.csv(modDgeData()$dat, file)
    }
  )
  
  observe({
    req(input$Vcomparison)
    updateSelectizeInput(session, "Include", "Include:", choices = input$Vcomparison, selected = input$Vcomparison[1:length(input$Vcomparison)])
  })
  
  observe({
    req(input$Include)
    comparisons <- setdiff(input$Vcomparison, input$Include)
    updateSelectizeInput(session,"Exclude", "Exclude:", choices = comparisons, selected = comparisons[1])
  })

  venndata <- reactive({
    input$Vcomparison
    venn.data <- list()
    for(i in 1:length(input$Vcomparison)){
      venn.data[[i]] <- callModule(filterOpts, "dgeVenn", reactive(values$results.file), reactive(input$Vcomparison[i]), "genes")
      venn.data[[i]] <- as.character(venn.data[[i]][,1])
    }
    names(venn.data) <- input$Vcomparison
    if(length(venn.data) > 5){
    venn.data <- venn.data[1:5]
    }
    return(venn.data)
  })

  Venn.intersection <- reactive({
    req(venndata())
    venndata1 <- venndata()[which(names(venndata()) %in% input$Include)]
    intersections <- Reduce(intersect, venndata1)
    n = length(input$Exclude)
    if(n > 0){
      venndata2 <- venndata()[which(names(venndata()) %in% input$Exclude)]
      venndata2 = unlist(venndata2)
      excl <- list(intersections, venndata2)
      intersections <- Reduce(setdiff, excl)
    }
    return(intersections)
  })
  
  output$testing <- renderText({
    length(venndata())
  })

  Venn.union <- reactive({
    req(venndata())
    venndata1 <- venndata()[which(names(venndata()) %in% input$Include)]
    unions <- Reduce(union, venndata1)
    n = length(input$Exclude)
    if(n > 0){
      venndata2 <- venndata()[which(names(venndata()) %in% input$Exclude)]
      venndata2 <- unlist(venndata2)
      excl <- list(unions, venndata2)
      unions <- Reduce(setdiff, excl)
    }
    return(unions)
  })

  genelist2 <- reactive({
    n = length(input$Include)
    comp <- input$Include
    fcomp <- list()
    for(i in 1:n){
      fcomp[[i]] <- c(paste0("Estimate.", input$Include[i]), paste0("Test.statistic.", input$Include[i]),
                      paste0("P.Value.", input$Include[i]))
    }
    if(input$UorI == 1){
      y <- which(values$results.file[,1] %in% Venn.intersection())
      cols <- list()
      for(i in 1:n){
        cols[[i]] <- match(fcomp[[i]], colnames(values$results.file))
      }
    }
    if(input$UorI == 2){
      y <- which(values$results.file[,1] %in% Venn.union())
      cols <- list()
      for(i in 1:n){
        #cols[[i]] <- grep(fcomp[i], colnames(values$results.file))
        cols[[i]] <- match(fcomp[[i]], colnames(values$results.file))
      }
    }
    x = list()
    pcols = list()
    pvals_fdr = list()
    pvals_bonf = list()
    pvals <- list()
    estimates <- list()
    tstats <- list()
    for(i in 1:n){
      if(i == 1){
        x[[i]] <- values$results.file[y,c(1,2,cols[[i]])]
      }
      if(i > 1){
        x[[i]] <- values$results.file[y,cols[[i]]]
      }
      pvals[[i]] <- grep("^P.Value", colnames(x[[i]]))
      estimates[[i]] <- grep("^Estimate", colnames(x[[i]]))
      tstats[[i]] <- grep("^Test.statistic", colnames(x[[i]]))
      if(i == 1){
        x[[i]] <- x[[i]][,c(1,2,estimates[[i]],pvals[[i]])]
      }
      if(i > 1){
        x[[i]] <- x[[i]][,c(estimates[[i]],pvals[[i]])]
      }
    }
    for(i in 1:n){
      if(i == 1){
        colnames(x[[i]]) = c("Transcript.ID", "Gene.Symbol", "Log2FC","P.Value")
      }
      if(i > 1){
        colnames(x[[i]]) = c("Log2FC", "P.Value")
      }
    }
    if(n > 1){
      for(i in 1:n){
        if(i == 1){
          colnames(x[[i]]) = c("Transcript.ID", "Gene.Symbol", paste("Log2FC for", comp[i], sep = " "), paste("P.Value.", comp[i], sep =" "))
        }
        if(i > 1){
          colnames(x[[i]]) = c(paste("Log2FC for", comp[i], sep = " "), paste("P.Value.", comp[i], sep = " "))
        }
      }
    }
    x = do.call("cbind", x)
    if(nrow(x) == 0){
      x <- data.frame("No Genes Present")
      names(x) <- "Transcript.ID"
    }
    return(x)
  })

  output$vennIntersection <- renderDT({
    y <- genelist2()
    y$Gene.Symbol <- paste("<a href=http://www.genecards.org/cgi-bin/carddisp.pl?gene=",y$Gene.Symbol," target = '_blank'",'>',y$Gene.Symbol,"</a>",sep='')
    datatable(y,
              class = "table-condensed",
              style = "bootstrap4",
              rownames = FALSE,
              escape = FALSE,
              options = list(autowidth = TRUE, scrollX = TRUE)
    )
  })

  output$vennDiagram <- renderPlot({
    dummy <- data.frame("No Genes Present")
    names(dummy) <- "Transcript.ID"
    if(identical(genelist2(),dummy)){return(NULL)}
    color.choices = c("blue", "green", "red", "orange","purple")
    color.choices <- color.choices[1:length(vendata())]
    grid.draw(venn.diagram(venndata(), filename = NULL, lwd = 1, col = color.choices, cat.cex = .9, fil = color.choices, margin = .12, ext.text = FALSE,
                           euler.d = TRUE))
  }, height = 400)

  output$downloadVennData <- downloadHandler(
    filename = function() {paste0("VennDiagram_", input$UorI, '.csv')  },
    content = function(file) {
      dummy <- data.frame("No Genes Present")
      names(dummy) <- "Transcript.ID"
      if(identical(genelist2(),dummy)){return(NULL)}
      write.csv(genelist2(), file,row.names = FALSE)
    }
  )

  output$downloadVennPic <- downloadHandler(
    filename = function() {paste0("VennDiagram_", input$UorI, '.png')  },
    content = function(file) {
      dummy <- data.frame("No Genes Present")
      names(dummy) <- "Transcript.ID"
      if(identical(genelist2(),dummy)){return(NULL)}
      color.choices = c("blue", "green", "red", "orange","purple")
      color.choices <- color.choices[1:length(vendata())]
      png(file)
      grid.draw(venn.diagram(venndata(), filename = NULL, lwd = 1, col = color.choices, cat.cex = .9, fil = color.choices, margin = .12, ext.text = FALSE,euler.d = TRUE))
      dev.off()
    }
  )


  output$downloadLMMModMap2 <- downloadHandler(
    filename = function() {paste('LMMModMap','.png', sep = '')},
    content = function(file){
      png(file, width = (resolution1()/72)*LMMmap1(), height = (resolution1()/72)*LMMmap2(), res = resolution1())
      if(min(lmm_module()$dat) >= 0 & max(lmm_module()$dat) > 0){
        print(aheatmap_circle(lmm_module()$dat, Rowv = lmm_module()$rowv, circleSize = input$LMMradius, fontsize = input$LMMFont, Colv = NA , color = lmm_module()$color_palette[500:1000]))
      }
      if(max(lmm_module()$dat) <= 0 & min(lmm_module()$dat) < 0){
        print(aheatmap_circle(lmm_module()$dat, Rowv = lmm_module()$rowv, circleSize = input$LMMradius, fontsize = input$LMMFont, Colv = NA , color = lmm_module()$color_palette[1:500]))
      }
      if(min(lmm_module()$dat) < 0 & max(lmm_module()$dat) >0){
        print(aheatmap_circle(lmm_module()$dat, Rowv = lmm_module()$rowv, circleSize = input$LMMradius, fontsize = input$LMMFont, Colv = NA , color = lmm_module()$color_palette))
      }
      if(all(lmm_module()$dat == 0)){
        return(NULL)
      }
      dev.off()
    }
  )

############## QUSAGE #####################

  output$qusage <- renderMenu({
    if(is.null(values$qusage.results)){
      return(strong(""))
    }
    if(is.null(values$qusage.results) == FALSE){
      menuItem("Q-Gen", icon = icon("filter"), tabName = "qusage")
    }
  })

  qgenRes <- reactive({
    req(values$qusage.results)
    qgenRes <- values$qusage.results
    colnames(qgenRes)[which(colnames(qgenRes) == "log.fold.change")] <- "Log2FC"
    bonf <- list()
    n <- length(unique(qgenRes$Comparison))
    for(i in 1:n){
      bonf[[i]] <- p.adjust(qgenRes$p.Value[which(qgenRes$Comparison == unique(qgenRes$Comparison)[i])], method = "bonferroni")
    }
    qgenRes$Bonf <- unlist(bonf)
    if(!is.null(values$annots)){
     n <- length(values$annots[,1])
     qgenRes$Modulev2_Annotation <- ""
     for(i in 1:n){
      qgenRes$Modulev2_Annotation[which(qgenRes$pathway.name %in% values$annots[,1][i])] <- as.character(values$annots[,2][i])
     }
    }
    qgenRes$Comparison <- as.factor(qgenRes$Comparison)
    return(qgenRes)
  })
  
  observe({
    req(qgenRes())
    updateSelectizeInput(session,"qgenComp", "Comparison(s) selection:", choices = unique(qgenRes()$Comparison), selected = unique(qgenRes()$Comparison)[1])
    updateSelectizeInput(session,"qgenSingleComp", "Comparison selection:", choices = unique(qgenRes()$Comparison), selected = unique(qgenRes()$Comparison)[1])
    updateSelectizeInput(session,"geneSetSelect", "Module selection:", choices = gtools::mixedsort(unique(as.character(qgenRes()$pathway.name))), 
                         selected = gtools::mixedsort(unique(as.character(qgenRes()$pathway.name)))[1])
  })

  output$PaloOrFirst1  <- renderUI({
    if(length(input$qgenComp) <= 1){return(NULL)}
    else{
      selectizeInput("PaloOrFirst", "Significant gene sets:", choices = c("In at least one comparison chosen" = 1, "In first comparison chosen" = 2), selected = 1)
    }
  })

  output$ColorChoice <- renderUI({
    if(length(input$qgenComp) <= 1){return(NULL)}
    else{
      cols = c("blue", "red","green","lightskyblue","purple","hotpink","brown","gold")
      selectizeInput("ColorChoice1", "Choose line colors:", choices = cols, selected = cols[1:length(input$qgenComp)], multiple = T)
    }
  })
  
  qgenResFilt <- reactive({
    input$qgenComp
    if(isTruthy(input$PaloOrFirst)){
      qgenRes <- callModule(filterOpts, "qusFcPlot", data = qgenRes, comparison = reactive(input$qgenComp), data.type = "geneSets", 
                            paloOrFirst = reactive(input$PaloOrFirst)) 
    } else{
      qgenRes <- callModule(filterOpts, "qusFcPlot", data = qgenRes, comparison = reactive(input$qgenComp), data.type = "geneSets",
                            paloOrFirst = reactive(NULL))
    }
    return(qgenRes)
  })

  venn.data <- reactive({
    req(qgenResFilt())
    dat <- qgenResFilt()
    comparisons <- as.character(unique(dat$Comparison))
    n <- length(comparisons)
    dat_sub <- list()
    for(i in 1:n){
      dat_sub[[i]] <- dat %>% filter(Comparison == comparisons[i], SIG == "Significant") %>% pull(pathway.name) 
    }
    names(dat_sub) <- comparisons
    return(dat_sub)
  })

  compOverview <- reactive({
    req(qgenRes())
    resOverview <- qgenRes() %>% 
      group_by(Comparison) %>%
      summarise(Raw = case_when(input$foldchange.q == "+" ~ sum(p.Value <= input$sigLevel & Log2FC > 0),
                                input$foldchange.q == "-" ~ sum(p.Value <= input$sigLevel & Log2FC < 0),
                                TRUE ~ sum(p.Value <= input$sigLevel)),
                FDR = case_when(input$foldchange.q == "+" ~ sum(FDR <= input$sigLevel & Log2FC > 0),
                                input$foldchange.q == "-" ~ sum(FDR <= input$sigLevel & Log2FC < 0),
                                TRUE ~ sum(FDR <= input$sigLevel)),
                Bonf = case_when(input$foldchange.q == "+" ~ sum(Bonf <= input$sigLevel & Log2FC > 0),
                                 input$foldchange.q == "-" ~ sum(Bonf <= input$sigLevel & Log2FC < 0),
                                 TRUE ~ sum(Bonf <= input$sigLevel))) %>%
      arrange(desc(Raw))
    return(resOverview)
  })

  singleGeneSetDat <- reactive({
    req(qgenRes())
    input$qgenSingleComp
    qgenRes <- qgenRes() %>% filter(Comparison == input$qgenSingleComp & pathway.name == input$geneSetSelect)
    lowerCI <- values$lower.ci
    upperCI <- values$upper.ci
    rownames(lowerCI) <- rownames(upperCI) <- as.character(values$results.file$Transcript.ID)
    geneLevelRes <- callModule(filterOpts, "qusIndFcPlot", data = reactive(values$results.file), comparison = reactive(input$qgenSingleComp), data.type = "genes")
    geneLevelRes <- geneLevelRes %>% filter(Transcript.ID %in% values$gene.sets[[input$geneSetSelect]])
    if(nrow(geneLevelRes) == 0){return(NULL)}
    lowerCI <- lowerCI[match(geneLevelRes$Transcript.ID, rownames(lowerCI), nomatch = 0), input$qgenSingleComp]
    upperCI <- upperCI[match(geneLevelRes$Transcript.ID, rownames(upperCI), nomatch = 0), input$qgenSingleComp]
    gene.symbols <- as.character(geneLevelRes %>% pull(Gene.Symbol))
    gene.ids <- as.character(geneLevelRes %>% pull(Transcript.ID))
    missing_symb <- unlist(sapply(gene.symbols, function(x) identical("",x)))
    if(sum(missing_symb) > 0){
      gene.symbols[missing_symb == TRUE] <- gene.ids[missing_symb == TRUE]
    }
    duplicates <- duplicated(gene.symbols)
    if(sum(duplicates) > 0){
      gene.symbols[duplicates == TRUE] <- gene.ids[duplicates == TRUE]
    }
    geneLevelRes$Gene.Symbol <- gene.symbols
    geneLevelRes <- cbind(geneLevelRes, lowerCI, upperCI)
    for(i in 1:ncol(geneLevelRes)){
      if(is.numeric(geneLevelRes[,i])){
        geneLevelRes[,i] <- as.numeric(formatC(geneLevelRes[,i], digits = 4))
      }
    }
    geneLevelRes <- geneLevelRes[order(geneLevelRes$Log2FC, decreasing = T),]
    geneLevelRes <- cbind(Index = c(1:nrow(geneLevelRes)), geneLevelRes)
    return(list(geneLevelRes = geneLevelRes, qgenRes = qgenRes))
  })
  
  geneSetLogFcPlot <- reactive({
    req(qgenResFilt())
    if(input$only_annotated){
      qgenResFilt <- qgenResFilt() %>% filter(Modulev2_Annotation != "")
    } else{
      qgenResFilt <- qgenResFilt()
    }
    qgenResFilt$ylow <- round(-max(c(abs(qgenResFilt$low), abs(qgenResFilt$up))),1)
    qgenResFilt$yhigh <- round(max(c(abs(qgenResFilt$low), abs(qgenResFilt$up))),1)
    if(input$showAnnotated){
      qgenResFilt$pathway.name <- factor(paste(qgenResFilt$pathway.name,qgenResFilt$Modulev2_Annotation, sep = " "), 
                                         levels = unique(paste(qgenResFilt$pathway.name,qgenResFilt$Modulev2_Annotation, sep = " ")))
    }
    if(length(input$qgenComp) > 1){
      n <- length(input$qgenComp)
      if(input$graphics6 == FALSE){
        line_colors = c("blue", "red","green","lightskyblue","purple","hotpink","brown","gold")
      } else{
        line_colors = input$ColorChoice1
      }
      qusage_plot <- ggplot(data = qgenResFilt, aes(x = pathway.name, y = Log2FC, ymin = ylow, ymax = yhigh)) + 
        scale_x_discrete(labels = as.character(unique(qgenResFilt$pathway.name))) + xlab("Modules") + ylab("Pathway Activity") + 
        scale_colour_manual(values = line_colors) + geom_line(aes(colour = Comparison, group = Comparison), size = 1) + 
        geom_hline(yintercept = 0, colour = "black", size = .5) + geom_point(aes(shape = SIG), size = 2) + scale_shape_manual(values = c(16,8))
    } else{
      qusage_plot <- ggplot(data = qgenResFilt, aes(x = pathway.name, y = Log2FC, ymin = ylow, ymax = yhigh)) + 
        ylab("Pathway Activity") + xlab("Modules") + geom_point() + geom_errorbar(aes(x = pathway.name, y = Log2FC, ymin = low, ymax = up)) + 
        geom_hline(yintercept = 0, colour = "red", size = 1)
    }
    qusage_plot <- qusage_plot + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
    return(qusage_plot)
  })

  singleGeneSetPlot <- reactive({
    req(singleGeneSetDat())
    geneLevelRes <- singleGeneSetDat()$geneLevelRes
    qgenRes <- singleGeneSetDat()$qgenRes
    geneLevelRes$ylow = round(-max(c(abs(min(geneLevelRes$lowerCI)), abs(max(geneLevelRes$upperCI)))), 1)
    geneLevelRes$yhigh = round(max(c(abs(min(geneLevelRes$lowerCI)), abs(max(geneLevelRes$upperCI)))), 1)
    geneLevelRes$Pathway.Activity = qgenRes$Log2FC
    
    geneSetPlot <- ggplot(data = geneLevelRes, aes(x = Index, y = Log2FC, ymin = ylow, ymax = yhigh)) + xlim(0, length(geneLevelRes$Index)) + 
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = qgenRes$low, ymax = qgenRes$up, alpha = .3, fill = "lightblue") + 
      scale_x_discrete(limits = as.character(geneLevelRes$Gene.Symbol[geneLevelRes$Index])) + geom_point() + 
      geom_errorbar(aes(x = Index, y = Log2FC, ymin = lowerCI, ymax = upperCI)) + xlab("Gene.Symbol ID") + ylab("Probe Level Log2 FC") + 
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_hline(yintercept = 0, colour = "red", size = 1) + 
      geom_hline(aes(yintercept = Pathway.Activity), linetype = "dashed",colour = "black") + theme(plot.margin = unit(c(1,10,0,0), "lines")) + 
      annotation_custom(textGrob(label = "Pathway Activity", hjust = 0), xmin = length(geneLevelRes$Index) + 1.08, 
                        xmax = length(geneLevelRes$Index) + 1.08, ymin = qgenRes$Log2FC, ymax = qgenRes$Log2FC)
    gt <- ggplot_gtable(ggplot_build(geneSetPlot))
    gt$layout$clip[gt$layout$name == "panel"] <- "off"
    gt
  })

  output$geneSetTable <- renderDT({
    req(singleGeneSetDat())
    dat <- singleGeneSetDat()$geneLevelRes
    dat$Gene.Symbol <- paste("<a href=http://www.genecards.org/cgi-bin/carddisp.pl?gene=",dat$Gene.Symbol," target = '_blank'",'>',dat$Gene.Symbol,"</a>",sep='')
    dat <- dat[,-1]
    datatable(dat,
              class = "table-condensed",
              style = "bootstrap4",
              rownames = FALSE,
              escape = FALSE,
              options = list(autowidth = TRUE, scrollX = TRUE)
    )
  })

  output$geneSetPlot <- renderPlot({
    req(singleGeneSetPlot())
    grid.draw(singleGeneSetPlot())
  },width = function(){input$PlotWidth2}, height = 400)

  output$compOverview <- renderDT({
    req(qgenRes())
    datatable(compOverview(),
              class = "table-condensed",
              style = "bootstrap4",
              rownames = FALSE,
              options = list(autowidth = TRUE, scrollX = TRUE)
    )
  })

  output$logFcPlot <- renderPlot({
    req(qgenRes())
    geneSetLogFcPlot()
  }, width = function(){input$PlotWidth1}, height = 400)

  output$venn <- renderUI({
    req(qgenRes())
    req(qgenResFilt())
    if(length(venn.data()) == 1){return(NULL)}
    else{
      plotOutput("venndiagram")
    }
  })

  output$venn.download <- renderUI({
    req(venn.data())
    if(length(venn.data()) == 1){return(NULL)}
    else{
      downloadButton("downloadVenn2", "Download Figure")
    }
  })

  vennPlot <- reactive({
    req(venn.data())
    color.choices <- c("blue", "red","green","lightskyblue","purple","hotpink","brown","gold")
    color.choices <- color.choices[1:length(venn.data())]
    if(input$graphics6){
      color.choices <- input$ColorChoice1
    }
    return(venn.diagram(venn.data(), filename = NULL, lwd = 1, col = color.choices, 
                        cat.cex = .9, fil = color.choices, margin = .12, ext.text = FALSE,euler.d = TRUE))
  })

  output$venndiagram <- renderPlot({
    req(qgenResFilt())
    grid.draw(vennPlot())
  }, height = 400)

  output$multipleCompTab <- renderDT({
    req(qgenResFilt())
    qgenResFilt <- qgenResFilt() %>% dplyr::rename(P.Value = p.Value)
    if(!is.null(values$annots)){
      if(input$only_annotated){
        qgenResFilt <- qgenResFilt %>% filter(Modulev2_Annotation != "")
      } 
    }
    qgenResFilt <- qgenResFilt[,-c(1,11)]
    for(i in 1:ncol(qgenResFilt)){
      if(is.numeric(qgenResFilt[,i])){
        qgenResFilt[,i] <- as.numeric(formatC(qgenResFilt[,i], digits = 4))
      }
    }
    datatable(qgenResFilt[,c(1,9,2,3,4,8,5,6,7)],
              class = "table-condensed",
              style = "bootstrap4",
              rownames = FALSE,
              options = list(autowidth = TRUE, scrollX = TRUE)
    )
  })

  output$downloadSigComps <- downloadHandler(
    filename = function() {paste(values$project.name,'_','Qusage_Sig_Comps_Table','.csv', sep='')  },
    content = function(file) {
      write.csv(compOverview(), file,row.names = FALSE)
    }
  )

  output$downloadPlot2 <- downloadHandler(
    filename = function() {paste(values$project.name,'_','Multi_Comparisons_Plot','.png', sep = '')},
    content = function(file){
      png(file, width = (plotres1()/72)*plotWidth1(), height = (plotres1()/72)*480, res = plotres1())
      print(geneSetLogFcPlot())
      dev.off()
    }
  )

  output$downloadVenn2 <- downloadHandler(
    filename = function() {paste(values$project.name,'_','Multi_Comparisons_Venn', '.png', sep = '')},
    content = function(file){
      png(file)
      grid.draw(vennPlot())
      dev.off()
    }
  )

  output$downloadTable2 <- downloadHandler(
    filename = function() {paste(values$project.name,'_','Qusage_Data_Table_Multi_Comp','.csv', sep='')  },
    content = function(file) {
      qgenResFilt <- qgenResFilt()[,-c(1,11)]
      qgenResFilt <- qgenResFilt[,c(1,9,2,3,4,8,5,6,7)]
      write.csv(qgenResFilt, file,row.names = FALSE)
    }
  )

  output$downloadPlot3 <- downloadHandler(
    filename = function() {paste(values$project.name,'_', 'Individual_GeneSet_Plot','.png', sep = '')},
    content = function(file){
      png(file, width = (plotres2()/72)*plotWidth2(), height = (plotres2()/72)*480, res = plotres2())
      grid.draw(singleGeneSetPlot())
      dev.off()
    }
  )

  output$downloadTable3 <- downloadHandler(
    filename = function() {paste(values$project.name,'_','Qusage_Individual_GeneSet_Data_Table','.csv', sep='')  },
    content = function(file) {
      write.csv(singleGeneSetDat()$geneLevelRes[,-1], file,row.names = FALSE)
    }
  )

  ################### ROAST ############################

  output$roast <- renderMenu({
    if(is.null(values$roast.results)){
      return(strong(""))
    }
    if(is.null(values$roast.results) == FALSE){
      menuItem("Roast", icon = icon("th-list"), tabName = "roast")
    }
  })

  observe({
    req(values$roast.results)
    updateSelectizeInput(session,"setStat", "Select gene set statistic:", choices = names(values$roast.results), selected = names(values$roast.results)[1])
  })

  roast.overview <- reactive({
    results <- values$roast.results[[which(names(values$roast.results) %in% input$setStat)]]
    Comparison <- names(results)
    Raw <- FDR <- Bonf <- c()
    for(i in 1:length(results)){
      Bonferroni <- p.adjust(results[[i]]$PValue, method = "bonferroni")
      Raw[i] <- length(which(results[[i]]$PValue <= input$sigLevelR))
      FDR[i] <- length(which(results[[i]]$FDR <= input$sigLevelR))
      Bonf[i] <- length(which(Bonferroni <= input$sigLevelR))
    }
    overview <- data.frame(Comparison = Comparison, Raw = Raw, FDR = FDR, Bonf = Bonf)
    overview
  })

  roast.results <- reactive({
    results <- values$roast.results[[which(names(values$roast.results) %in% input$setStat)]]
    for(i in 2:length(results)){
      results[[i]] <- results[[i]][match(rownames(results[[1]]), rownames(results[[i]]), nomatch = 0),]
    }
    Gene.set <- rownames(results[[1]])
    results <- do.call("cbind", results)
    results <- cbind(Gene.set, results)
    results
  })

  output$compOverviewR <- renderDT({
    datatable(roast.overview(),
              class = "table-condensed",
              style = "bootstrap4",
              rownames = FALSE,
              options = list(autowidth = TRUE, scrollX = TRUE)
    )
  })

  ###############Correlations########################
  
  output$correlations <- renderMenu({
    if(is.null(values$corrs)){
      return(strong(""))
    }
    else{
      return(menuItem("Correlations", icon = icon("bolt"), tabName = "corr"))
    }
  })
  
  output$TypeVariable <- renderUI({
    selectizeInput("TypeVariable1", "Choose Correlation Data:", choices = values$corr.names, selected = values$corr.names[1])
  })
  
  output$TypeVariable2 <- renderUI({
    selectizeInput("TypeVariable3", "Choose Correlation Data:", choices = values$corr.names, selected = values$corr.names[1])
  })
  
  output$TypeVariable4 <- renderUI({
    selectizeInput("TypeVariable5", "Choose Correlation Data:", choices = values$corr.names, selected = values$corr.names[1])
  })
  
  output$TypeVariable6 <- renderUI({
    selectizeInput("TypeVariable7", "Choose Correlation Data:", choices = values$corr.names, selected = values$corr.names[1])
  })
  
  output$var_switch <- renderUI({
    checkboxInput("var_switch", strong(paste0("Change Figure Rows to ", values$x.var[which(values$corr.names == input$TypeVariable1)])), value = FALSE)
  })
  
  corrResultsOverview <- reactive({
    correlations <- results <- Raw <- FDR <- Bonf <- list()
    for(i in 1:length(values$corrs)){
      corrs <- data.table::as.data.table(values$corrs[[i]][,c(2,1,3:ncol(values$corrs[[i]]))])
      corrs$FDR <- 1
      corrs$Bonf <- 1
      corrs[,FDR := p.adjust(Raw.P.Value, method = "fdr"), by = c(names(corrs)[c(2,3)])]
      corrs[,Bonf := p.adjust(Raw.P.Value, method = "bonferroni"), by = c(names(corrs)[c(2,3)])]
      corrs <- as.data.frame(corrs)
      corrs <- corrs[,c(2,1,3:ncol(corrs))]
      comps <- raw <- fdr <- bonf <- c()
      for(j in 1:length(unique(values$corrs[[i]][,3]))){
        comps[j] <- paste0(values$corr.names[i],"_",as.character(unique(values$corrs[[i]][,3])[j]))
        corrs.sub <- corrs[corrs[,3] == unique(corrs[,3])[j],]
        if(input$overviewCorrVal){
          if(input$selectCorrSign == "+"){
            corrs.sub <- corrs.sub[corrs.sub[,4] > input$selectCorrVal,]
          } else if(input$selectCorrSign == "-"){
            corrs.sub <- corrs.sub[corrs.sub[,4] < -input$selectCorrVal,]
          } else{
            corrs.sub <- corrs.sub[corrs.sub[,4] < -input$selectCorrVal | corrs.sub[,4] > input$selectCorrVal,]
          }
        }
        raw[j] <- length(which(corrs.sub$Raw.P.Value <= input$alphaCorr))
        fdr[j] <- length(which(corrs.sub$FDR <= input$alphaCorr))
        bonf[j] <- length(which(corrs.sub$Bonf <= input$alphaCorr))
      }
      results[[i]] <- data.frame(Correlation_Comparisons = comps, Raw = raw, FDR = fdr, Bonf = bonf)
      correlations[[i]] <- corrs
    }
    results <- do.call("rbind", results)
    z <- list(correlations = correlations, results = results)
  })
  
  output$corrOverview <- renderDT({
    withProgress(message = 'Making Results Table',
                 detail = 'This may take a while...', value = .1,{
                   datatable(corrResultsOverview()$results[order(corrResultsOverview()$results$Raw, decreasing = TRUE),],
                             class = "table-condensed",
                             style = "bootstrap4",
                             rownames = FALSE,
                             options = list(autowidth = TRUE, scrollX = TRUE)
                   )
                 }
    )
  })
  
  output$downloadOverviewCorrTable <- downloadHandler(
    filename = function() {paste(values$project.name, "_", "Correlation_Results_Overview_Table.csv")},
    content = function(file) {
      write.csv(corrResultsOverview()$results[order(corrResultsOverview()$results$Raw, decreasing = TRUE),], file, row.names = FALSE)
    }
  )
  
  output$subsetcorr <- renderUI({
    correlations <- values$corrs[[which(values$corr.names == input$TypeVariable1)]]
    if(input$var_switch == TRUE){
      correlations <- correlations[,c(2,1,3:ncol(correlations))]
      colnames(correlations) <- colnames(correlations)[c(2,1,3:ncol(correlations))]
    }
    visit.name = colnames(correlations)[3]
    visit <- unique(correlations[,3])
    visit <- gtools::mixedsort(as.character(visit))
    selectizeInput("subsetcorr1", "", choices = visit, selected = visit, multiple = TRUE)
  })
  
  output$WithVariable1 <- renderUI({
    correlations <- values$corrs[[which(values$corr.names == input$TypeVariable1)]]
    if(input$var_switch == TRUE){
      correlations <- correlations[,c(2,1,3:ncol(correlations))]
      colnames(correlations) <- colnames(correlations)[c(2,1,3:ncol(correlations))]
    }
    selectizeInput("WithVariable", "Choose Correlation Variable(s):", choices = as.character(gtools::mixedsort(unique(correlations$With))), selected = gtools::mixedsort(unique(correlations$With))[1], multiple = TRUE)
  })
  
  output$WithVariable2 <- renderUI({
    correlations <- values$corrs[[which(values$corr.names == input$TypeVariable3)]]
    if(input$var_switch2 == TRUE){
      correlations <- correlations[,c(2,1,3:ncol(correlations))]
      colnames(correlations) <- colnames(correlations)[c(2,1,3:ncol(correlations))]
    }
    selectizeInput("WithVariable3", "Choose Correlation Variable(s):", choices = as.character(gtools::mixedsort(unique(correlations$With))), selected = as.character(gtools::mixedsort(unique(correlations$With)))[1])
  })
  
  output$corr_Variable <- renderUI({
    correlations <- values$corrs[[which(values$corr.names == input$TypeVariable7)]]
    selectizeInput("corr_Variable2", paste("Choose", values$y.var[which(values$corr.names == input$TypeVariable7)], "Variable:", sep = " "), choices = as.character(gtools::mixedsort(unique(correlations$Variable))), selected = as.character(gtools::mixedsort(unique(correlations$With)))[1])
  })
  
  output$WithVariable4 <- renderUI({
    correlations <- values$corrs[[which(values$corr.names == input$TypeVariable7)]]
    selectizeInput("WithVariable5", paste("Choose ", values$x.var[which(values$corr.names == input$TypeVariable7)], " Variable:", sep = ""), choices = as.character(gtools::mixedsort(unique(correlations$With))), selected = as.character(gtools::mixedsort(unique(correlations$With)))[1])
  })
  
  output$uploadmod1 <- renderUI({
    if(input$var_switch == FALSE){
      return(checkboxInput("uploadmod", strong(paste("Upload ",values$y.var[which(values$corr.names == input$TypeVariable1)], " List:", sep = ""), style = "color:blue"), FALSE))
    } else{
      return(checkboxInput("uploadmod", strong(paste("Upload ",values$x.var[which(values$corr.names == input$TypeVariable1)], " List:", sep = ""), style = "color:blue"), FALSE))
    }
  })
  
  output$fileupload1 <- renderUI({
    if(input$uploadmod == TRUE){
      return(fileInput('modselect', '', accept = ".csv"))
    }
    else{
      return(NULL)
    }
  })
  
  output$uploadmod3 <- renderUI({
    if(input$var_switch3 == FALSE){
      return(checkboxInput("uploadmod2", strong(paste("Upload ",values$y.var[which(values$corr.names == input$TypeVariable5)], " List:", sep = ""), style = "color:blue"), FALSE))
    } else{
      return(checkboxInput("uploadmod2", strong(paste("Upload ",values$x.var[which(values$corr.names == input$TypeVariable5)], " List:", sep = ""), style = "color:blue"), FALSE))
    }
  })
  
  output$fileupload2 <- renderUI({
    if(input$uploadmod2 == TRUE){
      return(fileInput('modselect2', '', accept = ".csv"))
    }
    else{
      return(NULL)
    }
  })
  
  heatmap.data <- reactive({
    correlations <- corrResultsOverview()$correlations[[which(values$corr.names == input$TypeVariable1)]]
    if(input$var_switch == TRUE){
      correlations <- correlations[,c(2,1,3:ncol(correlations))]
      colnames(correlations) <- colnames(correlations)[c(2,1,3:ncol(correlations))]
    }
    if(input$uploadmod == TRUE){
      modnames.corr <- read.csv(input$modselect$datapath, header = TRUE)
      modnames.corr <- as.character(modnames.corr[,1])
      modnames.corr <- gsub(".", "_", modnames.corr,fixed = TRUE)
      correlations <- correlations[which(correlations$Variable %in% modnames.corr),]
      
      col_anno <- list()
      col_anno2 <- list()
      dat <- list()
      dat2 <- list()
      keep_names <- list()
      
      names <- gtools::mixedsort(unique(as.character(correlations[,3])))
      
      if(input$subsetModcorr == TRUE){
        names <- input$subsetcorr1
      }
      
      for(j in 1:length(input$WithVariable)){
        correlations1 <- correlations[which(correlations$With %in% input$WithVariable[j]),]
        for(i in 1:length(names)){
          dat[[i]] <- correlations1[which(correlations1[,3] == names[i]),]
          dat[[i]] <- dat[[i]][which(dat[[i]]$Variable %in% modnames.corr),]
          dat[[i]] <- as.character(dat[[i]]$Variable)
        }
        
        keep_names[[j]] <- Reduce(intersect, dat)
      }
      
      keep_names <- unlist(keep_names)
      
      for(j in 1:length(input$WithVariable)){
        correlations1 <- correlations[which(correlations$With %in% input$WithVariable[j]),]
        for(i in 1:length(names)){
          dat[[i]] <- correlations1[which(correlations1[,3] == names[i]),]
          dat[[i]] <- dat[[i]][which(dat[[i]]$Variable %in% keep_names),]
          dat[[i]] <- dat[[i]][order(dat[[i]]$Variable),]
          dat[[i]] <- dat[[i]][,4]
          col_anno[[i]] <- input$WithVariable[j]
        }
        
        col_anno2[[j]] <- unlist(col_anno)
        dat2[[j]] <- do.call("cbind", dat)
        rownames(dat2[[j]]) <- sort(unique(keep_names))
        colnames(dat2[[j]]) <- names
      }
      
      col_anno <- data.frame(unlist(col_anno2))
      colnames(col_anno)[1] <- paste(values$x.var[which(values$corr.names == input$TypeVariable1)], "variable", sep = " ")
      if(input$var_switch){
        colnames(col_anno)[1] <- paste(values$y.var[which(values$corr.names == input$TypeVariable1)], "variable", sep = " ")
      }
      dat <- do.call("cbind", dat2)
      z <- list(col_anno = col_anno, dat = dat)
    } else{
      dat3 <- list()
      colnames3 <- list()
      rownames3 <- list()
      
      for(j in 1:length(input$WithVariable)){
        correlations1 <- correlations[which(correlations$With %in% input$WithVariable[j]),]
        
        dat <- list()
        sig <- list()
        cor_vals <- list()
        names <- gtools::mixedsort(unique(as.character(correlations1[,3])))
        
        if(input$subsetModcorr == TRUE){
          names <- input$subsetcorr1
        }
        
        if(input$corrval2 == TRUE){
          if(input$corrsign1 == "+"){
            for(i in 1:length(names)){
              dat[[i]] <- correlations1[which(correlations1[,3] == names[i]),]
              cor_vals[[i]] <- dat[[i]][which(dat[[i]][,4] >= input$corrval3),]
              cor_vals[[i]] <- as.character(cor_vals[[i]]$Variable)
            }
          }
          if(input$corrsign1 == "-"){
            for(i in 1:length(names)){
              dat[[i]] <- correlations1[which(correlations1[,3] == names[i]),]
              cor_vals[[i]] <- dat[[i]][which(dat[[i]][,4] <= -input$corrval3),]
              cor_vals[[i]] <- as.character(cor_vals[[i]]$Variable)
            }
          }
          if(input$corrsign1 == "Both"){
            for(i in 1:length(names)){
              dat[[i]] <- correlations1[which(correlations1[,3] == names[i]),]
              cor_vals[[i]] <- dat[[i]][c(which(dat[[i]][,4] <= -input$corrval3),which(dat[[i]][,4] >= input$corrval3)),]
              cor_vals[[i]] <- as.character(cor_vals[[i]]$Variable)
            }
          }
          
          cor_vals <- Reduce(union,cor_vals)
          correlations1 <- correlations1[which(correlations1$Variable %in% cor_vals),]
        }
        
        for(i in 1:length(names)){
          dat[[i]] <- correlations1[which(correlations1[,3] == names[i]),]
          
          if(input$correction_method.corr == "Raw"){
            sig[[i]] <- dat[[i]][which(dat[[i]]$Raw.P.Value <= input$Alpha1),]
          }
          
          if(input$correction_method.corr == "FDR"){
            sig[[i]] <- dat[[i]][which(dat[[i]]$FDR <= input$Alpha1),]
          }
          
          if(input$correction_method.corr == "Bonferroni"){
            sig[[i]] <- dat[[i]][which(dat[[i]]$Bonf <= input$Alpha1),]
          }
          
          sig[[i]] <- as.character(sig[[i]]$Variable)
        }
        
        sigvars <- Reduce(union,sig)
        
        correlations1 <- correlations1[which(correlations1$Variable %in% sigvars),]
        
        dat2 <- list()
        for(i in 1:length(names)){
          dat2[[i]] <- correlations1[which(correlations1[,3] == names[i]),]
          dat2[[i]] <- dat2[[i]][which(dat2[[i]]$Variable %in% sigvars),]$Variable
        }
        
        dat2 <- Reduce(intersect, dat2)
        
        for(i in 1:length(names)){
          dat[[i]] <- correlations1[which(correlations1[,3] == names[i]),]
          dat[[i]] <- dat[[i]][which(dat[[i]]$Variable %in% dat2),]
          dat[[i]] <- dat[[i]][order(dat[[i]]$Variable),]
          dat[[i]] <- dat[[i]][,4]
        }
        
        dat3[[j]] <- do.call("cbind", dat)
        colnames(dat3[[j]]) <- names
        rownames(dat3[[j]]) <- sort(dat2)
        colnames3[[j]] <- names
        rownames3[[j]] <- sort(dat2)
      }
      rownames3 <- Reduce(union, rownames3)
      dat2 <- list()
      col_anno <- list()
      col_anno2 <- list()
      
      for(j in 1:length(input$WithVariable)){
        correlations1 <- correlations[which(correlations$With %in% input$WithVariable[j]),]
        dat <- list()
        for(i in 1:length(names)){
          dat[[i]] <- correlations1[which(correlations1[,3] == names[i]),]
          dat[[i]] <- dat[[i]][which(dat[[i]]$Variable %in% rownames3),]
          dat[[i]] <- dat[[i]][order(dat[[i]]$Variable),]
          dat[[i]] <- dat[[i]][,4]
          col_anno[[i]] <- input$WithVariable[j]
        }
        
        col_anno2[[j]] <- unlist(col_anno)
        dat2[[j]] <- do.call("cbind",dat)
        rownames(dat2[[j]]) <- sort(rownames3)
        colnames(dat2[[j]]) <- names
      }
      dat <- do.call("cbind", dat2)
      col_anno <- as.data.frame(unlist(col_anno2))
      colnames(col_anno)[1] <- paste(values$x.var[which(values$corr.names == input$TypeVariable1)], "variable", sep = " ")
      if(input$var_switch){
        colnames(col_anno)[1] <- paste(values$y.var[which(values$corr.names == input$TypeVariable1)], "variable", sep = " ")
      }
      z <- list(col_anno = col_anno, dat = dat)
    }
    
    return(z)
  })
  
  output$download_heatmap_data <- downloadHandler(
    filename = function() {paste(values$project.name, "_", "Correlation_heatmap_data.csv",sep = "")},
    content = function(file) {
      write.csv(heatmap.data()$dat, file, row.names = TRUE)
    }
  )
  
  row_clust.corr <- reactive({
    if(input$rowclustcorr == TRUE){return(TRUE)}
    else{
      return(NA)
    }
  })
  
  col_clust.corr <- reactive({
    if(input$colclustcorr == TRUE){return(TRUE)}
    else{
      return(NA)
    }
  })
  
  rows.plot <- reactive({
    if(nrow(heatmap.data()$dat) > 300){return(NA)}
    else{return(NULL)}
  })
  
  color.heatmap <- reactive({
    color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
    color_palette[c(48:52)] = "#FFFFFF"
    color_palette
  })
  
  plotresolution.corr <- reactive({
    input$PlotRes4
  })
  
  output$correlations_plotOverview <- renderPlot({
    if(ncol(heatmap.data()$dat) > 1){
      if(input$colorRange){
        return(aheatmap(heatmap.data()$dat,Colv = col_clust.corr(),annCol = heatmap.data()$col_anno,labRow = rows.plot(),treeheight = input$corrTreeHeight, fontsize = input$corrFontSize, 
                         color = color.heatmap(), border_color = "grey60", Rowv = row_clust.corr(), breaks = seq(-1,1,by=.04)))
      }
      return(aheatmap(heatmap.data()$dat,Colv = col_clust.corr(),annCol = heatmap.data()$col_anno,labRow = rows.plot(),treeheight = input$corrTreeHeight, fontsize = input$corrFontSize, color = color.heatmap(), 
                      border_color = "grey60", Rowv = row_clust.corr(), breaks = 0))
    }
    else{
      if(input$colorRange){
        return(aheatmap(heatmap.data()$dat,Colv = col_clust.corr(),labRow = rows.plot(), treeheight = input$corrTreeHeight, fontsize = input$corrFontSize, color = color.heatmap(), border_color = "grey60", 
                        Rowv = row_clust.corr(), breaks = seq(-1,1,by=.04)))
      }
      return(aheatmap(heatmap.data()$dat,Colv = col_clust.corr(),labRow = rows.plot(),treeheight = input$corrTreeHeight, fontsize = input$corrFontSize, color = color.heatmap(), border_color = "grey60", 
                      Rowv = row_clust.corr(), breaks = c(-1,0,1)))
    }
  }, height = function(){input$PlotHeight4}, width = function(){input$PlotWidth4})
  
  output$download_corr_heatmap <- downloadHandler(
    filename = function() {paste(values$project.name,"_","Heatmap.png",sep = "")},
    content = function(file){
      png(file, width = (plotresolution.corr()/72)*width1(), height = (plotresolution.corr()/72)*height1(), res = plotresolution.corr())
      if(input$colorRange){
        print(aheatmap(heatmap.data()$dat,Colv = col_clust.corr(),annCol = heatmap.data()$col_anno,labRow = rows.plot(),treeheight = input$corrTreeHeight, fontsize = input$corrFontSize,  
                        color = color.heatmap(), border_color = "grey60", Rowv = row_clust.corr(), breaks = seq(-1,1,by=.04)))
      } else{
        print(aheatmap(heatmap.data()$dat,Colv = col_clust.corr(),annCol = heatmap.data()$col_anno,labRow = rows.plot(),treeheight = input$corrTreeHeight, fontsize = input$corrFontSize, color = color.heatmap(), 
                       border_color = "grey60", Rowv = row_clust.corr(), breaks = 0))
      }
      dev.off()
    }
  )
  
  output$Visit <- renderUI({
    
    correlations <- values$corrs[[which(values$corr.names == input$TypeVariable3)]]
    if(input$var_switch2 == TRUE){
      correlations <- correlations[,c(2,1,3:ncol(correlations))]
      colnames(correlations) <- colnames(correlations)[c(2,1,3:ncol(correlations))]
    }
    visit.name = colnames(correlations)[3]
    visit <- unique(correlations[,3])
    visit <- gtools::mixedsort(as.character(visit))
    selectizeInput("visit", paste("Correlations by"," ",visit.name, ":", sep = ""), choices = visit, selected = visit[1])
  })
  
  sub.dat <- reactive({
    correlations <- corrResultsOverview()$correlations[[which(values$corr.names == input$TypeVariable3)]]
    if(input$var_switch2){
      correlations <- correlations[,c(2,1,3:ncol(correlations))]
      colnames(correlations) <- colnames(correlations)[c(2,1,3:ncol(correlations))]
    }
    correlations <- correlations[which(correlations$With == input$WithVariable3),]
    correlations <- correlations[which(correlations[,3] == input$visit),]
    z <- list(correlations = correlations)
    return(z)
  })
  
  subset_correlations <- reactive({
    correlations <- sub.dat()$correlations
    
    if(input$corrval == TRUE){
      if(input$corrsign == "+"){
        correlations <- correlations[which(correlations[,4] >= input$corrval1),]
      }
      if(input$corrsign == "-"){
        correlations <- correlations[which(correlations[,4] <= -input$corrval1),]
      }
      if(input$corrsign == "Both"){
        correlations <- correlations[c(which(correlations[,4] <= -input$corrval1), which(correlations[,4] >= input$corrval1)),]
      }
    }
    
    if(input$correction_method.corr2 == "Raw"){
      correlations <- correlations[which(correlations$Raw.P.Value <= input$Alpha),]
    }
    
    if(input$correction_method.corr2 == "FDR"){
      correlations <- correlations[which(correlations$FDR <= input$Alpha),]
    }
    
    if(input$correction_method.corr2 == "Bonferroni"){
      correlations <- correlations[which(correlations$Bonf <= input$Alpha),]
    }
    
    correlations <- correlations[,-c(which(colnames(correlations) == "Base_subtracted"), which(colnames(correlations) == "NObs"),
                                     which(colnames(correlations) == "Sign_NegLog10_p"))]
    colnames(correlations)[1] <- values$y.var[which(values$corr.names == input$TypeVariable3)]
    colnames(correlations)[2] <- values$x.var[which(values$corr.names == input$TypeVariable3)]
    if(input$var_switch2){
      colnames(correlations)[1] <- values$x.var[which(values$corr.names == input$TypeVariable3)]
      colnames(correlations)[2] <- values$y.var[which(values$corr.names == input$TypeVariable3)]
    }
    return(correlations)
  })
  
  output$download_data <- downloadHandler(
    filename = function() {paste(values$project.name, "_", "Correlations_subset.csv")},
    content = function(file) {
      write.csv(subset_correlations(), file, row.names = FALSE)
    }
  )
  
  output$var_switch2 <- renderUI({
    checkboxInput("var_switch2", strong(paste("Swap ",values$y.var[which(values$corr.names == input$TypeVariable3)], " and ", values$x.var[which(values$corr.names == input$TypeVariable3)], " columns", sep = "")),
                  value = FALSE)
  })
  
  output$correlation_table <- renderDT({
    datatable(subset_correlations(),
              class = "table-condensed",
              style = "bootstrap4",
              rownames = FALSE,
              options = list(autowidth = TRUE, scrollX = TRUE)
    )
  })
  
  output$Visit2 <- renderUI({
    correlations <- values$corrs[[which(values$corr.names == input$TypeVariable5)]]
    if(input$var_switch3 == TRUE){
      correlations <- correlations[,c(2,1,3:ncol(correlations))]
      colnames(correlations) <- colnames(correlations)[c(2,1,3:ncol(correlations))]
    }
    visit.name = colnames(correlations)[3]
    visit <- unique(correlations[,3])
    visit <- gtools::mixedsort(as.character(visit))
    selectizeInput("visit2", paste("Correlations by", visit.name, ":", sep = " "), choices = visit, selected = visit[1])
  })
  
  output$subsetcorr2 <- renderUI({
    correlations <- values$corrs[[which(values$corr.names == input$TypeVariable5)]]
    if(input$var_switch3 == TRUE){
      correlations <- correlations[,c(2,1,3:ncol(correlations))]
      colnames(correlations) <- colnames(correlations)[c(2,1,3:ncol(correlations))]
    }
    var <- as.character(sort(unique(correlations$With)))
    selectizeInput("subsetcorr3", paste("Subset by 'with' variable:",sep = ""), choices = var, selected = var[1:5], multiple = TRUE)
  })
  
  correlation.data <- reactive({
    #correlations <- values$corrs[[which(values$corr.names == input$TypeVariable5)]]
    correlations <- corrResultsOverview()$correlations[[which(values$corr.names == input$TypeVariable5)]]
    if(input$var_switch3 == TRUE){
      correlations <- correlations[,c(2,1,3:ncol(correlations))]
      colnames(correlations) <- colnames(correlations)[c(2,1,3:ncol(correlations))]
    }
    correlations <- correlations[which(correlations[,3] == input$visit2),]
    
    if(input$subsetModcorr2 == TRUE){
      correlations <- correlations[which(correlations$With %in% input$subsetcorr3),]
    }
    
    if(input$var_switch3 == TRUE){
      y <- values$x.var[which(values$corr.names == input$TypeVariable5)]
      x <- values$y.var[which(values$corr.names == input$TypeVariable5)]
    }
    else{
      y <- values$y.var[which(values$corr.names == input$TypeVariable5)]
      x <- values$x.var[which(values$corr.names == input$TypeVariable5)]
    }
    dat <- list()
    cor_vals <- list()
    sig <- list()
    names <- as.character(unique(correlations$With))
    
    if(input$corrval4 == TRUE){
      if(input$corrsign2 == "+"){
        for(i in 1:length(names)){
          dat[[i]] <- correlations[which(correlations$With == names[i]),]
          cor_vals[[i]] <- dat[[i]][which(dat[[i]][,4] >= input$corrval5),]
          cor_vals[[i]] <- as.character(cor_vals[[i]]$Variable)
        }
      }
      if(input$corrsign2 == "-"){
        for(i in 1:length(names)){
          dat[[i]] <- correlations[which(correlations$With == names[i]),]
          cor_vals[[i]] <- dat[[i]][which(dat[[i]][,4] <= -input$corrval5),]
          cor_vals[[i]] <- as.character(cor_vals[[i]]$Variable)
        }
      }
      if(input$corrsign2 == "Both"){
        for(i in 1:length(names)){
          dat[[i]] <- correlations[which(correlations$With == names[i]),]
          cor_vals[[i]] <- dat[[i]][c(which(dat[[i]][,4] <= -input$corrval5),which(dat[[i]][,4] >= input$corrval5)),]
          cor_vals[[i]] <- as.character(cor_vals[[i]]$Variable)
        }
      }
      
      cor_vals <- Reduce(union,cor_vals)
      correlations <- correlations[which(correlations$Variable %in% cor_vals),]
    }
    
    
    for(i in 1:length(names)){
      dat[[i]] <- correlations[which(correlations$With == names[i]),]
      if(input$correction_method.corr3 == "Raw"){
        sig[[i]] <- dat[[i]][which(dat[[i]]$Raw.P.Value <= input$Alpha2),]
      }
      if(input$correction_method.corr3 == "FDR"){
        sig[[i]] <- dat[[i]][which(dat[[i]]$FDR <= input$Alpha2),]
      }
      if(input$correction_method.corr3 == "Bonferroni"){
        sig[[i]] <- dat[[i]][which(dat[[i]]$Bonf <= input$Alpha2),]
      }
      sig[[i]] <- as.character(sig[[i]]$Variable)
    }
    
    sigvars <- Reduce(union,sig)
    
    correlations <- correlations[which(correlations$Variable %in% sigvars),]
    
    correlation.type <- colnames(correlations)[4]
    
    if(input$uploadmod2 == TRUE){
      modnames.corr <- read.csv(input$modselect2$datapath, header = TRUE)
      modnames.corr <- as.character(modnames.corr[,1])
      modnames.corr <- gsub(".", "_", modnames.corr,fixed = TRUE)
      correlations <- correlations[which(correlations$Variable %in% modnames.corr),]
    }
    
    n = length(unique(correlations$With))
    cor <- list()
    cor_sub <- list()
    names <- unique(correlations$With)
    
    for(i in 1:n){
      cor[[i]] <- correlations[which(correlations$With == names[i]),]
      
      if(i > 1){
        cor[[i]] <- cor[[i]][match(cor[[i-1]]$Variable, cor[[i]]$Variable, nomatch = 0),]
      }
      
      if(input$correction_method.corr3 == "Raw"){
        cor_sub[[i]] <- cor[[i]][,which(colnames(cor[[i]]) %in% c(correlation.type, "Raw.P.Value"))]
      }
      if(input$correction_method.corr3 == "FDR"){
        cor_sub[[i]] <- cor[[i]][,which(colnames(cor[[i]]) %in% c(correlation.type, "FDR"))]
      }
      if(input$correction_method.corr3 == "Bonferroni"){
        cor_sub[[i]] <- cor[[i]][,which(colnames(cor[[i]]) %in% c(correlation.type, "Bonf"))]
      }
      colnames(cor_sub[[i]]) <- c(paste(names[i], "r.val", sep = "."), paste(names[i], "pVal", sep = "."))
    }
    
    mat1 <- do.call("cbind", cor_sub)
    mat1 <- as.matrix(mat1)
    rownames(mat1) <- cor[[1]]$Variable
    mat1 <- as.data.frame(mat1)
    
    t = mat1 #mat1 is what you created with the correlation script.
    t$modName<-rownames(t)
    t.pval=reshape2::melt(t[c(length(t),grep("pVal",colnames(t), fixed = TRUE))]) #combine all columns with "pVal" into one loooooong object.
    t.rval=reshape2::melt(t[c(length(t),grep("r.val",colnames(t), fixed = TRUE))]) #combine all columns with "r.val" into one loooooong object.
    GraphFrame=cbind(t.pval,t.rval[,3])  #concaternate them into one matrix
    if(input$correction_method.corr3 == "Raw"){
      colnames(GraphFrame)=c(y,"Correlation","Raw.pVal","rVal")  #change column names
    }
    if(input$correction_method.corr3 == "FDR"){
      colnames(GraphFrame)=c(y,"Correlation","FDR.pVal","rVal")  #change column names
    }
    if(input$correction_method.corr3 == "Bonferroni"){
      colnames(GraphFrame)=c(y,"Correlation","Bonf.pVal","rVal")  #change column names
    }
    
    GraphFrame[,which(colnames(GraphFrame) == y)] = factor(GraphFrame[,which(colnames(GraphFrame) == y)],levels=unique(GraphFrame[,which(colnames(GraphFrame) == y)],order=T))  #this is super important as it makes sure the order on the axis is correct!
    GraphFrame$Correlation <- gsub(".pVal", "",GraphFrame$Correlation, fixed = TRUE)
    z <- list(GraphFrame = GraphFrame, y = y, x = x)
    return(z)
  })
  
  output$download_plot_data <- downloadHandler(
    filename = function() {paste(values$project.name, "_", "correlations_data.csv")},
    content = function(file) {
      write.csv(correlation.data()$GraphFrame, file, row.names = FALSE)
    }
  )
  
  plotresolution.corr1 <- reactive({
    input$PlotRes3
  })
  
  correlations_makePlot <- reactive({
    y <- correlation.data()$y
    x <- correlation.data()$x
    GraphFrame <- correlation.data()$GraphFrame
    if(input$correction_method.corr3 == "Raw"){
      GraphFrame <- GraphFrame[which(GraphFrame$Raw.pVal <= input$Alpha2),]
      p <- GraphFrame$Raw.pVal
    }
    if(input$correction_method.corr3 == "FDR"){
      GraphFrame <- GraphFrame[which(GraphFrame$FDR.pVal <= input$Alpha2),]
      p <- GraphFrame$FDR.pVal
    }
    if(input$correction_method.corr3 == "Bonferroni"){
      GraphFrame <- GraphFrame[which(GraphFrame$Bonf.pVal <= input$Alpha2),]
      p <- GraphFrame$Bonf.pVal
    }
    p[p < .0001] <- .0001 #set log-10 transformed p-values <4 (=0.0001) to 4
    p[p > .1] <- .1
    RV <- ifelse(GraphFrame$rVal>0,GraphFrame$rVal^4,(GraphFrame$rVal^4)*-1)
    
    if(length(p)/length(unique(GraphFrame$Correlation)) <= 260){
      if(input$Alpha2 > .05 & input$Alpha2 <= .075){
        graph <- ggplot(data = GraphFrame, aes(Correlation, GraphFrame[,1],colour=RV, size = p), environment = environment()) + labs(y = y, x = x)  + 
          scale_size_continuous("P-value",range=c(12,4),limits = c(min(p), .1), breaks = c(0.0001,0.025,0.05,0.075), labels = c(".0001", "0.025", "0.05", "0.075")) + 
          geom_point(colour="black",aes(size = p),shape=21,alpha=I(1)) + geom_point(alpha=I(.75))  + theme_minimal()  +
          theme(panel.background=element_rect(fill="white"),
                panel.grid = element_line(),panel.grid.major.y=element_blank(),
                panel.grid.major.x = element_blank(),
                axis.text.y=element_text(hjust=0),
                axis.text.x=element_text(angle=45,vjust=1,hjust=1)) +
          scale_x_discrete(expand=c(0,1))
      }
      if(input$Alpha2 > .075 & input$Alpha2 <= .1){
        graph <- ggplot(data = GraphFrame, aes(Correlation, GraphFrame[,1],colour=RV, size = p), environment = environment()) + labs(y = y, x = x)  + 
          scale_size_continuous("P-value",range=c(12,4),limits = c(min(p), .1), breaks = c(0.0001,0.025,0.05,0.075,0.1),
                                labels = c(".0001", "0.025", "0.05", "0.075","0.1")) + geom_point(colour="black",aes(size = p),shape=21,alpha=I(1)) + 
          geom_point(alpha=I(.75))  + theme_minimal()  +
          theme(panel.background=element_rect(fill="white"),
                panel.grid = element_line(),panel.grid.major.y=element_blank(),
                panel.grid.major.x = element_blank(),
                axis.text.y=element_text(hjust=0),
                axis.text.x=element_text(angle=45,vjust=1,hjust=1)) +
          scale_x_discrete(expand=c(0,1))
      }
      if(input$Alpha2 > .1){
        graph <- ggplot(data = GraphFrame, aes(Correlation, GraphFrame[,1],colour=RV, size = p), environment = environment()) + labs(y = y, x = x)  + 
          scale_size_continuous("P-value",range=c(12,4),limits = c(min(p), .1), breaks = c(0.0001,0.025,0.05,0.075,0.1),
                                labels = c(".0001", "0.025", "0.05", "0.075","> 0.1")) + geom_point(colour="black",aes(size = p),shape=21,alpha=I(1)) + 
          geom_point(alpha=I(.75))  + theme_minimal()  +
          theme(panel.background=element_rect(fill="white"),
                panel.grid = element_line(),panel.grid.major.y=element_blank(),
                panel.grid.major.x = element_blank(),
                axis.text.y=element_text(hjust=0),
                axis.text.x=element_text(angle=45,vjust=1,hjust=1)) +
          scale_x_discrete(expand=c(0,1))
      }
      if(input$Alpha2 <= .05){
        graph <- ggplot(data = GraphFrame, aes(Correlation, GraphFrame[,1],colour=RV, size = p), environment = environment()) + 
          labs(y = y, x = x)  + scale_size_continuous("P-value",range=c(12,4),limits = c(min(p), .1), breaks = c(0.0001,0.025,0.05),
                                                      labels = c(".0001","0.025", "0.05")) + geom_point(colour="black",aes(size = p),shape=21,alpha=I(1)) + 
          geom_point(alpha=I(.75))  + theme_minimal()  +
          theme(panel.background=element_rect(fill="white"),
                panel.grid = element_line(),panel.grid.major.y=element_blank(),
                panel.grid.major.x = element_blank(),
                axis.text.y=element_text(hjust=0),
                axis.text.x=element_text(angle=45,vjust=1,hjust=1)) +
          scale_x_discrete(expand=c(0,1))
      }
    }
    if(length(p)/length(unique(GraphFrame$Correlation)) > 260){
      graph <- ggplot(data = GraphFrame, aes(Correlation, GraphFrame[,1],colour=RV, size = p), environment = environment()) + labs(y = y, x = x)  + scale_size_continuous("P-value",range=c(12,4),limits = c(min(p), .06), breaks = c(round(min(p),4),round((min(p)+.05)/2,4),.05)) + geom_point(colour="black",aes(size = p),shape=21,alpha=I(1)) + geom_point(alpha=I(.75))  + theme_minimal()  +
        theme(panel.background=element_rect(fill="white"),
              panel.grid = element_line(),panel.grid.major.y=element_blank(),
              panel.grid.major.x = element_blank(),
              axis.text.y=element_blank(),
              axis.text.x=element_text(angle=45,vjust=1,hjust=1)) +
        scale_x_discrete(expand=c(0,1)) + scale_y_discrete(breaks = NULL)
    }
    if(input$bubbleColorRange){
      graph <- graph + scale_color_gradient2(paste(values$corr.method,"'s r",sep = ""),guide="colorbar",low="navy", high="red", midpoint=0, limits = c(-1,1))
    } else{
      graph <- graph + scale_color_gradient2(paste(values$corr.method,"'s r",sep = ""),guide="colorbar",low="navy", high="red", midpoint=0) 
    }
    return(graph)
  })
  
  output$correlations_plot <- renderPlot({
    correlations_makePlot()
  }, height = function(){input$PlotHeight3}, width = function(){input$PlotWidth3})
  
  output$download_corr_plot <- downloadHandler(
    filename = function() {paste(values$project.name, "_","significant_correlations_plot.png")},
    content = function(file){
      png(file, width = (plotresolution.corr1()/72)*width(), height = (plotresolution.corr1()/72)*height(), res = plotresolution.corr1())
      print(correlations_makePlot())
      dev.off()
    }
  )
  
  output$Visit3 <- renderUI({
    correlations <- values$corrs[[which(values$corr.names == input$TypeVariable7)]]
    visit.name = colnames(correlations)[3]
    visit <- unique(correlations[,3])
    visit <- gtools::mixedsort(as.character(visit))
    selectizeInput("visit3", paste("Correlations by"," ",visit.name, ":", sep = ""), choices = visit, selected = visit[1])
  })
  
  
  scatter_plot <- reactive({
    correlations <- corrResultsOverview()$correlations[[which(values$corr.names == input$TypeVariable7)]]
    correlation_file <- values$corr.files[[which(values$corr.names == input$TypeVariable7)]]
    x <- correlation_file[,input$WithVariable5][which(correlation_file[,colnames(correlations)[3]] == input$visit3)]
    y <- correlation_file[,input$corr_Variable2][which(correlation_file[,colnames(correlations)[3]] == input$visit3)]
    z <- list(x = x, y = y)
    return(z)
  })
  
  plot_res5 <- reactive({
    input$PlotRes5
  })
  
  
  axis_text_size <- reactive({
    input$axis_text_size
  })
  
  axis_label_size <- reactive({
    input$axis_label_size
  })
  
  scatter_makePlot <- reactive({
    correlations <- corrResultsOverview()$correlations[[which(values$corr.names == input$TypeVariable7)]]
    if(input$log_scale == TRUE){
      dat <- data.frame(x = log2(scatter_plot()$x + 1), y = log2(scatter_plot()$y + 1))
    }
    else{
      dat <- data.frame(x = scatter_plot()$x, y = scatter_plot()$y)
    }
    pearson.r <- cor.test(dat$x, dat$y)$estimate
    p.val <- cor.test(dat$x, dat$y)$p.value
    scat.plot <-  ggplot(data = dat, aes(x = x, y = y)) + geom_point(size = input$Point_size) + theme_bw()+
      labs(x = paste(input$WithVariable5, colnames(correlations)[3], input$visit3, sep = " "), y = paste(input$corr_Variable2, colnames(correlations)[3], input$visit3, sep = " ")) +
      theme(axis.text = element_text(size = axis_text_size()), axis.title = element_text(size = axis_label_size()))
    
    if(input$plot_reg == TRUE){
      scat.plot <- scat.plot + geom_smooth(method = "lm", se = FALSE)
    }
    if(input$plot_loess == TRUE){
      if(input$plot_reg == TRUE){
        scat.plot <- scat.plot + geom_smooth(method = "lm", se = FALSE)
      }
      else{
        scat.plot <- scat.plot + geom_smooth(method = "loess", se = FALSE, span = input$span)
      }
    }
    grob <- grid::grobTree(textGrob(paste("Pearson r = ", round(pearson.r, 4),"\n P-Value = ", round(p.val, 4), sep = ""), x = .99, y=.965, just=1,
                                    gp=gpar(col="blue", fontsize=14)))
    scat.plot <- scat.plot + annotation_custom(grob)
    z <- list(scatplot = scat.plot, dat = dat)
    return(z)
  })
  
  
  output$correlations_scatter_plot <- renderPlot({
    scatter_makePlot()$scatplot
  }, width = function(){input$PlotWidth5}, height = function(){input$PlotHeight5})
  
  
  output$download_scatter_data <- downloadHandler(
    filename = function() {correlations <- values$corrs[[which(values$corr.names == input$TypeVariable7)]]
    paste(values$project.name, "_",input$corr_Variable2, "vs", input$WithVariable5, "_", colnames(correlations)[3], input$visit3,"_data.csv")},
    content = function(file) {
      dat <- scatter_makePlot()$dat
      colnames(dat) <- c(input$WithVariable5, input$corr_Variable2)
      write.csv(dat, file, row.names = FALSE)
    }
  )
  
  output$download_scatter <- downloadHandler(
    filename = function() {correlations <- values$corrs[[which(values$corr.names == input$TypeVariable7)]]
    paste(values$project.name, "_",input$corr_Variable2, "vs", input$WithVariable5, "_", colnames(correlations)[3], input$visit3,".png")},
    content = function(file){
      png(file, width = (plot_res5()/72)*width(), height = (plot_res5()/72)*height(), res = plot_res5())
      print(scatter_makePlot()$scatplot)
      dev.off()
    }
  )
  
  ################End of Correlations Part###################

})
