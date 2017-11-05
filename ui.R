#User Interface for PALO FILTER Script!
#options(shiny.deprecation.messages = FALSE)
require(shiny)
#library(rmarkdown)
#library(RColorBrewer)
#library(fastcluster)
#library(NMF)
#library(grid)
#library(clValid)
#library(qusage)
#library(VennDiagram)
#library(gtools)
#library(scales)
#library(reshape2)
#library(data.table)
#library(pca3d)
#library(shinyjs)
#library(stringr)
library(shinydashboard)
library(ggplot2)

function(request){

header <- dashboardHeader(title = "Biostatistical Analysis Reporting Tool (BART)",
                          dropdownMenuOutput("messageMenu"),
                          titleWidth = 475)

sidebar <- dashboardSidebar(
  sidebarMenu(
    #menuItem("About BART", icon = icon("info-circle"), tabName = "intro"),
    menuItem("Upload", icon = icon("cloud-upload"), tabName = "upload"),
    menuItemOutput("SumStat"),
    menuItemOutput("Unsupervised"),
    menuItemOutput("diffge"),
    menuItemOutput("qusage"),
    menuItemOutput("roast"),
    menuItemOutput("flow.data"),
    menuItemOutput("metab.data"),
    menuItemOutput("correlations")
  )
)

body <-  dashboardBody(
  tags$style(HTML(".skin-blue .main-header .navbar {background-color: #2c343a}
              .skin-blue .left-side, .skin-blue .main-sidebar, .skin-blue .wrapper {background-color: #363f47}
              .skin-blue .main-header .logo {background-color: #2c343a}
              .skin-blue .main-header .logo:hover {background-color: #456dae}
              .skin-blue .main-header .navbar .sidebar-toggle:hover {background-color: #456dae}
              .skin-blue .sidebar-menu>li.active>a, .skin-blue .sidebar-menu>li:hover>a {border-left-color: #456dae}
              .box.box-primary {border-top-color: #18c4be}
              .nav-tabs-custom>.nav-tabs>li.active {border-top-color: #456dae}
              .box-body {width: 100%}
              #designDataTable.shiny-datatable-output.shiny-bound-output {overflow-x: auto}
              #sigcomptable.shiny-datatable-output.shiny-bound-output {overflow-x: auto}
              #genelisttable.shiny-datatable-output.shiny-bound-output {overflow-x: auto}
              #vennIntersection.shiny-datatable-output.shiny-bound-output {overflow-x: auto}
              #MultipleCompTab.shiny-datatable-output.shiny-bound-output {overflow-x: auto}
              #GeneSetTab.shiny-datatable-output.shiny-bound-output {overflow-x: auto}
              #summary1.shiny-html-output.shiny-bound-output {overflow-x: auto}
              #summary0.shiny-html-output.shiny-bound-output {overflow-x: auto}
              #modHeatmap.shiny-plot-output.shiny-bound-output {height: auto !important; overflow-x: auto}
              #heatmap.shiny-plot-output.shiny-bound-output {height: auto !important; overflow-x: auto}
              #modMap.shiny-plot-output.shiny-bound-output {height: auto !important}
              #heatmap1.shiny-plot-output.shiny-bound-output {height: auto !important; overflow-x: auto}
              #modDgeMap.shiny-plot-output.shiny-bound-output {height: auto !important; overflow-x: auto}
              #indFcPlot .col-sm-12 {padding-left: 0px; padding-right: 0px}
              #GeneSetPlot.shiny-plot-output.shiny-bound-output {height: auto !important; overflow-x: auto}
              #foldchange1.shiny-plot-output.shiny-bound-output {height: auto !important; overflow-x: auto}
              #correlations_plotOverview.shiny-plot-output.shiny-bound-output {height: auto !important; overflow-x: auto}
              #correlations_plot.shiny-plot-output.shiny-bound-output {height: auto !important; overflow-x: auto}
              #correlations_scatter_plot.shiny-plot-output.shiny-bound-output {height: auto !important; overflow-x: auto}")),
  shinyjs::useShinyjs(),
  tabItems(
    #tabItem(tabName = "intro",
            #htmltools::includeMarkdown("bart-vignette.md")),
    tabItem(tabName = "upload", 
            fluidRow(
              box(title = "Upload Results", width = 4, status = "primary", solidHeader = FALSE,
                  verbatimTextOutput('version'),
                  hr(),
                  fileInput('file1', 'Upload BART Files and any QC reports for the Project:', multiple = TRUE,
                            accept=c(".Rdata",".png",".html",".csv")),
                  bookmarkButton(),
                  tags$script('
                              Shiny.addCustomMessageHandler("resetFileInputHandler", function(x) {      
                              var id = "#" + x + "_progress";
                              var idBar = id + " .bar";
                              $(id).css("visibility", "hidden");
                              $(idBar).css("width", "0%");
                              });
                              ')
              ),
              box(title = "Summary", width = 8, status = "primary", solidHeader = FALSE,
                  helpText('Name of Uploaded Project: '),
                  textOutput('path'),
                  helpText(" "),
                  helpText("The tabs on the left contain reports of various aspects of the project. Be aware that
                           the app utilizes 'lazy' loads, meaning that it only uses the necessary pieces of information to produce the image at hand.
                           So when you arrive at each tab for the first time there will be some small red error messages that are residuals of the lazy loading.
                           Give the app a few seconds (maybe a minute for the heatmaps) and the errors will go away and the reports generated.  Thank you.")
              )
            )
    ),
    tabItem(tabName = "design",
            fluidRow(
              box(title = "Sample Annotation File", width = 12, status = "primary", solidHeader = FALSE,
                  tags$head(
                    tags$style(type="text/css", "tfoot {display: table-header-group}") # Move filters to top of datatable
                  ),
                  downloadButton('downloadDesign', "Download Table"),
                  dataTableOutput('designDataTable')
              )
            )
    ),
    
    tabItem(tabName = "summary",
            fluidRow(
              box(title = "Summary Options", width = 4, status = "primary", solidHeader = FALSE,
                  uiOutput('summaryName'),
                  uiOutput('respVar'),   
                  selectInput("digits", 
                              label = "Number of decimal places:",
                              choices = c(0,1,2,3,4),
                              selected = 1)
              ),
              box(title = "Table 1", width = 8, status = "primary", solidHeader = FALSE,
                  downloadButton('downloadSummary0', 'Download Table'),
                  br(),
                  div(style = "display:inline-block", uiOutput("summary0Text")),
                  div(style = "display:inline-block", helpPopup("Summary Table 1", "Given the summary variable selected by the user, statistics will be provided for 
                                                                each category of the users BY variable selection.  In longitudinal settings, this table will display 
                                                                a warning if the summary variable is changing over time.  Refer to Table 2.", placement = "right",
                                                                trigger = "click")),
                  tableOutput('summary0')
              )
                  ),
            fluidRow(
              box(title = "Table 2", width = 6, status = "primary", solidHeader = FALSE,
                  downloadButton('downloadSummary1', 'Download Table'),
                  br(),
                  div(style = "display:inline-block", uiOutput("summary1Text")),
                  div(style = "display:inline-block", helpPopup("Summary Table 2", "Given the summary variable selected by the user, statistics will be provided for 
                                                                each category of the users BY variable selection and for each time point in longitudinal settings.", placement = "right",
                                                                trigger = "click")),
                  tableOutput('summary1')
              ),
              box(title = "Table 3", width = 6, status = "primary", solidHeader = FALSE,
                  downloadButton('downloadSummary2', 'Download Table'),
                  br(),
                  div(style = "display:inline-block", uiOutput("summary2Text")),
                  div(style = "display:inline-block", helpPopup("Summary Table 3", "For longitudinal data, this table provides the counts of subjects categorized by 
                                                                their observed longitudinal profiles and by the specified BY variable.", placement = "right",
                                                                trigger = "click")),
                  tableOutput('summary2')
              )
                  )
              ),
    
    tabItem(tabName = "pvca",
            helpText("Data Filtering and Batch Assesment Summary:"),
            textOutput("PVCAtext"), #reactive text with the new content
            imageOutput('PVCA1'),  
            imageOutput('PVCA2')
    ),
    
    tabItem(tabName = "fastqc",
            fluidRow(
              box(title = "FastQC Samples", width = 3, status = "primary",
                  uiOutput('qc_select'),
                  uiOutput('qc_select2')
              ),
              box(title = "FastQC Report", width = 9, status = "primary",
                  uiOutput("fqc"),
                  plotOutput("ex"),
                  plotOutput("ex2"),
                  plotOutput("ex3")
              )
            )
    ),
    
    tabItem(tabName = "qcupload",
            sidebarLayout(
              sidebarPanel(
                uiOutput("dropdown"),
                uiOutput("dropdowncolumns"),
                uiOutput("coordflip"),
                br(),
                div(style = "display:inline-block", checkboxInput("graphicsqc", strong("Graphing options", style = "color:#456dae"), TRUE)),
                div(style = "display:inline-block", helpPopup("Graphing options", "These options allow the user to make adjustments to the plot (i.e. width and height).")),
                conditionalPanel(condition = "input.graphicsqc",
                                 sliderInput('qcplotsizew', "Grid width", min = 500, max = 2000, value = 800, step = 25),
                                 sliderInput('qcplotsizeh', "Grid height", min = 500, max = 2000, value = 600, step = 25),
                                 numericInput("axis_text_sizeqc", "Axis text size:", min = 10, value = 10, step = 1),
                                 numericInput("axis_label_sizeqc", "Axis label size:", min = 12, value = 12, step = 1),
                                 numericInput('PlotResolutionqc', "Plot resolution (for downloaded plot):", min = 72, max = 400, value = 72, step = 1))),
              mainPanel(
                downloadButton('downloadQC', 'Download Figure'),
                uiOutput("QCbox")
              )
            )
    ),
   
    tabItem(tabName = "probeheatmap",
            fluidRow(
              tabBox(title = "Heatmap Tool", width = 12, id = "unsupervisedHeatTool",
                     tabPanel("Gene Level Heatmaps",
                              fluidRow(
                                box(title = "Heatmap Options", width = 4, status = "primary", solidHeader = FALSE,
                                    div(style = "display:inline-block", actionButton("go", "Plot")),
                                    div(style = "display:inline-block", infoPopup("Plot", 'The heatmap will not update until the "Plot" button is clicked.
                                                                                  This allows the user to make multiple adjustments at once, without having to wait for 
                                                                                  each individual adjustment to update on the heatmap.', placement = "right", trigger = "click")),
                                    uiOutput('test'),
                                    br(),
                                    uploadVarsUI("transcripts", varType = "probes"),
                                    maxValuesUI("unsupervisedMaxValues"),
                                    br(),
                                    subsetAndOrderUI("unsupervisedSubOrder"),
                                    br(),
                                    colClusterUI("unsupervisedCluster"),
                                    br(),
                                    graphOptionsUI("unsupervisedGraphOptions")
                                ),
                                box(title = "Heatmap", width = 8, status = "primary", solidHeader = FALSE, 
                                    div(style = "display:inline-block", strong("Help Message")),
                                    div(style = "display:inline-block", infoPopup("Gene Level Heat Maps", 'The heat map below is constructed on individual samples for a number of scenarios, 
                                                                                  specifically for baseline samples only (cross sectional) or all samples at all time points. For baseline samples, 
                                                                                  heatmaps can be generated based on normalizing the expression data to the mean, or to a control group if applicable.  
                                                                                  When all samples at all time points are to be plotted, heatmaps can be generated by normalizing to the mean, 
                                                                                  a control group, or each subjects own baseline value. Only samples that have a corresponding baseline sample are 
                                                                                  included in the map. The initial graph of the heat map may not be very appealing depending on the number of samples. 
                                                                                  The inputs on the left have a wide variety of user options that range from the type of normalized data, clustering rows 
                                                                                  and columns, subsetting samples and probes, etc. These options are consistent across all the unsupervised analysis plots.  
                                                                                  One addition, unique to the probe level heat map is the "Max value on color key" option that specifies the coloring legend 
                                                                                  for the heatmap index. The default is set to +/- 2. If the user would like the index to have the most extreme red and blue 
                                                                                  to be set to a value of +/- 4 then user simply needs to enter 4. The user can also enter 0 which will make the limits based 
                                                                                  on the max and min of the entire expression file.',
                                                                                  placement = "bottom", trigger = "click")),
                                    helpText(""),
                                    textOutput('OptimalNumber'),
                                    downloadButton('downloadHeatmapData', 'Download Data'),
                                    downloadButton("downloadHeatmap", "Download Figure"),
                                    helpText(""),
                                    #div(style = "display:inline-block", checkboxInput("modulemeans", strong("Download module scores:", style = "color:#456dae"), FALSE)),
                                    #div(style = "display:inline-block", infoPopup("Download module scores", "The downloaded data is calculated by averaging all the probes within each module.",
                                                                                  #placement = "right", trigger = "click")),
                                    plotOutput('heatmap')
                                )
                              )
                     ),
                     tabPanel("Cluster Number Diagnostics",
                              div(class = "container-fluid",
                                  div(class = "row",
                                      div(class = "col-md-4",
                                        h4(strong("Dunn Index")),
                                        p("The Dunn Index is a metric used to assess the quality of a given set of clusters. When the number of clusters (n) is not known beforehand (such is the case with hierarchical clustering),
                                           the Dunn Index is computed for n = 2 to n = half the sample size. The optimal number of clusters is chosen by the highest Dunn Index. While in some sense, metrics such as the Dunn Index
                                           provide an objective means by which to select the optimal number of clusters, choosing cluster groupings should always entail looking at the heat maps and dendrograms to see if the number
                                           suggested by the Dunn Index makes sense for the specific data at hand."),
                                        p(strong("The bar plot on the right is based on the samples clustered in the previous tab."))
                                      ),
                                      div(class = "col-md-8",
                                          downloadButton('downloadClusterPlot3', 'Download Figure'),
                                          plotOutput('clusterPlot')
                                      )
                                  )
                              )
                     ),
                     tabPanel("Cluster Association Analysis",
                              fluidRow(
                                box(title = "Options", width = 4, status = "primary", solidHeader = FALSE,
                                    clusterAssociationUI("unsupervisedAssociation")
                                ),
                                box(title = "Association Tests", width = 8, status = "primary", solidHeader = FALSE,
                                    div(style = "overflow-x: auto",
                                        tableOutput('cluster_output4'),
                                        helpText("The Chi-square test below is used to test for association between the cluster groups and user specified variable. 
                                                 The Chi-square test is ideal when expected cell counts are large (expected values greater than 5). 
                                                 Because the test only makes sense for 2x2 contingency tables and above, test statistics are not shown 
                                                 when the data is subsetted on only one value."),
                                        helpText(""),
                                        helpText("Fisher's exact test is also shown. This test is ideal when cell counts are small. Like the Chi-square test, 
                                                 when the data is subsetted on only one value, test statistics are not shown."),
                                        helpText(""),
                                        textOutput('explanation3'),
                                        tableOutput('cluster_tab3'),
                                        textOutput('chisquare_test4'),
                                        textOutput('fisher4')
                                    )
                                )
                              )
                     )
              )
            )
    ),
                              
    tabItem(tabName = "moduleMap",
            fluidRow(
              tabBox(title = "Longitudinal Modular Maps", width = 12, id = "modMap",
                     tabPanel("Modular Map",
                              fluidRow(
                                box(title = "Heatmap Options", width = 4, status = "primary", solidHeader = FALSE,
                                    div(style = "display:inline-block", actionButton("goMod", "Plot")),
                                    uiOutput("baseOrCtrl"),
                                    uiOutput("modSelection"),
                                    uploadVarsUI("mod", varType = "modules"),
                                    br(),
                                    subsetAndOrderUI("modSubOrder"),
                                    br(),
                                    colClusterUI("modCluster"),
                                    br(),
                                    graphOptionsUI("modGraphOptions", varType = "modules")
                                ),
                                box(title = "Module Map", width = 8, status = "primary", solidHeader = FALSE,
                                    div(style = "display:inline-block", strong("Help Message")),
                                    div(style = "display:inline-block", infoPopup("Individual Longitudinal Module Map", "The module map below is constructed on individual samples for all time points. Healthy Control samples are used to determine an upper and lower threshold (mean healthy controls +/- 2 sd). The module proportion for each sample is then calculated based on the percentage of probes within a module that are above or below this threshold.  If there are no healthy controls in the study, the baseline samples will be treated as a healthy control and will serve as the threshold.
                                                                                  The initial graph of the modules may not be very appealing depending on the number of samples and modules. The inputs on the left are designed to let the user modify the map to their liking. The first 3 inputs allow the user to order the samples how they would like.
                                                                                  The default setting is to order the samples by responder status and then by donor id. The ordering variables are labeled across the top and included in the legend.
                                                                                  The next box allows for additional labeling of the sample across the top but do not contribute to the ordering of the samples.
                                                                                  The remaining options are either graphical parameters or subsetting paramaters. The user can specify to look at all the Baylor modules at once, only the first 6 rounds, or only annotated ones.
                                                                                  If there are too many samples the user can subset the data (Subset the module map checkbox) by selecting up to two variables two subset on. This allows the user to only look at certain subgroups like just the male, high responders.
                                                                                  Use the graphing inputs to make the heatmap pretty. It is suggested to start with the size and aspect ratio and only edit the circle sizes if needed after.",
                                                                                  placement = "bottom", trigger = "click")),
                                    helpText(""),
                                    textOutput('modOptimalNumber'),
                                    downloadButton('downloadModMap3Other', 'Download Data'),
                                    downloadButton("downloadModPlot3Other", "Download Figure"),
                                    plotOutput('modHeatmap')
                                )
                              )
                     ),
                     tabPanel("Cluster Number Diagnostics",
                              downloadButton('downloadClusterPlot2Other', 'Download Figure'),
                              plotOutput('modClusterPlot')
                     ),
                     tabPanel("Cluster Association Analysis",
                              fluidRow(
                                box(title = "Options", width = 4, status = "primary", solidHeader = FALSE,
                                    clusterAssociationUI("modAssociation")
                                ),
                                box(title = "Association Tests", width = 8, status = "primary", solidHeader = FALSE,
                                    div(style = "overflow-x: auto",
                                        tableOutput('cluster_output3Other'),
                                        helpText("The Chi-square test below is used to test for association between the cluster groups and user specified variable. 
                                                 The Chi-square test is ideal when expected cell counts are large (expected values greater than 5). 
                                                 Because the test only makes sense for 2x2 contingency tables and above, test statistics are not shown 
                                                 when the data is subsetted on only one value."),
                                        helpText(""),
                                        helpText("Fisher's exact test is also shown. This test is ideal when cell counts are small. Like the Chi-square test, 
                                                 when the data is subsetted on only one value, test statistics are not shown."),
                                        helpText(""),
                                        textOutput('explanation2Other'),
                                        tableOutput('cluster_tab2Other'),
                                        textOutput('chisquare_test3Other'),
                                        textOutput('fisher3Other')
                                    )
                                )
                              )
                     )
              )
            )
    ),
    
    tabItem(tabName = "overview",
            fluidRow(
              box(title = "Filtering Options", width = 4, status = "primary", solidHeader = FALSE,
                  numericInput('alphalevel2',"Significance threshold:", min = 0, max = 1, value = 0.05, step = 0.025),  
                  checkboxInput("overviewFc",label = strong("Filter on fold change:"),value = FALSE),
                  conditionalPanel(condition = "input.overviewFc",
                                   numericInput("selectFcValue", label = "Fold change threshold:",min = 0, value = 0, step = 0.1),
                                   selectInput("selectFcSign","Fold change sign:", c("+", "-", "Both"), c("Both"))),
                  #selectInput("sigsign", "Fold change sign:", c("All (include 0)", "+", "-"), selected = "all (include 0)"),
                  tags$hr() # add a horizontal line
              ),
              box(title = "Results Overview", width = 8, status = "primary", solidHeader = FALSE,
                  helpText("Cell value represents the number of significant probes under the selected significance level."),
                  downloadButton('downloadSC', 'Download Table'),
                  dataTableOutput("sigcomptable")
              )
            )
    ),
    
    tabItem(tabName = "genelistmaker",
            fluidRow(
              uiOutput("dgeTabs")
            )
    ),
    
    tabItem(tabName = "genesearch",
            fluidRow(
              box(title = "Gene Selection", width = 4, status = "primary", solidHeader = FALSE,
                  uiOutput("specgene1"),
                  uiOutput("probeid1"),
                  uiOutput("genTable")
              ),
              box(title = "Results for Specific Genes", width = 8, status = "primary", solidHeader = FALSE,
                  helpText("The processing speed depends on the number of gene symbols and your computer's performance. It may take up to 30 seconds to proceed."),
                  downloadButton('downloadCG','Download Table'),
                  dataTableOutput("specgenetable")
              )
            )
    ),
    
    tabItem(tabName = "qusage",
            fluidRow(
              tabBox(title = "Qusage Results", width = 12, id = "qusResults",
                     tabPanel("Overview",
                              fluidRow(
                                box(title = "Filtering Options", width = 4, status = "primary", solidHeader = FALSE,
                                    numericInput("SigLevel", "Significance threshold:", min = 0, max = 1, step = .025, value = .05),
                                    selectInput("foldchange.q", strong("Fold change sign:"), choices = c("All (include 0)", "+","-"), selected = "All (include 0)")
                                ),
                                box(title = "Results Overview", width = 8, status = "primary", solidHeader = FALSE,
                                    downloadButton("downloadSigComps", "Download Table"),
                                    dataTableOutput('CompOverview')
                                )
                              )
                     ),
                     tabPanel("Fold Change Plot and Venn Diagram",
                              fluidRow(
                                box(title = "Filtering/Plotting Options", width = 4, status = "primary", solidHeader = FALSE,
                                    div(style = "display:inline-block", uiOutput("Comparisons1")),
                                    div(style = "display:inline-block", infoPopup("Comparisons", "The user may select up to five comparisons.", placement = "right",
                                    trigger = "click")),
                                    uiOutput("PaloOrFirst1"),
                                    checkboxInput("only_annotated", strong("Plot only annotated modules"), FALSE),
                                    checkboxInput("showAnnotated", strong("Show annotations on plot"), FALSE),
                                    filterOptsUI("qusFcPlot"),
                                    br(),
                                    div(style = "display:inline-block", checkboxInput("graphics6", strong("Graphing options", style = "color:#456dae"), FALSE)),
                                    div(style = "display:inline-block", helpPopup("Graphing options", "These options allow the user to make adjustments to the plot (i.e. width and height).")),
                                    conditionalPanel(condition = "input.graphics6",
                                                     div(style = "display:inline-block",uiOutput("ColorChoice")), 
                                                     div(style = "display:inline-block",helpPopup("Help Message", "The order of the colors corresponds to the order in the legend, not the order selected.",
                                                                                                  placement = "right", trigger = "click")),
                                                     sliderInput("PlotWidth1", "Plot Width", min = 500, max = 5000, value = 800, step = 100),
                                                     numericInput("PlotRes1", "Plot resolution (for downloaded plot):", min = 72, max = 400, value = 72, step = 1))
                                ),
                                div(id = "fcPlot", style = "overflow-x: auto", 
                                    box(title = "FC Plot", width = 12, status = "primary", solidHeader = FALSE,
                                        downloadButton("downloadPlot2", "Download Figure"),
                                        plotOutput("foldchange1")
                                    )
                                )
                              ),
                              fluidRow(
                                box(title = "Venn Diagram", width = 5, status = "primary", solidHeader = FALSE,
                                    uiOutput("venn.download"),
                                    uiOutput('venn')
                                ),
                                box(title = "Results Table", width = 7, status = "primary", solidHeader = FALSE,
                                    downloadButton("downloadTable2", "Download Table"),
                                    dataTableOutput('MultipleCompTab')
                                )
                              )
                     ),
                     tabPanel("Individual Gene Set Plot and Table",
                              fluidRow(
                                box(title = "Filtering Options", width = 4, status = "primary", solidHeader = FALSE,
                                    uiOutput("Module_Select"),
                                    uiOutput("Comparisons2"),
                                    br(),
                                    filterOptsUI("qusIndFcPlot"),
                                    br(),
                                    div(style = "display:inline-block", checkboxInput("graphics7", strong("Graphing options", style = "color:#456dae"), FALSE)),
                                    div(style = "display:inline-block", helpPopup("Graphing options", "These options allow the user to make adjustments to the plot (i.e. width and height).")),
                                    conditionalPanel(condition = "input.graphics7",
                                                     sliderInput("PlotWidth2", "Plot Width", min = 500, max = 5000, value = 800, step = 100),
                                                     numericInput("PlotRes2", "Plot resolution (for downloaded plot):", min = 72, max = 400, value = 72, step = 1))
                                ),
                                column(width = 8, class = "container-fluid",
                                       div(id = "indFcPlot", style = "overflow-x: auto",
                                           box(title = "FC Plot", width = 12, status = "primary", solidHeader = FALSE,
                                               downloadButton("downloadPlot3", "Download Figure"),
                                               plotOutput('GeneSetPlot')
                                           )
                                       ),
                                       box(title = "Results Table", width = NULL, status = "primary", solidHeader = FALSE,
                                           downloadButton("downloadTable3", "Download Table"),
                                           dataTableOutput('GeneSetTab')
                                       )
                                )
                              )
                     )
              )
            )
    ),
    
      tabItem(tabName = "roast",
              fluidRow(
                tabBox(title = "Roast Results", width = 12, id = "roastResults",
                       tabPanel("Overview",
                                fluidRow(
                                  box(title = "Filtering Options", width = 4, status = "primary", solidHeader = FALSE,
                                      uiOutput("setStat"),
                                      numericInput('SigLevelR', "Significance threshold:", min = 0, max = 1, value = .05, step = .025)
                                  ),
                                  box(title = "Results Overview", 
                                      downloadButton("downloadSigCompsR", "Download Table"),
                                      dataTableOutput('CompOverviewR')
                                  )
                                )
                       )
                )
              )
      ),
		
      tabItem(tabName = "flow",
              fluidRow(
                tabBox(title = "Flow Results", width = 12, id = "flowResults",
                       tabPanel("Overview",
                                fluidRow(
                                  box(title = "Filtering Options", width = 4, status = "primary", solidHeader = FALSE,
                                      numericInput('alphaFlow_1',"Significance threshold:", min = 0, max = 1, value = 0.05, step = 0.025)
                                  ),
                                  box(title = "Results Overview", width = 8, status = "primary", solidHeader = FALSE,
                                      downloadButton('downloadFC', 'Download Table'),
                                      dataTableOutput("FlowOverview")
                                  )
                                )
                       ),
                       tabPanel("Significant Variables List",
                                fluidRow(
                                  box(title = "Filtering Options", width = 4, status = "primary", solidHeader = FALSE,
                                      uiOutput("CompF"),
                                      filterOptsUI('flowResultsTable')
                                  ),
                                  box(title = "Flow List", width = 8, status = "primary", solidHeader = FALSE,
                                      helpText("Right click on hyperlinks to open in new window"),
                                      downloadButton('downloadFL', 'Download Table'),
                                      dataTableOutput("flowlisttable")
                                  )
                                )
                       ),
                       tabPanel("Flow Figures",
                                fluidRow(
                                  box(title = "Options", width = 4, status = "primary", solidHeader = FALSE,
                                      uiOutput("FlowDataNames"),
                                      uiOutput("flowsub"),               
                                      checkboxInput("FlowTransform","lg2 Transform",FALSE),
                                      helpText("Options for 1st Plot:"),
                                      uiOutput("flowmax"),
                                      uiOutput("flowmin"),
                                      checkboxInput("FlowSamples","Individual curves",FALSE),
                                      helpText("Options for 2nd Plot:"),
                                      checkboxInput("flowbox","Box Plot View",FALSE),
                                      downloadButton('downloadFlowData', 'Download Data')
                                  ),
                                  box(title = "Figures", width = 8, status = "primary", solidHeader = FALSE,
                                      downloadButton('downloadFlowPlot', 'Download Figure'),
                                      plotOutput("FlowPlot"),
                                      downloadButton('downloadFlowPlot2', 'Download Figure'),
                                      plotOutput("FlowPlot2"),
                                      downloadButton('downloadFlowSummaries','Download Table'),
                                      dataTableOutput("FlowPlotSummary")
                                  )
                                )
                       )
                )
              )
      ),
    
    tabItem(tabName = "pcaM",
            fluidRow(
              box(title = "PCA Options", width = 4, status = "primary", solidHeader = FALSE,
                  uiOutput("pcaAnnotM"),
                  uiOutput("pcaAnnotValsM"),
                  uiOutput("pcaBlockingM"),
                  selectInput("PCSnumM", "Select number of PCs to plot:", choices = c(2,3), selected = 2),
                  div(style = "display:inline-block", checkboxInput("graphicsPCAM", strong("Graphing options", style = "color:#456dae"), FALSE)),
                  div(style = "display:inline-block", helpPopup("Graphing options", "These options allow the user to make adjustments to the plot (i.e. width and height).")),
                  conditionalPanel(condition = "input.graphicsPCAM",
                                   sliderInput("PlotWidthPCAM", "Plot Width", min = 500, max = 2000, value = 800, step = 50),
                                   sliderInput("PlotHeightPCAM", "Plot Height", min = 500, max = 2000, value = 500, step = 50),
                                   numericInput("CircleSizePCAM","Circle size",min=3,max=15,value=7,step=1),
                                   numericInput("PlotResPCAM", "Plot resolution (for downloaded plot):", min = 72, max = 400, value = 72, step = 1))
              ),
              box(title = "PCA Plot", width = 8, status = "primary", solidHeader = FALSE,
                  actionButton("showM", "View/hide screeplot", style = "color: white; background-color: green"),
                  uiOutput("screePlotM"),
                  br(),
                  plotOutput("PCAplotM"),
                  br(),
                  uiOutput("PCAplot3dM")
              )
            )
    ),
    
    tabItem(tabName = "metab",
            fluidRow(
              tabBox(title = "Metabolomic Results",  width = 12, id = "metabResults",
                     tabPanel("Overview",
                              fluidRow(
                                box(title = "Filtering Options", width = 4, status = "primary", solidHeader = FALSE,
                                    numericInput('alphaMetab_1',"Significance threshold:", min = 0, max = 1, value = 0.05, step = 0.025)
                                ),
                                box(title = "Results Overview", width = 8, status = "primary", solidHeader = FALSE,
                                    downloadButton('downloadFC2', 'Download Table'),
                                    dataTableOutput("MetabOverview")
                                )
                              )
                     ),
                     tabPanel("Significant Variables List",
                              fluidRow(
                                box(title = "Filtering Options", width = 4, status = "primary", solidHeader = FALSE,
                                    uiOutput("CompM"),
                                    selectInput("methodM", "Multiple testing correction:",
                                                c("FDR", "Bonferroni", "Raw"), "Raw"),
                                    numericInput('alphaMetab_2',"Significance threshold:", min = 0, max = 1, value = 0.05, step = 0.025)
                                ),
                                box(title = "Metabolite List", width = 8, status = "primary", solidHeader = FALSE,
                                    helpText("Right click on hyperlinks to open in new window"),
                                    downloadButton('downloadME', 'Download Table'),
                                    dataTableOutput("metablisttable")
                                )
                              )
                     ),
                     tabPanel("Metabolite Figures",
                              fluidRow(
                                box(title = "Options", width = 4, status = "primary", solidHeader = FALSE,
                                    uiOutput("MetabDataNames"),
                                    uiOutput("metabsub"),               
                                    checkboxInput("MetabTransform","lg2 Transform",FALSE),
                                    helpText("Options for 1st Plot:"),
                                    uiOutput("metabmax"),
                                    uiOutput("metabmin"),
                                    checkboxInput("MetabSamples","Individual curves",FALSE),
                                    helpText("Options for 2nd Plot:"),
                                    checkboxInput("metabbox","Box Plot View",FALSE),
                                    downloadButton('downloadMetabData', 'Download Data')
                                ),
                                box(title = "Figures", width = 8, status = "primary", solidHeader = FALSE,
                                    downloadButton('downloadMetabPlot', 'Download Figure'),
                                    plotOutput("MetabPlot"),
                                    downloadButton('downloadMetabPlot2', 'Download Figure'),
                                    plotOutput("MetabPlot2"),
                                    downloadButton('downloadMetabSummaries','Download Table'),
                                    dataTableOutput("MetabPlotSummary")
                                )
                              )
                     ),
                     tabPanel("Significant Variables Heatmap",
                              fluidRow(
                                box(title = "Options", width = 4, status = "primary", solidHeader = FALSE,
                                    div(style = "display:inline-block", actionButton("go4", "Plot")),
                                    div(style = "display:inline-block", infoPopup("Plot", 'The heatmap will not update until the "Plot" button is clicked.
                                                                                  This allows the user to make multiple adjustments at once, without having to wait for 
                                                                                  each individual adjustment to update on the heatmap.', placement = "right", trigger = "click")),
                                    uiOutput("CompM2"),
                                    filterOptsUI("metabHeatmap"),
                                    uiOutput("normMetab"),
                                    div(style = "display:inline-block", checkboxInput("metabRowCluster", strong("Row Cluster"), value = FALSE)),
                                    div(style = "display:inline-block", infoPopup("Row Cluster", "Based on the PC performance, it might not be recommended to cluster the rows if the number of rows exceed 7000",
                                                                                  placement = "right", trigger = "click")),
                                    br(),
                                    maxValuesUI("metabMaxValues"),
                                    br(),
                                    subsetAndOrderUI("metabSubOrder"),
                                    br(),
                                    colClusterUI("metabCluster"),
                                    br(),
                                    graphOptionsUI("metabGraphOptions")
                                ),
                                box(title = "Heatmap", width = 8, status = "primary", solidHeader = FALSE,
                                    downloadButton("downloadHeatmapm"),
                                    plotOutput("heatmapm")
                                )
                              )
                     )
              )
            )
    ),
    
    tabItem(tabName = "corr",
            fluidRow(
              tabBox(title = "Correlation Analysis", width = 12, id = "corrResults",
                     tabPanel("Significance Overview Table",
                              fluidRow(
                                box(title = "Filtering Options", width = 4, status = "primary", solidHeader = FALSE,
                                    numericInput('alphaCorr',"Significance threshold:", min = 0, max = 1, value = 0.05, step = 0.025),
                                    checkboxInput("overviewCorrVal",label = strong("Filter on correlation value:"),value = FALSE),
                                    conditionalPanel(condition = "input.overviewCorrVal",
                                                     numericInput("selectCorrVal", label = "Correlation value cut off:",
                                                                  min = 0, max = 1, value = 0, step = 0.1),
                                                     selectInput("selectCorrSign","Correlation value sign:", c("+", "-", "Both"), c("Both")))
                                ),
                                box(title = "Results Overview", width = 8, status = "primary", solidHeader = FALSE,
                                    downloadButton('downloadOverviewCorrTable', 'Download Table'),
                                    dataTableOutput("corrOverview")
                                )
                              )
                     ),
                     tabPanel("Overview Heatmap",
                              fluidRow(
                                box(title = "Filtering/Plotting Options", width = 4, status = "primary", solidHeader = FALSE,
                                    uiOutput("TypeVariable"),
                                    uiOutput("var_switch"),
                                    uiOutput("WithVariable1"),
                                    div(style = "display:inline-block",uiOutput("uploadmod1")),
                                    div(style = "display:inline-block", helpPopup("Upload Variable List", "Allows the user to provide their own list of variables to plot. The user should create a CSV file that contains 
                                                                                  a single column of variable names to upload. The column should be given a name of the user's choosing.")),
                                    uiOutput("fileupload1"),
                                    div(style = "display:inline-block", checkboxInput("TSoptions2", strong("P-Value and Correlation Value Subsetting Options:", style = "color:#456dae"), FALSE)),
                                    conditionalPanel(condition = "input.TSoptions2",
                                                     selectInput("correction_method.corr","Multiple testing correction:",c("FDR","Bonferroni","Raw"),"Raw"),
                                                     numericInput("Alpha1", "Significance threshold:", min = 0, max = 1, step = .01, value = 1),
                                                     checkboxInput("corrval2",label = strong("Filter on correlation value:"),value = FALSE),
                                                     conditionalPanel(condition = "input.corrval2 == true",
                                                                      numericInput("corrval3", label = "Correlation value cut off:",
                                                                                   min = 0, max = 1, value = 0, step = 0.1),
                                                                      selectInput("corrsign1","Correlation value sign:", c("+", "-", "Both"), c("Both")))
                                    ),
                                    checkboxInput("subsetModcorr",strong("Subset or Order Columns:", style = "color:#456dae"),FALSE),
                                    conditionalPanel(condition="input.subsetModcorr",
                                                     uiOutput("subsetcorr")),
                                    div(style = "display:inline-block", checkboxInput("rowclustcorr", strong("Row cluster"),FALSE)),
                                    div(style = "display:inline-block",infoPopup("Row Cluster", "Based on the PC performance, it might not be recommended to cluster the rows if the number of rows exceed 7000",
                                                                                 placement = "right", trigger = "click")),
                                    checkboxInput("colclustcorr", strong("Column cluster"), FALSE),
                                    checkboxInput("colorRange", strong("Force Color Range From -1 to 1"), FALSE),
                                    div(style = "display:inline-block", checkboxInput("graphics9", strong("Graphing options", style = "color:#456dae"), TRUE)),
                                    div(style = "display:inline-block", helpPopup("Graphing options:", "These options allow the user to make adjustments to the plot (i.e. width and height).")),
                                    conditionalPanel(condition = "input.graphics9",
                                                     sliderInput("PlotWidth4", "Plot Width", min = 500, max = 5000, value = 800, step = 100),
                                                     sliderInput("PlotHeight4", "Plot Height", min = 500, max = 5000, value = 600, step = 100),
                                                     numericInput("corrFontSize", "Font size:", min = 10, value = 10, step = 2),
                                                     numericInput("corrLegendSize", "Legend size:", min = 1, max = 5, value = 1, step = .50),
                                                     numericInput("corrTreeHeight", "Tree height:", min = 0, value = 50, step = 5),
                                                     numericInput("PlotRes4", "Plot resolution (for downloaded plot):", min = 72, max = 400, value = 72, step = 1))
                                ),
                                box(title = "Heatmap", width = 8, status = "primary", solidHeader = FALSE,
                                    downloadButton('download_heatmap_data', "Download Data"),
                                    downloadButton('download_corr_heatmap', "Download Figure"),
                                    plotOutput("correlations_plotOverview")
                                )
                              )
                     ),
                     tabPanel("Correlation Bubble Plot",
                              fluidRow(
                                box(title = "Filtering/Plotting Options", width = 4, status = "primary", solidHeader = FALSE,
                                    uiOutput("TypeVariable4"),
                                    checkboxInput("var_switch3", strong("Swap Figure Rows and Columns"), value = FALSE),
                                    div(style = "display:inline-block",uiOutput("uploadmod3")),
                                    div(style = "display:inline-block", helpPopup("Upload Variable List", "Allows the user to provide their own list of variables to plot. The user should create a CSV file that contains 
                                                                                  a single column of variable names to upload. The column should be given a name of the user's choosing.")),
                                    uiOutput("fileupload2"),
                                    uiOutput("Visit2"),
                                    div(style = "display:inline-block", checkboxInput("TSoptions4", strong("P-Value and Correlation Value Subsetting Options:", style = "color:#456dae"), FALSE)),
                                    conditionalPanel(condition = "input.TSoptions4",
                                                     selectInput("correction_method.corr3","Multiple testing correction:",c("FDR","Bonferroni","Raw"),"Raw"),
                                                     numericInput("Alpha2", "Significance threshold:", min = 0, max = 1, step = .01, value = .05),
                                                     checkboxInput("corrval4",label = strong("Filter on correlation value:"),value = FALSE),
                                                     conditionalPanel(condition = "input.corrval4 == true",
                                                                      numericInput("corrval5", label = "Correlation value cut off:",
                                                                                   min = 0, max = 1, value = 0, step = 0.1),
                                                                      selectInput("corrsign2","Correlation value sign:", c("+", "-", "Both"), c("Both")))
                                    ),
                                    checkboxInput("subsetModcorr2",strong("Subset Column Options:", style = "color:#456dae"),FALSE),
                                    conditionalPanel(condition="input.subsetModcorr2",
                                                     uiOutput("subsetcorr2")),
                                    checkboxInput("bubbleColorRange", strong("Force Color Range from -1 to 1"), FALSE),
                                    div(style = "display:inline-block", checkboxInput("graphics8", strong("Graphing options", style = "color:#456dae"), TRUE)),
                                    div(style = "display:inline-block", helpPopup("Graphing options:", "These options allow the user to make adjustments to the plot (i.e. width and height).")),
                                    conditionalPanel(condition = "input.graphics8",
                                                     sliderInput("PlotWidth3", "Plot Width", min = 500, max = 5000, value = 800, step = 100),
                                                     sliderInput("PlotHeight3", "Plot Height", min = 500, max = 5000, value = 600, step = 100),
                                                     numericInput("PlotRes3", "Plot resolution (for downloaded plot):", min = 72, max = 400, value = 72, step = 1))
                                ),
                                box(title = "Plot", width = 8, status = "primary", solidHeader = FALSE,
                                    downloadButton('download_plot_data', "Download Data"),
                                    downloadButton('download_corr_plot', "Download Figure"),
                                    plotOutput("correlations_plot")
                                )
                              )
                     ),
                     tabPanel("Result Tables",
                              fluidRow(
                                box(title = "Filtering Options", width = 4, status = "primary", solidHeader = FALSE,
                                    uiOutput("TypeVariable2"),
                                    uiOutput("var_switch2"),
                                    uiOutput("WithVariable2"),
                                    div(style = "display:inline-block", checkboxInput("TSoptions3", strong("P-Value and Correlation Value Subsetting Options:", style = "color:#456dae"), FALSE)),
                                    conditionalPanel(condition = "input.TSoptions3",
                                                     selectInput("correction_method.corr2","Multiple testing correction:",c("FDR","Bonferroni","Raw"),"Raw"),
                                                     numericInput("Alpha", "Significance threshold:", min = 0, max = 1, step = .01, value = .05),
                                                     checkboxInput("corrval",label = strong("Filter on correlation value:"),value = FALSE),
                                                     conditionalPanel(condition = "input.corrval == true",
                                                                      numericInput("corrval1", label = "Correlation value cut off:",
                                                                                   min = 0, max = 1, value = 0, step = 0.1),
                                                                      selectInput("corrsign","Correlation value sign:", c("+", "-", "Both"), c("Both")))
                                    ),
                                    uiOutput("Visit")
                                ),
                                box(title = "Results Table", width = 8, status = "primary", solidHeader = FALSE,
                                    downloadButton('download_data', "Download Table"),
                                    dataTableOutput("correlation_table")
                                )
                              )
                     ),
                     tabPanel("Pairwise Scatterplots",
                              fluidRow(
                                box(title = "Plotting Options", width = 4, status = "primary", solidHeader = FALSE, 
                                    uiOutput("TypeVariable6"),
                                    uiOutput("Visit3"),
                                    uiOutput("corr_Variable"),
                                    uiOutput("WithVariable4"),
                                    checkboxInput("log_scale", "lg2 Transform", value = FALSE),
                                    checkboxInput("plot_reg", "Plot Regression Line", value = FALSE),
                                    checkboxInput("plot_loess", "Fit loess curve", value = FALSE),
                                    conditionalPanel(condition = "input.plot_loess",
                                                     numericInput("span", "Span (0 to 1):", min = 0, max = 1, value = 1, step = .05)),
                                    div(style = "display:inline-block", checkboxInput("graphics10", strong("Graphing options", style = "color:#456dae"), TRUE)),
                                    div(style = "display:inline-block", helpPopup("Graphing options:", "These options allow the user to make adjustments to the plot (i.e. width and height).")),
                                    conditionalPanel(condition = "input.graphics10",
                                                     numericInput("Point_size", "Circle size (1 to 5):",min=1,max=5,value=2,step=.5),
                                                     numericInput("axis_text_size", "Axis text size:", min = 10, value = 10, step = 1),
                                                     numericInput("axis_label_size", "Axis label size:", min = 12, value = 12, step = 1),
                                                     sliderInput("PlotWidth5", "Plot Width", min = 500, max = 5000, value = 800, step = 100),
                                                     sliderInput("PlotHeight5", "Plot Height", min = 500, max = 5000, value = 600, step = 100),
                                                     numericInput("PlotRes5", "Plot resolution (for downloaded plot):", min = 72, max = 400, value = 72, step = 1))
                                ),
                                box(title = "Plot", width = 8, status = "primary", solidHeader = FALSE,
                                    downloadButton('download_scatter_data', "Download Data"),
                                    downloadButton('download_scatter', "Download Figure"),
                                    plotOutput("correlations_scatter_plot")
                                )
                              )
                     )
              )
            )
    )
  )
)

shinyUI(dashboardPage(skin="blue",
  header,
  sidebar,
  body
))

}


