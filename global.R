enableBookmarking(store = "server")
load("moduleinfo.rda")
load("moduleinfo_rna.rda")
moduleinfo2<-qusage::read.gmt("BaylorModules.gmt")
x1<-unlist(moduleinfo2)
x2<-rep(names(moduleinfo2),times=lapply(moduleinfo2,length))
x2<-gsub("_",".",x2)
moduleinfo2<-data.frame(SYMBOL=x1,Module=x2)
names(moduleinfo)[6] <- "affy"
names(moduleinfo)[8] <- "Module_V3"
module_annotations <- read.csv("v2_annotated_module_list.csv", header = T)
moduleinfo$Modulev2_Annotation <- module_annotations[match(moduleinfo$Module, module_annotations$Module), 2]
moduleinfo_rna$Modulev2_Annotation <- module_annotations[match(moduleinfo_rna$Module, module_annotations$Module), 2]

get_all_tables <- function(dat, my_i){
  nam <- names(dat)
  index1 <- grep("p.Value",nam,fixed=T)
  index2 <- grep("log10",nam,fixed=T)
  p.nams <- nam[setdiff(index1,index2)]
  nams <- gsub("p.Value.for.Estimate.of.", "", p.nams)
  sc_n <- paste0("sc_", my_i)
  assign(sc_n, data.frame(dat[, c(1, 2)], dat[, grep(nams[my_i], names(dat))]))
  write.csv(get(paste0("sc_", my_i)), file = paste0(nams[my_i], ".csv"))
}

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

source("HeatMap2.txt")

mycorrection<-function(dat,alpha,correction){
  
  k<-1:length(dat)
  index<-alpha*k/length(dat)
  sort.dat<-dat[order(dat)]
  
  if(correction=="FDR") maxk<-ifelse( length(which(sort.dat<=index))==0,0,max(which(sort.dat<=index)))
  if(correction=="BONF") maxk<-ifelse( length(which(sort.dat<=(alpha/length(sort.dat))))==0,0,max(which(sort.dat<=(alpha/length(sort.dat)))))
  if(correction=="RAW") maxk<-max(which(sort.dat<=alpha))
  
  percentage<-maxk/length(k)
  return(percentage)
}

flowlistmaker <- function(dat, var_name, alpha = 0.05, method = "FDR"){
  
  var_name_1 <- paste0("P.Value for ", var_name)
  pvalues <- dat[, var_name_1]
  sorted_p <- pvalues[order(pvalues)]
  index3 <- grep(var_name, names(dat), fixed = TRUE)
  
  #Fixing bug if contrast are written very similarly where
  #grep can't tell the difference        
  allpvalnames_ind <- grep("P.Value", names(dat)[index3])
  pvalnames <- names(dat)[index3][allpvalnames_ind]
  #log10pval_index <- grep("Sign_neg_log10", names(dat)[index3])
  #ind<-setdiff(allpvalnames_ind,log10pval_index)
  #pvalnames <- names(dat)[index3][ind]
  n_pvalnames <- length(pvalnames)
  
  ############################
  newdat <- data.frame(dat[,1], dat[, index3], p.adjust(pvalues,"fdr"), p.adjust(pvalues, "bonferroni"))
  names(newdat) <- c("Flow_Variable", names(dat)[index3], paste0(var_name,"FDR"), paste0(var_name,"BONF"))
  
  sort.dat <- newdat[order(dat[, var_name_1]),]
  if(method=="FDR"){    #k<-1:(dim(dat)[1])
                        #index<-alpha*k/dim(dat)[1]
                        maxk<-ifelse( length(which(p.adjust(sorted_p, "fdr")<=alpha))==0,0,max(which(p.adjust(sorted_p, "fdr")<=alpha)))
                        if(maxk!=0){flow_list<-sort.dat[1:maxk,]}
                        if(maxk==0){flow_list<- "No Variables Present"}       
  }
  if(method=="Raw"){       maxp <- ifelse( length(which(sorted_p<=alpha))==0,0,max(which(sorted_p<=alpha)))
                           if(maxp!=0){flow_list<-sort.dat[1:maxp,]}
                           if(maxp==0){flow_list<- "No Variables Present"}
  }
  if(method=="Bonferroni"){
                              maxp<-ifelse( length(which(p.adjust(sorted_p, "bonferroni")<=alpha))==0,0,max(which(p.adjust(sorted_p, "bonferroni")<=alpha)))
                              if(maxp!=0){flow_list<-sort.dat[1:maxp,]}
                              if(maxp==0){flow_list<- "No Variables Present"}
  }
  return(data.frame(Flow_Variable = flow_list))
}

movetolast <- function(dat, move) {
  dat[c(setdiff(names(dat), move), move)]
}

genelistmaker<-function(dat,var_name,alpha=.05,method="FDR",module_merge=FALSE){
  pvalues<-dat[,var_name]
  sorted_p<-pvalues[order(pvalues)]
  index3<-grep(substring(var_name,24),names(dat),fixed=T)
  
  #Fixing bug if contrast are written very similarly where
  #grep can't tell the difference        
  allpvalnames_ind<-grep("p.Val",names(dat)[index3])
  allpvalnames<-names(dat)[index3][allpvalnames_ind]
  log10pval_index<-grep("log10",names(dat)[index3])
  ind<-setdiff(allpvalnames_ind,log10pval_index)
  pvalnames<-names(dat)[index3][ind]
  n_pvalnames<-length(pvalnames)
  
  if(n_pvalnames>1){div<-length(index3)/n_pvalnames
  indicator1<-nchar(pvalnames)  
  trueval<-nchar(var_name)
  trueindex<-which(indicator1==trueval)
  newindex<-c()
  for(i in 1:div){
    newindex[i]<-index3[(trueindex+n_pvalnames*(i-1))]     }
  index3<-newindex  
  }
  ############################
  newdat<-cbind(dat$PROBE_ID,dat$SYMBOL,dat[,index3],p.adjust(pvalues,"fdr"),p.adjust(pvalues,"bonferroni"))
  names(newdat)<-c("PROBE_ID","SYMBOL",names(dat)[index3],paste(var_name,"FDR",sep=""),paste(var_name,"BONF",sep=""))
  
  sort.dat<-newdat[order(dat[,var_name]),]
  if(method=="FDR"){       #k<-1:(dim(dat)[1])
                           #index<-alpha*k/dim(dat)[1]
                           maxk<-ifelse( length(which(p.adjust(sorted_p,"fdr")<=alpha))==0,0,max(which(p.adjust(sorted_p,"fdr")<=alpha)))
                           if(maxk!=0){PROBE_ID<-sort.dat[1:maxk,]}
                           if(maxk==0){PROBE_ID<- "No Genes Present"}       
  }
  if(method=="Raw"){       maxp<- ifelse(length(which(sorted_p<=alpha)) == 0,0,max(which(sorted_p<=alpha)))
                           if(maxp != 0){PROBE_ID<-sort.dat[1:maxp,]}
                           if(maxp == 0){PROBE_ID<- "No Genes Present"}
  }
  if(method=="Bonferroni"){     maxp<-ifelse( length(which(p.adjust(sorted_p,"bonferroni")<=alpha))==0,0,max(which(p.adjust(sorted_p,"bonferroni")<=alpha)))
                                if(maxp!=0){PROBE_ID<-sort.dat[1:maxp,]}
                                if(maxp==0){PROBE_ID<- "No Genes Present"}
  }
  if(module_merge==TRUE){ PROBE_ID<-merge(PROBE_ID,moduleinfo,by="PROBE_ID",all.x = TRUE)
  index4<- which( !(names(PROBE_ID)%in%c("PROBE_ID","SYMBOL","Module","Module_V3","Modulev2_Annotation","Modulev3_Annotation")))   
  name2<-names(PROBE_ID)[index4]
  PROBE_ID<-cbind(PROBE_ID[,c("PROBE_ID","SYMBOL","Module","Modulev2_Annotation")],PROBE_ID[,index4])    
  names(PROBE_ID)<-c("PROBE_ID","SYMBOL","Module","Modulev2_Annotation",name2)
  PROBE_ID<-PROBE_ID[order(PROBE_ID[,var_name]),]
  }
  PROBE_ID <- data.frame(PROBE_ID)
  return(PROBE_ID)
}


setPalettes = function(n) {
  pal = list("Spectral", "PuOr", "PiYG","RdGy","PRGn","RdBu","RdYlBu")
  sel = pal[1:n]
  lapply(sel, getColor)
}

#BrBG PiYG PRGn PuOr RdBu RdGy RdYlBu RdYlGn Spectral


getLevels = function(x) {
  as.character(levels(x))
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

dataManipulate <- function(y, x, colname, ref_var, ref_level, long = FALSE,
                           subjects, keep = TRUE, format = "Probes",
                           allsamples = TRUE) {
  index_samples <- which(names(y) %in% x[, colname])
  if (format == "Probes") {
    index_ps <- which(names(y) %in% c("PROBE_ID", "SYMBOL"))
  }
  if (format == "Modules") {
    index_ps <- which(names(y) %in% c("Module"))
  }
  y <- y[, c(index_ps, index_samples)]
  if (allsamples == FALSE) {
    if (long == TRUE) {
      index_refvar <- which(names(x) == ref_var)
      lev <- x[which(x[, index_refvar] == ref_level), ]
      lev_rm <- x[which(x[, index_refvar] != ref_level), ]
      index_subid <- which(names(x) == subjects)
      lev_subids <- unique(lev[, index_subid])
      lev_rm_subids <- unique(lev_rm[, index_subid])
      index_include <- lev_subids %in% lev_rm_subids
      complete_ids <- lev_subids[index_include]
      lev_norm <- vector("list", length(complete_ids))
      for (i in 1:(length(lev_norm))) {
        lev_norm[[i]] <- getNorm(y = y, x = x, colname = colname,
                                 id = complete_ids[i], index_sid = index_subid,
                                 mynames = c(subjects, "columnname", ref_var),
                                 index_refvar = index_refvar,
                                 ref_level = ref_level, keep = keep)
      }
      lev_norm <- do.call("cbind", lev_norm)
      #design_norm <- x[(x[, index_subid] %in% complete_ids), ]
      design_norm <- x[match(colnames(lev_norm), x$columnname, nomatch = 0),]
      #if (keep == FALSE) {
        #design_norm <- design_norm[(design_norm[, index_refvar] != ref_level), ]
      #}
    }
    if (long == FALSE) {
      index_refvar <- which(names(x) == ref_var)
      lev <- x[which(x[, index_refvar] == ref_level), ]
      lev_rm <- x[which(x[, index_refvar] != ref_level), ]
      lev_names <- lev[, colname]
      lev_rm_names <- lev_rm[, colname]
      lev_norm <- getNorm2(y = y, lev_names, lev_rm_names, keep = keep)
      design_norm <- x[match(colnames(lev_norm), x$columnname, nomatch = 0),]
      #if (keep == FALSE) {
        #design_norm <- lev_rm
      #}
      #if (keep == TRUE) {
        #design_norm <- x
      #}
    }
  } else {
    means_y <- as.vector(apply(y[, -(1:2)], 1, mean))
    lev_norm <- y[, -c(1:2)] - means_y
    design_norm <- x
  }
  final_norm_shuff <- as.matrix(lev_norm)
  colnames(final_norm_shuff) <- names(lev_norm)
  if (format == "Probes") {
    rownames(final_norm_shuff) <- y$SYMBOL
  }
  if (format == "Modules") {
    rownames(final_norm_shuff) <- y$Module
  }
  return(list(heatexp = final_norm_shuff, heatdes = design_norm))
}


data.manipulate<-function(exp,des,basevariable,baselevel,longitudinal=FALSE,subjects,lg2=FALSE,keepbase=TRUE,format="Probes",allsamples=T){
  index.samples<-which(names(exp) %in% des$columnname)

  if(format=="Probes"){ 
            index.ps<-which(names(exp) %in% c("PROBE_ID","SYMBOL"))
                               }
  if(format=="Modules"){
            index.ps<-which(names(exp) %in% c("Module"))
                                }
exp<-exp[,c(index.ps,index.samples)]

if(allsamples==FALSE){
if (longitudinal==TRUE){
            index.basevar<-which(names(des)==basevariable)
            base<-des[which(des[,index.basevar]==baselevel),]
		all_time_base_rm<-des[which(des[,index.basevar]!=baselevel),]

            index.subjectid<-which(names(des)==subjects)
		base_subject_ids<-unique(base[,index.subjectid])
		
            all_time_base_rm_subject_ids<-unique(all_time_base_rm[,index.subjectid])

		index.include<-base_subject_ids%in%all_time_base_rm_subject_ids

		#"complete_ids" is a list of subject ids that have a "basevariable" of "baselevel" as well as other "baselevel" observations

		complete_ids<-base_subject_ids[index.include]
         
		columnname<-names(exp)[-index.ps]

get.norm<-function(id,mynames,index.sid,index.basevar,baselevel,keep=TRUE){ data1<-des[which(des[,index.sid]==id),mynames]
				base_name<-data1[which(data1[,3]==baselevel),"columnname"]
                        if(keep==TRUE){
                        base_norm1<-exp[,as.character(data1$columnname)]-exp[,as.character(base_name)]
                        }
                        if(keep==FALSE){
                        
                        base_norm1<-exp[,setdiff(as.character(data1$columnname),as.character(base_name))]-exp[,as.character(base_name)]
                        if(length(setdiff(as.character(data1$columnname),as.character(base_name)))==1){
                                                  base_norm1=data.frame(base_norm1)
                                                  names(base_norm1)<-as.character(data1$columnname[which(as.character(data1$columnname)!=as.character(base_name))])}
				}
				return(base_norm1)
			     }


            base_norm<-get.norm(id=complete_ids[1],mynames=c(subjects,"columnname",basevariable),index.sid=index.subjectid,index.basevar=index.basevar,baselevel=baselevel,keep=keepbase)
		for(i in 2:(length(complete_ids))){base_norm<-cbind(base_norm,get.norm(complete_ids[i],c(subjects,"columnname",basevariable),index.subjectid,index.basevar,baselevel,keep=keepbase))}

		design_base_norm<-des[(des[,index.subjectid]%in%complete_ids),]
            if(keepbase==FALSE){design_base_norm<-design_base_norm[which(design_base_norm[,index.basevar]!=baselevel),]}
                        }



if (longitudinal==FALSE){
		index.basevar<-which(names(des)==basevariable)
            base<-des[which(des[,index.basevar]==baselevel),]
            base_rm<-des[which(des[,index.basevar]!=baselevel),]
            base_names<-base[,"columnname"]
            base_rm_names<-base_rm[,"columnname"]
             
get.norm2<-function(baseline,rest,keep=TRUE){index.base<-which(names(exp)%in%baseline)
                        index.rest<-which(names(exp) %in%rest )
				if(keep==FALSE){
				normdat<-exp[,index.rest]-apply(exp[,index.base],1,mean)
				#names(normdat)<-names(exp)[index.rest]
				}
                        if(keep==TRUE){
				#normdat<-exp[,-c(1,2)]-apply(exp[,index.base],1,mean)
                        normdat<-exp[,c(index.base,index.rest)]-apply(exp[,index.base],1,mean)
                        #names(normadat)<-names(exp[,-c(1,2)])  
                        }
                        return(normdat)
		           }
		base_norm<-get.norm2(base_names,base_rm_names,keep=keepbase)
            
            if(keepbase==FALSE){design_base_norm<-base_rm}
            if(keepbase==TRUE){design_base_norm<-des}
			}	
    }

if(allsamples==TRUE){
    myfactor<-as.vector(apply(exp[,-(1:2)],1,mean))
    base_norm<-exp[,-(1:2)]-myfactor
    design_base_norm<-des
    }
if(lg2==TRUE){base_norm=log(base_norm,2)}

#Data is normalized with an updated design file.  Normalized exp file still doesn't have
#Probe and Symbol attached will do in a moment.

#Now we rearrange the columns based on "ordernames".
# if(length(ordernames)>0){
#            data1<-design_base_norm[,c("Columnname",ordernames)]
#            inside<-paste("data1$",ordernames,sep="")
#            column_index<-eval(parse(text=paste("order(",paste(inside,collapse=","),")",sep="")))
#	      sorted.data<-base_norm[,column_index]}
#        if(length(ordernames)==0){sorted.data=base_norm}
		
#final_norm_shuff<-cbind(exp[,index.ps],sorted.data)
#names(final_norm_shuff)<-c(names(exp[,index.ps]),names(sorted.data))

#final_norm_shuff<-as.matrix(sorted.data)
 final_norm_shuff<-as.matrix(base_norm)
 colnames(final_norm_shuff)<-names(base_norm)
 if(format=="Probes"){rownames(final_norm_shuff)<-exp$SYMBOL}
 if(format=="Modules"){rownames(final_norm_shuff)<-exp$Module}


return(list(heatexp=final_norm_shuff,heatdes=design_base_norm))
}


moduleinfo1<-read.table("module_v2only_complete.txt",header=T,sep=",")
modnames<-unique(moduleinfo1$Module)
y<-unlist(strsplit(as.character(modnames),".",fixed=TRUE))
y<-matrix(y,byrow=T,nrow=260,ncol=2)
sortmodnames<-as.character(modnames)[order(y[,1],y[,2])]
Modulelist<-list(as.character(moduleinfo1$PROBE_ID[which(moduleinfo1$Module==sortmodnames[1])]),as.character(moduleinfo1$PROBE_ID[which(moduleinfo1$Module==sortmodnames[2])]))
for(i in 3:length(sortmodnames)){Modulelist<-c(Modulelist,list(as.character(moduleinfo1$PROBE_ID[which(moduleinfo1$Module==sortmodnames[i])])))}
names(Modulelist)<-sortmodnames

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  #library(grid)
  
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


## Shiny Modules

# Filtering Modules
filterOptsUI <- function(id){
  # Create a namespace function using the provided id
  ns <- NS(id)
  condCall1 <- paste0("input['",ns("TSoptions"),"']")
  condCall2 <- paste0("input['",ns("showfc"),"'] == true")
  
  tagList(
    div(style = "display:inline-block", checkboxInput(ns("TSoptions"), strong("Testing and fold change subsetting options", style = "color:#456dae"), FALSE)),
    div(style = "display:inline-block", helpPopup("Testing and fold change subsetting options", "These options allow the user to filter on signifance threshold and fold change.")),
    conditionalPanel(condition = condCall1,
                     selectInput(ns("correction_method"),"Multiple testing correction:",
                                 c("FDR","Bonferroni","Raw"),"Raw"),
                     numericInput(ns('alphalevel'),"Significance threshold:", min = 0, max = 1, value = 0.05, step = 0.025),
                     checkboxInput(inputId = ns("showfc"),label = strong("Filter on fold change"),value = FALSE),
                     conditionalPanel(condition = condCall2,
                                      numericInput(inputId = ns("fcval"), label = "FC cut off:",
                                                   min = 0, max = 10, value = 0, step = 0.25),
                                      selectInput(ns("sign"),"Fold change sign:", c("+", "-", "Both"), c("Both")))
    )
  )
}

filterOpts <- function(input, output, session, data, comparison = NULL, data.type = NULL, paloOrFirst = NULL, geneList = NULL){
  results <- data()
  if(data.type %in% c("genes", "flow", "metab")){
    if(data.type == "genes"){
      results <- results[,c(1,2,grep(comparison(), colnames(results), fixed = TRUE))]
    } else{
      results <- results[,c(1,grep(comparison(), colnames(results), fixed = TRUE))]
    }
    fdr <- p.adjust(results[,5], method = "fdr")
    bonf <- p.adjust(results[,5], method = "bonferroni")
    results <- do.call("cbind", list(results, FDR = fdr, Bonf = bonf))
    if(data.type == "genes"){
      colnames(results) <- c("PROBE_ID", "SYMBOL", "Log2FC", "Test.Statistic", "P.Value", "FDR", "Bonferroni")
    } else if(data.type == "flow"){
      results[,2] <- NULL
      colnames(results) <- c("Flow.Variables","Log2FC", "Test.Statistic", "P.Value", "FDR", "Bonferroni")
    } else if(data.type == "metab"){
      results[,2] <- NULL
      colnames(results) <- c("Metab.Variables","Log2FC", "Test.Statistic", "P.Value", "FDR", "Bonferroni")
    }
    results <- data.frame(results)
    if(input$correction_method == "Raw"){
      results <- results[which(results$P.Value <= input$alphalevel),]
      results <- results[order(results$P.Value, decreasing = FALSE),]
    }
    if(input$correction_method == "FDR"){
      results <- results[which(results$FDR <= input$alphalevel),]
      results <- results[order(results$FDR, decreasing = FALSE),]
    }
    if(input$correction_method == "Bonferroni"){
      results <- results[which(results$Bonferroni <= input$alphalevel),]
      results <- results[order(results$Bonferroni, decreasing = FALSE),]
    }
    if(input$showfc){
      if(input$sign=="Both"){results <- results[which(abs(results$Log2FC) > input$fcval),]}
      if(input$sign=="+"){results <- results[which(results$Log2FC > input$fcval),]}
      if(input$sign=="-"){results <- results[which(results$Log2FC < -input$fcval),]}
    }
    if(nrow(results) == 0){
      if(data.type == "genes"){
        results <- data.frame("No Genes Present")
        names(results) <- "PROBE_ID"
      } else if(data.type == "flow"){
        results <- data.frame("No Flow Variables Present")
        names(results) <- "Flow.Variables"
      } else if(data.type == "metab"){
        results <- data.frame("No Metabolites Present")
        names(results) <- "Metab.Variables"
      }
    }
  } 
  if(data.type == "geneSets") {
    if(is.null(results)){return(NULL)}
    if(length(comparison()) < 1){return(NULL)}
    results <- results[which(results$Comparison %in% comparison()),]
    first_comp = results[which(results$Comparison %in% comparison()[1]),]
    if(input$showfc){
      if(input$sign == "Both"){
        first_cut <- unique(first_comp[which(abs(first_comp$Log2FC) >= input$fcval), ]$pathway.name)
        palo_cut <- unique(results[which(abs(results$Log2FC) >= input$fcval), ]$pathway.name)
      } else if(input$sign == "+"){
        first_cut <- unique(first_comp[which(first_comp$Log2FC >= input$fcval), ]$pathway.name)
        palo_cut <- unique(results[which(results$Log2FC >= input$fcval), ]$pathway.name)
      } else if(input$sign == "-"){
        first_cut <- unique(first_comp[which(first_comp$Log2FC <= -input$fcval), ]$pathway.name)
        palo_cut <- unique(results[which(results$Log2FC <= -input$fcval), ]$pathway.name)
      }
      if(length(comparison()) == 1 || paloOrFirst() == 1){
        results <- results[which(results$pathway.name %in% palo_cut), ]
      }
      else{
        results <- results[which(results$pathway.name %in% first_cut), ]
      }
    }
    if(input$correction_method == "Raw"){
      pval <- "p.Value"
    } else if(input$correction_method == "FDR"){
      pval <- "FDR"
    } else if(input$correction_method == "Bonerroni"){
      pval <- "Bonf"
    }
    if(paloOrFirst() == 1 || length(comparison()) == 1){
      keep <- unique(results[which(results[[pval]] <= input$alphalevel),]$pathway.name)
    } else {
      results1 <- results[which(results$Comparison %in% comparison()[1]),]
      keep <- unique(results1[which(results1[[pval]] <= input$alphalevel),]$pathway.name)
    }
    if(length(keep) == 0){return(NULL)}
    results <- results[which(results$pathway.name %in% keep),]
    results.temp <- list()
    for(i in 1:length(comparison())){
      results.temp[[i]] <- results[which(results$Comparison %in% comparison()[i]),]
      if(i == 1){
        results.temp[[i]] <- results.temp[[i]][order(results.temp[[i]]$Log2FC, decreasing = TRUE),]
      } else {
        results.temp[[i]] <- results.temp[[i]][match(results.temp[[1]]$pathway.name, results.temp[[i]]$pathway.name, nomatch = 0),]
      }
      Index <- 1:nrow(results.temp[[i]])
      results.temp[[i]] <- cbind(Index, results.temp[[i]])
    }
    results <- do.call("rbind", results.temp)
    results[,2] = factor(results[,2], levels = unique(results[,2]))
    sig = which(results[[pval]] <= input$alphalevel)
    results$SIG = "Not Significant"
    results$SIG[sig] = "Significant"
  }
  if(data.type == "percents"){
    dataset <- data()
    nam <- colnames(dataset)
    index <- grep("P.Value",nam,fixed=T)
    p.names <- gsub("P.Value for ","",nam[index])
    y <- as.matrix(dataset[,index])
    if(input$correction_method == "Raw"){
      dat<-y
    } else if(input$correction_method == "FDR"){
      dat<-apply(y,2,p.adjust,method="fdr")
    } else if(input$correction_method == "Bonferroni"){
      dat<-apply(y,2,p.adjust,method="bonferroni")
    }
    dat[dat <= input$alphalevel] <- 2
    dat[dat != 2] <- 0
    dat[dat == 2] <- 1
    
    fcindex <- grep("Estimate",colnames(dataset),fixed=T)
    dataset.fcindex <- dataset[,fcindex]
    if(input$showfc == TRUE){
      dataset.fcindex[dataset.fcindex < input$fcval & dataset.fcindex > -input$fcval] <- 0
      sign.dat <- sign(dataset.fcindex)
      dat <- sign.dat*dat
      if(input$sign == "+"){
        dat[dat == -1] <- 0
      }
      if(input$sign == "-"){
        dat[dat == 1] <- 0
      }
    } else {
      sign.dat <- sign(dataset.fcindex)
      dat <- sign.dat*dat
    }
    dat <- data.frame(cbind(PROBE_ID = as.character(dataset$PROBE_ID),dat))
    dat <- merge(dat, geneList(), by = "PROBE_ID")
    dat$PROBE_ID <- NULL
    count_matrix <- list()
    for(i in 1:length(index)){
      count_matrix[[i]] <- aggregate(as.formula(paste(names(dat)[i],"~","Module",sep="")),data=dat,sum)
      rownames(count_matrix[[i]]) <- as.character(count_matrix[[i]]$Module)
      count_matrix[[i]]$Module <- NULL
    }
    count_matrix <- do.call(cbind, count_matrix)
    freq <- data.frame(table(geneList()$Module))$Freq
    names(freq) <- data.frame(table(geneList()$Module))$Var1
    freq.matched <- freq[match(rownames(count_matrix), names(freq), nomatch = 0)]
    prop_matrix <- count_matrix/freq.matched
    missingMods <- names(freq)[which(!names(freq) %in% names(freq.matched))]
    missingData <- data.frame(matrix(0,nrow = length(missingMods), ncol = ncol(prop_matrix)))
    rownames(missingData) <- missingMods
    colnames(missingData) <- colnames(prop_matrix)
    colnames(prop_matrix) <- p.names
    prop_matrix2 <- prop_matrix
    prop_matrix2[prop_matrix < .1 & prop_matrix > -.1] <- 0
    results <- list(prop_matrix = prop_matrix, prop_matrix2 = prop_matrix2)
  }
  return(results)
}




## Create data to be plotted in heatmap

# Subset and Ordering Modules
subsetAndOrderUI <- function(id){
  ns <- NS(id)
  condCall1 <- paste0("input['",ns("subsetCols"),"']")
  condCall2 <- paste0("input['",ns("orderCols"),"'] == true")
  
  tagList(
    div(style = "display:inline-block", checkboxInput(ns("subsetCols"),strong("Subset column options", style = "color:#456dae"),FALSE)),
    div(style = "display:inline-block", helpPopup("Subset column options", "These options allow the user to plot a subset of the samples within the expression data set.")),
    conditionalPanel(condition = condCall1,
                     uiOutput(ns('subsetVar1')),
                     uiOutput(ns('subsetVal1')),
                     uiOutput(ns('subsetVar2')),
                     uiOutput(ns('subsetVal2'))
    ),
    br(),
    div(style = "display:inline-block", checkboxInput(ns('orderCols'), strong("Ordering column options", style = "color:#456dae"), FALSE)),
    div(style = "display:inline-block", helpPopup("Ordering column options", "These options allow the user to specify column ordering by column variable.")),
    conditionalPanel(condition = condCall2,
                     uiOutput(ns('orderingVars'))
    )
  )
}

subsetAndOrderRenderUI <- function(input, output, session, des){
  return(
    tagList(
      output$subsetVar1 <- renderUI({
        ns <- session$ns
        selectizeInput(ns("subsetVar1"),"Variable 1 to subset heatmap:",c(colnames(des())), selected = NULL,multiple = TRUE,options = list(maxItems = 1))
      }),
      output$subsetVal1 <- renderUI({
        ns <- session$ns
        selectizeInput(ns("subsetVal1"),"Value(s) of Variable to subset heatmap:", unique(des()[,input$subsetVar1]), multiple=TRUE)
      }),
      output$subsetVar2 <- renderUI({
        ns <- session$ns
        selectizeInput(ns("subsetVar2"),"Variable 2 to subset heatmap:",c(colnames(des())),selected = NULL, multiple = TRUE,options = list(maxItems = 1))
      }),
      output$subsetVal2 <- renderUI({
        ns <- session$ns
        selectizeInput(ns("subsetVal2"),"Value(s) of Variable to subset heatmap:", unique(des()[,input$subsetVar2]), multiple=TRUE)
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
    if(!is.null(input$subsetVar1) & !is.null(input$subsetVar2)){
      x.sub <- x[,which(des()[,input$subsetVar1] %in% input$subsetVal1 & des()[,input$subsetVar2] %in% input$subsetVal2)]
      design <- des()[which(des()[,input$subsetVar1] %in% input$subsetVal1 & des()[,input$subsetVar2] %in% input$subsetVal2),]
      design[,input$subsetVar1] <- factor(as.character(design[,input$subsetVar1]), levels = input$subsetVal1)
      design[,input$subsetVar2] <- factor(as.character(design[,input$subsetVar2]), levels = input$subsetVal2)
    } else if(!is.null(input$subsetVar1) & is.null(input$subsetVar2)){
      x.sub <- x[,which(des()[,input$subsetVar1] %in% input$subsetVal1)]
      design <- des()[which(des()[,input$subsetVar1] %in% input$subsetVal1),]
      design[,input$subsetVar1] <- factor(as.character(design[,input$subsetVar1]), levels = input$subsetVal1)
    } else if(is.null(input$subsetVar1) & is.null(input$subsetVar2)){
      x.sub <- x
      design <- des()
    }
    if(input$orderCols){
      order.vars <- do.call(order, design[,input$orderingVars, drop = FALSE])
      x.ord <- x.sub[,order.vars]
      colAnnot <- design[order.vars,input$orderingVars,drop = FALSE]
      design <- design[order.vars,]
    } else {
      orderingVars <- c(input$subsetVar1, input$subsetVar2)
      order.vars <- do.call(order, design[,orderingVars, drop = FALSE])
      x.ord <- x.sub[,order.vars]
      colAnnot <- design[order.vars,orderingVars, drop = FALSE]
      if(is.null(input$subsetVar1) & is.null(input$subsetVar2)){
        orderingVars <- colnames(des())[1]
        x.ord <- x.sub[,order(design[,orderingVars])]
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
      colAnnot[,i] <- as.numeric(as.factor(as.character(colAnnot[,i])))
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


# Column Clustering Modules
colClusterUI <- function(id){
  ns <- NS(id)
  condCall1 <- paste0("input['",ns("clusterOptions"),"']")
  condCall2 <- paste0("input['",ns("clusterGroups"),"'] == true")
  
  tagList(
    div(style = "display:inline-block", checkboxInput(ns('clusterOptions'), strong("Clustering column options", style = "color:#456dae"), FALSE)),
    div(style = "display:inline-block", helpPopup("Cluster column options", "These options allow the user to cluster the samples and label the groups they are clustered in.")),
    conditionalPanel(condition = condCall1,
                     checkboxInput(ns("colCluster"),"Cluster samples (columns)",FALSE),
                     checkboxInput(ns('clusterGroups'), "Show cluster groups", FALSE),
                     conditionalPanel(condition = condCall2,
                                      div(style = "display:inline-block", uiOutput(ns('clusterCuts')))),
                                      div(style = "display:inline-block", infoPopup("Cluster Cuts", "Optimal number calculated using Dunn's Index.",
                                                                                    placement = "right", trigger = "click")))
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
  uploadDescription <- paste0("Allows the user to provide their own list of ", varType, " (CSV) to plot. The CSV file should contain a single column named 'PROBE_ID' or 'SYMBOL', 
                               depending on whether the list provided is the PROBE ID's or gene symbols.") 
  
  tagList(
    div(style = "display:inline-block", checkboxInput(ns("uploadVars"), strong(uploadType, style = "color:#456dae"), FALSE)),
    div(style = "display:inline-block", infoPopup(uploadType, uploadDescription, placement = "right", trigger = "click")),
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
    labelRows <- NULL
    #x <- data()[which(rownames(data()) %in% varnames[,1]),]
    x <- data()[match(varnames[,1], rownames(data()), nomatch = 0),]
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
    ddm <- NA
    labelRows = NA
  }
  return(list(ddm = ddm, x = x, labelRows = labelRows))
}

# Heatmap Graphing Options

graphOptionsUI <- function(id, varType = "transcripts"){
  ns <- NS(id)
  conCall1 <- condCall1 <- paste0("input['",ns("graphOptions"),"']")
  
  if(varType == "modules"){
    return(tagList(
      div(style = "display:inline-block", checkboxInput(ns("graphOptions"), strong("Graphing options", style = "color:#456dae"), FALSE)),
      div(style = "display:inline-block", helpPopup(ns("Graphing options"), "These options allow the user to make adjustments to the plot (i.e. width and height).")),
      conditionalPanel(condition = condCall1,
                       sliderInput(ns('graphWidth'), "Plot width", min = 400, max = 2000, value = 750, step = 25),
                       sliderInput(ns('graphHeight'), "Plot height", min = 400, max = 2000, value = 800, step = 25),
                       numericInput(ns('circleSize'),"Circle size (.25 to 4):",min=.25,max=4,value=1,step=.25),
                       numericInput(ns('fontSize'), "Font size:", min = 10, value = 10, step = 2),
                       numericInput(ns('legendSize'), "Legend size:", min = 1, max = 5, value = 1, step = .50),
                       numericInput(ns('treeHeight'), "Tree height:", min = 0, value = 50, step = 5),
                       numericInput(ns('plotResolution'), "Plot resolution (for downloaded plot):", min = 72, max = 400, value = 72, step = 1))
    ))
  }
  return(tagList(
    div(style = "display:inline-block", checkboxInput(ns("graphOptions"), strong("Graphing options", style = "color:#456dae"), FALSE)),
    div(style = "display:inline-block", helpPopup(ns("Graphing options"), "These options allow the user to make adjustments to the plot (i.e. width and height).")),
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


