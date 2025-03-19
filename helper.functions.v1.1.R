require(tidyverse)
require(cmapR)
require(factoextra)
require(FactoMineR)
require(ICC)
require(qvalue)
require(RColorBrewer)
require(gridExtra)
require(UpSetR)
require(gprofiler2)
require(cowplot)
require(ComplexHeatmap)
require(circlize)


#Color specification
sex_cols = c('M' = '#5555ff',
             'Male' = '#5555ff',
             'male' = '#5555ff',
             'F' = '#ff6eff',
             'Female' = '#ff6eff',
             'female' = '#ff6eff')

###1. GCT Parsing, Normalization, and Filtering function ------------------------------------
sm_to_gct_lfq <- function(x,md,col.suffix="totalIntensity",join.by="sample"){
  #Function to convert a spectrum mill LFQ ratio report into a gct file
  #Args:
  #	x: a spectrum mill lfq report stored as data frame 
  #	species: species to be selected from the table. If NA, all species are allowed
  #	col.suffix: the suffix of the columns that contain the data. Spetrum Mill
  #	provides total intensities (defualt), peptide counts, and spectra counts.
  #	join.by: Specifies the column that is common between the 
  #Returns:
  #	A gct object
  
  mat <- x %>% select(ends_with(col.suffix)) %>% as.matrix
  cn <- data.frame(col = colnames(mat)) %>% separate("col",c("sample",NA)," ") 
  
  rownames(mat) <- x$id
  rdesc <- select(x,-ends_with(col.suffix)) 
  cdesc <- left_join(cn,md,by=join.by) %>% rename(id = sample)
  
  colnames(mat) <- cdesc$id
  if(ncol(mat)!=nrow(cdesc)|nrow(mat)!=nrow(rdesc)){
    stop("Column or row annotations do not match 
		     the matrix size.")
  }
  
  return(new("GCT", mat=mat, rdesc=rdesc, cdesc=cdesc))
}

sm_to_gct <- function(x,md=NA,md.join = "sample"){
  #Function to convert a spectrum mill tmt ratio report into a gct file
  #Args:
  #	x: a spectrum mill tmt ratio report stored as data frame 
  #	md: Column metadata
  # md.join: Variable to join the metadata table
  #Returns:
  #	A gct object
  
  mat <- select(x,contains(":")) %>% as.matrix
  cdesc <- parse_sm_colnames_ratio(colnames(mat))
  rdesc <- select(x,-contains(" ")) 
  
  #If the current ids are not unique, use plex_tmtlabels
  if(sum(duplicated(cdesc$sample))>0){
    cdesc <- cdesc %>%
      mutate(id = paste0(plex,"_",tmtlabel))
  } else {
    cdesc <- cdesc %>% mutate(id = sample)
  } 
  
  #If there is sample metadata, add to cdesc
  if(!is.na(md)){
    cdesc <- left_join(cdesc,md,by=md.join) 
  } 
  
  if(ncol(mat)!=nrow(cdesc)|nrow(mat)!=nrow(rdesc)){
    stop("Column or row annotations do not match 
		     the matrix size.")
  }
  
  rownames(mat) <- x$id
  colnames(mat) <- cdesc$id
  
  return(new("GCT", mat=mat, rdesc=rdesc, cdesc=cdesc))
}

sm_to_gct_tmtintensity <- function(x,md=NA,md.join = "id"){
  #Function to convert a spectrum mill tmt ratio report into a gct file
  #Args:
  #	x: a spectrum mill tmt ratio report stored as data frame 
  #	md: Column metadata
  #Returns:
  #	A gct object
  
  mat <- select(x,contains(", ")) %>% as.matrix
  cdesc <- parse_sm_colnames_intensity(colnames(mat))
  rdesc <- select(x,-contains(" ")) 
  

  #If the current ids are not unique, use plex_tmtlabels
  if(sum(duplicated(cdesc$sample))>0){
    cdesc <- cdesc %>%
      mutate(id = paste0(plex,"_",tmtlabel))
  } else {
    cdesc <- cdesc %>% mutate(id = sample)
  } 
  
  #If there is sample metadata, add to cdesc
  if(!is.na(md)){
    cdesc <- left_join(cdesc,md,by=md.join) 
  } 
  
  
  if(ncol(mat)!=nrow(cdesc)|nrow(mat)!=nrow(rdesc)){
    stop("Column or row annotations do not match 
		     the matrix size.")
  }
  
  rownames(mat) <- x$id
  colnames(mat) <- cdesc$id
  
  return(new("GCT", mat=mat, rdesc=rdesc, cdesc=cdesc))
}



parse_sm_colnames_ratio <- function(x){
  #Function to parse the column names for TMT ratio datasets obtained from the process report of spectrum mill
  #Args:
  #	x: the column names for tmt ratios obtained by spectrum mill process report
  #Returns:
  #	A data frame with the plex, tmt label used, and sample name as columns
  
  data.frame(col = x) %>% separate("col",c("plex","tmtlabel","sample"),",") %>% 
    transmute(plex = plex,
              tmtlabel = str_match(tmtlabel," (.*):")[,2],
              sample = str_match(sample," (.*):")[,2])
  
}

parse_sm_colnames_intensity <- function(x){
  #Function to parse the column names for TMT ratio datasets obtained from the process report of spectrum mill
  #Args:
  #	x: the column names for tmt ratios obtained by spectrum mill process report
  #Returns:
  #	A data frame with the plex, tmt label used, and sample name as columns
  
  data.frame(col = x) %>% separate("col",c("plex","tmtlabel","sample"),", ")
}

gct_ratio <- function(x,denom){
  #Convert data to log2 ratios by dividing to the median of samples selected as denominators.
  #Args:
  #	x: a GCT object
  #	denom: Character vector indicating the sample ids to be used as denominators
  
  if(length(denom)>1){
    d <- apply(x@mat[,denom],1,mean,na.rm=T)
  } else {
    d <- x@mat[,denom]
  }
  
  
  #x <- subset_gct(x,cid = setdiff(colnames(x@mat),denom)) 
  x@mat <- log2(x@mat/d)
  return(x)
}

#' Convert intensities to ratios for multiplex studies
#' Convert data to log2 ratios by dividing to the median of samples selected as denominators.
#' 
#' @param x a GCT object
#' @param denom_sample Character vector indicating the id or ids to be used as 
#' denominators as found in the 'cdesc' column identified by `denom_column`
#' @param denom_column A string indicating the column in 'cdesc' which 
#' identifies the denominator or denominators to be used per plex
#' @param plex_column A string indicating the column in 'cdesc' which separate the datset by plex
#'
#' @return
#' @export
#'
#' @examples
gct_ratio_multiplex <- function(x,denom_samples, denom_column,plex_column = "plex"){
  x_full <- x
  
  for(i in unique(x@cdesc[,plex_column])  ){
    plex_samples <-  x@cdesc[,plex_column] == i 
    
    x <- subset_gct(x_full,cid = plex_samples)
    
    denom <- x@cid[x@cdesc[,denom_column] %in% denom_samples]
    
    if(length(denom)>1){
      d <- apply(x@mat[,denom],1,mean,na.rm=T)
    } else {
      d <- x@mat[,denom]
    }

    x@mat <- log2(x@mat/d)
    
    x_full@mat[,plex_samples] <- x@mat
  }
}


median_mad_norm <- function(x,mad=T){
  #Perform median normalization
  #Args:
  #	x: a GCT object
  #	mad: logic indicating to use mad normalization
  #Returns:
  #	A normalized GCT object
  if(mad){
    scale.factor <- mean(apply(x@mat,2,mad,na.rm=T))
    x@mat <- scale(x@mat,center = apply(x@mat,2,median,na.rm=T),
                   scale = apply(x@mat,2,mad,na.rm=T))
    x@mat <- x@mat*scale.factor
  } else{
    x@mat <- scale(x@mat,center = apply(x@mat,2,median,na.rm=T),
                   scale = F)
  }
  return(x)
}



remove_all_na <- function(x){
  #Remove proteins that have all missing values
  #Args:
  #	x: a GCT object
  #Returns:
  #	A GCT object without proteins that have all missing values
  x <- subset_gct(x,rid = which(rowSums(is.na(x@mat)) < ncol(x@mat))) 
  return(x)
}

remove_na <- function(x,pct){
  #Remove proteins that have more than the set percentage of missing values
  #Args:
  #	x: a GCT object
  #	pct: the percent cutoff
  #Returns:
  #	A GCT object without proteins that have any missing values
  
  x <- subset_gct(x,rid = which(rowSums(is.na(x@mat)) <= pct*ncol(x@mat))) 
  return(x)
}

filter_gct <- function(x,column, value){
  #Filter GCT object based on row metadata
  #Args:
  #	x: a GCT object
  #	column: name of metadata column to use for filtering
  #	value: The value that should be matched for the filtering
  #Returns:
  #	A GCT object with species filtered
  x <- subset_gct(x,rid = which(x@rdesc[,column] == value)) 
  return(x)
}


###2. Format and manipulate data-------------------------------------------

make_results <- function(x,contrast.matrix,fit.eb,verbose=F,logFC_se = F,storey_q = T){
  #Calculate q values and generates a result table
  #Args:
  #	x: a GCT object
  #	contrast.matrix: contrast matrix used to calculate contrasts
  #	fit.eb: The fitted linear model with empirical bayes calculation for adj. pvalue
  #Returns:
  #	A dataframe with the original expression values, pvalues, and qvalues
  res <- as.data.frame(x@mat) %>% rownames_to_column(var="id")
  
  for(i in colnames(contrast.matrix)){
    
    stats <- topTable(fit.eb, number = nrow(x@mat), sort.by = "none", coef = i) %>%
      select(logFC, P.Value, adj.P.Val,t) %>%
      mutate_at(vars(), signif, 4)
    
    if(storey_q){
      stats$Q.Value <- qvalue(stats$P.Value, fdr.level=0.05, pi0.method="bootstrap")$qvalues
    } else {
      stats$Q.Value <- stats$adj.P.Val
    }
    
    
    if(logFC_se){
      stats <- stats %>%
        mutate(logFC_se = (sqrt(fit.eb$s2.post) * fit.eb$stdev.unscaled)[,i])
    }
    
    colnames(stats) <- paste0(i,".",colnames(stats))
    res <- cbind(res,stats)
    

    
    if(verbose){
      print(i)
      summary(qvalue(fit.eb$p.value[,i]))
    }
    
  }
  
  return(res)
}

make_nested_results <- function(x, contrast.matrix, fit.eb,verbose=F,logFC_se = F,storey_q=T){
  #Calculate q values and generates a result table in nested format
  #Args:
  #	x: a GCT object
  #	contrast.matrix: contrast matrix used to calculate contrasts
  #	fit.eb: The fitted linear model with empirical bayes calculation for adj. pvalue
  #Returns:
  #	A nested dataframe for each contrast with the original expression values, 
  #       pvalues, and qvalues
  
  res <- as.data.frame(x@mat) %>% rownames_to_column(var="id")
  
  #Holder for list of data to be place in the nested dataframe
  data <- list()
  
  for(i in colnames(contrast.matrix)){
    
    stats <- topTable(fit.eb, number = nrow(res), sort.by = "none", coef = i) %>%
      select(logFC, P.Value, adj.P.Val,t) %>%
      mutate_at(vars(), signif, 4)
    
    if(storey_q){
      stats$Q.Value <- qvalue(stats$P.Value, fdr.level=0.05, pi0.method="bootstrap")$qvalues
    } else {
      stats$Q.Value <- stats$adj.P.Val
    }
    
    
    
    
    if(logFC_se){
      stats <- stats %>%
        mutate(logFC_se = (sqrt(fit.eb$s2.post) * fit.eb$stdev.unscaled)[,i])
    }
    
    data <- c(data,list(cbind(res,stats)))

    
    if(verbose){
      print(i)
      summary(qvalue(fit.eb$p.value[,i]))
    }
  }
  return(tibble(contrasts = colnames(contrast.matrix),
                data = data))
}

make_results_ftest <- function(x,fit.eb,coef=NULL,verbose=F,...){
  #Calculate q values and generates a result table
  #Args:
  #	x: a GCT object
  #	coef: lsit of coefficients to perform the F-test (or empty for all)
  #	fit.eb: The fitted linear model with empirical bayes calculation for adj. pvalue
  #Returns:
  #	A dataframe with the original expression values, pvalues, and qvalues
  res <- as.data.frame(x@mat) %>% rownames_to_column(var="id")
  
  stats <- topTable(fit.eb, number = nrow(x@mat), sort.by = "none", coef = coef,...) %>%
    mutate(Q.Value = qvalue(P.Value, fdr.level=0.05, pi0.method="bootstrap")$qvalues) %>%
    mutate_at(vars(), signif, 4)
  
  
  
  res <- cbind(res,stats)
  if(verbose){
    summary(qvalue(stats$P.Value))
  }
  
  
  
  return(res)
}

get_significant_summary <- function(results,qtreshold=0.05,logFC = 0){
  #Creates a table that summarizes the number of significantly regulated genes
  #Args:
  #	results: A nested result table as produced by make_nested_results
  #	qtreshold: q-value treshold to count a significant differential abundance
  output <- results %>% mutate(
    total = map_dbl(data,function(x){sum(!is.na(x$Q.Value<qtreshold))}),
    significant = map_dbl(data,function(x){sum(x$Q.Value<qtreshold & abs(x$logFC)>logFC,na.rm=T)}),
    up = map_dbl(data,function(x){sum(x$logFC>logFC & x$Q.Value<qtreshold,na.rm=T)}),
    down = map_dbl(data,function(x){sum(x$logFC< (-1*logFC) & x$Q.Value<qtreshold,na.rm=T)})
  ) %>% select(-data)
  return(output)
}


###3. Produce plots-----------------------------------------------------------------

pca_plots <- function(pca,cdesc,axes = c(1,2),title="PCA",addEllipses = T,
                      geom= "point", ...){
  #Generate pca plots colored with specific metadata
  #Args:
  #	pca: a pca object as generated by PCA
  #	cdesc: a data.frame indicating the features to color
  #	axes: the pca dimensions to plot
  #	title: the desired title for the plot
  #	...: Other arguments passed to fviz_pca
  #Returns:
  #	A list of pca plots
  output <- list()
  for(i in colnames(cdesc)){
    
    g <- fviz_pca_ind(pca,col.ind = cdesc[,i],
                      mean.point=FALSE,geom=geom,axes=axes,
                      addEllipses = addEllipses,...)+
      labs(title=title,subtitle = paste0("Labeled by ",i))+
      theme_bw()
    output <- c(output,list(g))
  }
  names(output) <- colnames(cdesc)
  return(output)
}

label_pca_plots <- function(pca,cdesc,axes = c(1,2),labels= FALSE,title="PCA",addEllipses = T, ...){
  #Generate pca plots colored with specific metadata
  #Args:
  #	pca: a pca object as generated by PCA
  #	cdesc: a data.frame indicating the features to color
  #	axes: the pca dimensions to plot
  #	title: the desired title for the plot
  #	...: Other arguments passed to fviz_pca
  #Returns:
  #	A list of pca plots
  output <- list()
  for(i in colnames(cdesc %>% select(sex,group,plex))){
    if (labels == TRUE) {
      g <- fviz_pca_ind(pca,col.ind = as.factor(cdesc[,i]),
                        mean.point=FALSE,geom="point",axes=axes,
                        addEllipses = addEllipses,...)+
        labs(title=title,subtitle = paste0("Labeled by ",i))+
        theme_bw()
      g1 <- g+geom_text_repel(aes(label = cdesc$bid))
    }else{
      g1 <- fviz_pca_ind(pca,col.ind = as.factor(cdesc[,i]),
                         mean.point=FALSE,geom="point",axes=axes,
                         addEllipses = addEllipses,...)+
        labs(title=title,subtitle = paste0("Labeled by ",i))+
        theme_bw()
    }
    output <- c(output,list(g1))
  }
  names(output) <- colnames(cdesc%>% select(sex,group,plex))
  return(output)
}

pca_correlations <- function(pca,cdesc,components = c(1:20)){
  #Calculates the intragroup correlation coefficient between each principal
  #component and the sample annotation
  #Args:
  #	pca: a pca object as generated by PCA
  #	cdesc: a data.frame indicating the features to color
  #	components: the components to be correlated. defaults to the first 10
  #
  
  #Obtain the principal component coordinates
  p <- data.frame(pca$ind$coord[,components])
  data <- data.frame(dims=character(),
                     correlation=numeric(),
                     experimental.factor = character() )
  
  for(i in colnames(cdesc)){
    #Create the class column to be used for correlation with the PC coordinates
    class = data.frame(matrix(rep(as.factor(cdesc[,i]),ncol(p)),ncol = ncol(p)))
    
    #Calculate intragroup correlations
    y <- map2(class,p,function(class,p){
      return(ICCest(class,p)$ICC)
    }) %>% unlist
    
    dimensions <- as.numeric(gsub("Dim\\.","",colnames(p)))
    
    data <- rbind(data,
                  data.frame(dims = dimensions,
                             correlation = y,
                             experimental.factor=rep(i,ncol(p))))
    
  }
  ggplot(data=data,aes(as.factor(dims),correlation,group=experimental.factor,color=experimental.factor)) +
    geom_line() + geom_point() + labs(x="Component (% Variance Explained)") + theme_bw()+
    scale_x_discrete(labels = paste0("PC",data$dims, " (",round(pca$eig[,2],1),"%)"))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

pca_variance_explained <- function(pca,cdesc,components=c(1:10)){
  #Calculates the percent of variance explained from sample metadata for each principal
  #component. The calculation is performed by fitting a linear model individually for 
  #each of the metadata variables
  #
  #Args:
  #	pca: a pca object as generated by PCA
  #	cdesc: a data.frame indicating the features to color
  #	components: the components to be correlated. defaults to the first 10
  #
  
  #Obtain the principal component coordinates
  p <- data.frame(pca$ind$coord[,components])
  
  #Intitialize the result data frame
  data <- data.frame(dims=character(),
                     pct.exp=numeric(),
                     experimental.factor = character() )
  
  #Loop through all component-metadata combinations 
  for(i in colnames(cdesc)){
    
    for(j in colnames(p)){
      
      #Check if the current metadata vector is valid
      if( (sum(is.na(cdesc[,i])) == length(cdesc[,i])) |
          (length(levels(as.factor(cdesc[,i])))<2)       ){
        next
      }
      
      #Fit a linear model between the principal component and metadata variable
      fit <- lm(p[,j] ~ cdesc[,i])
      af <- anova(fit)
      afss <- af$"Sum Sq"
      
      dimensions <- as.numeric(gsub("Dim\\.","",j))
      
      data <- rbind(data,
                    data.frame(dims = dimensions,
                               pct.exp = afss[1]/sum(afss)*100,
                               experimental.factor=i))
    }
  }
  
  
  g <- ggplot(data=data,aes(as.factor(dims),pct.exp,group=experimental.factor,color=experimental.factor)) +
    geom_line() + geom_point() + labs(x="Component (% Total Variance Explained)",y="% Variance Explaines within Component") + theme_bw()+
    scale_x_discrete(labels = paste0("PC",data$dims, " (",round(pca$eig[,2],1),"%)"))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(g)
  
}

pca_variance_explained_heatmap <- function(pca,cdesc,components=c(1:10)){
  #Calculates the percent of variance explained from sample metadata for each principal
  #component. The calculation is performed by fitting a linear model individually for 
  #each of the metadata variables
  #
  #Args:
  #	pca: a pca object as generated by PCA
  #	cdesc: a data.frame indicating the features to color
  #	components: the components to be correlated. defaults to the first 10
  #
  
  #Obtain the principal component coordinates
  p <- data.frame(pca$ind$coord[,components])
  
  #Intitialize the result data frame
  data <- data.frame(dims=character(),
                     pct.exp=numeric(),
                     experimental.factor = character() )
  
  #Loop through all component-metadata combinations 
  for(i in colnames(cdesc)){
    
    for(j in colnames(p)){
      
      #Check if the current metadata vector is valid
      if( (sum(is.na(cdesc[,i])) == length(cdesc[,i])) |
          (length(levels(as.factor(cdesc[,i])))<2)       ){
        next
      }
      
      #Fit a linear model between the principal component and metadata variable
      fit <- lm(p[,j] ~ cdesc[,i])
      af <- anova(fit)
      afss <- af$"Sum Sq"
      
      dimensions <- as.numeric(gsub("Dim\\.","",j))
      
      data <- rbind(data,
                    data.frame(dims = dimensions,
                               pct.exp = afss[1]/sum(afss)*100,
                               experimental.factor=i))
    }
  }
  
  dat <- 
    data %>%
    mutate(dims = paste0("dim",dims)) %>%
    pivot_wider(names_from = dims,values_from = pct.exp) %>%
    column_to_rownames("experimental.factor") %>%
    as.matrix()
  
  g <- Heatmap(dat,cluster_columns = F,
               col = colorRamp2((0:100), viridis::viridis(101)))
  
  return(g)
  
}

plot_protein <- function(x,id,cdesc.group="treatment",cdesc.fill=NA,cdesc.color=NA){
  #Plot the values of a gct object for a specific proteins
  #Args:
  #	x: a GCT object
  #	id: id of the protein to be plotted
  #	cdesc.value: The cdesc column to be used for the grouping
  #	cdesc.color: The fill color to use for additioanl separation
  #	             If NA, no separation occurs. 
  if(class(x)[1] != "GCT"){
    data <- x %>% filter(id.x==id)
  }else{
    data <- subset_gct(x,rid = id) %>% melt_gct(x) #%>% filter(id.x==id) 		
  }
  
  g<- ggplot(data, aes_string(x=cdesc.group,y="value"))+labs(title=id)+theme_bw()
  if(!is.na(cdesc.fill) & !is.na(cdesc.color)){
    g<- g+geom_boxplot(aes_string(fill=cdesc.fill),outlier.shape=NA)+
      geom_point(position=position_dodge(width=0.75),
                 aes_string(group=cdesc.fill,color=cdesc.color))	
  } else if(!is.na(cdesc.fill) ){
    g<- g+geom_boxplot(aes_string(fill=cdesc.fill),outlier.shape=NA)+
      geom_point(position=position_dodge(width=0.75),aes_string(group=cdesc.fill))
  } else if(!is.na(cdesc.color) ){
    g<- g+geom_boxplot(aes_string(color=cdesc.color),outlier.shape=NA)+
      geom_point(position=position_dodge(width=0.75),aes_string(group=cdesc.color))
  }
  
  else{
    g <- g<- g+geom_boxplot(outlier.shape=NA)+
      geom_point(position=position_dodge(width=0.75))
  }
  return(g)
}

plot_protein_dea <- function(x,id){
  #Plot the values of a MoTrPAC DEA table for a specific protein
  #Args:
  #	x: training_dea table
  #	id: id of the protein to be plotted
  for (i in unique(x$tissue)){
    x.sub <- x %>% filter(feature_ID == id & tissue == i) 
    if(nrow(x.sub)>0){
      g <-  x.sub %>%
        ggplot(aes(x= comparison_group, y=logFC, color = sex,group=sex)) + 
        geom_point()+
        geom_line()+
        geom_errorbar(aes(ymin = logFC - logFC_se,
                          ymax = logFC + logFC_se),
                      width=0.2)+
        theme_bw()+
        scale_color_manual(values=sex_cols)+
        geom_hline(yintercept = 0)+
        labs(x = "Sex and time point", y = "log2 fold-change vs sedentary", title = id, subtitle = i) 
      plot(g)
    }
    
  }
  return()
}

plot_metab_dea <- function(x,id){
  #Plot the values of a MoTrPAC DEA table for a specific metabolite
  #Args:
  #	x: training_dea table
  #	id: id of the metaboliteto be plotted
  for (i in unique(x$tissue)){
    x.sub <- x %>% filter(feature_ID == id & tissue == i) 
    if(nrow(x.sub)>0){
      g <-  x.sub %>%
        ggplot(aes(x= comparison_group, y=effect_size, color = sex,group=sex)) + 
        geom_point()+
        geom_line()+
        geom_errorbar(aes(ymin = effect_size - effect_size_se,
                          ymax = effect_size + effect_size_se),
                      width=0.2)+
        theme_bw()+
        scale_color_manual(values=sex_cols)+
        geom_hline(yintercept = 0)+
        labs(x = "Sex and time point", y = "log2 fold-change vs sedentary", title = id, subtitle = i) 
      plot(g)
    }
    
  }
  return()
}

plot_volcano <- function(contrasts,results,qtreshold=0.05,xlims = NA, ylims = NA){
  #Generates a volcano plot using p-values and logFC ratios
  #Args:
  #	contrasts: Name of the contrast being plotted
  #	results: a results table with logFC, q values, and p values
  #	qtreshold: The q-value treshold to be used. Default is 0.05.
  results <- results[!is.na(results$Q.Value),]
  
  g <- ggplot(results,aes(logFC,-log10(P.Value),color = Q.Value < qtreshold)) + 
    geom_point(size=2,alpha=0.5) + 
    scale_color_manual(values= c("Grey","#0000C9"),#brewer.pal(n=3,name = "Set2")[1]),
                       labels=c("Non-significant","Significant"),
                       name = paste0("q value < ",qtreshold))+
    labs(x= paste0("Log2(",regmatches(contrasts,regexpr("^[[:alnum:]]+",contrasts))," over ",regmatches(contrasts,regexpr("[[:alnum:]]+$",contrasts)),")"), 
         y="-Log10(p-value)",
         title = paste0("Comparison ",contrasts),
         subtitle = paste0(sum(results$Q.Value<qtreshold)," significant out of ",nrow(results),
                           ". ",sum(results$logFC>0&results$Q.Value<qtreshold), " increasing and ",
                           sum(results$logFC<0&results$Q.Value<qtreshold)," decreasing"))+
    theme_bw()
  
  
  if(min(results$Q.Value) < qtreshold){
    p.tresh <- (results %>% mutate(min = abs(qtreshold-results$Q.Value)) %>% 
                  arrange(min))[[1,"P.Value"]]
    
    g<- g+	annotate("text",x=max(results$logFC),y=-log10(p.tresh),
                    label=paste0("q value < ",qtreshold),vjust = -1,hjust="inward")+
      geom_hline(yintercept = -log10(p.tresh))
  }
  
  if(!is.na(xlims[1])){
    g <- g+scale_x_continuous(limits=xlims)
  }
  
  if(!is.na(ylims[1])){
    g <- g+scale_y_continuous(limits=ylims)
  }
  
  
  return(g)
}

plot_volcano_label <- function(contrasts,results,qtreshold=0.05,xlims = NA, ylims = NA){
  #Generates a volcano plot using p-values and logFC ratios
  #Args:
  #	contrasts: Name of the contrast being plotted
  #	results: a results table with logFC, q values, and p values
  #	qtreshold: The q-value treshold to be used. Default is 0.05.
  results <- results[!is.na(results$Q.Value),]
  
  g <- ggplot(results,aes(logFC,-log10(P.Value),color = Q.Value < qtreshold,
                          protein = id)) + 
    geom_point(size=2,alpha=0.7) + 
    scale_color_manual(values= c("Grey",brewer.pal(n=3,name = "Set2")[1]),
                       labels=c("Non-significant","Significant"),
                       name = paste0("q value < ",qtreshold))+
    labs(x= paste0("Log2(",regmatches(contrasts,regexpr("^[[:alnum:]]+",contrasts))," over ",regmatches(contrasts,regexpr("[[:alnum:]]+$",contrasts)),")"), 
         y="-Log10(p-value)",
         title = paste0("Comparison ",contrasts),
         subtitle = paste0(sum(results$Q.Value<qtreshold)," significant out of ",nrow(results),
                           ". ",sum(results$logFC>0&results$Q.Value<qtreshold), " increasing and ",
                           sum(results$logFC<0&results$Q.Value<qtreshold)," decreasing"))+
    theme_bw()
  
  
  if(min(results$Q.Value) < qtreshold){
    p.tresh <- (results %>% mutate(min = abs(qtreshold-results$Q.Value)) %>% 
                  arrange(min))[[1,"P.Value"]]
    
    g<- g+	annotate("text",x=max(results$logFC),y=-log10(p.tresh),
                    label=paste0("q value < ",qtreshold),vjust = -1,hjust="inward")+
      geom_hline(yintercept = -log10(p.tresh))
  }
  
  if(!is.na(xlims[1])){
    g <- g+scale_x_continuous(limits=xlims)
  }
  
  if(!is.na(ylims[1])){
    g <- g+scale_y_continuous(limits=ylims)
  }
  
  
  return(toWebGL(g))
}

plot_volcano_color <- function(contrasts,results,color.list=NA,qtreshold=0.05,color.name="",
                               colors=c("Grey",brewer.pal(n=3,name = "Set2")[2]),
                               xlims = NA, ylims = NA,labels=character(),alpha.list = 0.7){
  #Generates a volcano plot using p-values and logFC ratios
  #Args:
  #	contrasts: Name of the contrast being plotted
  #	results: a results table with logFC, q values, and p values
  #	color.list: A vector with color categories to be plotted
  #	qtreshold: The q-value treshold to be used. Default is 0.05.
  to.keep <- !is.na(results$Q.Value)
  results <- results[to.keep,]
  p.tresh <- (results %>% mutate(min = abs(qtreshold-results$Q.Value)) %>% arrange(min))[[1,"P.Value"]]
  
  g <- ggplot(results,aes(logFC,-log10(P.Value),color = id %in% color.list)) + 
    geom_point(size=2,alpha=0.3) + 
    scale_color_manual(values=colors,labels=c("False","True"), name = color.name)+
    geom_hline(yintercept = -log10(p.tresh))+
    labs(x= paste0("Log2(",regmatches(contrasts,regexpr("^[[:alnum:]]+",contrasts)),
                   " over ",regmatches(contrasts
                                       ,regexpr("[[:alnum:]]+$",contrasts)),")"), 
         y="-Log10(p-value)",
         title = paste0("Comparison ",contrasts),
         subtitle = paste0(sum(results$Q.Value<qtreshold)," significant proteins out of ",nrow(results)))+
    theme_bw()+
    annotate("text",x=max(results$logFC),y=-log10(p.tresh),
             label=paste0("q value < ",qtreshold),vjust = -1,hjust="inward")	
  
  
  if(!is.na(xlims[1])){
    g <- g+scale_x_continuous(limits=xlims)
  }
  
  if(!is.na(ylims[1])){
    g <- g+scale_y_continuous(limits=ylims)
  }
  
  if(length(labels)>0){
    labels <- labels[to.keep]
    g<-g+geom_text_repel(aes(label= as.character(labels)),nudge_y = 1)
  }
  
  
  return(g)
}

plot_volcano_color2 <- function(contrasts,results,color.list="Gray",qtreshold=0.05,color.name="",
                                colors=c("Grey",brewer.pal(n=3,name = "Set2")[2]),
                                labels=character(), label.size = 5,
                                alpha.list = 0.7,
                                shape.list="", shape.name="",shapes=c(16,4,15,17),
                                xlims = NA, ylims = NA){
  #Generates a volcano plot using p-values and logFC ratios
  #Args:
  #	contrasts: Name of the contrast being plotted
  #	results: a results table with logFC, q values, and p values
  #	color.list: A vector with color categories to be plotted
  #	qtreshold: The q-value treshold to be used. Default is 0.05.
  #	labels: A character vector with equal rows as results 
  #results <- results[!is.na(results$adj.P.Val),]
  p.tresh <- (results %>% mutate(min = abs(qtreshold-results$Q.Value)) %>% arrange(min))[[1,"P.Value"]]
  
  g <- ggplot(results,aes(logFC,-log10(P.Value),color = color.list)) + 
    geom_point(size=2,aes(shape=shape.list,alpha = alpha.list)) + 
    scale_shape_manual(values=shapes,name=shape.name)+
    scale_color_manual(values=colors, name = color.name)+
    geom_hline(yintercept = -log10(p.tresh))+
    labs(x= paste0("Log2(",regmatches(contrasts,regexpr("^[[:alnum:]]+",contrasts)),
                   " over ",regmatches(contrasts
                                       ,regexpr("[[:alnum:]]+$",contrasts)),")"), 
         y="-Log10(p-value)",
         title = paste0("Comparison ",contrasts),
         subtitle = paste0(sum(results$Q.Value<qtreshold,na.rm=T)," significant proteins out of ",nrow(results)))+
    theme_bw()+
    guides(alpha=F)+
    annotate("text",x=max(results$logFC),y=-log10(p.tresh),
             label=paste0("q value < ",qtreshold),vjust = -1,hjust="inward")
  if(length(labels)>0){
    g<-g+geom_text_repel(aes(label= as.character(labels)),size=label.size)
  }
  if(!is.na(xlims[1])){
    g <- g+scale_x_continuous(limits=xlims)
  }
  
  if(!is.na(ylims[1])){
    g <- g+scale_y_continuous(limits=ylims)
  }
  return(g)
}

plot_time_clust <- function(results,cont, clust, k, y.limits = NULL,
                            alpha = 0.1, average.line = T){
  #Plot temporal time traces clustered according to an hclust object
  # and a set number of clusters (k)
  #Args:
  #	results: A results data frame
  #	cont: The column names for the contrasts to be plotted
  #	clust: A vector indicating clusters (as produced by cuttree or pheatmap)
  #	k: THe number of clusters to make
  
  results <- cbind(results,cluster = clust)
  
  
  line_plot <- function(cluster,data,cont,y.limits){
    data.l <- gather(data,"timepoint","value",cont)
    data.l$timepoint <- factor(data.l$timepoint,levels = cont)
    g <- ggplot(data.l,aes(x=timepoint,y=value,group=id))+geom_line(alpha=alpha)+
      theme_bw() +labs(title=paste0("Cluster ",cluster),
                       x="Time point",
                       y="log2(Ratio to control)",
                       subtitle= paste0(nrow(data)," proteins or peptides"))+
      scale_y_continuous(limits=y.limits)+
      geom_hline(aes(yintercept=0),color="red")
    
    if(average.line){
      data.mean <- data.l %>% group_by(timepoint) %>%
        summarise(med = median(value)) %>% cbind(data.frame(id=rep("m",nrow(.))))
      g <- g + geom_line(aes(x=timepoint,y=med,color="brickred"),
                         data=data.mean,color="blue",size=1)
    }
    
    return(g)
  }
  
  if(length(y.limits) != 2){
    #y.limits <- c( results %>% select(cont) %>% unlist() %>% min() ,
    #               results %>% select(cont) %>% unlist() %>% max() )
    y.limits <- c( results %>% select(cont) %>% unlist() %>% quantile(0.01) ,
                   results %>% select(cont) %>% unlist() %>% quantile(0.99) )
    
  }
  
  gl <- results %>% group_by(cluster) %>% nest %>%
    pmap(line_plot,cont,y.limits)
  
  
  #grid.arrange(grobs=gl,ncol=2)
  return(gl)
}

plot_upset_plot <- function(results,qtreshold=0.05,type="up",text.scale=2){
  #Creates a table that summarizes the number of significantly regulated genes
  #Args:
  #	results: A nested result table as produced by make_nested_results
  #	qtreshold: q-value treshold to count a significant differential abundance
  upset.data.up <- data.frame(id = results$data[[1]]$id)
  upset.data.down <- data.frame(id = results$data[[1]]$id)
  for(i in 1:nrow(results)){
    x <- results$data[[i]]
    
    upset.data.up <- cbind(upset.data.up,
                           as.numeric(x$logFC>0 & x$Q.Value<qtreshold))
    
    upset.data.down <- cbind(upset.data.down,
                             as.numeric(x$logFC<0 & x$Q.Value<qtreshold))
  }
  
  colnames(upset.data.up) <- c("id",results$contrasts)
  colnames(upset.data.down) <- c("id",results$contrasts)
  
  upset.data.up[which(is.na(upset.data.up),arr.ind=T)] <- 0
  upset.data.down[which(is.na(upset.data.down),arr.ind=T)] <- 0
  
  if(type=="up"){
    upset(upset.data.up,sets= rev(results$contrasts),keep.order=T,text.scale=text.scale)
  }else{
    
    upset(upset.data.down,sets= rev(results$contrasts),keep.order=T,text.scale=text.scale)
  }
}


###4. Gprofiler helper functions---------------------------------------------------------------------------

make_gem <- function(x){
  #Formats gprofiler:gost results into Generich Enrichment Map format
  #Args:
  #	x: The results table of the output of gprofiler2::gost function.
  
  output <-
    transmute(x,
              GO.ID = term_id,
              Description = term_name,
              p.Val = p_value,
              FDR = p_value,
              Phenotype = rep(1,nrow(output)),
              Genes = intersection
    )
  return(output)
}

find_goterm_intersect_refseq <- function(goterm, id.list,
                                         organism = "rnorvegicus",
                                         return.genes = T){
  #Returns the intersect of proteins in the goterm and the id.list
  #Args:
  #	goterm: A character value for the term to be searched
  #	id.list: A list of ids to be searched for
  #	organism: The name of the organism
  output  <- rbind(
    gconvert(goterm, organism = organism, target = "REFSEQ_PEPTIDE_PREDICTED_ACC"),
    gconvert(goterm, organism = organism, target = "REFSEQ_PEPTIDE_ACC")
  ) %>% 
    filter(target %in% id.list)
  
  if(return.genes){
    return(output$name)
  }else{
    return(output$target)
  }
}

gconvert_unique <- function(query,organism,target,...){
  #Helper function to perform ID mapping with gconvert
  #while preserving only the first entry. This is useful
  #when trying you want a vector of the same length
  #as your query
  #
  #Args: same arguments as gconvert
  #
  
  x <- gconvert(query,organism,target,...)
  x <- x %>% distinct(input,.keep_all=T)
  return(x)
  
}


###5. Google Cloud helper functions --------------------------------------------------------------------------------

dl_read_gcp <- function(gcp_path,sep='\t',tmpdir='/srv/tmp',GSUTIL_PATH='~/google-cloud-sdk/bin/gsutil'){
  system(sprintf('mkdir -p %s',tmpdir))
  # download
  new_path = sprintf('%s/%s',tmpdir,basename(gcp_path))
  if(!file.exists(new_path)){
    cmd = sprintf('%s cp %s %s', GSUTIL_PATH, gcp_path, tmpdir)
    system(cmd,ignore.stdout = T,ignore.stderr = T)
  }
  # read in the data
  dt <- fread(new_path,sep=sep,header=T)
  return(dt)
}


###6. MoTrPAC MAWG visualization functions -------------------------------------------------------------------------
plot_dea <- function(dea, id,flip_sex = F){
  if(!exists("sex_cols")){devtools::source_url("https://raw.githubusercontent.com/MoTrPAC/motrpac-mawg/master/pass1b-06/figures/motrpac_colors_abbr.R?token=ADAK7J5SXWUO4S43IQV467K72DXXM")
}
  dea <- dea %>% filter(prot_id == id)
  if(flip_sex) dea$sex <- factor(dea$sex,levels = c("Male","Female")) 
  results <- list()
  for(i in unique(dea$feature_id)){
    g<- dea %>% filter(feature_id == i) %>% 
      ggplot(aes(x = paste0(sex,"_",comparison_group), y = logfc, fill = sex))+
      geom_col(color="gray25")+
      geom_errorbar(aes(ymin=logfc-logfc_se,ymax=logfc+logfc_se,width=.2))+
      geom_hline(yintercept = 0,color="gray25")+
      scale_fill_manual(values=sex_cols)+
      facet_wrap(~sex,scales="free_x")+
      labs(x="",y="logFC",title=i)+
      theme_cowplot()+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
    results[[i]] <- g
    
  }
  return(results)
}




###7. Statistical analysis------------------------------------------------------

limma_res_extract_se<-function(limma_res,e_fit,
                               effect_col="logFC",t_col="t",colname=NULL){
  #Returns the standard error estimate
  #Args:
  #	limma_res: Limma results table as extracted by topTable
  #	e_fit: Empirical bayes fit 
  #	effect_col: The name of the column with the effect size
  # t_col: The name of the column with the T statistics
  # colname: Not used
  
  # method 1
  effects = limma_res[[effect_col]]
  ts = limma_res[[t_col]]
  ses1 = effects/ts
  if(is.null(colname)){
    return(ses1)
  }
  # method 2
  # ses2 = sqrt(e_fit$s2.post) * e_fit$stdev.unscaled
  # ts1 = effects/ses1
  # ts2 = effects/ses2[,colname]
  # if(max(abs(ts-ts1),na.rm = T) > 1e-10 || 
  #    max(abs(ts-ts2),na.rm = T) > 1e-10 || 
  #    max(abs(ts1-ts2),na.rm = T)>1e-10){
  #   print("Warning: t-stats from computed ses are incompatible with limma's output")
  # }
  # return(ses2[,colname])
}


pfizer_colors <- function(color = NA,secondary = F){
  if(!secondary){
    colors <- c(
      primary1 = "#0000C9",
      primary2 = "#0095FF"
    )
  } else{
    colors <- c(
      secondary1 = "#0DBDBA",
      secondary2 = "#67BB6E",
      secondary3 = "#9D73F7",
      secondary4 = "#D95776",
      secondary5 = "#F49C34",
      secondary6 = "#F8DF5A",
      secondary7 = "#F4DDBA"
    )
  }
  
  if(is.na(color)){
    return(unname(colors))
  } else{
    return(colors[[color]])
  }
}




#' Find kinase motifs within the protein by querying the scansite4.mit.edu proteinscan tool.
#'
#' @param id Protein identifier to be queried. This should match the database 
#' defined in data_source.The default is 'swissprot' which requires unirprot 
#' entry id (e.g., ZO3_MOUSE)
#' @param data_source Database used for the query. Defaults to 'swissprot' 
#' which queries the uniprot database.
#' @param motif_class Organism class used for the motif query. Currently available options are 'MAMMALIAN' and 'YEAST'.
#' @param motif_name Specific motifs to be scanned for. An empty string 
#' (default) will query all available motifs 
#' @param stringency Character specifying motif score stringency. Three levels are available [High|Medium|Low]. These would return motifs with scores in the top 0.002, 0.01, and 0.05 percentile respectively.
#'
#' @return A table showing motif name, motif name abbreviation, score 
#' percentile, site position, and site flanking sequence (+/- 7 amino acids)
#' @export
#'
#' @examples
motif_scan <- function(id, data_source = "swissprot", 
                       motif_class = "MAMMALIAN", motif_name = "",
                       stringency = "High"){
  query_string <- paste0(
    "http://scansite4.mit.edu/webservice/proteinscan/identifier=",
    id,
    "/dsshortname=",
    data_source,
    "/motifclass=",
    motif_class,
    "/motifshortnames=",
    motif_name,
    "/stringency=",
    stringency
  )
  
  request_timeout <- T
  
  tryCatch({
    r <- GET(query_string)
    request_timeout <- F
  }, error = function(e){
    print(e)
  })
  
  #Check for success of query
  if(status_code(r) != 200 | request_timeout){
    print(paste0("Query failed for entry ",id,"! Please double check spelling of input parameters.",
                 "HTML error = ", status_code(r)))
    return(data.frame())
  }
  
  silent_name_repair <- function(x){vctrs::vec_as_names(names=x, repair = "unique", quiet = T)}
  
  #Extract query content
  return_table <- r %>%
    content(encoding= "UTF-8") %>% 
    as_list
  
  #Remove information not related to the predicted site
  return_table$proteinScanResult <- 
    return_table$proteinScanResult[names(return_table$proteinScanResult) == "predictedSite"]
  
  #Transform to table format
  return_table <- return_table %>%  
    as_tibble %>%
    unnest_wider(proteinScanResult, names_repair = ~ vctrs::vec_as_names(..., quiet = TRUE)) %>%
    select(-starts_with(".")) %>%
    drop_na %>%
    mutate_all(unlist) %>%
    as_data_frame()
  
  return(return_table)
  
}


#' Make a contrast list using all possible nums and denoms pairs
#'
#' @param nums A character vector of variables to use as numerators
#' @param denoms A character vector of variables to use as denominator
#' @param contrast_table Contrast table to add values to
#'
#' @return A contrast table
#' @export
#'
#' @examples
define_contrasts <- function(nums,denoms,contrast_table = NULL){
  
  #Initialize contrast table if not one available
  if(is.null(contrast_table)){
    contrast_table <- data.frame(nums = c(),denoms =c())
  }
  
  #Add all combinations of nums and denoms
  contrast_table <- 
    rbind(
      contrast_table,
      data.frame(
        nums = rep(nums,length(denoms)),
        denoms = rep(denoms,each = length(nums))
      )
    ) %>%
    
    #Check for duplicates and same nums and denoms
    dplyr::distinct() %>%
    dplyr::filter(nums != denoms) %>%
    
    #Sort columns
    dplyr::arrange(nums, denoms)
  
  return(contrast_table)
}

###8. Sequence Processing------------------------------------------------------
extract_flanking_sequence <- function(prot_id,sites,fasta,n=7){
  #Extracts the flanking sequence of a site in a protein
  #Args:
  #	prot_id: Character - The protein name
  #	sites: Numeric - The sites to be extracted
  #	fasta: A Biostring MultipleAASequence object
  #	n: The number of amino acids to extract
  #
  
  if(prot_id %in% names(fasta)){
    motifs <- character(length(sites))
    for(i in 1:length(sites)){
      if(sites[i] < n+1){
        motifs[i] <- as.character(fasta[[prot_id]][1:(sites[i]+n)])
        #Total character length must be 15
        motifs[i] <- paste0(str_c(rep("_",(n*2+1)-str_length(motifs[i])),collapse=""),motifs[i])
      } else if(sites[i] > (length(fasta[[prot_id]])-n)){
        motifs[i] <- as.character(fasta[[prot_id]][(sites[i]-n):length(fasta[[prot_id]])])
        motifs[i] <- paste0(motifs[i],str_c(rep("_",(n*2+1)-str_length(motifs[i])),collapse=""))
      } else{
        motifs[i] <- as.character(fasta[[prot_id]][(sites[i]-n):(sites[i]+n)])
      }
    }
    return(motifs)
  } else{
    warning(paste0("Protein ",prot_id," not found in fasta file"))
    return(NULL)
  }
}









