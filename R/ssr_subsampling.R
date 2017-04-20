require(adegenet)
require(poppr)
require(readr)
require(dplyr)
require(magrittr)

#' Subsampling SSR Data function
#'
#' This function subsample a portion of ssr data and calculate the mean and sd of the subsampled set
#' @param pct Percentage of data to subsample
#' @param times Number of run
#' @param outputname Prefix of output files
#' @return 4 csv files: 1) full_locustable: A locustable of the full dataset; 2) subsample_locustable: locustable of all subsample dataset; 3) subsample_list; 4)subsample_locusallelmeansd: allel mean and sd of all subsample dataset
#' Genind object of the whole dataset
#' list object of all subsample data
#' Matrix descriping meand and sd of the subsample dataset
#' @keywords subsample ssr
#' @export
#' @examples
#' subsample_countmean("ssr_input.csv",pct=0.2,times=20,outputname="speciesA")

subsample_countmean <- function(filepath,pct=0.1,times=20,outputname="sample"){
  
  ssrtable <- read_csv(filepath)
  colnames(ssrtable)[1] <- ""
  
  #list for holding the locus table data
  subsamplingls <- list()
  locustablels <- list()
  
  #Subsamling
  for(i in seq_along(1:times)){
    sub_ssrtable <- sample_frac(ssrtable,pct)
    
    #Concantanate each codominant locus in one column seperated by /
    col_number <- ncol(sub_ssrtable)
    target_index <- seq(2,col_number,2)
    codominant_index <- seq(3,col_number,2)
    
    start_index <- 2
    
    for(n in seq_along(target_index)){
      col_name <- paste0(colnames(sub_ssrtable)[start_index],"_c")
      sub_ssrtable[,col_name] <- apply(sub_ssrtable[,start_index:(start_index+1)],1,paste,collapse="/")
      start_index <- start_index + 2
    }
    
    sub_ssrtable_c <- sub_ssrtable[,(col_number+1):(col_number+length(target_index))]
    sub_ssrtable_c %<>% as.data.frame()
    rownames(sub_ssrtable_c) <- as.list(sub_ssrtable[,1]) %>% unlist()
    
    subsamplingls[[i]] <- sub_ssrtable_c
    
    
    #Convert to genind then genclone object
    ogenind <- df2genind(sub_ssrtable_c,sep="/")
    ogenclone <- as.genclone(ogenind)
    to_append_locustb <- locus_table(ogenclone) %>% unclass()
    
    #locus table
    locustablels[[i]] <- to_append_locustb
  }
  
  binded_locustable <- do.call(rbind,locustablels)
  binded_locustable <- cbind(binded_locustable,locus=rownames(binded_locustable))
  
  binded_locustable %<>% tbl_df()
  binded_locustable <- mutate(binded_locustable,allele=as.numeric(allele))
  binded_locustable <- mutate(binded_locustable,`1-D`=as.numeric(`1-D`))
  binded_locustable <- mutate(binded_locustable,Hexp=as.numeric(Hexp))
  binded_locustable <- mutate(binded_locustable,Evenness=as.numeric(Evenness))
  
  rowperreplicate <- length(target_index) + 1
  binded_locustable$Subsample_group <- rep(1:times,each=rowperreplicate)

  write.csv(binded_locustable,paste0(outputname,"_subsample_locustable.csv"),row.names = FALSE)
  
  binded_subsamplingls <- do.call(rbind,subsamplingls)
  #nrow in each subsample
  subsample_row <- subsamplingls[1] %>% as.data.frame() %>% nrow()
  binded_subsamplingls <- cbind(binded_subsamplingls,Subsample_group=rep(1:times,each=subsample_row))

  assign("subsamplels",binded_subsamplingls,envir = .GlobalEnv)
  write.csv(binded_subsamplingls,paste0(outputname,"_subsample_list.csv"))
  
  
  #allele mean and SD for each locus (subsampled)
  locus_name_ls <- binded_locustable$locus[1:length(target_index)]  #list of locus
  mean_ls <- list()
  sd_ls <- list()
  for(name in locus_name_ls){
    subsample_locus_mean <- binded_locustable[binded_locustable$locus==name,"allele"] %>% colMeans()
    subsample_locus_sd <- binded_locustable[binded_locustable$locus==name,"allele"] %>% as.matrix() %>% as.list() %>% unlist() %>% sd()
    mean_ls[name] <- subsample_locus_mean
    sd_ls[name] <- subsample_locus_sd
  }
  
  subsample_locus_meansd_total <- cbind(mean_ls,sd_ls)
  assign("subsample_mean_sd",subsample_locus_meansd_total,envir = .GlobalEnv)
  
  write.csv(subsample_locus_meansd_total,paste0(outputname,"_subsample_locusallelemeansd.csv"))

  
  ###info_table for the whole dataset
  start_index <- 2
  
  for(n in seq_along(target_index)){
    col_name <- paste0(colnames(ssrtable)[start_index],"_c")
    ssrtable[,col_name] <- apply(ssrtable[,start_index:(start_index+1)],1,paste,collapse="/")
    start_index <- start_index + 2
  }
  
  ssrtable_c <- ssrtable[,(col_number+1):(col_number+length(target_index))]
  ssrtable_c %<>% as.data.frame()
  rownames(ssrtable_c) <- as.list(ssrtable[,1]) %>% unlist()
  
  #Convert to genind then genclone object
  ogenind <- df2genind(ssrtable_c,sep="/")
  assign("fullsample_genind",ogenind,envir = .GlobalEnv)
  ogenclone <- as.genclone(ogenind)
  full_locustable <- locus_table(ogenclone) %>% unclass()
  write.csv(full_locustable,paste0(outputname,"_full_locustable.csv"))
}







