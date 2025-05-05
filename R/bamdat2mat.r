#' Write BAM dat file values to a set of matrices and data frames
#'
#' The script converts objects defined in the DATA SECTION of BAM into a set of
#' data objects that are more easily read by humans.
#' @param CommonName Common name of species modeled in BAM files. Only used when accessing dat, tpl, and cxx character vectors named as e.g. dat_CommonName
#' @param bam Object returned from \code{bam2r}.
#' @param writeFiles File extension of files to write output data to. Current options include "csv".
#' Set to any other value to run the function without writing output.
#' @param writeFiles_dir Name of directory to write files to.
#' @param prompt_me Turn on/off prompts that ask for user input before deleting files
#' @param nm_out_exclude Names of output objects to exclude. By default, excludes age and length composition matrices converted to numbers of fish which are computed in this function. BAM uses comp matrices in proportions.
#' @keywords bam stock assessment fisheries
#' @export
#' @returns Silently returns a \code{list} of all data objects,
#' @examples
#' \dontrun{
#' # Read in any of the current BAM models
#' dm_AtMe <- bamdat2mat("AtlanticMenhaden")
#' dm_BlSB <- bamdat2mat("BlackSeaBass")
#' dm_BlTi <- bamdat2mat("BluelineTilefish")
#' dm_Cobi <- bamdat2mat("Cobia")
#' dm_GagG <- bamdat2mat("GagGrouper")
#' dm_GrTr <- bamdat2mat("GrayTriggerfish")
#' dm_GrAm <- bamdat2mat("GreaterAmberjack")
#' dm_ReGr <- bamdat2mat("RedGrouper")
#' dm_RePo <- bamdat2mat("RedPorgy")
#' dm_ReSn <- bamdat2mat("RedSnapper")
#' dm_SnGr <- bamdat2mat("SnowyGrouper")
#' dm_Tile <- bamdat2mat("Tilefish")
#' dm_VeSn <- bamdat2mat("VermilionSnapper")
#'
#' # Write data objects to csv
#' dm_AtMe <- bamdat2mat("AtlanticMenhaden",writeFiles="csv")
#'
#' # Or don't write them and just do stuff with the data objects
#' dm_BlSB <- bamdat2mat("BlackSeaBass",writeFiles = "")
#' matplot(rownames(dm_BlSB$obs_L),dm_BlSB$obs_L,type="l")
#'
#' agec <- dm_BlSB$obs_agec
#' par(mfrow=c(length(agec),1),mar=c(2,2,1,1),mgp=c(.8,.1,0),tck=-0.01)
#' lapply(seq_along(agec),function(i){
#' nmi <- names(agec)[i]
#' x <- agec[[i]]
#' image(as.numeric(rownames(x)),as.numeric(colnames(x)),x,xlab="year",ylab="age",main=nmi)
#' })
#' }

bamdat2mat <- function(CommonName=NULL,bam=NULL,
                      writeFiles="csv",
                      writeFiles_dir=paste("bamdat2mat",writeFiles,sep="_"),
                      prompt_me=TRUE,
                      nm_out_exclude=c("obs_agec_nfish","obs_lenc_nfish")
){
  if(!is.null(CommonName)){
    bam <- bam2r(CommonName)
    bam <- standardize_bam(
      dat_obj=bam$dat,
      tpl_obj=bam$tpl,
      cxx_obj=bam$cxx
    )
  }
  if(is.list(bam)){
    init <- bam$init
  }

  inm <- names(init)
  inm_return <- c() # Names of objects accounted for in returned object
  out_list_misc <- list() # Extra stuff to add to return from function
  lh_misc <- list()       # Extra life history stuff to add to out_list_misc

  # Data input
  # General
  styr <- init$styr
  endyr <- init$endyr
  yr <- styr:endyr

  agebins_types <- inm[grepl("^agebins",inm)]
  agebins_nbins <- unlist(lapply(init[agebins_types],length))
  agebins_agec <- NULL
  for(i in agebins_types){
    assign(i,as.character(as.numeric(init[[i]])))
  }
  # Check to see if any agebins_types have the same number of agebins
  if(length(agebins_types)>1){
    apairs <- combn(1:length(agebins_nbins),2)
    abmatches <- list()
    for(i in 1:dim(apairs)[2]){
      x1_i <- apairs[1,i]
      x2_i <- apairs[2,i]
      nm1_i <- agebins_nbins[x1_i]
      nm2_i <- agebins_nbins[x2_i]
      if(nm1_i==nm2_i){
        abmatches[[length(abmatches)+1]] <- c(names(nm1_i),names(nm2_i))
      }
    }

    if(length(abmatches)>0){
      # Don't report a warning if it's just agebins and agebins_agec that match
      if(!"agebins"%in%abmatches[[1]]&"agebins_agec"%in%abmatches[[1]]){
        warning(paste0("At least two agebins objects have the same number of ages. Check bam files to assure that obs_agec data are labeled with the correct ages."))
        print(agebins_nbins)
      }
    }
  }

  lenbins_types <- inm[grepl("^lenbins?.*[^width]$",inm)]
  lenbins_nbins <- unlist(lapply(init[lenbins_types],length))
  for(i in lenbins_types){
    assign(i,init[[i]])
  }
  # Check to see if any lenbins_types have the same number of lenbins
  if(length(lenbins_types)>1){
    apairs <- combn(1:length(lenbins_nbins),2)
    abmatches <- list()
    for(i in 1:dim(apairs)[2]){
      x1_i <- apairs[1,i]
      x2_i <- apairs[2,i]
      nm1_i <- lenbins_nbins[x1_i]
      nm2_i <- lenbins_nbins[x2_i]
      if(nm1_i==nm2_i){
        abmatches[[length(abmatches)+1]] <- c(names(nm1_i),names(nm2_i))
      }
    }

    if(length(abmatches)>0){
        warning(paste0("At least two lenbins objects have the same number of lenbins. Check bam files to assure that obs_lenc data are labeled with the correct lenbins."))
        print(lenbins_nbins)
    }
  }

  ## Landings
  nm_obs_L <- inm[grepl("^obs_L",inm)]

  # Assemble obs_L and obs_cv_L data
  abb_L <- gsub("obs_L_","",nm_obs_L) # Fleet abbreviations
  obs_L <- matrix(NA,nrow=length(yr),ncol=length(abb_L),dimnames=list("year"=yr,"fleet"=abb_L))
  obs_cv_L <- obs_L
  for(i in abb_L){
    nm_obs_L_i <- paste0("obs_L_",i)
    nm_obs_cv_L_i <- paste0("obs_cv_L_",i)
    nm_styr_L_i <- paste0("styr_L_",i)
    nm_endyr_L_i <- paste0("endyr_L_",i)

    obs_L_i <- as.numeric(init[[nm_obs_L_i]])
    obs_cv_L_i <- as.numeric(init[[nm_obs_cv_L_i]])
    styr_L_i <- init[[nm_styr_L_i]]
    endyr_L_i <- init[[paste0("endyr_L_",i)]]
    yr_L_i <- styr_L_i:endyr_L_i

    obs_L[paste(yr_L_i),i] <- obs_L_i
    obs_cv_L[paste(yr_L_i),i] <- obs_cv_L_i

    inm_return <- c(inm_return,nm_obs_L_i,nm_obs_cv_L_i,nm_styr_L_i,nm_endyr_L_i)
  }

  ## Discards
  nm_obs_released <- inm[grepl("^obs_released",inm)]
  # Assemble obs_released and obs_cv_D data
  abb_released <- gsub("obs_released_","",nm_obs_released) # Fleet abbreviations
  obs_released <- matrix(NA,nrow=length(yr),ncol=length(abb_released),dimnames=list("year"=yr,"fleet"=abb_released))
  obs_cv_D <- obs_released
  for(i in abb_released){
    nm_obs_released_i <- paste0("obs_released_",i)
    nm_obs_cv_D_i <- paste0("obs_cv_D_",i)
    nm_styr_D_i <- paste0("styr_D_",i)
    nm_endyr_D_i <- paste0("endyr_D_",i)

    obs_released_i <- as.numeric(init[[nm_obs_released_i]])
    obs_cv_D_i <- as.numeric(init[[nm_obs_cv_D_i]])
    styr_D_i <- init[[nm_styr_D_i]]
    endyr_D_i <- init[[nm_endyr_D_i]]
    yr_D_i <- styr_D_i:endyr_D_i
    if(length(obs_released_i)>0){
      obs_released[paste(yr_D_i),i] <- obs_released_i
    }
    if(length(obs_cv_D_i)>0){
      obs_cv_D[paste(yr_D_i),i] <- obs_cv_D_i
    }
    inm_return <- c(inm_return,nm_obs_released_i,nm_obs_cv_D_i,nm_styr_D_i,nm_endyr_D_i)
  }

  ## cpue
  nm_obs_cpue <- inm[grepl("^obs_cpue",inm)]
  # Assemble obs_cpue and obs_cv_cpue data
  abb_cpue <- gsub("obs_cpue_","",nm_obs_cpue) # Fleet abbreviations
  obs_cpue <- matrix(NA,nrow=length(yr),ncol=length(abb_cpue),dimnames=list("year"=yr,"fleet"=abb_cpue))
  obs_cv_cpue <- obs_cpue
  for(i in abb_cpue){
    nm_obs_cpue_i <- paste0("obs_cpue_",i)
    nm_obs_cv_cpue_i <- paste0("obs_cv_cpue_",i)
    nm_styr_cpue_i <- paste0("styr_cpue_",i)
    nm_endyr_cpue_i <- paste0("endyr_cpue_",i)
    nm_yrs_cpue_i <- paste0("yrs_cpue_",i)

    obs_cpue_i <- as.numeric(init[[paste0("obs_cpue_",i)]])
    obs_cv_cpue_i <- as.numeric(init[[paste0("obs_cv_cpue_",i)]])
    styr_cpue_i <- init[[paste0("styr_cpue_",i)]]
    endyr_cpue_i <- init[[paste0("endyr_cpue_",i)]]
    if(is.null(init[[nm_yrs_cpue_i]])){
      yrs_cpue_i <- styr_cpue_i:endyr_cpue_i
    }else{
      yrs_cpue_i <- init[[nm_yrs_cpue_i]]
    }

    obs_cpue[paste(yrs_cpue_i),i] <- obs_cpue_i
    obs_cv_cpue[paste(yrs_cpue_i),i] <- obs_cv_cpue_i

    inm_return <- c(inm_return,nm_obs_cpue_i,nm_obs_cv_cpue_i,nm_styr_cpue_i,nm_endyr_cpue_i,nm_yrs_cpue_i)
    rm(obs_cpue_i,obs_cv_cpue_i,styr_cpue_i,endyr_cpue_i,yrs_cpue_i,
       nm_obs_cpue_i,nm_obs_cv_cpue_i,nm_styr_cpue_i,nm_endyr_cpue_i,nm_yrs_cpue_i)
  }

  ## lenc
  nm_obs_lenc <- inm[grepl("^obs_lenc",inm)]
  # Assemble obs_lenc data
  abb_lenc <- gsub("obs_lenc_","",nm_obs_lenc) # Fleet abbreviations
  obs_lenc <- obs_lenc_nfish <- list()
  nfish_lenc <- matrix(NA,nrow=length(yr),ncol=length(abb_lenc),dimnames=list("year"=yr,"fleet"=abb_lenc))
  nsamp_lenc <- nfish_lenc
  for(i in abb_lenc){
    nm_obs_lenc_i <- paste0("obs_lenc_",i)
    nm_yrs_lenc_i <- paste0("yrs_lenc_",i)
    nm_nyr_lenc_i <- paste0("nyr_lenc_",i)
    nm_nfish_lenc_i <- paste0("nfish_lenc_",i)
    nm_nsamp_lenc_i <- paste0("nsamp_lenc_",i)

    obs_lenc_matrix_i <- local({
      a_i <- init[[nm_obs_lenc_i]]
      b_i <- apply(a_i,2,as.numeric)
      matrix(b_i,nrow=nrow(a_i))
    })
    yrs_lenc_i <-   init[[nm_yrs_lenc_i]]
    nfish_lenc_i <- init[[nm_nfish_lenc_i]]
    nsamp_lenc_i <- init[[nm_nsamp_lenc_i]]
    lenbins_nbins_i <- dim(obs_lenc_matrix_i)[2]

    # Two step approach to trying to identify lenbins
    if(!is.null(lenbins)&length(lenbins)==lenbins_nbins_i){ # Use lenbins first if it's the right number of lenbins..
      lenbins_i <- lenbins
    }else{ #.. otherwise use whatever option has the same number of lenbins.
        lenbins_i <- get(names(lenbins_nbins[match(lenbins_nbins_i,lenbins_nbins)]))
      }
    lenbins_i <- as.character(as.numeric(lenbins_i)) # Should remove ".0" if length bins are entered that way (e.g. "100.0, 105.0" etc.)

    obs_lenc_i <- matrix(NA,nrow=length(yr),ncol=length(lenbins_i),dimnames=list("year"=yr,"len"=lenbins_i))
    obs_lenc_i[yrs_lenc_i,lenbins_i] <- obs_lenc_matrix_i
    obs_lenc[[i]] <- obs_lenc_i
    nfish_lenc[paste(yrs_lenc_i),i] <- as.numeric(nfish_lenc_i)
    nsamp_lenc[paste(yrs_lenc_i),i] <- as.numeric(nsamp_lenc_i)

    # Convert lenc to numbers
    obs_lenc_nfish_i <- local({
      a_i <- round(obs_lenc_i*nfish_lenc[,i])
      a_i[a_i==0] <- 1 # Don't let nfish be less than 1
      a_i
    })
    obs_lenc_nfish[[i]] <- obs_lenc_nfish_i

    inm_return <- c(inm_return,nm_obs_lenc_i,nm_yrs_lenc_i,nm_nyr_lenc_i,nm_nfish_lenc_i,nm_nsamp_lenc_i)
    rm(obs_lenc_matrix_i,obs_lenc_i,obs_lenc_nfish_i,yrs_lenc_i,lenbins_i,nfish_lenc_i,nsamp_lenc_i,
       nm_obs_lenc_i,nm_yrs_lenc_i,nm_nyr_lenc_i,nm_nfish_lenc_i,nm_nsamp_lenc_i)
  }

  ## agec
  nm_obs_agec <- inm[grepl("^obs_agec",inm)]
  # Assemble obs_agec data
  abb_agec <- gsub("obs_agec_","",nm_obs_agec) # Fleet abbreviations
  obs_agec <- obs_agec_nfish <- list()
  nfish_agec <- matrix(NA,nrow=length(yr),ncol=length(abb_agec),dimnames=list("year"=yr,"fleet"=abb_agec))
  nsamp_agec <- nfish_agec
  for(i in abb_agec){
    nm_obs_agec_i <- paste0("obs_agec_",i)
    nm_yrs_agec_i <- paste0("yrs_agec_",i)
    nm_nyr_agec_i <- paste0("nyr_agec_",i)
    nm_nfish_agec_i <- paste0("nfish_agec_",i)
    nm_nsamp_agec_i <- paste0("nsamp_agec_",i)

    obs_agec_matrix_i <- local({
      a_i <- init[[nm_obs_agec_i]]
      b_i <- apply(a_i,2,as.numeric)
      matrix(b_i,nrow=nrow(a_i))
    })
    yrs_agec_i <- init[[nm_yrs_agec_i]]
    nfish_agec_i <- init[[nm_nfish_agec_i]]
    nsamp_agec_i <- init[[nm_nsamp_agec_i]]
    agebins_nbins_i <- dim(obs_agec_matrix_i)[2]

    # Two step approach to trying to identify agebins
    if(!is.null(agebins_agec)&length(agebins_agec)==agebins_nbins_i){ # Use agebins_agec first if it's the right number of agebins..
      agebins_i <- agebins_agec
    }else{ #.. or agebins  if it's the right number of agebins..
      if(!is.null(agebins)&length(agebins)==agebins_nbins_i){
        agebins_i <- agebins
        }else{ #.. otherwise see use whatever option has the same number of agebins.
        agebins_i <- get(names(agebins_nbins[match(agebins_nbins_i,agebins_nbins)]))
      }

    }

    obs_agec_i <- matrix(NA,nrow=length(yr),ncol=length(agebins_i),dimnames=list("year"=yr,"age"=agebins_i))
    obs_agec_i[yrs_agec_i,agebins_i] <- obs_agec_matrix_i
    obs_agec[[i]] <- obs_agec_i
    nfish_agec[paste(yrs_agec_i),i] <- as.numeric(nfish_agec_i)
    nsamp_agec[paste(yrs_agec_i),i] <- as.numeric(nsamp_agec_i)

    # Convert agec to numbers
    obs_agec_nfish_i <- local({
      a_i <- round(obs_agec_i*nfish_agec[,i])
      a_i[a_i==0] <- 1 # Don't let nfish be less than 1
      a_i
    })
    obs_agec_nfish[[i]] <- obs_agec_nfish_i

    inm_return <- c(inm_return,nm_obs_agec_i,nm_yrs_agec_i,nm_nyr_agec_i,nm_nfish_agec_i,nm_nsamp_agec_i)
    rm(obs_agec_matrix_i,obs_agec_i,obs_agec_nfish_i,yrs_agec_i,agebins_i,nfish_agec_i,nsamp_agec_i,
       nm_obs_agec_i,nm_yrs_agec_i,nm_nyr_agec_i,nm_nfish_agec_i,nm_nsamp_agec_i)
  }

  # set parameters (not devs)
  nm_set <- inm[grepl("^set",inm,perl=TRUE)]
  init_set <- init[nm_set]
  init_set_length <- unlist(lapply(init_set,length))

  nm_set_len1 <- names(init_set_length)[init_set_length==1]
  nm_set_len3 <- names(init_set_length)[init_set_length==3]
  nm_set_len7 <- names(init_set_length)[init_set_length==7]
  nm_set_len7 <- nm_set_len7[grepl("^(?!.*dev).*",nm_set_len7,perl=TRUE)] # Remove any dev vectors
  nm_set_yr <- nm_set[grepl("^set_log_dev_vals_(F|rec)",nm_set)] # Set vectors with a year dimension
  nm_set_age <- c("set_log_dev_vals_Nage","set_M") # Set vectors with an age dimension

  nm_set_len1_Dmort <- sort(nm_set_len1[grepl("^set_Dmort_",nm_set_len1)])
  nm_set_len1_q <- sort(nm_set_len1[grepl("^set_.*q_.*",nm_set_len1)])
  nm_set_len1_w <- sort(nm_set_len1[grepl("^set_w_",nm_set_len1)])
  nm_set_len1_misc <- nm_set_len1[!nm_set_len1%in%c(nm_set_len1_Dmort,nm_set_len1_q,nm_set_len1_w)]

  set_len1_Dmort <- unlist(lapply(init[nm_set_len1_Dmort],as.numeric))
  if(!is.null(set_len1_Dmort)){
  set_len1_Dmort <- matrix(set_len1_Dmort,ncol=1,dimnames=list(names(set_len1_Dmort),"value"))
  }

  set_len1_q <- unlist(lapply(init[nm_set_len1_q],as.numeric))
  set_len1_q <- matrix(set_len1_q,ncol=1,dimnames=list(names(set_len1_q),"value"))

  set_len1_w <- unlist(lapply(init[nm_set_len1_w],as.numeric))
  set_len1_w <- matrix(set_len1_w,ncol=1,dimnames=list(names(set_len1_w),"value"))

  set_len1_misc <- unlist(lapply(init[nm_set_len1_misc],as.numeric))
  set_len1_misc <- matrix(set_len1_misc,ncol=1,dimnames=list(names(set_len1_misc),"value"  ))

  set_len3 <- do.call(rbind,lapply(init[nm_set_len3],as.numeric))
  set_len7 <- do.call(rbind,lapply(init[nm_set_len7],as.numeric))
  colnames(set_len3) <- c("lower","upper","phase")
  colnames(set_len7) <- c("init","lower","upper","phase","prior_mean","prior_var","prior_pdf")

  set_age <- cbind(c(NA,as.numeric(init[["set_log_dev_vals_Nage"]])),as.numeric(init[["set_M"]]))
  dimnames(set_age) <- list("agebins"=agebins,nm_set_age)

  # Possible names of life history vectors
  nm_lh_age <- c("obs_maturity_f","obs_maturity_m","obs_prop_f","obs_prop_m","fecpar_batches","set_fecpar_batches",
                 "fecundity","wgt_spawn","wgt_middle")
  nm_lh_age_avail <- nm_lh_age[nm_lh_age%in%inm]
  lh_age <- local({
    out <- matrix(NA,ncol=length(nm_lh_age_avail),nrow=length(agebins),dimnames = list("agebins"=agebins,nm_lh_age_avail))
    for(i in nm_lh_age_avail){
      obj_i <- init[[i]]
      if(is.matrix(obj_i)){
        lh_misc[[i]] <- apply(obj_i,2,as.numeric,simplify=TRUE)
        out_i <- as.numeric(obj_i[nrow(obj_i),])
        message(paste(i,"is a matrix. lh_age will contain the last row. The full matrix is added to end of the output."))
      }else{
        out_i <- as.numeric(obj_i)
      }
      out[,i] <- out_i
    }
    list(lh_age=out,lh_misc=lh_misc)
  })

  nm_set_log_F_dev_vals_L <- nm_set_yr[grepl("^set_log_dev_vals_F_L",nm_set_yr)]
  set_yr <- local({
    out <- matrix(NA,ncol=length(nm_set_yr),nrow=length(yr),dimnames = list("year"=yr,nm_set_yr))
    # L

    L_abb <- gsub("^set_log_dev_vals_F_L_","",nm_set_log_F_dev_vals_L)
    obs_L_yr <- apply(obs_L,2,function(x){names(which(!is.na(x)))},simplify = FALSE)
    for(abb_i in L_abb){
      dev_i <- paste0("set_log_dev_vals_F_L_",abb_i)
      out[obs_L_yr[[abb_i]],dev_i] <- as.numeric(init[[dev_i]])
    }

    # released (live discards)
    nm_set_log_F_dev_vals_D <- nm_set_yr[grepl("^set_log_dev_vals_F_D",nm_set_yr)]
    D_abb <- gsub("^set_log_dev_vals_F_D_","",nm_set_log_F_dev_vals_D)
    obs_released_yr <- apply(obs_released,2,function(x){names(which(!is.na(x)))},simplify=FALSE)
    for(abb_i in D_abb){
      dev_i <- paste0("set_log_dev_vals_F_D_",abb_i)
      out[obs_released_yr[[abb_i]],dev_i] <- as.numeric(init[[dev_i]])
    }

    # deviations of recruitment
    rec_yr <- paste(as.numeric(init[["styr_rec_dev"]])+0:(length(init[["set_log_dev_vals_rec"]])-1))
    out[rec_yr,"set_log_dev_vals_rec"] <- as.numeric(init[["set_log_dev_vals_rec"]])

    out
  })

  # age_error
  age_error <- matrix(as.numeric(init$age_error),nrow=length(agebins),ncol=length(agebins),dimnames=list("agebins"=agebins,"agebins"=agebins))


  inm_return <- c(inm_return,
                  nm_set_len1_Dmort, nm_set_len1_q, nm_set_len1_w, nm_set_len1_misc,
                  nm_set_len3, nm_set_len7, nm_set_age, nm_set_yr,
                  nm_lh_age_avail,
                  "age_error")
  inm_misc <- inm[!inm%in%inm_return]

  inm_misc_length <- unlist(lapply(init[inm_misc],length))

  nm_misc_1 <- names(inm_misc_length)[inm_misc_length==1]
  nm_misc_2 <- inm_misc[!inm_misc%in%nm_misc_1]

  # Bag everything else up in miscellaneous categories
  misc_1 <- unlist(lapply(init[nm_misc_1],as.numeric))
  misc_1 <- matrix(misc_1,ncol=1,dimnames=list(names(misc_1),"value"))
  misc_2 <- lapply(init[nm_misc_2],function(x){
    if(is.matrix(x)){
      apply(x,2,as.numeric)
    }else{
      as.numeric(x)
    }
  }
  )

  # Identify miscellaneous vectors with an age dimension
  nm_misc_age <- local({
    a <- lapply(1:length(misc_2),function(x){
    x_obj <- misc_2[[x]]
    x_nm <- nm_misc_2[x]
    if(length(x_obj)==length(agebins) & grepl("agebins|sel",x_nm)){
      TRUE
    }else{
      FALSE
    }
  })
    nm_misc_2[unlist(a)]
  })

  misc_age <- do.call(cbind,misc_2[nm_misc_age])

  # Modify misc_2
    nm_misc_2 <- nm_misc_2[!nm_misc_2%in%nm_misc_age]
    misc_2 <- misc_2[nm_misc_2]

  # Build by_age
    by_age <- cbind(lh_age$lh_age,set_age,misc_age)

  # Identify any two-dimensional objects in misc_2
    nm_misc_2_mat <- names(misc_2)[unlist(lapply(misc_2,is.matrix))] # Identify matrices
    nm_misc_2_vec <- names(misc_2)[!names(misc_2)%in%nm_misc_2_mat]  # Identify vectors

  # Separate misc_2 into vectors and matrices
    misc_2_mat <- misc_2[nm_misc_2_mat]
    misc_2_vec <- misc_2[nm_misc_2_vec]

    # Convert misc_2_vec to matrix
    misc_2_vec_mat <- local({
      nr <- max(unlist(lapply(misc_2_vec,length)))
      m <- matrix(NA,nrow=nr,ncol=length(misc_2_vec),dimnames=list(1:nr,names(misc_2_vec)))
      for(nmi in names(misc_2_vec)){
        i <- misc_2_vec[[nmi]]
        m[1:length(i),nmi] <- misc_2_vec[[nmi]]
      }
      m
    })

    out_list_misc <- c(out_list_misc,lh_age$lh_misc)

  # Output
  out <- c(list(
    init=init,
    points=misc_1,
    vectors=misc_2_vec_mat,
    obs_L=obs_L,
    obs_cv_L=obs_cv_L,
    obs_released=obs_released,
    obs_cv_D=obs_cv_D,
    obs_cpue=obs_cpue,
    obs_cv_cpue=obs_cv_cpue,
    obs_lenc=obs_lenc,
    obs_agec=obs_agec,
    obs_lenc_nfish=obs_lenc_nfish,
    obs_agec_nfish=obs_agec_nfish,
    nfish_lenc=nfish_lenc,
    nsamp_lenc=nsamp_lenc,
    nfish_agec=nfish_agec,
    nsamp_agec=nsamp_agec,
    set_Dmort=set_len1_Dmort,
    set_q=set_len1_q,
    set_w=set_len1_w,
    set_misc=set_len1_misc,
    set_par_point=set_len7,
    set_par_dev=set_len3,
    # set_age=set_age,
    # lh_age=lh_age,
    by_age=by_age,
    set_yr=set_yr,
    age_error=age_error
  ),
  out_list_misc)

  out <- out[!names(out)%in%nm_out_exclude]

  if(length(misc_2_mat)>0){
    txt1 <- ifelse(length(misc_2_mat)==1,"matrix","matrices")
    message(paste("dat contains non-typical",txt1,". The following",length(misc_2_mat),"objects will be added to the end of the output:\n",paste(names(misc_2_mat),collapse=", ")))

    out <- c(out,misc_2_mat)
    }

  # Write files
  if(writeFiles=="csv"){
       if(dir.exists(writeFiles_dir)){
        if(prompt_me){
          delete_files <- readline(prompt=paste0("Are you sure you want to delete all files in ",writeFiles_dir,"? (Enter TRUE or FALSE)"))
        }else{
          delete_files <- TRUE
        }
        if(delete_files){
          unlink(paste0(writeFiles_dir,"/*"))
          message(paste("The contents of ",writeFiles_dir,"have been deleted."))
        }else{
          message(paste("If you do not want to delete the contents of ",writeFiles_dir,", please specify a new value of writeFiles_dir and rerun bamdat2mat."))
        }
      }else{
        dir.create(writeFiles_dir,recursive=TRUE)
        message(paste("Created",writeFiles_dir))
      }

      # Write all output objects to csv files
      for(nmi in names(out)[names(out)!="init"]){
        i <- out[[nmi]]
        if(is.matrix(i)){
        write.csv(i,file=file.path(writeFiles_dir,paste0(nmi,".csv")))
        }
        if(is.list(i)){
          for(nmj in names(i)){
            j <- i[[nmj]]
            write.csv(j,file=file.path(writeFiles_dir,paste0(paste(nmi,nmj,sep="_"),".csv")))
          }
        }
      }
    }

  # Return output
  invisible(out)

  }
