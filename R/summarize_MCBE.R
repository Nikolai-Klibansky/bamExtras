#' Summarize results from MCBE uncertainty analysis
#'
#' @param dir_bam_sim Name of directory where run_MCBE results files are stored, relative to the working directory.
#' @param dir_bam_base Name of directory where base (i.e. reference) model results are stored, relative to the working directory.
#' @param fileName Name given to BAM base run files, not including file extensions
#' @param nm_model character string indicating the model name as part of the MCBE file names.
#' For example, "bam" is the \code{nm_model} in "1-bam.rdat". This optional argument should be specified if there
#' are .rdat files in \code{dir_bam_sim} other than the MCBE output files. Otherwise, this function
#' will try to read all .rdat files.
#' @param obj_collect Names of objects in MCBE rdat files to collect and summarize. Possible values include "info","parms","like","spr.brps","mse", "a.series","t.series","parm.cons", "N.age", "N.age.mdyr", "N.age.spawn","Z.age", and "sel.age"
#' @param obj_labels list of names for each value of each parameter in parm.cons
#' @param nm_Fref Name of F alternative F reference point to use for calculating Fref.65, Fref.75, and Fref.85. There is surely a better way to do those calculations, but this is the current method.
#' @param ncores number of cores to use for parallel processing
#' @details This function also runs the \code{bamExtras::check_bounds} function
#' on the rdat for each simulation and adds the logical value \code{nearbound} to
#' the \code{rdat$parms} indicating if any parameters were near the bounds.
#' @returns list of summarized outputs, most of which are named similar to typical BAM output (rdat) files.
#' This list can be passed to \code{plot_MCBE} to make a variety of useful plots
#' @keywords bam MCBE stock assessment fisheries
#' @export
#' @examples
#' \dontrun{
#' # Run MCBE, writing files to dir_bam_sim
#' run_MCBE("GrayTriggerfish", dir_bam_sim="GrTr_sim", dir_bam_base="GrTr_base")
#' # Summarize results of MCBE and assign to object
#' MCBE_GrTr <- summarize_MCBE(dir_bam_sim="sim_GrTr")
#' }



summarize_MCBE <-
  function(dir_bam_sim="sim",
                           dir_bam_base="base",
                           fileName="bam",
                           nm_model="",
                           obj_collect = c("info","parms","like","spr.brps","mse",
                                           "a.series","t.series","parm.cons",
                                           "N.age","Z.age","sel.age"),
                           obj_labels = list(parm.cons = c("init", "lower", "upper", "phase",
                                                           "prior_mean","prior_var", "prior_pdf",
                                                           "result")
                           ),
                           nm_Fref = c("F40","Fref"),
                           ncores = NULL
                     ){

#### parallel setup
if(is.null(ncores)){
coresAvail <- parallel::detectCores()
ncores <- coresAvail-1
}
ncores  <- min(c(ncores,coresAvail))
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

### Do stuff with base run
rdat <- dget(file.path(dir_bam_base,paste(fileName,'rdat', sep=".")))


# nm_sim <- sprintf(paste("%0",nchar(nsim),".0f",sep=""),1:nsim)
#
#
# ## Get base model stuff
#
# # Identify working directory
#   wd <- getwd()
#   message(paste("working directory:",wd))
#
simFiles <- list.files(dir_bam_sim)
pattern_rdat <- paste0(nm_model,".rdat$")
simFilesRdat <- simFiles[grepl(x=simFiles,pattern=pattern_rdat)]
modRunName <- local({
  a <- gsub(x=simFilesRdat,pattern=".rdat",replacement = "",fixed=TRUE)
  num.parts <- unique(unlist(lapply(strsplit(a,split="-"),length))) # number of parts to file names separated by hyphens
  b <- matrix(unlist(strsplit(a,split="-")),ncol=num.parts,byrow=TRUE)
  c <- b[,1]
  if(num.parts==3){ # Give longer names, when sets of runs exist (e.g. profiling)
    c <- paste(b[,1],b[,3])
  }
  return(c)})

a.series <- rdat$a.series
t.series <- rdat$t.series
a.series.names <- names(a.series)
t.series.names <- names(t.series)
parm.cons.names <- names(rdat$parm.cons)
ages <- a.series$age
years <- t.series$year

# Retrieve values from sim rdat files
simSummary <- foreach::foreach(iter.j=seq_along(simFilesRdat),
                      .packages="bamExtras",
                       .multicombine=TRUE
) %dopar% {
  # dget each rdat file and pull stuff out of it
  rdat.j.name <- simFilesRdat[iter.j]
  modRunName.j <- modRunName[iter.j]#unlist(strsplit(rdat.j.name, "\\-|\\."))[1]
  modRunGroup.j <- unlist(strsplit(modRunName.j,split=" "))[2]
  rdat.j <- dget(paste(dir_bam_sim,rdat.j.name,sep="/"))

  ## Calculate some values needed in the benchmarks table
  eq <- rdat.j$eq.series
  Fmsy <- ifelse("Fmsy"%in%names(rdat.j$parms),rdat.j$parms$Fmsy,NA)
  Fref <- if(any(nm_Fref%in%names(rdat.j$parms))){
    a <- nm_Fref[nm_Fref%in%names(rdat.j$parms)]
    if(length(a)>1){
      nm_Fref <- a[1]
      message(paste0(paste(a,collapse=", ")," found in parms. The first match (",out,") will be used."))
      a[1]
    }else{
      nm_Fref <- a
    }
    rdat.j$parms[[nm_Fref]]
  }else{
   NA
  }

  F.eq <- eq$F.eq
  # Revisit at a later time 2023-12-19 NPK
  # Revisiting 2024-08-20 NPK
  nm_Leq <- local({
    nm_eq <- names(eq)
    a <- nm_eq[grepl("^L\\.eq\\.[A-Za-z]*klb",nm_eq)]
    if(length(a)>1){
    out <- a[1]
      message(paste0("More than one variable in eq.series has a name that starts with L.eq and ends with klb. The first match (",out,") will be used."))
    a[1]
    }else{
    out <- a
    }
    out
  })
  nm_Deq <- local({
    nm_eq <- names(eq)
    a <- nm_eq[grepl("^D\\.eq\\.[A-Za-z]*knum",nm_eq)]
    if(length(a)>1){
      out <- a[1]
      message(paste0("More than one variable in eq.series has a name that starts with D.eq and ends with knum. The first match (",out,") will be used."))
    }else if(length(a)==0){
      out <- NA
      message(paste0("No variables in eq.series have a name that starts with D.eq. Are there discards in this model?"))
      }else{
      out <- a
    }
    out
  })

  L.eq <- eq[,nm_Leq]
  D.eq <- if(is.na(nm_Deq)){
    NA}else{
    eq[,nm_Deq]
  }


  if(!is.na(Fmsy)){
  Fmsy.65 <- 0.65*Fmsy
  Fmsy.75 <- 0.75*Fmsy
  Fmsy.85 <- 0.85*Fmsy

  Lmsy.65.klb  <- L.eq[which.min(abs(F.eq-Fmsy.65))]
  Lmsy.75.klb  <- L.eq[which.min(abs(F.eq-Fmsy.75))]
  Lmsy.85.klb  <- L.eq[which.min(abs(F.eq-Fmsy.85))]

  Dmsy.65.knum <- D.eq[which.min(abs(F.eq-Fmsy.65))]
  Dmsy.75.knum <- D.eq[which.min(abs(F.eq-Fmsy.75))]
  Dmsy.85.knum <- D.eq[which.min(abs(F.eq-Fmsy.85))]
  }else{
    Fmsy.65 <- Fmsy.75 <- Fmsy.85 <- NA
    Lmsy.65.klb <- Lmsy.75.klb <- Lmsy.85.klb <- NA
    Dmsy.65.knum <- Dmsy.75.knum <- Dmsy.85.knum <- NA
  }

  if(!is.na(Fref)){
    Fref.65 <- 0.65*Fref
    Fref.75 <- 0.75*Fref
    Fref.85 <- 0.85*Fref

    Lref.65.klb  <- L.eq[which.min(abs(F.eq-Fref.65))]
    Lref.75.klb  <- L.eq[which.min(abs(F.eq-Fref.75))]
    Lref.85.klb  <- L.eq[which.min(abs(F.eq-Fref.85))]
    Dref.65.knum <- D.eq[which.min(abs(F.eq-Fref.65))]
    Dref.75.knum <- D.eq[which.min(abs(F.eq-Fref.75))]
    Dref.85.knum <- D.eq[which.min(abs(F.eq-Fref.85))]
  }else{
    Fref.65 <- Fref.75 <- Fref.85 <- NA
    Lref.65.klb <- Lref.75.klb <- Lref.85.klb <- NA
    Dref.65.knum <- Dref.75.knum <- Dref.85.knum <- NA
  }

  rdat.j$parms$Lmsy.65.klb <- Lmsy.65.klb
  rdat.j$parms$Lmsy.75.klb <- Lmsy.75.klb
  rdat.j$parms$Lmsy.85.klb <- Lmsy.85.klb
  rdat.j$parms$Dmsy.65.knum <- Dmsy.65.knum
  rdat.j$parms$Dmsy.75.knum <- Dmsy.75.knum
  rdat.j$parms$Dmsy.85.knum <- Dmsy.85.knum

  rdat.j$parms$Lref.65.klb <- Lref.65.klb
  rdat.j$parms$Lref.75.klb <- Lref.75.klb
  rdat.j$parms$Lref.85.klb <- Lref.85.klb
  rdat.j$parms$Dref.65.knum <- Dref.65.knum
  rdat.j$parms$Dref.75.knum <- Dref.75.knum
  rdat.j$parms$Dref.85.knum <- Dref.85.knum

  # Compute mean squared errors associated with data fitted
  rdat.j$mse <- mse_calc(rdat.j)

  # Check for bounding issues
  rdat.j$parms$nearbound <-  any(check_bounds(rdat.j)$nearbound)

  # Create output list from each iteration of foreach
  # NOTE: This output will be combined with the .combine function
  out <- list("modRunName"=modRunName.j, "modRunGroup"=modRunGroup.j)
  out <- c(out,rdat.j[obj_collect])

  return(out)
}##### END FOREACH LOOP ######

# Redefine objects for collecting values from each sim run
names(simSummary) <- unlist(lapply(simSummary,FUN=function(x){x[["modRunName"]]}))

## Vectors

if("info"%in%obj_collect){
L.sim.info <- lapply(simSummary,FUN=function(x){x[["info"]]})    # Generate lists with potentially uneven levels
L.sim.info.itemnames <- unique(unlist(lapply(L.sim.info,names))) # Identify the set of item names among all levels in L.sim.info
# Fill in all item names in L.sim.info so that all levels have the same items in the same order
L.sim.info <- lapply(L.sim.info,FUN=function(x){
  missing.itemnames <- L.sim.info.itemnames[!L.sim.info.itemnames%in%names(x)]
  missing.items <- rep(NA,length(missing.itemnames))
  names(missing.items) <- missing.itemnames
  a <- c(x,missing.items)
  a <- a[L.sim.info.itemnames]
  a
})
L.sim.info <- lapply(1:length(L.sim.info),FUN=function(x){c("modRunName"=names(L.sim.info)[x],L.sim.info[[x]])})
D.sim.info <- do.call(rbind.data.frame, L.sim.info)
rownames(D.sim.info) <- D.sim.info$modRunName
out.info <- D.sim.info
}

if("parms"%in%obj_collect){
L.sim.parms <- lapply(simSummary,FUN=function(x){x[["parms"]]})
L.sim.parms.itemnames <- unique(unlist(lapply(L.sim.parms,names)))
L.sim.parms <- lapply(L.sim.parms,FUN=function(x){
  missing.itemnames <- L.sim.parms.itemnames[!L.sim.parms.itemnames%in%names(x)]
  missing.items <- rep(NA,length(missing.itemnames))
  names(missing.items) <- missing.itemnames
  a <- c(x,missing.items)
  a <- a[L.sim.parms.itemnames]
  a
})
L.sim.parms <- lapply(1:length(L.sim.parms),FUN=function(x){c("modRunName"=names(L.sim.parms)[x],L.sim.parms[[x]])})
D.sim.parms <- do.call(rbind.data.frame, L.sim.parms)
rownames(D.sim.parms) <- D.sim.parms$modRunName
out.parms <- D.sim.parms
}

if("like"%in%obj_collect){
L.sim.like <- lapply(simSummary,FUN=function(x){as.list(x[["like"]])})
L.sim.like.itemnames <- unique(unlist(lapply(L.sim.like,names)))
L.sim.like <- lapply(L.sim.like,FUN=function(x){
  missing.itemnames <- L.sim.like.itemnames[!L.sim.like.itemnames%in%names(x)]
  missing.items <- rep(NA,length(missing.itemnames))
  names(missing.items) <- missing.itemnames
  a <- c(x,missing.items)
  a <- a[L.sim.like.itemnames]
  a
})
L.sim.like <- lapply(1:length(L.sim.like),FUN=function(x){c("modRunName"=names(L.sim.like)[x],L.sim.like[[x]])})
D.sim.like <- do.call(rbind.data.frame, L.sim.like)
rownames(D.sim.like) <- D.sim.like$modRunName
out.like <- D.sim.like
}

if("spr.brps"%in%obj_collect){
L.sim.spr.brps <- lapply(simSummary,FUN=function(x){x[["spr.brps"]]})
L.sim.spr.brps.itemnames <- unique(unlist(lapply(L.sim.spr.brps,names)))
L.sim.spr.brps <- lapply(L.sim.spr.brps,FUN=function(x){
  missing.itemnames <- L.sim.spr.brps.itemnames[!L.sim.spr.brps.itemnames%in%names(x)]
  missing.items <- rep(NA,length(missing.itemnames))
  names(missing.items) <- missing.itemnames
  a <- c(x,missing.items)
  a <- a[L.sim.spr.brps.itemnames]
  a
})
L.sim.spr.brps <- lapply(1:length(L.sim.spr.brps),FUN=function(x){c("modRunName"=names(L.sim.spr.brps)[x],L.sim.spr.brps[[x]])})
D.sim.spr.brps <- do.call(rbind.data.frame, L.sim.spr.brps)
rownames(D.sim.spr.brps) <- D.sim.spr.brps$modRunName
out.spr.brps <- D.sim.spr.brps
}

if("mse"%in%obj_collect){
L.sim.mse <- lapply(simSummary,FUN=function(x){as.list(x[["mse"]])})
L.sim.mse.itemnames <- unique(unlist(lapply(L.sim.mse,names)))
L.sim.mse <- lapply(L.sim.mse,FUN=function(x){
  missing.itemnames <- L.sim.mse.itemnames[!L.sim.mse.itemnames%in%names(x)]
  missing.items <- rep(NA,length(missing.itemnames))
  names(missing.items) <- missing.itemnames
  a <- c(x,missing.items)
  a <- a[L.sim.mse.itemnames]
  a
})
L.sim.mse <- lapply(1:length(L.sim.mse),FUN=function(x){c("modRunName"=names(L.sim.mse)[x],L.sim.mse[[x]])})
D.sim.mse <- do.call(rbind.data.frame, L.sim.mse)
rownames(D.sim.mse) <- D.sim.mse$modRunName
out.mse <- D.sim.mse
}

## Matrices
if("a.series"%in%obj_collect){
# a.series
aser.vec <- setNames(rep(NA,length(ages)),ages)
L.sim.a.series <- foreach::foreach(i=seq_along(a.series.names),
                      .multicombine=TRUE
) %dopar% {
  a <- lapply(simSummary,FUN=function(x){
    age_sim <- paste(x[["a.series"]][["age"]])
    out <- aser.vec
    out[age_sim] <- x[["a.series"]][[a.series.names[i]]]
    out
    })
  b <- do.call(rbind.data.frame,a)
  dimnames(b) <- list(names(simSummary),ages) # Rename columns
  return(b)
}
names(L.sim.a.series) <- a.series.names
out.a.series <- L.sim.a.series
}

if("t.series"%in%obj_collect){
# t.series
tser.vec <- setNames(rep(NA,length(years)),years)
L.sim.t.series <- foreach::foreach(i=seq_along(t.series.names),
                          .multicombine=TRUE
) %dopar% {
  # Define
  #year_sim <- lapply(simSummary,FUN=function(x){x[["t.series"]][["year"]]}) # years in each sim run

  a <- lapply(simSummary,FUN=function(x){
    year_sim <- paste(x[["t.series"]][["year"]])
    out <- tser.vec
    out[year_sim] <- x[["t.series"]][[t.series.names[i]]]
    out
    })
  b <- do.call(rbind.data.frame,a)
  dimnames(b) <- list(names(simSummary),years) # Rename columns
  return(b)
}
names(L.sim.t.series) <- t.series.names
out.t.series <- L.sim.t.series
}

if("N.age"%in%obj_collect){
# N.age
N.age.dn1 <- dimnames(rdat$N.age)[[1]] # ages
N.age.dn2 <- dimnames(rdat$N.age)[[2]] # years
N.age.vec <- setNames(rep(NA,length(N.age.dn1)),N.age.dn1)
N.year.vec <- setNames(rep(NA,length(N.age.dn2)),N.age.dn2)

L.sim.N.age <- foreach::foreach(i=seq_along(N.age.dn1),
                          .multicombine=TRUE
) %dopar% {
  a <- lapply(simSummary,FUN=function(x){
    age_sim <- dimnames(x[["N.age"]])[[2]]
    out <- N.year.vec
    out[age_sim] <- x[["N.age"]][N.age.dn1[i],]
    out
    })
  b <- do.call(rbind.data.frame,a)
  dimnames(b) <- list(names(simSummary),N.age.dn2)
  return(b)
}
names(L.sim.N.age) <- N.age.dn1
out.N.age <- L.sim.N.age
}

if("N.age.spawn"%in%obj_collect){
  # N.age.spawn
  N.age.spawn.dn1 <- dimnames(rdat$N.age.spawn)[[1]] # ages
  N.age.spawn.dn2 <- dimnames(rdat$N.age.spawn)[[2]] # years
  N.age.spawn.vec <- setNames(rep(NA,length(N.age.spawn.dn1)),N.age.spawn.dn1)
  N.year.vec <- setNames(rep(NA,length(N.age.spawn.dn2)),N.age.spawn.dn2)

  L.sim.N.age.spawn <- foreach::foreach(i=seq_along(N.age.spawn.dn1),
                                  .multicombine=TRUE
  ) %dopar% {
    a <- lapply(simSummary,FUN=function(x){
      age_sim <- dimnames(x[["N.age.spawn"]])[[2]]
      out <- N.year.vec
      out[age_sim] <- x[["N.age.spawn"]][N.age.spawn.dn1[i],]
      out
    })
    b <- do.call(rbind.data.frame,a)
    dimnames(b) <- list(names(simSummary),N.age.spawn.dn2)
    return(b)
  }
  names(L.sim.N.age.spawn) <- N.age.spawn.dn1
  out.N.age.spawn <- L.sim.N.age.spawn
}

if("N.age.mdyr"%in%obj_collect){
  # N.age.mdyr
  N.age.mdyr.dn1 <- dimnames(rdat$N.age.mdyr)[[1]] # ages
  N.age.mdyr.dn2 <- dimnames(rdat$N.age.mdyr)[[2]] # years
  N.age.mdyr.vec <- setNames(rep(NA,length(N.age.mdyr.dn1)),N.age.mdyr.dn1)
  N.year.vec <- setNames(rep(NA,length(N.age.mdyr.dn2)),N.age.mdyr.dn2)

  L.sim.N.age.mdyr <- foreach::foreach(i=seq_along(N.age.mdyr.dn1),
                                        .multicombine=TRUE
  ) %dopar% {
    a <- lapply(simSummary,FUN=function(x){
      age_sim <- dimnames(x[["N.age.mdyr"]])[[2]]
      out <- N.year.vec
      out[age_sim] <- x[["N.age.mdyr"]][N.age.mdyr.dn1[i],]
      out
    })
    b <- do.call(rbind.data.frame,a)
    dimnames(b) <- list(names(simSummary),N.age.mdyr.dn2)
    return(b)
  }
  names(L.sim.N.age.mdyr) <- N.age.mdyr.dn1
  out.N.age.mdyr <- L.sim.N.age.mdyr
}
# if("Z.age"%in%obj_collect){
#   # Z.age
#   Z.age.dn1 <- dimnames(rdat$Z.age)[[1]]
#   Z.age.dn2 <- dimnames(rdat$Z.age)[[2]]
#   Z.age.vec <- setNames(rep(NA,length(Z.age.dn1)),Z.age.dn1)
#   L.sim.Z.age <- foreach::foreach(i=seq_along(Z.age.dn2),
#                                   .multicombine=TRUE
#   ) %dopar% {
#     a <- lapply(simSummary,FUN=function(x){
#       year_sim <- dimnames(x[["Z.age"]])[[1]]
#       out <- Z.age.vec
#       out[year_sim] <- x[["Z.age"]][,Z.age.dn2[i]]
#       out
#     })
#     b <- do.call(rbind.data.frame,a)
#     dimnames(b) <- list(names(simSummary),Z.age.dn1)
#     return(b)
#   }
#   names(L.sim.Z.age) <- Z.age.dn2
#   out.Z.age <- L.sim.Z.age
# }

if("Z.age"%in%obj_collect){
  # Z.age
  Z.age.dn1 <- dimnames(rdat$Z.age)[[1]] # ages
  Z.age.dn2 <- dimnames(rdat$Z.age)[[2]] # years
  Z.age.vec <- setNames(rep(NA,length(Z.age.dn1)),Z.age.dn1)
  N.year.vec <- setNames(rep(NA,length(Z.age.dn2)),Z.age.dn2)

  L.sim.Z.age <- foreach::foreach(i=seq_along(Z.age.dn1),
                                  .multicombine=TRUE
  ) %dopar% {
    a <- lapply(simSummary,FUN=function(x){
      age_sim <- dimnames(x[["Z.age"]])[[2]]
      out <- N.year.vec
      out[age_sim] <- x[["Z.age"]][Z.age.dn1[i],]
      out
    })
    b <- do.call(rbind.data.frame,a)
    dimnames(b) <- list(names(simSummary),Z.age.dn2)
    return(b)
  }
  names(L.sim.Z.age) <- Z.age.dn1
  out.Z.age <- L.sim.Z.age
}


if("parm.cons"%in%obj_collect){
  # parm.cons
  L.sim.parm.cons <- foreach::foreach(i=seq_along(parm.cons.names),
                            .multicombine=TRUE
  ) %dopar% {
    a <- lapply(simSummary,FUN=function(x){x[["parm.cons"]][[parm.cons.names[i]]]})
    b <- do.call(rbind.data.frame,a)
    names(b) <- obj_labels$parm.cons # Rename columns
    return(b)
  }
  names(L.sim.parm.cons) <- parm.cons.names
  out.parm.cons <- L.sim.parm.cons
}

## Lists
if("sel.age"%in%obj_collect){
  # sel.age
  sel.age.names <- names(rdat$sel.age)
  L.sim.sel.age <- list()

  for(i in seq_along(sel.age.names)){
    if(is.null(dim(rdat$sel.age[[i]]))){
      # vectors
      a <- lapply(simSummary,FUN=function(x){
        age_sim <- paste(x[["a.series"]][["age"]])
        out <- aser.vec
        out[age_sim] <- x[["sel.age"]][[sel.age.names[i]]]
        out
        })
      b <- do.call(rbind,a)
    }else{
      # matrices
      sel.age.dn1 <- dimnames(rdat$sel.age[[i]])[[1]]
      sel.age.dn2 <- dimnames(rdat$sel.age[[i]])[[2]]
      sel.age.vec <- setNames(rep(NA,length(sel.age.dn1)),sel.age.dn1)
      b <- foreach::foreach(j=seq_along(sel.age.dn2),
                             .multicombine=TRUE
      ) %dopar% {
        a <- lapply(simSummary,FUN=function(x){
          year_sim <- dimnames(x[["sel.age"]][[i]])[[1]]
          out <- sel.age.vec
          out[year_sim] <- x[["sel.age"]][[i]][,sel.age.dn2[j]]
          out
          })
        b <- do.call(rbind,a)
        return(b)
      }
      names(b) <- sel.age.dn2
    }
    L.sim.sel.age[[i]] <- b

  }
  names(L.sim.sel.age) <- sel.age.names
  out.sel.age <- L.sim.sel.age
}


# Assign out to object
out <- mget(paste0("out.",obj_collect))
names(out) <- gsub("^out.","",names(out))


# parallel shut down
  parallel::stopCluster(cl)

invisible(out)
}
