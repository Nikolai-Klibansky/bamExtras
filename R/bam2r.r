#' BAM inputs to R object conversion function
#'
#' This script reads in BAM dat, tpl, and cxx files, and converts them to R objects (character vectors of each line of the code). It also builds an R list on the tpl init objects, from the tpl and dat files.
#' @param CommonName Common name of species modeled in BAM files. Only used when accessing dat, tpl, and cxx character vectors named as e.g. dat_CommonName
#' @param dat_file dat file path
#' @param tpl_file tpl file path
#' @param cxx_file cxx file path
#' @param dat_obj dat file read in as a character vector with readLines(con=dat_file)
#' @param tpl_obj tpl file read in as a character vector with readLines(con=tpl_file)
#' @param cxx_obj cxx file read in as a character vector with readLines(con=cxx_file)
#' @param init user supplied list of tpl "init" object names defined in tpl Data Section. If missing, the function builds one from values in the dat and tpl files.
#' @param object_labels list of labels for bounded points and vectors
#' @param is_bam Are you using this function on an actual bam model? If FALSE,
#' code that is very specific to bam will not be run. The function was designed
#' to read bam model files into R, but in many cases it should work to read other
#' dat and tpl files.
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' # Read in any of the current BAM models
#' bam_AtMe <- bam2r("AtlanticMenhaden")
#' bam_BlSB <- bam2r("BlackSeaBass")
#' bam_BlTi <- bam2r("BluelineTilefish")
#' bam_Cobi <- bam2r("Cobia")
#' bam_GagG <- bam2r("GagGrouper")
#' bam_GrTr <- bam2r("GrayTriggerfish")
#' bam_GrAm <- bam2r("GreaterAmberjack")
#' bam_ReGr <- bam2r("RedGrouper")
#' bam_RePo <- bam2r("RedPorgy")
#' bam_ReSn <- bam2r("RedSnapper")
#' bam_ScGr <- bam2r("ScampGrouper")
#' bam_SnGr <- bam2r("SnowyGrouper")
#' bam_Tile <- bam2r("Tilefish")
#' bam_VeSn <- bam2r("VermilionSnapper")
#'
#' # Run a bam model and assign rdat output to object
#' rdat_RePo <- run_bam(bam=bam_RePo,fileName="RePo")
#'
#' # Modify data input from a BAM model and incorporate it back into the dat file object
#' init2 <- bam_RePo$init
#' init2$set_steep[c(1,5)] <- paste(0.6) # Change steepness
#' init2$set_steep[4] <- paste(-abs(as.numeric(init2$set_steep[4]))) # Fix steepness
#' bam_RePo2 <- bam2r("RedPorgy",init=init2)
#' rdat_RePo2 <- run_bam(bam=bam_RePo2,fileName="RePo2")
#'
#' # Compare models
#' plot(rdat_RePo$t.series$year, rdat_RePo$t.series$SSB.msst,type="o", ylim=range(c(rdat_RePo$t.series$SSB.msst,rdat_RePo2$t.series$SSB.msst),na.rm=TRUE), xlab="year",ylab="SSB.msst")
#' points(rdat_RePo2$t.series$year,rdat_RePo2$t.series$SSB.msst,type="o",col="blue")
#' legend("topright",legend=c("RePo","RePo2"),pch=1,lty=1,col=c("black","blue"),bty="n")
#' }



bam2r <- function(CommonName=NULL,init=NULL,
                  dat_file=NULL,tpl_file=NULL,cxx_file=NULL,
                  dat_obj=NULL, tpl_obj=NULL,cxx_obj=NULL,
                  object_labels=list("bounded_point"=c("init","lower","upper","phase","prior_mean","prior_var","prior_pdf"),
                                     "bounded_vector"=c("lower","upper","phase")
                                     ),
                  is_bam=TRUE

){
  fcn_args <- as.list(environment()) # Gets function arguments
  env_tmp <- new.env() # Create new empty environment

  if(!is.null(CommonName)){
    dat <- get(paste0("dat_",CommonName))
    tpl <- get(paste0("tpl_",CommonName))
    cxx <- get(paste0("cxx_",CommonName))
  }

  if(!is.null(dat_obj) & !is.null(tpl_obj) & !is.null(cxx_obj)){
    dat <- dat_obj
    tpl <- tpl_obj
    cxx <- cxx_obj
  }

  # Read in dat, tpl, and cxx files
  if(!is.null(dat_file) & !is.null(tpl_file) & !is.null(cxx_file)){
    dat <- readLines(con=dat_file)
    tpl <- readLines(con=tpl_file)
    cxx <- readLines(con=cxx_file)
  }

  dat <- trimws(dat) # Remove leading and/or trailing whitespace from character strings

  ### tpl ###

  # Isolate tpl DATA_SECTION
  tpl_trimws <- trimws(x=tpl,which="left")
  tplDATA <- tpl_trimws[which(substr(tpl_trimws,start=1,stop=nchar("DATA_SECTION"))=="DATA_SECTION"):
                          (which(substr(tpl_trimws,start=1,stop=nchar("PARAMETER_SECTION"))=="PARAMETER_SECTION")-1)]

  # Identify which items in tplDATA correspond to initialization of variables corresponding to values in the dat file
  tplDATA_items <- tplDATA[grepl("^init_",tplDATA)]

  # Create a data frame that can be used to parse the dat file and create an R list
  D_tplVar <- local({
    # Split items into meaningful parts and assemble a data frame
    sp1 <- strsplit(tplDATA_items,split=";")

    sp2 <- lapply(sp1,FUN=function(x){x[[1]]})

    sp2ab <- strsplit(unlist(sp2),split=" ",fixed=TRUE)
    type <- stringr::str_extract(tplDATA_items,"^[a-z_0-9]*")
    nameArg <- stringr::str_extract(tplDATA_items,"(?<= )[A-Za-z_0-9,\\(\\)]*")
    name <- stringr::str_extract(tplDATA_items,"(?<= )[A-Za-z_0-9]*")
    arg <- stringr::str_extract(nameArg,"(?<=\\().+?(?=\\))")
    comment <- stringr::str_extract(tplDATA_items,"(?<=//).*")

    # Create a data frame that can be used to parse the dat file and create an R list
    return(
      data.frame("type"=type,"name"=name,"arg"=arg,"comment"=comment,stringsAsFactors=FALSE)
    )
  })

  #### LnonInitDdatnum ####

  # Identify what type of information is in each row of the dat file (numeric or non-numeric)
  dat_num_rows_logical <- grepl("^[0-9\\.-]",dat)
  dat_num_rows <- which(dat_num_rows_logical)
  dat_nonNum_rows <- which(!dat_num_rows_logical) # Rows of dat file with non-numeric row values (i.e. comments and whitespace)

  # List of non-numeric row values from dat file
  L_nonInit <- local({
    a <- as.list(dat[dat_nonNum_rows])
    names(a) <- paste("datRow",sprintf("%05.0f", dat_nonNum_rows))
    a
  })

  # Create a data frame that includes only numeric values from the dat file that the tpl reads, and adjacent comments
  D_dat_num <- local({
    dat_num_row <- dat[dat_num_rows] # dat rows that begin with a numeric value
    dat_value <- trimws(stringr::str_extract(dat_num_row,"^[^#]*")) # start at the beginning and get everything until you hit an octothorp
    dat_comment <- trimws(stringr::str_extract(dat_num_row,"(?<=#).*")) # Get everything after the first octothorp
    data.frame("value"=dat_value,"datRow"=dat_num_rows,"comment"=dat_comment,stringsAsFactors=FALSE)
  })

  # Create a vector all (uncommented) numeric values from the dat file
  dat_vals <- unlist(D_dat_num$value) %>% paste(collapse=" ") %>% gsub(pattern="\t",replacement = " ") %>% strsplit(split="\\s+") %>% unlist()


  #### Build init ####
  # Build R list with names equal to BAM tpl init object names and values equal to init objects from dat file
  initNames <- D_tplVar[,"name"]
  init <- vector("list", length(initNames))
  names(init) <- initNames
  # Row of dat file associated with the first row of each init object
  datRowInit <- local({
    a <- rep(NA,length(init))
    names(a) <- names(init)
    return(a)
  })

  for(i in 1:nrow(D_tplVar)){
    if(i==1){
      datRow_i <- 0 # Initialize dat row counter
      dat_vals_read_i <- 0 # Initialize value of the last dat value to be read in and assigned to a variable
    }

    type_i <- D_tplVar[i,"type"]
    name_i <- D_tplVar[i,"name"]
    arg_i <- D_tplVar[i,"arg"]

    datRow_first <- datRow_i + 1 # dat row counter associated with the first row of init object i
    datRowInit[name_i] <- D_dat_num[datRow_first,"datRow"]

    if(type_i%in%c("init_int","init_number")){ # all should be length 1
      datRow_i <- datRow_first

      varValueLength_i <- 1
      dat_vals_i <- dat_vals_read_i+(1:varValueLength_i)
      varValue_i <- dat_vals[dat_vals_i]
      dat_vals_read_i <- tail(dat_vals_i,1) # update value

      if(type_i=="init_int"){
        varValue_i <- as.integer(varValue_i)
      }
      if(type_i=="init_number"){
        varValue_i <- as.numeric(varValue_i)
      }
    }

    if(type_i%in%c("init_vector","init_ivector")){
      datRow_i <- datRow_first
      # Initialize vector
      varValueNames_i <- local({
        lim1 <- strsplit(arg_i,split=",",fixed=TRUE)[[1]][1]
        lim2 <- strsplit(arg_i,split=",",fixed=TRUE)[[1]][2]
        lim1 <- type.convert(lim1,as.is=TRUE) # Converts character vector to numeric if it can. Otherwise leaves it as character
        lim2 <- type.convert(lim2,as.is=TRUE) # Converts character vector to numeric if it can. Otherwise leaves it as character
        if(is_bam & type_i=="init_vector" & grepl("^set_(?!log_dev_vals)",name_i,perl=TRUE) & lim1==1 & lim2%in%c(3,7)){
          if(lim2==7){
            out <- object_labels$bounded_point
          } else if(lim2==3){
            out <- object_labels$bounded_vector
          }
        } else if (is_bam & grepl("^nyr",lim2) & exists(gsub("nyr","yrs",lim2),envir=env_tmp,inherits=FALSE)){
        out <- get(gsub("nyr","yrs",lim2), envir=env_tmp)
          } else if (is_bam & grepl("^nages",lim2) & exists(paste0("agebins",gsub("nages","",lim2)),envir=env_tmp,inherits=FALSE)){
            out <- get(paste0("agebins",gsub("nages","",lim2)),envir=env_tmp,inherits=FALSE)[lim1:get(lim2,envir=env_tmp,inherits=FALSE)]
            } else {
          if(is.character(lim1)){lim1 <- get(lim1,envir=env_tmp,inherits=FALSE)}
          if(is.character(lim2)){lim2 <- get(lim2,envir=env_tmp,inherits=FALSE)}
        out <- lim1:lim2
        }
        return(out)
      })

      varValueLength_i <-  as.numeric(length(varValueNames_i))

      if(type_i=="init_vector"){
        varValue_i <- vector(mode = "numeric", length = varValueLength_i)
      }
      if(type_i=="init_ivector"){
        varValue_i <- vector(mode = "integer", length = varValueLength_i)
      }

      # Assign values to vector
      dat_vals_i <- dat_vals_read_i+(1:varValueLength_i)
      varValue_i <- setNames(dat_vals[dat_vals_i],varValueNames_i)
      dat_vals_read_i <- tail(dat_vals_i,1) # update value
    }

    if(type_i%in%c("init_matrix","init_imatrix")){
      datRow_first <- datRow_first # First row of matrix
      matrix_rownames <- local({
        rowlim1 <- strsplit(arg_i,split=",",fixed=TRUE)[[1]][1]
        rowlim1 <- type.convert(rowlim1,as.is=TRUE) # Converts character vector to numeric if it can. Otherwise leaves it as character
        if(is.character(rowlim1)){rowlim1 <- get(rowlim1,envir=env_tmp,inherits=FALSE)}

        rowlim2 <- strsplit(arg_i,split=",",fixed=TRUE)[[1]][2]
        rowlim2 <- type.convert(rowlim2,as.is=TRUE) # Converts character vector to numeric if it can. Otherwise leaves it as character
        if(is_bam & grepl("^nyr",rowlim2)){
          out <-  get(gsub("nyr","yrs",rowlim2),envir=env_tmp,inherits=FALSE)
        } else if (is_bam & grepl("^nages",rowlim2) & exists(paste0("agebins",gsub("nages","",rowlim2)),envir=env_tmp,inherits=FALSE)){
          out <- get(paste0("agebins",gsub("nages","",rowlim2)),envir=env_tmp,inherits=FALSE)[rowlim1:get(rowlim2,envir=env_tmp,inherits=FALSE)]
        } else if(is.character(rowlim2)){
          rowlim2 <- get(rowlim2,envir=env_tmp,inherits=FALSE)
          out <- rowlim1:rowlim2
          }
      return(out)
    })

      matrix_nrow <- length(matrix_rownames)

      matrix_colnames <- local({
        collim1 <- strsplit(arg_i,split=",",fixed=TRUE)[[1]][3]
        collim1 <- type.convert(collim1,as.is=TRUE) # Converts character vector to numeric if it can. Otherwise leaves it as character
        if(is.character(collim1)){collim1 <- get(collim1,envir=env_tmp,inherits=FALSE)}

        collim2 <- strsplit(arg_i,split=",",fixed=TRUE)[[1]][4]
        collim2 <- type.convert(collim2,as.is=TRUE) # Converts character vector to numeric if it can. Otherwise leaves it as character
        if(is_bam & grepl("^nlenbins",collim2)){
          out <- get(gsub("^nlenbins","lenbins",collim2),envir=env_tmp,inherits=FALSE)
        } else if (is_bam & grepl("^nages",collim2) & exists(paste0("agebins",gsub("nages","",collim2)),envir=env_tmp,inherits=FALSE)){
          out <- get(paste0("agebins",gsub("nages","",collim2)),envir=env_tmp,inherits=FALSE)[collim1:get(collim2,envir=env_tmp,inherits=FALSE)]
        } else if(is_bam & is.character(collim2) & !grepl("^nlenbins",collim2)){
          collim2 <- get(collim2,envir=env_tmp,inherits=FALSE)
          out <- collim1:collim2
        } else{
          if(collim2%in%ls(env_tmp)){
            collim2 <- get(collim2,envir=env_tmp,inherits=FALSE)
          }
          out <- collim1:collim2
        }
        return(out)
      })

      matrix_ncol <- length(matrix_colnames)

      datRow_last <- datRow_first+matrix_nrow-1 # last row of matrix

      varValueLength_i <- matrix_nrow*matrix_ncol
      dat_vals_i <- dat_vals_read_i+(1:varValueLength_i)
      varValue_i <- dat_vals[dat_vals_i]
      dat_vals_read_i <- tail(dat_vals_i,1) # update value

      # Restructure varValue_i as matrix
      varValue_i <- matrix(
        as.character(varValue_i),
        byrow=TRUE,
        nrow=matrix_nrow, ncol=matrix_ncol,
        dimnames=list(matrix_rownames,matrix_colnames))
      datRow_i <- datRow_last
    }
    varValue_i <- trimws(varValue_i) # Remove extra whitespace from values

    assign(name_i,varValue_i, envir = env_tmp) # Assign object to current environment (important for defining some later objects)
    init[[name_i]] <- varValue_i # Assign value to init
  }

  # If init was supplied in the function call, override init.
  if(!is.null(fcn_args$init)){
    init <- fcn_args$init
  }

  #### Create dat_list ####
  # Create new list of dat file elements incorporating init

  dat_list <- local({

    init2 <- init
    names(init2) <- paste("datRow",sprintf("%05.0f", datRowInit[names(init)]))

    dat_list <- c(init2,L_nonInit)

    dat_list <- dat_list[order(names(dat_list))]

    # convert all items to make one long character vector
    dat_list <- lapply(dat_list,FUN=function(x){
      if("matrix"%in%class(x)){
        out <- apply(x,1,FUN=function(x){paste(x,collapse="\t")})
      }else{
        out <- paste(x,collapse="\t")
      }
    })
    return(dat_list)
  })
  # Add available comments back to rows with numeric values
  dat_list <- local({
    datRowsWithComments <- D_dat_num[which(D_dat_num$comment!=""),"datRow"]
    datRowsWithCommentsNames <- paste("datRow",sprintf("%05.0f", datRowsWithComments)) # numeric dat rows with comments
    datRowsWithCommentsComments <- D_dat_num$comment[D_dat_num$datRow%in%datRowsWithComments]
    dat_list[datRowsWithCommentsNames] <- paste(dat_list[datRowsWithCommentsNames],datRowsWithCommentsComments,sep=" # ")

    dat_list
  })

  return(list("dat"=unlist(dat_list),"tpl"=tpl,"cxx"=cxx,"init"=init,"dat_list"=dat_list))
}

#' @rdname bam2r
#' @export
