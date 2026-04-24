#' Standardize object names in bam files (dat, tpl, and cxx)
#'
#' This script reads in BAM dat, tpl, and cxx files, and converts them to R objects (character vectors of each line of the code). It then runs a set of gsub() functions to
#' replace object names with preferred naming conventions and returns dat, tpl, and cxx character vectors to use to rewrite the corresponding files with writeLines)
#' @param CommonName Common name of species modeled in BAM files. Only used when accessing dat, tpl, and cxx character vectors named as e.g. dat_CommonName
#' @param bam Output of \code{bam2r}.
#' @param dat_file dat file path
#' @param tpl_file tpl file path
#' @param cxx_file cxx file path
#' @param dat_obj dat file read in as a character vector with readLines(con=dat_file)
#' @param tpl_obj tpl file read in as a character vector with readLines(con=tpl_file)
#' @param cxx_obj cxx file read in as a character vector with readLines(con=cxx_file)
#' @param nm_pattern Character vector of additional object names (i.e. patterns) to be replaced
#' @param nm_replace Character vector of of replacements for nm1
#' @param match_whole_words Should patterns replace whole words only (pattern is regex enclosed by \\b)
#' @param fleet_key key to fleet names as a names list of character vectors.
#' names of list items are the replacements while the character vectors are sets of old fleet names to replace
#' @keywords bam stock assessment fisheries
#' @export
#' @examples
#' \dontrun{
#' # Read in and standardize the current BAM models
#' bam_AtMe <- standardize_bam("AtlanticMenhaden")
#' bam_BlSB <- standardize_bam("BlackSeaBass")
#' bam_BlTi <- standardize_bam("BluelineTilefish")
#' bam_Cobi <- standardize_bam("Cobia")
#' bam_GagG <- standardize_bam("GagGrouper")
#' bam_GrTr <- standardize_bam("GrayTriggerfish")
#' bam_GrAm <- standardize_bam("GreaterAmberjack")
#' bam_ReGr <- standardize_bam("RedGrouper")
#' bam_RePo <- standardize_bam("RedPorgy")
#' bam_ReSn <- standardize_bam("RedSnapper")
#' bam_ScGr <- standardize_bam("ScampGrouper")
#' bam_SnGr <- standardize_bam("SnowyGrouper")
#' bam_Tile <- standardize_bam("Tilefish")
#' bam_VeSn <- standardize_bam("VermilionSnapper")
#'
#' # Run a bam model and assign rdat output to object
#' rdat_RePo <- run_bam(bam=bam_RePo,fileName="RePo")
#' }

standardize_bam <- function(CommonName=NULL,
                            bam = NULL,
                            dat_file=NULL,tpl_file=NULL,cxx_file=NULL,
                            dat_obj=NULL, tpl_obj=NULL,cxx_obj=NULL,
                            nm_pattern=NULL, nm_replace=NULL,match_whole_words=TRUE,
                            fleet_key=list("sCT"=c("CVT"),              # Chevron trap (could possibly include video)
                                           "sTV"=c("CVID","Mcvt","mcvt"),# Combined chevron trap/video data (sCT in Red Porgy; Mcvt used in Black Seabass bass for the combined index but also for the comps which are really only from the trap)
                                           "sVD"=c("VID"),              # Video data (from chevron trap survey)
                                           "sBT"=c("Mbft"),             # MARMAP blackfish trap (see Black Seabass)
                                           "sBL"=c("sM"),               # MARMAP bottom longline survey (see Tilefish)
                                           "sVL"=c("vll"),              # MARMAP vertical longline survey (see Snowy Grouper)
                                           "sFT"=c("FST"),              # MARMAP Florida Snapper Trap (see Vermilion Snapper)
                                           "sAN"=c("nad"),              # Northern Adult composite survey index (Menhaden)
                                           "sAM"=c("mad"),              # Middle Adult composite survey index (Menhaden)
                                           "sAS"=c("sad"),              # Southern Adult composite survey index (Menhaden)
                                           "sJA"=c("jai"),              # Juvenile Abundance composite survey index (Menhaden)
                                           "sME"=c("mareco"),           # MARMAP and ECOMON survey index (Menhaden)
                                           "cDV"=c("cD"),               # Commercial diving (see Gag)
                                           "cHL"=c("cH","cHl"),         # Commercial handlines (a.k.a. commercial lines)
                                           "cLL"=c("cL"),               # Commercial longlines (see Blueline Tilefish)
                                           "cOT"=c("cO"),               # Commercial other (see Red Grouper)
                                           "cPT"=c("cp","cP"),          # Commercial pots (see Black Seabass)
                                           "cTW"=c("cT","cTw", "cHTR"), # Commercial trawl (see Black Seabass, Red Porgy, Vermilion Snapper)
                                           "cGN"=c("comm"),             # Commercial all. general commercial (see Black Seabass "D.comm.ob")
                                           "cRN"=c("cRn"),              # Commercial Reduction fishery, North  (Menhaden)
                                           "cRS"=c("cRs"),              # Commercial Reduction fishery, South  (Menhaden)
                                           "cBN"=c("cBn"),              # Commercial Bait fishery, North  (Menhaden)
                                           "cBS"=c("cBs"),              # Commercial Reduction fishery, South  (Menhaden)
                                           "rHB"=c("HB","hb","rHb"),    # Recreational headboat
                                           "rHD"=c("hbd","HBD"),      # Recreational headboat discards (atypical abbreviation found in Black Sea Bass selectivity parameters)
                                           "rGN"=c("GR","mrip","rGe","rA")  # Recreational all (a.k.a. general recreational, i.e. not headboat)
                            ),
                            fleet_replace = TRUE
){
  if(!is.null(CommonName)){
    dat <- get(paste0("dat_",CommonName))
    tpl <- get(paste0("tpl_",CommonName))
    cxx <- get(paste0("cxx_",CommonName))
  }

  if(!is.null(bam)){
    dat <- bam$dat
    tpl <- bam$tpl
    cxx <- bam$cxx
  }

  if(!is.null(dat_obj)&!is.null(tpl_obj)&!is.null(cxx_obj)){
    dat <- dat_obj
    tpl <- tpl_obj
    cxx <- cxx_obj
  }

  # Read in dat, tpl, and cxx files
  if(!is.null(dat_file)&!is.null(tpl_file)&!is.null(cxx_file)){
    dat <- readLines(con=dat_file)
    tpl <- readLines(con=tpl_file)
    cxx <- readLines(con=cxx_file)
  }

  dat <- trimws(dat) # Remove leading and/or trailing whitespace from character strings

  # Run bam2r to generate init
  bamr <- bam2r(
    dat_obj = dat,
    tpl_obj = tpl,
    cxx_obj = cxx)

  dat <- bamr$dat
  dat <- trimws(dat) # Remove leading and/or trailing whitespace from character strings
  tpl <- bamr$tpl
  cxx <- bamr$cxx
  init <- bamr$init

  # Restructure variable names so that fleet name is at the end (this makes it easier to identify similar variables in search)
  inm <- names(init)

  ## Landings
  # styr_L
  nm1_styr_L <- inm[grepl("^styr_.*_L$",inm)]
  nm2_styr_L <- paste0("styr_L_",gsub("^styr_|_L$","",nm1_styr_L))

  # endyr_L
  nm1_endyr_L <- inm[grepl("^endyr_.*_L$",inm)]
  nm2_endyr_L <- paste0("endyr_L_",gsub("^endyr_|_L$","",nm1_endyr_L))

  # obs_L
  nm1_obs_L <- inm[grepl("^obs_.*_L$",inm)]
  nm2_obs_L <- paste0("obs_L_",gsub("^obs_|_L$","",nm1_obs_L))

  # obs_cv_L (landings cvs)
  nm1_obs_cv_L <- inm[grepl("_L_cv$",inm)]
  nm2_obs_cv_L <- paste0("obs_cv_L_",gsub("_L_cv$","",nm1_obs_cv_L))

  ## Discards
  # styr_D
  nm1_styr_D <- inm[grepl("^styr(?!.*(cpue|lenc|agec)).*_D$",inm,perl=TRUE)]
  nm2_styr_D <- paste0("styr_D_",gsub("^styr_|_D$","",nm1_styr_D))

  # endyr_D
  nm1_endyr_D <- inm[grepl("^endyr(?!.*(cpue|lenc|agec)).*_D$",inm,perl=TRUE)]
  nm2_endyr_D <- paste0("endyr_D_",gsub("^endyr_|_D$","",nm1_endyr_D))

  # obs_released (undead discards)
  nm1_obs_D <- inm[grepl("^obs_.*_released$",inm)]
  nm2_obs_D <- paste0("obs_released_",gsub("^obs_|_released$","",nm1_obs_D))

  # obs_cv_D (discards cvs)
  nm1_obs_cv_D <- inm[grepl("_D_cv",inm)]
  nm2_obs_cv_D <- paste0("obs_cv_D_",gsub("_D_cv","",nm1_obs_cv_D))

  ## CPUE
  # styr_cpue
  nm1_styr_cpue <- inm[grepl("^styr_.*_cpue$",inm)]
  nm2_styr_cpue <- paste0("styr_cpue_",gsub("^styr_|_cpue$","",nm1_styr_cpue))

  # endyr_cpue
  nm1_endyr_cpue <- inm[grepl("^endyr_.*_cpue$",inm)]
  nm2_endyr_cpue <- paste0("endyr_cpue_",gsub("^endyr_|_cpue$","",nm1_endyr_cpue))


  # nyr_cpue (not found in most assessments)
  nm1_nyr_cpue <- inm[grepl("^nyr_.*_cpue$",inm)]
  nm2_nyr_cpue <- paste0("nyr_cpue_",gsub("^nyr_|_cpue$","",nm1_nyr_cpue))

  # yrs_cpue (not found in most assessments)
  nm1_yrs_cpue <- inm[grepl("^yrs_.*_cpue$",inm)]
  nm2_yrs_cpue <- paste0("yrs_cpue_",gsub("^yrs_|_cpue$","",nm1_yrs_cpue))

  # obs_cpue
  nm1_obs_cpue <- inm[grepl("^obs_.*_cpue$",inm)]
  nm2_obs_cpue <- paste0("obs_cpue_",gsub("^obs_|_cpue$","",nm1_obs_cpue))

  # obs_cv_cpue (cpue cvs)
  nm1_obs_cv_cpue <- inm[grepl("_cpue_cv",inm)]
  nm2_obs_cv_cpue <- paste0("obs_cv_cpue_",gsub("_cpue_cv","",nm1_obs_cv_cpue))


  ## Length comps
  # nyr_lenc
  nm1_nyr_lenc <- inm[grepl("^nyr_.*_lenc$",inm)]
  nm2_nyr_lenc <- paste0("nyr_lenc_",gsub("^nyr_|_lenc$","",nm1_nyr_lenc))

  # yrs_lenc
  nm1_yrs_lenc <- inm[grepl("^yrs_.*_lenc$",inm)]
  nm2_yrs_lenc <- paste0("yrs_lenc_",gsub("^yrs_|_lenc$","",nm1_yrs_lenc))

  # nsamp_lenc
  nm1_nsamp_lenc <- inm[grepl("^nsamp_.*_lenc$",inm)]
  nm2_nsamp_lenc <- paste0("nsamp_lenc_",gsub("^nsamp_|_lenc$","",nm1_nsamp_lenc))

  # nfish_lenc
  nm1_nfish_lenc <- inm[grepl("^nfish_.*_lenc$",inm)]
  nm2_nfish_lenc <- paste0("nfish_lenc_",gsub("^nfish_|_lenc$","",nm1_nfish_lenc))

  # obs_lenc
  nm1_obs_lenc <- inm[grepl("^obs_.*_lenc$",inm)]
  nm2_obs_lenc <- paste0("obs_lenc_",gsub("^obs_|_lenc$","",nm1_obs_lenc))

  # pred_lenc (not in dat, but in tpl and cxx)
  nm1_pred_lenc <- inm[grepl("^pred_.*_lenc$",inm)]
  nm2_pred_lenc <- paste0("pred_lenc_",gsub("^pred_|_lenc$","",nm1_pred_lenc))


  # lenc (rearrange and replace remaining variable names)
  nm1_lenc <- paste(gsub("obs_","",nm1_obs_lenc))
  nm2_lenc <- paste("lenc",gsub("_lenc","",nm1_lenc),sep="_")

  ## Length comps pooled (not in most assessments)
  # nyr_lenc_pool
  nm1_nyr_lenc_pool <- inm[grepl("^nyr.*lenc.*pool",inm)]
  nm2_nyr_lenc_pool <- paste0("nyr_lenc_pool_",gsub("nyr|lenc|pool|_","",nm1_nyr_lenc_pool))

  # yrs_lenc_pool
  nm1_yrs_lenc_pool <- inm[grepl("^yrs.*lenc.*pool",inm)]
  nm2_yrs_lenc_pool <- paste0("yrs_lenc_pool_",gsub("yrs|lenc|pool|_","",nm1_yrs_lenc_pool))

  # nsamp_lenc_pool
  nm1_nsamp_lenc_pool <- inm[grepl("^nsamp.*lenc.*pool",inm)]
  nm2_nsamp_lenc_pool <- paste0("nsamp_lenc_pool_",gsub("nsamp|lenc|pool|_","",nm1_nsamp_lenc_pool))

  # nfish_lenc_pool
  nm1_nfish_lenc_pool <- inm[grepl("^nfish.*lenc.*pool",inm)]
  nm2_nfish_lenc_pool <- paste0("nfish_lenc_pool_",gsub("nfish|lenc|pool|_","",nm1_nfish_lenc_pool))

  # obs_lenc_pool
  nm1_obs_lenc_pool <- inm[grepl("^obs.*lenc.*pool",inm)]
  nm2_obs_lenc_pool <- paste0("obs_lenc_pool_",gsub("obs|lenc|pool|_","",nm1_obs_lenc_pool))

  # pred_lenc_pool (not in dat, but in tpl and cxx)
  nm1_pred_lenc_pool <- inm[grepl("^pred.*lenc.*pool",inm)]
  nm2_pred_lenc_pool <- paste0("pred_lenc_pool_",gsub("pred|lenc|pool|_","",nm1_pred_lenc_pool))

  ## Age comps
  # nyr_agec
  nm1_nyr_agec <- inm[grepl("^nyr_.*_agec$",inm)]
  nm2_nyr_agec <- paste0("nyr_agec_",gsub("^nyr_|_agec$","",nm1_nyr_agec))

  # yrs_agec
  nm1_yrs_agec <- inm[grepl("^yrs_.*_agec$",inm)]
  nm2_yrs_agec <- paste0("yrs_agec_",gsub("^yrs_|_agec$","",nm1_yrs_agec))

  # nsamp_agec
  nm1_nsamp_agec <- inm[grepl("^nsamp_.*_agec$",inm)]
  nm2_nsamp_agec <- paste0("nsamp_agec_",gsub("^nsamp_|_agec$","",nm1_nsamp_agec))

  # nfish_agec
  nm1_nfish_agec <- inm[grepl("^nfish_.*_agec$",inm)]
  nm2_nfish_agec <- paste0("nfish_agec_",gsub("^nfish_|_agec$","",nm1_nfish_agec))

  # obs_agec
  nm1_obs_agec <- inm[grepl("^obs_.*_agec$",inm)]
  nm2_obs_agec <- paste0("obs_agec_",gsub("^obs_|_agec$","",nm1_obs_agec))

  # pred_agec (not in dat, but in tpl and cxx)
  nm1_pred_agec <- inm[grepl("^pred_.*_agec$",inm)]
  nm2_pred_agec <- paste0("pred_agec_",gsub("^pred_|_agec$","",nm1_pred_agec))

  # agec (rearrange and replace remaining variable names)
  nm1_agec <- paste(gsub("obs_","",nm1_obs_agec))
  nm2_agec <- paste("agec",gsub("_agec","",nm1_agec),sep="_")

  ## Age comps pooled (not in most assessments)
  # nyr_agec_pool
  nm1_nyr_agec_pool <- inm[grepl("^nyr_.*_agec_pool$",inm)]
  nm2_nyr_agec_pool <- paste0("nyr_agec_pool_",gsub("^nyr_|_agec_pool$","",nm1_nyr_agec_pool))

  # yrs_agec_pool
  nm1_yrs_agec_pool <- inm[grepl("^yrs_.*_agec_pool$",inm)]
  nm2_yrs_agec_pool <- paste0("yrs_agec_pool_",gsub("^yrs_|_agec_pool$","",nm1_yrs_agec_pool))

  # nsamp_agec_pool
  nm1_nsamp_agec_pool <- inm[grepl("^nsamp_.*_agec_pool$",inm)]
  nm2_nsamp_agec_pool <- paste0("nsamp_agec_pool_",gsub("^nsamp_|_agec_pool$","",nm1_nsamp_agec_pool))

  # nfish_agec_pool
  nm1_nfish_agec_pool <- inm[grepl("^nfish_.*_agec_pool$",inm)]
  nm2_nfish_agec_pool <- paste0("nfish_agec_pool_",gsub("^nfish_|_agec_pool$","",nm1_nfish_agec_pool))

  # obs_agec_pool
  nm1_obs_agec_pool <- inm[grepl("^obs_.*_agec_pool$",inm)]
  nm2_obs_agec_pool <- paste0("obs_agec_pool_",gsub("^obs_|_agec_pool$","",nm1_obs_agec_pool))

  # pred_agec_pool (not in dat, but in tpl and cxx)
  nm1_pred_agec_pool <- inm[grepl("^pred_.*_agec_pool$",inm)]
  nm2_pred_agec_pool <- paste0("pred_agec_pool_",gsub("^pred_|_agec_pool$","",nm1_pred_agec_pool))

  ## minSS
  nm1_minSS_lenc <- inm[grepl("^minSS_.*_lenc$",inm)]
  nm2_minSS_lenc <- paste0("minSS_lenc_",gsub("minSS_|_lenc","",nm1_minSS_lenc))

  nm1_minSS_agec <- inm[grepl("^minSS_.*_agec$",inm)]
  nm2_minSS_agec <- paste0("minSS_agec_",gsub("minSS_|_agec","",nm1_minSS_agec))

  ## set_ parameters
  # set_log_dm_lenc
  nm1_set_log_dm_lenc <- inm[grepl("^set_log_dm_.*_lc$",inm)]
  nm2_set_log_dm_lenc <- paste0("set_log_dm_lenc_",gsub("^set_log_dm_|_lc$","",nm1_set_log_dm_lenc))

  # set_log_dm_agec
  nm1_set_log_dm_agec <- inm[grepl("^set_log_dm_.*_ac$",inm)]
  nm2_set_log_dm_agec <- paste0("set_log_dm_agec_",gsub("^set_log_dm_|_ac$","",nm1_set_log_dm_agec))

  # set_log_q_cpue
  nm1_set_log_q_cpue <- inm[grepl("^set_log_q_|set_logq_",inm)]
  nm2_set_log_q_cpue <- paste0("set_log_q_cpue_",gsub("^set_log_q_|set_logq_|cpue_","",nm1_set_log_q_cpue))

  # F pars for landings (L)
  # set_log_avg_F_L
  # nm1_set_log_avg_F_L <- inm[grepl("^set_log_avg_F_.*[^D]$",inm)]
  #nm1_set_log_avg_F_L <- inm[grepl("^set_log_avg_F_(?![LD]_).*[^_D]$",inm,perl = TRUE)]
  nm1_set_log_avg_F_L <- inm[grepl("^set_log_avg_F_(?![LD]_)(?!.*_D$)",inm,perl = TRUE)]
  nm2_set_log_avg_F_L <- paste0("set_log_avg_F_L_",gsub("^set_log_avg_F_","",nm1_set_log_avg_F_L))

  # set_log_dev_vals_F_L
  nm1_set_log_dev_vals_F_L <- inm[grepl("^set_log_F_dev(?!.*_D_)(?!.*_D$).*(?=.*vals)",inm,perl=TRUE)]
  nm2_set_log_dev_vals_F_L <- paste0("set_log_dev_vals_F_L_",gsub("^set_log_F_dev|[^a-zA-Z0-9][LD][^a-zA-Z0-9]|(vals)|_","",nm1_set_log_dev_vals_F_L,perl=TRUE))

  # set_log_dev_F_L
  # nm1_set_log_dev_F_L <- inm[grepl("^set_log_F_dev_.*[^_vals|_D]$",inm)]
  # nm1_set_log_dev_F_L <- inm[grepl("^set_log_F_dev_(?![LD]|vals_).*[^_vals|_D]$",inm,perl = TRUE)]
  # nm1_set_log_dev_F_L <- inm[grepl("^set_log_F_dev_(?!.*vals)(?!.*_D)(?!.*_L)",inm,perl = TRUE)]
  nm1_set_log_dev_F_L <- inm[grepl("^set_log_F_dev(?!.*vals)(?!.*_D_)(?!.*_D$)",inm,perl = TRUE)]
  nm2_set_log_dev_F_L <- paste0("set_log_dev_F_L_",gsub("^set_log_F_dev_|L_|D_","",nm1_set_log_dev_F_L))

  # F pars for discards (D)
  # set_log_avg_F_D
  # nm1_set_log_avg_F_D <- inm[grepl("^set_log_avg_F_.*_D$",inm)]
  nm1_set_log_avg_F_D <- inm[grepl("^set_log_avg_F_(?![LD]_).*_D$",inm,perl = TRUE)]
  nm2_set_log_avg_F_D <- paste0("set_log_avg_F_D_",gsub("^set_log_avg_F_|_D$","",nm1_set_log_avg_F_D))

  # set_log_dev_vals_F_D
  nm1_set_log_dev_vals_F_D <- inm[grepl("^set_log_F_dev(?=.*vals)(?=.*_D)",inm,perl = TRUE)]
  nm2_set_log_dev_vals_F_D <- paste0("set_log_dev_vals_F_D_",gsub("^set_log_F_dev|[^a-zA-Z0-9][LD][^a-zA-Z0-9]|(vals)|_","",nm1_set_log_dev_vals_F_D,perl=TRUE))

  # set_log_dev_F_D
  # nm1_set_log_dev_F_D <- inm[grepl("^set_log_F_dev_.*_D$",inm)]
  nm1_set_log_dev_F_D <- inm[grepl("^set_log_F_dev(?!.*vals)(?=.*_D)",inm,perl = TRUE)]
  nm2_set_log_dev_F_D <- paste0("set_log_dev_F_D_",gsub("^set_log_F_dev_|D_|_D$","",nm1_set_log_dev_F_D))

  # deviations from recruitment time series
  # set_log_dev_rec
  nm1_set_log_dev_rec <- "set_log_rec_dev"
  nm2_set_log_dev_rec <- "set_log_dev_rec"

  # set_log_dev_vals_rec
  nm1_set_log_dev_vals_rec <- "set_log_rec_dev_vals"
  nm2_set_log_dev_vals_rec <- "set_log_dev_vals_rec"

  # deviations from initial numbers at age
  # set_log_dev_Nage
  nm1_set_log_dev_Nage <- "set_log_Nage_dev"
  nm2_set_log_dev_Nage <- "set_log_dev_Nage"

  # set_log_dev_vals_Nage
  nm1_set_log_dev_vals_Nage <- "set_log_Nage_dev_vals"
  nm2_set_log_dev_vals_Nage <- "set_log_dev_vals_Nage"

  # deviations from random walk of catchability
  # set_log_dev_RWq
  nm1_set_log_dev_RWq <- "set_log_RWq_dev"
  nm2_set_log_dev_RWq <- "set_log_dev_RWq"

  # selectivity parameters in logit space (see Black Sea Bass assessment)
  nm1_set_selpar_logit <- inm[grepl("^set_selpar_.*_logit",inm)]
  nm2_set_selpar_logit <- gsub("^(set_selpar)_(Age[0-9])_(.*)_(logit)$","\\1_\\4_\\2_\\3",nm1_set_selpar_logit)

  ## weight parameters
  # set_w_cpue
  nm1_set_w_cpue <- inm[grepl("^set_w_I_",inm)]
  nm2_set_w_cpue <- paste0("set_w_cpue_",gsub("^set_w_I_","",nm1_set_w_cpue))

  # set_w_lenc
  nm1_set_w_lenc <- inm[grepl("^set_w_lc_",inm)]
  nm2_set_w_lenc <- paste0("set_w_lenc_",gsub("^set_w_lc_","",nm1_set_w_lenc))

  # set_p_lenc (not in most assessments. See Black Sea Bass)
  nm1_set_p_lenc <- inm[grepl("^set_p_lenc_",inm)]
  nm2_set_p_lenc <- paste0("set_p_lenc_",gsub("^set_p_lenc_","",nm1_set_p_lenc))

  # set_w_agec
  nm1_set_w_agec <- inm[grepl("^set_w_ac_",inm)]
  nm2_set_w_agec <- paste0("set_w_agec_",gsub("^set_w_ac_","",nm1_set_w_agec))

  # tv (time-varying) objects (found in AtlanticMenhaden)
  nm1_tv <- inm[grepl("tv_",inm)]
  nm2_tv <- paste0(gsub("tv_","",nm1_tv),"_tv")

  # sex and maturity
  nm1_obs_prop_f <- "prop_f_obs"
  nm2_obs_prop_f <- "obs_prop_f"

  nm1_obs_prop_m <- "prop_m_obs"
  nm2_obs_prop_m <- "obs_prop_m"

  # obs_maturity
  nm1_obs_maturity <- inm[grepl("^maturity_.*_obs$",inm)]
  nm2_obs_maturity <- gsub("^(maturity)_(.*)_(obs)$","\\3_\\1_\\2",nm1_obs_maturity)

  # obs_prop
  nm1_obs_prop <- inm[grepl("^prop_.*_obs$",inm)]
  nm2_obs_prop <- gsub("^(prop)_(.*)_(obs)$","\\3_\\1_\\2",nm1_obs_prop)

  # Combine all nm1 and nm2
  nm1 <- list(nm1_styr_L,nm1_endyr_L,nm1_obs_L,nm1_obs_cv_L,
              nm1_styr_D,nm1_endyr_D,nm1_obs_D,nm1_obs_cv_D,
              nm1_styr_cpue,nm1_endyr_cpue,nm1_nyr_cpue,nm1_yrs_cpue,nm1_obs_cpue,nm1_obs_cv_cpue,
              nm1_nyr_lenc, nm1_yrs_lenc ,nm1_nsamp_lenc, nm1_nfish_lenc, nm1_obs_lenc, nm1_pred_lenc, nm1_lenc,
              nm1_nyr_lenc_pool, nm1_yrs_lenc_pool ,nm1_nsamp_lenc_pool, nm1_nfish_lenc_pool, nm1_obs_lenc_pool, nm1_pred_lenc_pool,
              nm1_nyr_agec_pool, nm1_yrs_agec_pool ,nm1_nsamp_agec_pool, nm1_nfish_agec_pool, nm1_obs_agec_pool, nm1_pred_agec_pool,
              nm1_nyr_agec, nm1_yrs_agec ,nm1_nsamp_agec, nm1_nfish_agec, nm1_obs_agec, nm1_pred_agec, nm1_agec,
              nm1_minSS_lenc, nm1_minSS_agec,
              nm1_set_log_dm_lenc,nm1_set_log_dm_agec, nm1_set_log_q_cpue,
              nm1_set_log_avg_F_L, nm1_set_log_dev_vals_F_L, nm1_set_log_dev_F_L,
              nm1_set_log_avg_F_D, nm1_set_log_dev_vals_F_D, nm1_set_log_dev_F_D,
              nm1_set_log_dev_rec, nm1_set_log_dev_vals_rec,
              nm1_set_log_dev_Nage, nm1_set_log_dev_vals_Nage,
              nm1_set_log_dev_RWq,
              nm1_set_selpar_logit,
              nm1_set_w_cpue,nm1_set_w_lenc,nm1_set_p_lenc,nm1_set_w_agec,
              nm1_tv,
              nm1_obs_prop, nm1_obs_maturity
              )

  nm2 <- list(nm2_styr_L,nm2_endyr_L,nm2_obs_L,nm2_obs_cv_L,
              nm2_styr_D,nm2_endyr_D,nm2_obs_D,nm2_obs_cv_D,
              nm2_styr_cpue,nm2_endyr_cpue,nm2_nyr_cpue,nm2_yrs_cpue,nm2_obs_cpue,nm2_obs_cv_cpue,
              nm2_nyr_lenc, nm2_yrs_lenc ,nm2_nsamp_lenc, nm2_nfish_lenc, nm2_obs_lenc, nm2_pred_lenc, nm2_lenc,
              nm2_nyr_lenc_pool, nm2_yrs_lenc_pool ,nm2_nsamp_lenc_pool, nm2_nfish_lenc_pool, nm2_obs_lenc_pool, nm2_pred_lenc_pool,
              nm2_nyr_agec_pool, nm2_yrs_agec_pool ,nm2_nsamp_agec_pool, nm2_nfish_agec_pool, nm2_obs_agec_pool, nm2_pred_agec_pool,
              nm2_nyr_agec, nm2_yrs_agec ,nm2_nsamp_agec, nm2_nfish_agec, nm2_obs_agec, nm2_pred_agec, nm2_agec,
              nm2_minSS_lenc, nm2_minSS_agec,
              nm2_set_log_dm_lenc,nm2_set_log_dm_agec, nm2_set_log_q_cpue,
              nm2_set_log_avg_F_L, nm2_set_log_dev_vals_F_L, nm2_set_log_dev_F_L,
              nm2_set_log_avg_F_D, nm2_set_log_dev_vals_F_D, nm2_set_log_dev_F_D,
              nm2_set_log_dev_rec, nm2_set_log_dev_vals_rec,
              nm2_set_log_dev_Nage, nm2_set_log_dev_vals_Nage,
              nm2_set_log_dev_RWq,
              nm2_set_selpar_logit,
              nm2_set_w_cpue,nm2_set_w_lenc,nm2_set_p_lenc,nm2_set_w_agec,
              nm2_tv,
              nm2_obs_prop, nm2_obs_maturity
              )

  if(!is.null(nm_pattern)&!is.null(nm_replace)){
    if(length(nm_pattern)!=length(nm_replace)){
      stop("length(nm_pattern)!=length(nm_replace)")
    }
    nm1 <- c(nm1,nm_pattern)
    nm2 <- c(nm2,nm_replace)
  }

  # Identify which elements in nm1 were not found in inm and remove those elements from nm1 and nm2
  ef <- unlist(lapply(nm1,function(x){length(x)!=0}))
  nm1 <- unlist(nm1[ef])
  nm2 <- unlist(nm2[ef])

  # Identify root of set_ parameters, remove set_ and add roots to nm1
  nm1_set <- nm1[substring(nm1,1,4)=="set_"]
  nm1_wo_set <- gsub("set_","",nm1_set)
  nm1 <- c(nm1,nm1_wo_set)


  # Identify root of set_ parameters, remove set_ and add roots to nm2
  nm2_set <- nm2[substring(nm2,1,4)=="set_"]
  nm2_wo_set <- gsub("set_","",nm2_set)
  nm2 <- c(nm2,nm2_wo_set)

  # For all nm1, replace them with nm2
  for(i in seq_along(nm1)){
    if(match_whole_words){
    nm1_i <- paste0("\\b",nm1[i],"\\b") # match whole words only
    }else{
    nm1_i <- nm1[i]
    }
    nm2_i <- nm2[i]

    dat <- gsub(pattern=nm1_i,replacement=nm2_i,x=dat,fixed=FALSE)
    tpl <- gsub(pattern=nm1_i,replacement=nm2_i,x=tpl,fixed=FALSE)
    cxx <- gsub(pattern=nm1_i,replacement=nm2_i,x=cxx,fixed=FALSE)

  }

  if(fleet_replace){
    # apply general fleet key
    for(replacement_i in names(fleet_key)){
      # pattern_beg_i <- paste0("^(",paste(fleet_key[[i]],collapse="|"),")\\.")
      # pattern_mid_i <- paste0("\\.(",paste(fleet_key[[i]],collapse="|"),")\\.")
      # pattern_end_i <- paste0("\\.(",paste(fleet_key[[i]],collapse="|"),")([0-9]*)$")
      # pattern <- paste0("(?<=[^a-zA-Z])",paste(fleet_key[[i]],collapse="|"),"(?=[^a-zA-Z])")
      pattern_i <- fleet_key[[replacement_i]]

      for(pattern_ij in pattern_i){
          dat <- gsub(paste0("(?<=[^a-zA-Z]|^)",pattern_ij,"(?=[^a-zA-Z]|$)"),replacement_i,dat,perl=TRUE)
          tpl <- gsub(paste0("(?<=[^a-zA-Z]|^)",pattern_ij,"(?=[^a-zA-Z]|$)"),replacement_i,tpl,perl=TRUE)
          cxx <- gsub(paste0("(?<=[^a-zA-Z]|^)",pattern_ij,"(?=[^a-zA-Z]|$)"),replacement_i,cxx,perl=TRUE)
        }
      }
    }


  # Rerun bam2r to regenerate output
  bamr <- bam2r(
    dat_obj = dat,
    tpl_obj = tpl,
    cxx_obj = cxx)


  return(bamr)
}
