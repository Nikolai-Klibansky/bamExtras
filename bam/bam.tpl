//##  Author: NMFS, Beaufort Lab, Atlantic Fisheries Branch
//##  Analyst: Nikolai Klibansky
//##  Species: Black Sea Bass
//##  Region: US South Atlantic
//##  SEDAR: 76 update 2026
//##  Date: 2026-04-28

DATA_SECTION

!!cout << "Starting Beaufort Assessment Model" << endl;
!!cout << endl;
!!cout << "                BAM!" << endl;
!!cout << endl;

//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: set-up section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>

// Starting and ending year of the model (year data starts)
init_int styr;  // Start year of model
init_int endyr; // End year of model

init_int styr_rec_dev;  // Start year to estimate deviation from S-R curve
init_int endyr_rec_dev; // End year to estimate deviation from S-R curve

//Starting year to estimate arithmetic average of recruitment (styr_rec_mu through endyr_rec_dev); this mean recruitment applied after endyr_rec_dev
init_int styr_rec_mu;

// Periods of constraints on recruitment deviations
init_int endyr_rec_phase1; // end year of recruitment phase 1
init_int endyr_rec_phase2; // end year of recruitment phase 2

init_int nyr; // number of assessment years
init_int nyr_rec_dev; // number of rec dev years
init_int nyr_rec_mu; //number years used to compute mean recruitment for yrs after endyr_rec_dev

// Ending years for selectivity periods
init_int endyr_sel_b01;
init_int endyr_sel_b02;
init_int endyr_sel_b03;
init_int endyr_sel_b04;

init_int styr_closure_com; // Starting year for commercial closures. Really, just used for selectivity of commercial discards.

//size limits
init_number sizelim1;
init_number sizelim2;
init_number sizelim3;
init_number sizelim4;
init_number sizelim5;

//// Ages setup ////
init_int nages; // Total number of ages in population model
init_vector agebins(1,nages); // Vector of ages for age bins

//// Lengths setup ////
init_int nlenbins;          // Number of length bins used to match length comps
init_number lenbins_width;  // Width of length bins used to match length comps
// Vector of lengths for length bins
init_vector lenbins(1,nlenbins);
  
init_number max_F_spr_msy; // Max F used in spr and msy calcs

init_int n_iter_msy; // Number of iterations in msy calculations
init_int n_iter_spr; // Number of iterations in spr calculations


// computing arithmetic average recruitment for SPR-related values
init_int styr_rec_spr;  // Start year
init_int endyr_rec_spr; // End year
init_int nyr_rec_spr;  // Number of years
 
// Switch for recruitment calcs in SPR-related values
//	Switch=1 uses recruitment averaged over years specified directly above. Switch=1 preferable if SR function switch = 1 (B-H) or 2 (Ricker)
//	Switch=2 uses R0 (which may be bias-corrected). Switch=2 preferable if SR function switch = 3 (i.e., none). With this option, SPR-related years input above are ignored.
init_int SPR_rec_switch;

// SPR proxy (as a ratio): F30, F40, F50 computed, but most spr output will correspond to the selection below
init_number spr_proxy;

// Switch for spawner-recruit function (integer 1=Beverton-Holt, 2=Ricker, 3=None (constant average recruitment))
init_int SR_switch;

//Number years at end of time series over which to average sector F's, for weighted selectivities
init_int nyr_selpar_wgted;

// Bias correction (set to 1.0 for no bias correction or a negative value to compute from rec variance)
init_number set_BiasCor;
  
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: observed data section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>

// sBT
init_int styr_cpue_sBT;
init_int endyr_cpue_sBT;
init_vector obs_cpue_sBT(styr_cpue_sBT,endyr_cpue_sBT);
init_vector obs_cv_cpue_sBT(styr_cpue_sBT,endyr_cpue_sBT);

init_int nyr_lenc_sBT;
init_ivector yrs_lenc_sBT(1,nyr_lenc_sBT);
init_vector nsamp_lenc_sBT(1,nyr_lenc_sBT);
init_vector nfish_lenc_sBT(1,nyr_lenc_sBT);
init_matrix obs_lenc_sBT(1,nyr_lenc_sBT,1,nlenbins);

init_int nyr_agec_sBT;
init_ivector yrs_agec_sBT(1,nyr_agec_sBT);
init_vector nsamp_agec_sBT(1,nyr_agec_sBT);
init_vector nfish_agec_sBT(1,nyr_agec_sBT);
init_matrix obs_agec_sBT(1,nyr_agec_sBT,1,nages);

// sTV
init_int styr_cpue_sTV;
init_int endyr_cpue_sTV;
init_vector obs_cpue_sTV(styr_cpue_sTV,endyr_cpue_sTV);
init_vector obs_cv_cpue_sTV(styr_cpue_sTV,endyr_cpue_sTV);

init_int nyr_agec_sTV;
init_ivector yrs_agec_sTV(1,nyr_agec_sTV);
init_vector nsamp_agec_sTV(1,nyr_agec_sTV);
init_vector nfish_agec_sTV(1,nyr_agec_sTV);
init_matrix obs_agec_sTV(1,nyr_agec_sTV,1,nages);

// cHL
init_int styr_cpue_cHL;                                             
init_int endyr_cpue_cHL;                                            
init_vector obs_cpue_cHL(styr_cpue_cHL,endyr_cpue_cHL);
init_vector obs_cv_cpue_cHL(styr_cpue_cHL,endyr_cpue_cHL);
 
init_int styr_L_cHL;
init_int endyr_L_cHL;
init_vector obs_L_cHL(styr_L_cHL,endyr_L_cHL);
init_vector obs_cv_L_cHL(styr_L_cHL,endyr_L_cHL);

init_int nyr_lenc_cHL;
init_ivector yrs_lenc_cHL(1,nyr_lenc_cHL);
init_vector nsamp_lenc_cHL(1,nyr_lenc_cHL);
init_vector nfish_lenc_cHL(1,nyr_lenc_cHL);
init_matrix obs_lenc_cHL(1,nyr_lenc_cHL,1,nlenbins);

init_int nyr_agec_cHL;
init_ivector yrs_agec_cHL(1,nyr_agec_cHL);
init_vector nsamp_agec_cHL(1,nyr_agec_cHL);
init_vector nfish_agec_cHL(1,nyr_agec_cHL);
init_matrix obs_agec_cHL(1,nyr_agec_cHL,1,nages);

// // cHD all years, during the open season
// init_int styr_D_cHD;
// init_int endyr_D_cHD;
// init_vector obs_released_cHD(styr_D_cHD,endyr_D_cHD);
// init_vector obs_cv_D_cHD(styr_D_cHD,endyr_D_cHD);
// 
// // cHD during years were there was a close season, during the closed season
// init_int styr_D_cHDclosed;
// init_int endyr_D_cHDclosed;
// init_vector obs_released_cHDclosed(styr_D_cHDclosed,endyr_D_cHDclosed);

// cPT
init_int styr_L_cPT;
init_int endyr_L_cPT;
init_vector obs_L_cPT(styr_L_cPT,endyr_L_cPT);
init_vector obs_cv_L_cPT(styr_L_cPT,endyr_L_cPT);

init_int nyr_lenc_cPT;
init_ivector yrs_lenc_cPT(1,nyr_lenc_cPT);
init_vector nsamp_lenc_cPT(1,nyr_lenc_cPT);
init_vector nfish_lenc_cPT(1,nyr_lenc_cPT);
init_matrix obs_lenc_cPT(1,nyr_lenc_cPT,1,nlenbins);

init_int nyr_agec_cPT;
init_ivector yrs_agec_cPT(1,nyr_agec_cPT);
init_vector nsamp_agec_cPT(1,nyr_agec_cPT);
init_vector nfish_agec_cPT(1,nyr_agec_cPT);
init_matrix obs_agec_cPT(1,nyr_agec_cPT,1,nages);

// // cPD all years, during the open season
// init_int styr_D_cPD;
// init_int endyr_D_cPD;
// init_vector obs_released_cPD(styr_D_cPD,endyr_D_cPD);
// init_vector obs_cv_D_cPD(styr_D_cPD,endyr_D_cPD);
// 
// // cPD during years were there was a close season, during the closed season
// init_int styr_D_cPDclosed;
// init_int endyr_D_cPDclosed;
// init_vector obs_released_cPDclosed(styr_D_cPDclosed,endyr_D_cPDclosed);

// cGD
init_int styr_D_cGD;                                             
init_int endyr_D_cGD; 
init_vector obs_released_cGD(styr_D_cGD,endyr_D_cGD);
init_vector obs_cv_D_cGD(styr_D_cGD,endyr_D_cGD);

// cTW
init_int styr_L_cTW;
init_int endyr_L_cTW;
init_vector obs_L_cTW(styr_L_cTW,endyr_L_cTW);
init_vector obs_cv_L_cTW(styr_L_cTW,endyr_L_cTW);

// rHB
init_int styr_cpue_rHB;
init_int endyr_cpue_rHB;
init_vector obs_cpue_rHB(styr_cpue_rHB,endyr_cpue_rHB);
init_vector obs_cv_cpue_rHB(styr_cpue_rHB,endyr_cpue_rHB);

init_int styr_L_rHB;
init_int endyr_L_rHB;
init_vector obs_L_rHB(styr_L_rHB,endyr_L_rHB);
init_vector obs_cv_L_rHB(styr_L_rHB,endyr_L_rHB);

init_int nyr_lenc_rHB;
init_ivector yrs_lenc_rHB(1,nyr_lenc_rHB);
init_vector nsamp_lenc_rHB(1,nyr_lenc_rHB);
init_vector nfish_lenc_rHB(1,nyr_lenc_rHB);
init_matrix obs_lenc_rHB(1,nyr_lenc_rHB,1,nlenbins);

init_int nyr_agec_rHB;
init_ivector yrs_agec_rHB(1,nyr_agec_rHB);
init_vector nsamp_agec_rHB(1,nyr_agec_rHB);
init_vector nfish_agec_rHB(1,nyr_agec_rHB);
init_matrix obs_agec_rHB(1,nyr_agec_rHB,1,nages);

// rHD

init_int styr_D_rHD;
init_int endyr_D_rHD;
init_vector obs_released_rHD(styr_D_rHD,endyr_D_rHD);
init_vector obs_cv_D_rHD(styr_D_rHD,endyr_D_rHD);

init_int nyr_lenc_rHD;
init_ivector yrs_lenc_rHD(1,nyr_lenc_rHD);
init_vector nsamp_lenc_rHD(1,nyr_lenc_rHD);
init_vector nfish_lenc_rHD(1,nyr_lenc_rHD);
init_matrix obs_lenc_rHD(1,nyr_lenc_rHD,1,nlenbins);

// rGN
init_int styr_L_rGN;
init_int endyr_L_rGN;
init_vector obs_L_rGN(styr_L_rGN,endyr_L_rGN);
init_vector obs_cv_L_rGN(styr_L_rGN,endyr_L_rGN);

init_int nyr_lenc_rGN;
init_ivector yrs_lenc_rGN(1,nyr_lenc_rGN);
init_vector nsamp_lenc_rGN(1,nyr_lenc_rGN);
init_vector nfish_lenc_rGN(1,nyr_lenc_rGN);
init_matrix obs_lenc_rGN(1,nyr_lenc_rGN,1,nlenbins);

// rIN
init_int styr_L_rIN;
init_int endyr_L_rIN;
init_vector obs_L_rIN(styr_L_rIN,endyr_L_rIN);
init_vector obs_cv_L_rIN(styr_L_rIN,endyr_L_rIN);

// rGD
init_int styr_D_rGD;
init_int endyr_D_rGD;
init_vector obs_released_rGD(styr_D_rGD,endyr_D_rGD);
init_vector obs_cv_D_rGD(styr_D_rGD,endyr_D_rGD);

// rID
init_int styr_D_rID;
init_int endyr_D_rID;
init_vector obs_released_rID(styr_D_rID,endyr_D_rID);
init_vector obs_cv_D_rID(styr_D_rID,endyr_D_rID);

// Dmort
init_number set_Dmort_cGD;
// init_number set_Dmort_cPD1;
// init_number set_Dmort_cPD2;
// init_number set_Dmort_cHD;
init_number set_Dmort_rHD;
init_number set_Dmort_rGD;
init_number set_Dmort_rID;

//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--##
//-- BAM DATA_SECTION: parameter section
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--##

// ############ Parameter values and initial guesses #############
// ###############################################################
// ###prior PDF (1=none, 2=lognormal, 3=normal, 4=beta)          #
// ###############################################################
// ##initial # lower # upper #       #  prior  # prior   # prior #
// ## guess  # bound # bound # phase #  mean   # var/-CV #  PDF  #
// ##--------#-------#-------#-------#---------#---------#-------#
//// Von Bertalanffy parameters
init_vector set_Linf(1,7);
init_vector set_K(1,7);
init_vector set_t0(1,7);
init_vector set_len_cv(1,7);

////  Natural mortality
init_vector set_M_constant(1,7);

//// Spawner-recruit
init_vector set_steep(1,7);
init_vector set_log_R0(1,7);
init_vector set_R_autocorr(1,7);
init_vector set_rec_sigma(1,7);

//// Dirichlet-multinomial overdispersion parameters
init_vector set_log_dm_lenc_sBT(1,7);    
init_vector set_log_dm_lenc_cHL(1,7);    
init_vector set_log_dm_lenc_cPT(1,7);    
init_vector set_log_dm_lenc_rHB(1,7);  
init_vector set_log_dm_lenc_rHD(1,7);   
init_vector set_log_dm_lenc_rGN(1,7);   
init_vector set_log_dm_agec_sBT(1,7);    
init_vector set_log_dm_agec_sTV(1,7);    
init_vector set_log_dm_agec_cHL(1,7);    
init_vector set_log_dm_agec_cPT(1,7);    
init_vector set_log_dm_agec_rHB(1,7);   

//// Selectivity parameters
init_vector set_selpar_A50_sBT(1,7);
init_vector set_selpar_slope_sBT(1,7);
init_vector set_selpar_A502_sBT(1,7);
init_vector set_selpar_slope2_sBT(1,7);

init_vector set_selpar_A50_sTV(1,7);
init_vector set_selpar_slope_sTV(1,7);
init_vector set_selpar_A502_sTV(1,7);
init_vector set_selpar_slope2_sTV(1,7);

init_vector set_selpar_A50b02_cHL(1,7);
init_vector set_selpar_slopeb02_cHL(1,7);
init_vector set_selpar_A50b03_cHL(1,7);
init_vector set_selpar_slopeb03_cHL(1,7);
init_vector set_selpar_A50b04_cHL(1,7);
init_vector set_selpar_slopeb04_cHL(1,7);

init_vector set_selpar_A50b02_cPT(1,7);
init_vector set_selpar_slopeb02_cPT(1,7);
init_vector set_selpar_A50b03_cPT(1,7);
init_vector set_selpar_slopeb03_cPT(1,7);
init_vector set_selpar_A50b04_cPT(1,7);
init_vector set_selpar_slopeb04_cPT(1,7);

init_vector set_selpar_A50b01_rHB(1,7);
init_vector set_selpar_slopeb01_rHB(1,7);
init_vector set_selpar_A50b02_rHB(1,7);
init_vector set_selpar_slopeb02_rHB(1,7);
init_vector set_selpar_A50b03_rHB(1,7);
init_vector set_selpar_slopeb03_rHB(1,7);
init_vector set_selpar_A50b04_rHB(1,7);
init_vector set_selpar_slopeb04_rHB(1,7);
init_vector set_selpar_A50b05_rHB(1,7);
init_vector set_selpar_slopeb05_rHB(1,7);

init_vector set_selpar_logitAge0_rHD(1,7);
init_vector set_selpar_logitAge1_rHD(1,7);
init_vector set_selpar_logitAge2_rHD(1,7);

init_vector set_selpar_A50b04_rHD(1,7);
init_vector set_selpar_slopeb04_rHD(1,7);
init_vector set_selpar_A502b04_rHD(1,7);
init_vector set_selpar_slope2b04_rHD(1,7);
init_vector set_selpar_A50b05_rHD(1,7);
init_vector set_selpar_slopeb05_rHD(1,7);
init_vector set_selpar_A502b05_rHD(1,7);
init_vector set_selpar_slope2b05_rHD(1,7);

//// Initial F
init_vector set_F_init_ratio(1,7);

//// index catchability
init_vector set_log_q_cpue_sBT(1,7);
init_vector set_log_q_cpue_sTV(1,7);
init_vector set_log_q_cpue_cHL(1,7);
init_vector set_log_q_cpue_rHB(1,7);

// q rate increase
init_vector set_q_rate(1,7);  // negative value for phase turns it off

//// average F in log space
init_vector set_log_avg_F_L_cHL(1,7);
init_vector set_log_avg_F_L_cPT(1,7);
init_vector set_log_avg_F_L_cTW(1,7);
init_vector set_log_avg_F_L_rHB(1,7);
init_vector set_log_avg_F_L_rGN(1,7);
init_vector set_log_avg_F_L_rIN(1,7);
// init_vector set_log_avg_F_D_cHD(1,7);
// init_vector set_log_avg_F_D_cPD(1,7);
init_vector set_log_avg_F_D_cGD(1,7);
init_vector set_log_avg_F_D_rHD(1,7);
init_vector set_log_avg_F_D_rGD(1,7);
init_vector set_log_avg_F_D_rID(1,7);

// ## Dev vectors
// #########################
// # lower # upper #       #
// # bound # bound # phase #
// #-------#-------#-------#
init_vector set_log_dev_F_L_cHL(1,3);
init_vector set_log_dev_F_L_cPT(1,3);
init_vector set_log_dev_F_L_cTW(1,3);
init_vector set_log_dev_F_L_rHB(1,3);
init_vector set_log_dev_F_L_rGN(1,3);
init_vector set_log_dev_F_L_rIN(1,3);
// init_vector set_log_dev_F_D_cHD(1,3);
// init_vector set_log_dev_F_D_cPD(1,3);
init_vector set_log_dev_F_D_cGD(1,3);
init_vector set_log_dev_F_D_rHD(1,3);
init_vector set_log_dev_F_D_rGD(1,3);
init_vector set_log_dev_F_D_rID(1,3);

init_vector set_log_dev_RWq(1,3);
init_vector set_log_dev_rec(1,3);
init_vector set_log_dev_Nage(1,3);

// initial set_log_dev_vals
init_vector set_log_dev_vals_F_L_cHL(styr_L_cHL,endyr_L_cHL);
init_vector set_log_dev_vals_F_L_cPT(styr_L_cPT,endyr_L_cPT);
init_vector set_log_dev_vals_F_L_cTW(styr_L_cTW,endyr_L_cTW);
init_vector set_log_dev_vals_F_L_rHB(styr_L_rHB,endyr_L_rHB);
init_vector set_log_dev_vals_F_L_rGN(styr_L_rGN,endyr_L_rGN);
init_vector set_log_dev_vals_F_L_rIN(styr_L_rGN,endyr_L_rIN);
// init_vector set_log_dev_vals_F_D_cHD(styr_D_cHD,endyr_D_cHD);
// init_vector set_log_dev_vals_F_D_cPD(styr_D_cPD,endyr_D_cPD);
init_vector set_log_dev_vals_F_D_cGD(styr_D_cGD,endyr_D_cGD);
init_vector set_log_dev_vals_F_D_rHD(styr_D_rHD,endyr_D_rHD);
init_vector set_log_dev_vals_F_D_rGD(styr_D_rGD,endyr_D_rGD);
init_vector set_log_dev_vals_F_D_rID(styr_D_rGD,endyr_D_rID);
init_vector set_log_dev_vals_rec(styr_rec_dev,endyr_rec_dev);
init_vector set_log_dev_vals_Nage(2,nages);   

//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--##
//-- BAM DATA_SECTION: likelihood weights section
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--##
init_number set_w_L;
init_number set_w_D;

init_number set_w_lenc_sBT;
init_number set_w_lenc_cHL;
init_number set_w_lenc_cPT;
init_number set_w_lenc_rHB;
init_number set_w_lenc_rHD;
init_number set_w_lenc_rGN;

init_number set_w_agec_sBT;
init_number set_w_agec_sTV;
init_number set_w_agec_cHL;
init_number set_w_agec_cPT;
init_number set_w_agec_rHB;

init_number set_w_cpue_sBT;
init_number set_w_cpue_sTV;
init_number set_w_cpue_cHL;
init_number set_w_cpue_rHB;

init_number set_w_rec;             //for fitting S-R curve
init_number set_w_rec_early;       //additional constraint on early years recruitment
init_number set_w_rec_end;         //additional constraint on ending years recruitment 

//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--##
//-- BAM DATA_SECTION: miscellaneous stuff section
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--##

// TL(mm)-weight(whole weight in kg) relationship: W=aL^b
init_number wgtpar_a;
init_number wgtpar_b;

// weight(whole weight)- fecundity (units=eggs/batch) relationship: log(y)=a+bW
init_number fecpar_a;
init_number fecpar_b;
init_number fecpar_batches; // number of batches of eggs per year
init_number fecpar_scale;   // annual fecundity scaling factor
 
init_vector obs_maturity_f(1,nages); // proportion females mature at age
init_vector obs_prop_f(1,nages);     // proportion female at age (sex ratio)

init_vector set_M(1,nages);  // natural mortality at age

init_number spawn_time_frac; // time of year (as fraction) for spawning

// catchability (q) related parameters
// Estimation phase of random walk on q (negative value turns it off)
init_int set_q_DD_phase;
// Density dependent catchability exponent, value of zero is density independent, est range is (0.1,0.9)
init_number set_q_DD_beta;
// SE of density dependent catchability exponent (0.128 provides 95% CI in range 0.5)
init_number set_q_DD_beta_se;
// Minimum age to use when computing B0_q_DD (should be age near full exploitation)
init_int set_q_DD_stage;
// Estimation phase for random walk in catchability (negative value turns it off) 
init_int set_q_RW_phase;
// Variance (sd^2) of fishery dependent random walk catchabilities (0.03 is near the sd=0.17 of Wilberg and Bence)
init_number set_q_RW_var_cHL;
init_number set_q_RW_var_rHB;

// threshold sample sizes for including length comps, age comps, respectively (set to 99999.0 if sel is fixed)
init_number minSS_lenc;
init_number minSS_agec;

// maximum sample sizes for length comps, age comps, respectively [if N>val, N=val], set high (e.g., 99999.0) to turn off
init_number maxSS_lenc;
init_number maxSS_agec;

// Projection input
init_int endyr_proj;    // Projection end year, must be later than assessment endyr
init_int styr_regs;     // apply current F until styr_regs, then the projection F
init_int Fproj_switch;  // Value to use for Fproj: 1=current F, 2=Fmsy, 3=F30, 4=F40, 5=F50, 6=Fproxy
init_number Fproj_mult; // Multiplier "c" applied to compute projection F, for example Fproj=cFmsy
init_int Rproj_switch;  // Value to use for projected recruitment: 1=mean-unbiased spawner recruit curve, 2=mean recruitment defined for SPR calcs (rec_mean), 3=mean recruitment from end of assessment (rec_mu)
int styr_proj;
 LOCAL_CALCS
   styr_proj=endyr+1;
 END_CALCS
 
 // age_error: Ageing error matrix
 // 	 columns = true ages
 //   rows = ages from age comps
 //   columns should sum to one
init_matrix age_error(1,nages,1,nages);
 
// #######Indexing integers for year(iyear), age(iage),length(ilen) ###############
int iyear;
int iage;
int ilen;
int ff;

number sqrt2pi;
number g2mt;                    //conversion of grams to metric tons 
number g2kg;                    //conversion of grams to kg   
number g2klb;                   //conversion of grams to 1000 lb   
number mt2klb;                  //conversion of metric tons to 1000 lb
number mt2lb;                   //conversion of metric tons to lb
number dzero;                   //small additive constant to prevent division by zero
number huge_number;             //huge number, to avoid irregular parameter space
number onehalf;                 //0.5

init_number end_of_data_file;

//this section MUST BE INDENTED!!!
  LOCAL_CALCS
    if(end_of_data_file!=999)
    {
        cout << "*** WARNING: Data File NOT READ CORRECTLY ****" << endl;
        exit(0);  
    }
    else
    {cout << "Data File read correctly" << endl;} 
  END_CALCS     


//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
PARAMETER_SECTION
 LOCAL_CALCS
  const double Linf_LO=set_Linf(2); const double Linf_HI=set_Linf(3); const double Linf_PH=set_Linf(4);
  const double K_LO=set_K(2); const double K_HI=set_K(3); const double K_PH=set_K(4);
  const double t0_LO=set_t0(2); const double t0_HI=set_t0(3); const double t0_PH=set_t0(4);  
  const double len_cv_LO=set_len_cv(2); const double len_cv_HI=set_len_cv(3); const double len_cv_PH=set_len_cv(4); 
     
  const double M_constant_LO=set_M_constant(2); const double M_constant_HI=set_M_constant(3); const double M_constant_PH=set_M_constant(4);        
  const double steep_LO=set_steep(2); const double steep_HI=set_steep(3); const double steep_PH=set_steep(4);
  const double log_R0_LO=set_log_R0(2); const double log_R0_HI=set_log_R0(3); const double log_R0_PH=set_log_R0(4);
  const double R_autocorr_LO=set_R_autocorr(2); const double R_autocorr_HI=set_R_autocorr(3); const double R_autocorr_PH=set_R_autocorr(4);
  const double rec_sigma_LO=set_rec_sigma(2); const double rec_sigma_HI=set_rec_sigma(3); const double rec_sigma_PH=set_rec_sigma(4);
  
  const double log_dm_sBT_lc_LO=set_log_dm_lenc_sBT(2); const double log_dm_sBT_lc_HI=set_log_dm_lenc_sBT(3); const double log_dm_sBT_lc_PH=set_log_dm_lenc_sBT(4);
  // const double log_dm_sTV_lc_LO=set_log_dm_sTV_lc(2); const double log_dm_sTV_lc_HI=set_log_dm_sTV_lc(3); const double log_dm_sTV_lc_PH=set_log_dm_sTV_lc(4);
  const double log_dm_rHB_lc_LO=set_log_dm_lenc_rHB(2); const double log_dm_rHB_lc_HI=set_log_dm_lenc_rHB(3); const double log_dm_rHB_lc_PH=set_log_dm_lenc_rHB(4);
  const double log_dm_rGN_lc_LO=set_log_dm_lenc_rGN(2); const double log_dm_rGN_lc_HI=set_log_dm_lenc_rGN(3); const double log_dm_rGN_lc_PH=set_log_dm_lenc_rGN(4);
  const double log_dm_rHD_lc_LO=set_log_dm_lenc_rHD(2); const double log_dm_rHD_lc_HI=set_log_dm_lenc_rHD(3); const double log_dm_rHD_lc_PH=set_log_dm_lenc_rHD(4);
  const double log_dm_cHL_lc_LO=set_log_dm_lenc_cHL(2); const double log_dm_cHL_lc_HI=set_log_dm_lenc_cHL(3); const double log_dm_cHL_lc_PH=set_log_dm_lenc_cHL(4);
  const double log_dm_cPT_lc_LO=set_log_dm_lenc_cPT(2); const double log_dm_cPT_lc_HI=set_log_dm_lenc_cPT(3); const double log_dm_cPT_lc_PH=set_log_dm_lenc_cPT(4);
  const double log_dm_sBT_ac_LO=set_log_dm_agec_sBT(2); const double log_dm_sBT_ac_HI=set_log_dm_agec_sBT(3); const double log_dm_sBT_ac_PH=set_log_dm_agec_sBT(4);
  const double log_dm_sTV_ac_LO=set_log_dm_agec_sTV(2); const double log_dm_sTV_ac_HI=set_log_dm_agec_sTV(3); const double log_dm_sTV_ac_PH=set_log_dm_agec_sTV(4);
  const double log_dm_rHB_ac_LO=set_log_dm_agec_rHB(2); const double log_dm_rHB_ac_HI=set_log_dm_agec_rHB(3); const double log_dm_rHB_ac_PH=set_log_dm_agec_rHB(4);
  //const double log_dm_rGN_ac_LO=set_log_dm_rGN_ac(2); const double log_dm_rGN_ac_HI=set_log_dm_rGN_ac(3); const double log_dm_rGN_ac_PH=set_log_dm_rGN_ac(4);
  const double log_dm_cHL_ac_LO=set_log_dm_agec_cHL(2); const double log_dm_cHL_ac_HI=set_log_dm_agec_cHL(3); const double log_dm_cHL_ac_PH=set_log_dm_agec_cHL(4);
  const double log_dm_cPT_ac_LO=set_log_dm_agec_cPT(2); const double log_dm_cPT_ac_HI=set_log_dm_agec_cPT(3); const double log_dm_cPT_ac_PH=set_log_dm_agec_cPT(4);
  
  const double selpar_A50_sBT_LO=set_selpar_A50_sBT(2); const double selpar_A50_sBT_HI=set_selpar_A50_sBT(3); const double selpar_A50_sBT_PH=set_selpar_A50_sBT(4);
  const double selpar_slope_sBT_LO=set_selpar_slope_sBT(2); const double selpar_slope_sBT_HI=set_selpar_slope_sBT(3); const double selpar_slope_sBT_PH=set_selpar_slope_sBT(4);
  const double selpar_A502_sBT_LO=set_selpar_A502_sBT(2); const double selpar_A502_sBT_HI=set_selpar_A502_sBT(3); const double selpar_A502_sBT_PH=set_selpar_A502_sBT(4);
  const double selpar_slope2_sBT_LO=set_selpar_slope2_sBT(2); const double selpar_slope2_sBT_HI=set_selpar_slope2_sBT(3); const double selpar_slope2_sBT_PH=set_selpar_slope2_sBT(4);

  const double selpar_A50_sTV_LO=set_selpar_A50_sTV(2); const double selpar_A50_sTV_HI=set_selpar_A50_sTV(3); const double selpar_A50_sTV_PH=set_selpar_A50_sTV(4);
  const double selpar_slope_sTV_LO=set_selpar_slope_sTV(2); const double selpar_slope_sTV_HI=set_selpar_slope_sTV(3); const double selpar_slope_sTV_PH=set_selpar_slope_sTV(4);
  const double selpar_A502_sTV_LO=set_selpar_A502_sTV(2); const double selpar_A502_sTV_HI=set_selpar_A502_sTV(3); const double selpar_A502_sTV_PH=set_selpar_A502_sTV(4);
  const double selpar_slope2_sTV_LO=set_selpar_slope2_sTV(2); const double selpar_slope2_sTV_HI=set_selpar_slope2_sTV(3); const double selpar_slope2_sTV_PH=set_selpar_slope2_sTV(4);

  const double selpar_A50b02_cHL_LO=set_selpar_A50b02_cHL(2); const double selpar_A50b02_cHL_HI=set_selpar_A50b02_cHL(3); const double selpar_A50b02_cHL_PH=set_selpar_A50b02_cHL(4);
  const double selpar_slopeb02_cHL_LO=set_selpar_slopeb02_cHL(2); const double selpar_slopeb02_cHL_HI=set_selpar_slopeb02_cHL(3); const double selpar_slopeb02_cHL_PH=set_selpar_slopeb02_cHL(4);
  const double selpar_A50b03_cHL_LO=set_selpar_A50b03_cHL(2); const double selpar_A50b03_cHL_HI=set_selpar_A50b03_cHL(3); const double selpar_A50b03_cHL_PH=set_selpar_A50b03_cHL(4);
  const double selpar_slopeb03_cHL_LO=set_selpar_slopeb03_cHL(2); const double selpar_slopeb03_cHL_HI=set_selpar_slopeb03_cHL(3); const double selpar_slopeb03_cHL_PH=set_selpar_slopeb03_cHL(4);
  const double selpar_A50b04_cHL_LO=set_selpar_A50b04_cHL(2); const double selpar_A50b04_cHL_HI=set_selpar_A50b04_cHL(3); const double selpar_A50b04_cHL_PH=set_selpar_A50b04_cHL(4);
  const double selpar_slopeb04_cHL_LO=set_selpar_slopeb04_cHL(2); const double selpar_slopeb04_cHL_HI=set_selpar_slopeb04_cHL(3); const double selpar_slopeb04_cHL_PH=set_selpar_slopeb04_cHL(4);
  
  const double selpar_A50b02_cPT_LO=set_selpar_A50b02_cPT(2); const double selpar_A50b02_cPT_HI=set_selpar_A50b02_cPT(3); const double selpar_A50b02_cPT_PH=set_selpar_A50b02_cPT(4);
  const double selpar_slopeb02_cPT_LO=set_selpar_slopeb02_cPT(2); const double selpar_slopeb02_cPT_HI=set_selpar_slopeb02_cPT(3); const double selpar_slopeb02_cPT_PH=set_selpar_slopeb02_cPT(4);
  const double selpar_A50b03_cPT_LO=set_selpar_A50b03_cPT(2); const double selpar_A50b03_cPT_HI=set_selpar_A50b03_cPT(3); const double selpar_A50b03_cPT_PH=set_selpar_A50b03_cPT(4);
  const double selpar_slopeb03_cPT_LO=set_selpar_slopeb03_cPT(2); const double selpar_slopeb03_cPT_HI=set_selpar_slopeb03_cPT(3); const double selpar_slopeb03_cPT_PH=set_selpar_slopeb03_cPT(4);
  const double selpar_A50b04_cPT_LO=set_selpar_A50b04_cPT(2); const double selpar_A50b04_cPT_HI=set_selpar_A50b04_cPT(3); const double selpar_A50b04_cPT_PH=set_selpar_A50b04_cPT(4);
  const double selpar_slopeb04_cPT_LO=set_selpar_slopeb04_cPT(2); const double selpar_slopeb04_cPT_HI=set_selpar_slopeb04_cPT(3); const double selpar_slopeb04_cPT_PH=set_selpar_slopeb04_cPT(4);
  
  const double selpar_A50b01_rHB_LO=set_selpar_A50b01_rHB(2); const double selpar_A50b01_rHB_HI=set_selpar_A50b01_rHB(3); const double selpar_A50b01_rHB_PH=set_selpar_A50b01_rHB(4);
  const double selpar_slopeb01_rHB_LO=set_selpar_slopeb01_rHB(2); const double selpar_slopeb01_rHB_HI=set_selpar_slopeb01_rHB(3); const double selpar_slopeb01_rHB_PH=set_selpar_slopeb01_rHB(4);
  const double selpar_A50b02_rHB_LO=set_selpar_A50b02_rHB(2); const double selpar_A50b02_rHB_HI=set_selpar_A50b02_rHB(3); const double selpar_A50b02_rHB_PH=set_selpar_A50b02_rHB(4);
  const double selpar_slopeb02_rHB_LO=set_selpar_slopeb02_rHB(2); const double selpar_slopeb02_rHB_HI=set_selpar_slopeb02_rHB(3); const double selpar_slopeb02_rHB_PH=set_selpar_slopeb02_rHB(4);
  const double selpar_A50b03_rHB_LO=set_selpar_A50b03_rHB(2); const double selpar_A50b03_rHB_HI=set_selpar_A50b03_rHB(3); const double selpar_A50b03_rHB_PH=set_selpar_A50b03_rHB(4);
  const double selpar_slopeb03_rHB_LO=set_selpar_slopeb03_rHB(2); const double selpar_slopeb03_rHB_HI=set_selpar_slopeb03_rHB(3); const double selpar_slopeb03_rHB_PH=set_selpar_slopeb03_rHB(4);
  const double selpar_A50b04_rHB_LO=set_selpar_A50b04_rHB(2); const double selpar_A50b04_rHB_HI=set_selpar_A50b04_rHB(3); const double selpar_A50b04_rHB_PH=set_selpar_A50b04_rHB(4);
  const double selpar_slopeb04_rHB_LO=set_selpar_slopeb04_rHB(2); const double selpar_slopeb04_rHB_HI=set_selpar_slopeb04_rHB(3); const double selpar_slopeb04_rHB_PH=set_selpar_slopeb04_rHB(4);
  const double selpar_A50b05_rHB_LO=set_selpar_A50b05_rHB(2); const double selpar_A50b05_rHB_HI=set_selpar_A50b05_rHB(3); const double selpar_A50b05_rHB_PH=set_selpar_A50b05_rHB(4);
  const double selpar_slopeb05_rHB_LO=set_selpar_slopeb05_rHB(2); const double selpar_slopeb05_rHB_HI=set_selpar_slopeb05_rHB(3); const double selpar_slopeb05_rHB_PH=set_selpar_slopeb05_rHB(4);
  
  const double selpar_logitAge0_rHD_LO=set_selpar_logitAge0_rHD(2); const double selpar_logitAge0_rHD_HI=set_selpar_logitAge0_rHD(3); const double selpar_logitAge0_rHD_PH=set_selpar_logitAge0_rHD(4);
  const double selpar_logitAge1_rHD_LO=set_selpar_logitAge1_rHD(2); const double selpar_logitAge1_rHD_HI=set_selpar_logitAge1_rHD(3); const double selpar_logitAge1_rHD_PH=set_selpar_logitAge1_rHD(4);
  const double selpar_logitAge2_rHD_LO=set_selpar_logitAge2_rHD(2); const double selpar_logitAge2_rHD_HI=set_selpar_logitAge2_rHD(3); const double selpar_logitAge2_rHD_PH=set_selpar_logitAge2_rHD(4);
  
  const double selpar_A50b04_rHD_LO=set_selpar_A50b04_rHD(2); 
  const double selpar_A50b04_rHD_HI=set_selpar_A50b04_rHD(3); 
  const double selpar_A50b04_rHD_PH=set_selpar_A50b04_rHD(4);
  const double selpar_slopeb04_rHD_LO=set_selpar_slopeb04_rHD(2); 
  const double selpar_slopeb04_rHD_HI=set_selpar_slopeb04_rHD(3); 
  const double selpar_slopeb04_rHD_PH=set_selpar_slopeb04_rHD(4);
  const double selpar_A502b04_rHD_LO=set_selpar_A502b04_rHD(2); 
  const double selpar_A502b04_rHD_HI=set_selpar_A502b04_rHD(3); 
  const double selpar_A502b04_rHD_PH=set_selpar_A502b04_rHD(4);
  const double selpar_slope2b04_rHD_LO=set_selpar_slope2b04_rHD(2); 
  const double selpar_slope2b04_rHD_HI=set_selpar_slope2b04_rHD(3); 
  const double selpar_slope2b04_rHD_PH=set_selpar_slope2b04_rHD(4);
  const double selpar_A50b05_rHD_LO=set_selpar_A50b05_rHD(2); 
  const double selpar_A50b05_rHD_HI=set_selpar_A50b05_rHD(3); 
  const double selpar_A50b05_rHD_PH=set_selpar_A50b05_rHD(4);
  const double selpar_slopeb05_rHD_LO=set_selpar_slopeb05_rHD(2); 
  const double selpar_slopeb05_rHD_HI=set_selpar_slopeb05_rHD(3); 
  const double selpar_slopeb05_rHD_PH=set_selpar_slopeb05_rHD(4);
  const double selpar_A502b05_rHD_LO=set_selpar_A502b05_rHD(2); 
  const double selpar_A502b05_rHD_HI=set_selpar_A502b05_rHD(3); 
  const double selpar_A502b05_rHD_PH=set_selpar_A502b05_rHD(4);
  const double selpar_slope2b05_rHD_LO=set_selpar_slope2b05_rHD(2); 
  const double selpar_slope2b05_rHD_HI=set_selpar_slope2b05_rHD(3); 
  const double selpar_slope2b05_rHD_PH=set_selpar_slope2b05_rHD(4);
  
  const double log_q_sBT_LO=set_log_q_cpue_sBT(2); const double log_q_sBT_HI=set_log_q_cpue_sBT(3); const double log_q_sBT_PH=set_log_q_cpue_sBT(4);
  const double log_q_rHB_LO=set_log_q_cpue_rHB(2); const double log_q_rHB_HI=set_log_q_cpue_rHB(3); const double log_q_rHB_PH=set_log_q_cpue_rHB(4);
  const double log_q_sTV_LO=set_log_q_cpue_sTV(2); const double log_q_sTV_HI=set_log_q_cpue_sTV(3); const double log_q_sTV_PH=set_log_q_cpue_sTV(4);
  const double log_q_cHL_LO=set_log_q_cpue_cHL(2); const double log_q_cHL_HI=set_log_q_cpue_cHL(3); const double log_q_cHL_PH=set_log_q_cpue_cHL(4);
  //const double log_q_rHD_LO=set_logq_rHD(2); const double log_q_rHD_HI=set_logq_rHD(3); const double log_q_rHD_PH=set_logq_rHD(4);
  
  const double log_avg_F_cHL_LO=set_log_avg_F_L_cHL(2); const double log_avg_F_cHL_HI=set_log_avg_F_L_cHL(3); const double log_avg_F_cHL_PH=set_log_avg_F_L_cHL(4);
  const double log_avg_F_cPT_LO=set_log_avg_F_L_cPT(2); const double log_avg_F_cPT_HI=set_log_avg_F_L_cPT(3); const double log_avg_F_cPT_PH=set_log_avg_F_L_cPT(4);
  const double log_avg_F_cTW_LO=set_log_avg_F_L_cTW(2); const double log_avg_F_cTW_HI=set_log_avg_F_L_cTW(3); const double log_avg_F_cTW_PH=set_log_avg_F_L_cTW(4);
  const double log_avg_F_rHB_LO=set_log_avg_F_L_rHB(2); const double log_avg_F_rHB_HI=set_log_avg_F_L_rHB(3); const double log_avg_F_rHB_PH=set_log_avg_F_L_rHB(4); 
  const double log_avg_F_rGN_LO=set_log_avg_F_L_rGN(2); const double log_avg_F_rGN_HI=set_log_avg_F_L_rGN(3); const double log_avg_F_rGN_PH=set_log_avg_F_L_rGN(4);
  const double log_avg_F_rIN_LO=set_log_avg_F_L_rIN(2); const double log_avg_F_rIN_HI=set_log_avg_F_L_rIN(3); const double log_avg_F_rIN_PH=set_log_avg_F_L_rIN(4);
  
  // const double LO_log_avg_F_D_cHD=set_log_avg_F_D_cHD(2); const double HI_log_avg_F_D_cHD=set_log_avg_F_D_cHD(3); const double PH_log_avg_F_D_cHD=set_log_avg_F_D_cHD(4);
  // const double LO_log_avg_F_D_cPD=set_log_avg_F_D_cPD(2); const double HI_log_avg_F_D_cPD=set_log_avg_F_D_cPD(3); const double PH_log_avg_F_D_cPD=set_log_avg_F_D_cPD(4);
  const double LO_log_avg_F_D_cGD=set_log_avg_F_D_cGD(2); const double HI_log_avg_F_D_cGD=set_log_avg_F_D_cGD(3); const double PH_log_avg_F_D_cGD=set_log_avg_F_D_cGD(4);
  const double LO_log_avg_F_D_rHD=set_log_avg_F_D_rHD(2); const double HI_log_avg_F_D_rHD=set_log_avg_F_D_rHD(3); const double PH_log_avg_F_D_rHD=set_log_avg_F_D_rHD(4); 
  const double LO_log_avg_F_D_rGD=set_log_avg_F_D_rGD(2); const double HI_log_avg_F_D_rGD=set_log_avg_F_D_rGD(3); const double PH_log_avg_F_D_rGD=set_log_avg_F_D_rGD(4);
  const double LO_log_avg_F_D_rID=set_log_avg_F_D_rID(2); const double HI_log_avg_F_D_rID=set_log_avg_F_D_rID(3); const double PH_log_avg_F_D_rID=set_log_avg_F_D_rID(4);
  
  //-dev vectors-----------------------------------------------------------------------------------------------------------  
  const double log_dev_F_cHL_LO=set_log_dev_F_L_cHL(1); const double log_dev_F_cHL_HI=set_log_dev_F_L_cHL(2); const double log_dev_F_cHL_PH=set_log_dev_F_L_cHL(3);  
  const double log_dev_F_cPT_LO=set_log_dev_F_L_cPT(1); const double log_dev_F_cPT_HI=set_log_dev_F_L_cPT(2); const double log_dev_F_cPT_PH=set_log_dev_F_L_cPT(3);   
  const double log_dev_F_cTW_LO=set_log_dev_F_L_cTW(1); const double log_dev_F_cTW_HI=set_log_dev_F_L_cTW(2); const double log_dev_F_cTW_PH=set_log_dev_F_L_cTW(3);   
  const double log_dev_F_rHB_LO=set_log_dev_F_L_rHB(1); const double log_dev_F_rHB_HI=set_log_dev_F_L_rHB(2); const double log_dev_F_rHB_PH=set_log_dev_F_L_rHB(3);   
  const double log_dev_F_rGN_LO=set_log_dev_F_L_rGN(1); const double log_dev_F_rGN_HI=set_log_dev_F_L_rGN(2); const double log_dev_F_rGN_PH=set_log_dev_F_L_rGN(3);
  const double log_dev_F_rIN_LO=set_log_dev_F_L_rIN(1); const double log_dev_F_rIN_HI=set_log_dev_F_L_rIN(2); const double log_dev_F_rIN_PH=set_log_dev_F_L_rIN(3);   
  
  // const double log_dev_F_cHD_LO=set_log_dev_F_D_cHD(1); const double log_dev_F_cHD_HI=set_log_dev_F_D_cHD(2); const double log_dev_F_cHD_PH=set_log_dev_F_D_cHD(3);
  // const double log_dev_F_cPD_LO=set_log_dev_F_D_cPD(1); const double log_dev_F_cPD_HI=set_log_dev_F_D_cPD(2); const double log_dev_F_cPD_PH=set_log_dev_F_D_cPD(3);
  const double log_dev_F_cGD_LO=set_log_dev_F_D_cGD(1); const double log_dev_F_cGD_HI=set_log_dev_F_D_cGD(2); const double log_dev_F_cGD_PH=set_log_dev_F_D_cGD(3);   
  const double log_dev_F_rHD_LO=set_log_dev_F_D_rHD(1); const double log_dev_F_rHD_HI=set_log_dev_F_D_rHD(2); const double log_dev_F_rHD_PH=set_log_dev_F_D_rHD(3);   
  const double log_dev_F_rGD_LO=set_log_dev_F_D_rGD(1); const double log_dev_F_rGD_HI=set_log_dev_F_D_rGD(2); const double log_dev_F_rGD_PH=set_log_dev_F_D_rGD(3);
  const double log_dev_F_rID_LO=set_log_dev_F_D_rID(1); const double log_dev_F_rID_HI=set_log_dev_F_D_rID(2); const double log_dev_F_rID_PH=set_log_dev_F_D_rID(3);  
  const double F_init_ratio_LO=set_F_init_ratio(2); const double F_init_ratio_HI=set_F_init_ratio(3); const double F_init_ratio_PH=set_F_init_ratio(4); 

  const double log_RWq_LO=set_log_dev_RWq(1); const double log_RWq_HI=set_log_dev_RWq(2); const double log_RWq_PH=set_log_dev_RWq(3);  
  
  const double log_rec_dev_LO=set_log_dev_rec(1); const double log_rec_dev_HI=set_log_dev_rec(2); const double log_rec_dev_PH=set_log_dev_rec(3);          
  const double log_Nage_dev_LO=set_log_dev_Nage(1); const double log_Nage_dev_HI=set_log_dev_Nage(2); const double log_Nage_dev_PH=set_log_dev_Nage(3);

  const double q_rate_LO=set_q_rate(2); const double q_rate_HI=set_q_rate(3);  const double q_rate_PH=set_q_rate(4);
  
 END_CALCS
 
////--------------Growth---------------------------------------------------------------------------
 
  //Population growth parms and conversions
  init_bounded_number Linf(Linf_LO,Linf_HI,Linf_PH);
  init_bounded_number K(K_LO,K_HI,K_PH);
  init_bounded_number t0(t0_LO,t0_HI,t0_PH);
  init_bounded_number len_cv_val(len_cv_LO,len_cv_HI,len_cv_PH);  
  vector Linf_out(1,8);
  vector K_out(1,8);
  vector t0_out(1,8);
  vector len_cv_val_out(1,8);
  
  vector meanlen_TL(1,nages);   //mean Total length (mm) at age
  
  vector wgt_g(1,nages);        //whole wgt in g
  vector wgt_kg(1,nages);       //whole wgt in kg
  vector wgt_mt(1,nages);       //whole wgt in mt
  vector wgt_klb(1,nages);      //whole wgt in 1000 lb
  vector wgt_lb(1,nages);       //whole wgt in lb
  vector fecundity(1,nages);    //fecundity at age, perhaps scaled  

  //mean length and weight at age of landings and discards by fleet (may differ from population mean)
  matrix len_cHL_mm(styr,endyr,1,nages);  
  matrix wgt_cHL_klb(styr,endyr,1,nages);
  matrix len_cPT_mm(styr,endyr,1,nages);
  matrix wgt_cPT_klb(styr,endyr,1,nages);
  matrix len_cTW_mm(styr,endyr,1,nages);
  matrix wgt_cTW_klb(styr,endyr,1,nages);
  matrix len_rHB_mm(styr,endyr,1,nages);  
  matrix wgt_rHB_klb(styr,endyr,1,nages);
  matrix len_rGN_mm(styr,endyr,1,nages);   
  matrix wgt_rGN_klb(styr,endyr,1,nages);
  matrix len_rIN_mm(styr,endyr,1,nages);   
  matrix wgt_rIN_klb(styr,endyr,1,nages);

  // matrix len_cHD_mm(styr,endyr,1,nages);
  // matrix wgt_cHD_klb(styr,endyr,1,nages);  
  // matrix len_cPD_mm(styr,endyr,1,nages);
  // matrix wgt_cPD_klb(styr,endyr,1,nages);  
  matrix len_cGD_mm(styr,endyr,1,nages);
  matrix wgt_cGD_klb(styr,endyr,1,nages);
  matrix len_rHD_mm(styr,endyr,1,nages);
  matrix wgt_rHD_klb(styr,endyr,1,nages);
  matrix len_rGD_mm(styr,endyr,1,nages);
  matrix wgt_rGD_klb(styr,endyr,1,nages);
  matrix len_rID_mm(styr,endyr,1,nages);
  matrix wgt_rID_klb(styr,endyr,1,nages);  
 
  matrix lenprob(1,nages,1,nlenbins);           //distn of size at age (age-length key, 1 cm bins) in population
  
  number zscore_len;                            //standardized normal values used for computing lenprob
  vector cprob_lenvec(1,nlenbins);              //cumulative probabilities used for computing lenprob
  number zscore_lzero;                          //standardized normal values for length = 0
  number cprob_lzero;                           //length probability mass below zero, used for computing lenprob

  vector len_sd(1,nages);
  vector len_cv(1,nages); //for fishgraph 
  
//----Predicted length compositions
  matrix pred_lenc_sBT(1,nyr_lenc_sBT,1,nlenbins);
  matrix pred_lenc_cHL(1,nyr_lenc_cHL,1,nlenbins);
  matrix pred_lenc_cPT(1,nyr_lenc_cPT,1,nlenbins);
  matrix pred_lenc_rHB(1,nyr_lenc_rHB,1,nlenbins);
  matrix pred_lenc_rHD(1,nyr_lenc_rHD,1,nlenbins);
  matrix pred_lenc_rGN(1,nyr_lenc_rGN,1,nlenbins);  

// Predicted age compositions not including age error      
  matrix predNoErr_agec_sBT(1,nyr_agec_sBT,1,nages);
  matrix predNoErr_agec_sTV(1,nyr_agec_sTV,1,nages);  
  matrix predNoErr_agec_cHL(1,nyr_agec_cHL,1,nages);
  matrix predNoErr_agec_cPT(1,nyr_agec_cPT,1,nages);  
  matrix predNoErr_agec_rHB(1,nyr_agec_rHB,1,nages);
  
// Predicted age compositions including age error      
  matrix pred_agec_sBT(1,nyr_agec_sBT,1,nages);      
  matrix pred_agec_sTV(1,nyr_agec_sTV,1,nages);
  matrix pred_agec_cHL(1,nyr_agec_cHL,1,nages);
  matrix pred_agec_cPT(1,nyr_agec_cPT,1,nages);
  matrix pred_agec_rHB(1,nyr_agec_rHB,1,nages);

 
//effective sample size applied in multinomial distributions
  vector nsamp_allyr_lenc_sBT(styr,endyr);
  vector nsamp_allyr_lenc_cHL(styr,endyr);
  vector nsamp_allyr_lenc_cPT(styr,endyr);
  vector nsamp_allyr_lenc_rHB(styr,endyr);
  vector nsamp_allyr_lenc_rHD(styr,endyr);
  vector nsamp_allyr_lenc_rGN(styr,endyr);
  vector nsamp_allyr_agec_sBT(styr,endyr);
  vector nsamp_allyr_agec_sTV(styr,endyr);
  vector nsamp_allyr_agec_cHL(styr,endyr);
  vector nsamp_allyr_agec_cPT(styr,endyr);  
  vector nsamp_allyr_agec_rHB(styr,endyr);
  //vector nsamp_allyr_agec_rGN(styr,endyr);

//Nfish used in MCB analysis (not used in fitting)
  vector nfish_allyr_lenc_sBT(styr,endyr);
  vector nfish_allyr_lenc_cHL(styr,endyr);
  vector nfish_allyr_lenc_cPT(styr,endyr);
  vector nfish_allyr_lenc_rHB(styr,endyr);
  vector nfish_allyr_lenc_rHD(styr,endyr);
  vector nfish_allyr_lenc_rGN(styr,endyr);
  vector nfish_allyr_agec_sBT(styr,endyr);
  vector nfish_allyr_agec_sTV(styr,endyr);
  vector nfish_allyr_agec_cHL(styr,endyr);
  vector nfish_allyr_agec_cPT(styr,endyr);  
  vector nfish_allyr_agec_rHB(styr,endyr);
  //vector nfish_allyr_agec_rGN(styr,endyr);
  
//Computed effective sample size for output (not used in fitting)
  vector out_neff_allyr_lenc_sBT(styr,endyr);
  vector out_neff_allyr_lenc_cHL(styr,endyr);
  vector out_neff_allyr_lenc_cPT(styr,endyr);
  vector out_neff_allyr_lenc_rHB(styr,endyr);
  vector out_neff_allyr_lenc_rHD(styr,endyr);
  vector out_neff_allyr_lenc_rGN(styr,endyr);
  vector out_neff_allyr_agec_sBT(styr,endyr);
  vector out_neff_allyr_agec_sTV(styr,endyr);
  vector out_neff_allyr_agec_cHL(styr,endyr);
  vector out_neff_allyr_agec_cPT(styr,endyr);  
  vector out_neff_allyr_agec_rHB(styr,endyr);
  //vector out_neff_allyr_agec_rGN(styr,endyr);


//-----Population-----------------------------------------------------------------------------------
  matrix N(styr,endyr,1,nages);             //Population numbers by year and age at start of yr
  matrix N_mdyr(styr,endyr,1,nages);        //Population numbers by year and age at mdpt of yr: used for comps and cpue
  matrix N_spawn(styr,endyr,1,nages);       //Population numbers by year and age at peaking spawning: used for SSB  
  init_bounded_vector log_dev_Nage(2,nages,log_Nage_dev_LO,log_Nage_dev_HI,log_Nage_dev_PH);
  //vector log_dev_Nage(2,nages);
  vector log_Nage_dev_output(1,nages);      //used in output. equals zero for first age
  matrix B(styr,endyr,1,nages);             //Population biomass by year and age at start of yr
  vector totB(styr,endyr);                  //Total biomass by year
  vector totN(styr,endyr);                  //Total abundance by year
  vector SSB(styr,endyr);                   //Total spawning biomass by year (scaled popn fecundity)
  vector MatFemB(styr,endyr);               //Total spawning biomass by year (total mature female biomass)  
  vector rec(styr,endyr);                   //Recruits by year
  vector prop_f(1,nages);                   //Proportion female by age
  vector maturity_f(1,nages);               //Proportion of female mature at age
  vector reprod(1,nages);                   //vector used to compute spawning biomass (scaled popn fecundity)
  vector reprod2(1,nages);                  //vector used to compute mature female biomass 

////---Stock-Recruit Function (Beverton-Holt, steepness parameterization)----------
  init_bounded_number log_R0(log_R0_LO,log_R0_HI,log_R0_PH);        //log(virgin Recruitment)
  vector log_R0_out(1,8);
  number R0;                                  //virgin recruitment
  init_bounded_number steep(steep_LO,steep_HI,steep_PH); //steepness
  vector steep_out(1,8);
  init_bounded_number rec_sigma(rec_sigma_LO,rec_sigma_HI,rec_sigma_PH);  //sd recruitment residuals  
  vector rec_sigma_out(1,8);
  init_bounded_number R_autocorr(R_autocorr_LO,R_autocorr_HI,R_autocorr_PH);  //autocorrelation in SR  
  vector R_autocorr_out(1,8);
  
  number rec_sigma_sq;                        //square of rec_sigma      
  number rec_sigma_sqd2;                      //square of rec_sigma divided by two
  number rec_logL_add;                        //additive term in -logL term   
    
  init_bounded_dev_vector log_dev_rec(styr_rec_dev,endyr_rec_dev,log_rec_dev_LO,log_rec_dev_HI,log_rec_dev_PH);
  vector log_rec_dev_output(styr_rec_dev,endyr);             //used in output. equals zero except for yrs in log_dev_rec
  vector log_rec_dev_out(styr_rec_dev,endyr_rec_dev);  //used in output for bound checking

  number var_rec_dev;                                //variance of log recruitment deviations, from yrs with unconstrainted S-R(XXXX-XXXX)
  number sigma_rec_dev;                              //sample SD of log residuals (may not equal rec_sigma 
  number BiasCor;                               //Bias correction in equilibrium recruits
  //number R_autocorr;
  number S0;                                    //equal to spr_F0*R0 = virgin SSB
  number B0;                                    //equal to bpr_F0*R0 = virgin B  
  number R1;                                    //Recruits in styr
  number R_virgin;                              //unfished recruitment with bias correction
  vector SdS0(styr,endyr);                      //SSB / virgin SSB

  init_bounded_number log_dm_lenc_sBT(log_dm_sBT_lc_LO,log_dm_sBT_lc_HI,log_dm_sBT_lc_PH);
  // init_bounded_number log_dm_sTV_lc(log_dm_sTV_lc_LO,log_dm_sTV_lc_HI,log_dm_sTV_lc_PH);
  init_bounded_number log_dm_lenc_cHL(log_dm_cHL_lc_LO,log_dm_cHL_lc_HI,log_dm_cHL_lc_PH);
  init_bounded_number log_dm_lenc_cPT(log_dm_cPT_lc_LO,log_dm_cPT_lc_HI,log_dm_cPT_lc_PH);
  init_bounded_number log_dm_lenc_rHB(log_dm_rHB_lc_LO,log_dm_rHB_lc_HI,log_dm_rHB_lc_PH);
  init_bounded_number log_dm_lenc_rHD(log_dm_rHD_lc_LO,log_dm_rHD_lc_HI,log_dm_rHD_lc_PH);
  init_bounded_number log_dm_lenc_rGN(log_dm_rGN_lc_LO,log_dm_rGN_lc_HI,log_dm_rGN_lc_PH);
  init_bounded_number log_dm_agec_sBT(log_dm_sBT_ac_LO,log_dm_sBT_ac_HI,log_dm_sBT_ac_PH);
  init_bounded_number log_dm_agec_sTV(log_dm_sTV_ac_LO,log_dm_sTV_ac_HI,log_dm_sTV_ac_PH);
  init_bounded_number log_dm_agec_cHL(log_dm_cHL_ac_LO,log_dm_cHL_ac_HI,log_dm_cHL_ac_PH);
  init_bounded_number log_dm_agec_cPT(log_dm_cPT_ac_LO,log_dm_cPT_ac_HI,log_dm_cPT_ac_PH);
  init_bounded_number log_dm_agec_rHB(log_dm_rHB_ac_LO,log_dm_rHB_ac_HI,log_dm_rHB_ac_PH);
  //init_bounded_number log_dm_rGN_ac(log_dm_rGN_ac_LO,log_dm_rGN_ac_HI,log_dm_rGN_ac_PH);
  vector log_dm_sBT_lc_out(1,8);
  // vector log_dm_sTV_lc_out(1,8);
  vector log_dm_cHL_lc_out(1,8);
  vector log_dm_cPT_lc_out(1,8);
  vector log_dm_rHB_lc_out(1,8);
  vector log_dm_rGN_lc_out(1,8);
  vector log_dm_rHD_lc_out(1,8);
  vector log_dm_sBT_ac_out(1,8);
  vector log_dm_sTV_ac_out(1,8);
  vector log_dm_cHL_ac_out(1,8);
  vector log_dm_cPT_ac_out(1,8);
  vector log_dm_rHB_ac_out(1,8);
  //vector log_dm_rGN_ac_out(1,8);
  
//-------------------------------------------------------------------------------------------------------
////---Selectivity---------------------------------------------------- 

// sBT
  matrix sel_sBT(styr,endyr,1,nages);
  vector sel_sBT_b01(1,nages);
  vector sel_sBT_b02(1,nages);
  vector sel_sBT_b03(1,nages);
  vector sel_sBT_b04(1,nages);
  vector sel_sBT_b05(1,nages); 
  init_bounded_number selpar_A50_sBT(selpar_A50_sBT_LO,selpar_A50_sBT_HI,selpar_A50_sBT_PH);
  init_bounded_number selpar_slope_sBT(selpar_slope_sBT_LO,selpar_slope_sBT_HI,selpar_slope_sBT_PH);
  init_bounded_number selpar_A502_sBT(selpar_A502_sBT_LO,selpar_A502_sBT_HI,selpar_A502_sBT_PH);
  init_bounded_number selpar_slope2_sBT(selpar_slope2_sBT_LO,selpar_slope2_sBT_HI,selpar_slope2_sBT_PH);
  vector selpar_A50_sBT_out(1,8);
  vector selpar_slope_sBT_out(1,8);
  vector selpar_A502_sBT_out(1,8);
  vector selpar_slope2_sBT_out(1,8);
  
// sTV
  matrix sel_sTV(styr,endyr,1,nages);
  vector sel_sTV_b01(1,nages);
  vector sel_sTV_b02(1,nages);
  vector sel_sTV_b03(1,nages);
  vector sel_sTV_b04(1,nages);
  vector sel_sTV_b05(1,nages);  
  init_bounded_number selpar_A50_sTV(selpar_A50_sTV_LO,selpar_A50_sTV_HI,selpar_A50_sTV_PH);
  init_bounded_number selpar_slope_sTV(selpar_slope_sTV_LO,selpar_slope_sTV_HI,selpar_slope_sTV_PH);
  init_bounded_number selpar_A502_sTV(selpar_A502_sTV_LO,selpar_A502_sTV_HI,selpar_A502_sTV_PH);
  init_bounded_number selpar_slope2_sTV(selpar_slope2_sTV_LO,selpar_slope2_sTV_HI,selpar_slope2_sTV_PH);
  vector selpar_A50_sTV_out(1,8);
  vector selpar_slope_sTV_out(1,8);
  vector selpar_A502_sTV_out(1,8);
  vector selpar_slope2_sTV_out(1,8);
  
// cHL
  matrix sel_cHL(styr,endyr,1,nages);  
  vector sel_cHL_b01(1,nages);
  vector sel_cHL_b02(1,nages);
  vector sel_cHL_b03(1,nages);
  vector sel_cHL_b04(1,nages);
  vector sel_cHL_b05(1,nages); 
  
  init_bounded_number selpar_A50b02_cHL(selpar_A50b02_cHL_LO,selpar_A50b02_cHL_HI,selpar_A50b02_cHL_PH);
  init_bounded_number selpar_slopeb02_cHL(selpar_slopeb02_cHL_LO,selpar_slopeb02_cHL_HI,selpar_slopeb02_cHL_PH);
  init_bounded_number selpar_A50b03_cHL(selpar_A50b03_cHL_LO,selpar_A50b03_cHL_HI,selpar_A50b03_cHL_PH);
  init_bounded_number selpar_slopeb03_cHL(selpar_slopeb03_cHL_LO,selpar_slopeb03_cHL_HI,selpar_slopeb03_cHL_PH);
  init_bounded_number selpar_A50b04_cHL(selpar_A50b04_cHL_LO,selpar_A50b04_cHL_HI,selpar_A50b04_cHL_PH);
  init_bounded_number selpar_slopeb04_cHL(selpar_slopeb04_cHL_LO,selpar_slopeb04_cHL_HI,selpar_slopeb04_cHL_PH);
 
  vector selpar_A50b02_cHL_out(1,8);
  vector selpar_slopeb02_cHL_out(1,8);
  vector selpar_A50b03_cHL_out(1,8);
  vector selpar_slopeb03_cHL_out(1,8);
  vector selpar_A50b04_cHL_out(1,8);
  vector selpar_slopeb04_cHL_out(1,8);
     
// cGD
  matrix sel_cGD(styr,endyr,1,nages); 
  vector sel_cGD_b01(1,nages);
  vector sel_cGD_b02(1,nages);
  vector sel_cGD_b03(1,nages);
  vector sel_cGD_b04(1,nages);
  vector sel_cGD_b05(1,nages);   
 
//values used for weighting selex and avg weights of discards during yrs with quotas
  number Dopen_cHL4; number Dclosed_cHL4; number Lopen_cHL4;  
  number Dopen_cPT4; number Dclosed_cPT4; number Lopen_cPT4;
  number D_sum_cHLcPT4; 
  number Dprop_cGN_sel_D4; number Dprop_cGN_sel_cHL4; number Dprop_cGN_sel_cPT4;
  // number Dopen_cHL5; number Dclosed_cHL5; number Lopen_cHL5;  
  // number Dopen_cPT5; number Dclosed_cPT5; number Lopen_cPT5;
  // number D_sum_cHLcPT5; 
  // number Dprop_cGN_sel_D5; number Dprop_cGN_sel_cHL5; number Dprop_cGN_sel_cPT5;

// cPT        
  matrix sel_cPT(styr,endyr,1,nages); 
  vector sel_cPT_b01(1,nages);
  vector sel_cPT_b02(1,nages);
  vector sel_cPT_b03(1,nages);
  vector sel_cPT_b04(1,nages);
  vector sel_cPT_b05(1,nages);  
  
  init_bounded_number selpar_A50b02_cPT(selpar_A50b02_cPT_LO,selpar_A50b02_cPT_HI,selpar_A50b02_cPT_PH);
  init_bounded_number selpar_slopeb02_cPT(selpar_slopeb02_cPT_LO,selpar_slopeb02_cPT_HI,selpar_slopeb02_cPT_PH);
  init_bounded_number selpar_A50b03_cPT(selpar_A50b03_cPT_LO,selpar_A50b03_cPT_HI,selpar_A50b03_cPT_PH);
  init_bounded_number selpar_slopeb03_cPT(selpar_slopeb03_cPT_LO,selpar_slopeb03_cPT_HI,selpar_slopeb03_cPT_PH);
  init_bounded_number selpar_A50b04_cPT(selpar_A50b04_cPT_LO,selpar_A50b04_cPT_HI,selpar_A50b04_cPT_PH);
  init_bounded_number selpar_slopeb04_cPT(selpar_slopeb04_cPT_LO,selpar_slopeb04_cPT_HI,selpar_slopeb04_cPT_PH);
 
  vector selpar_A50b02_cPT_out(1,8);
  vector selpar_slopeb02_cPT_out(1,8);
  vector selpar_A50b03_cPT_out(1,8);
  vector selpar_slopeb03_cPT_out(1,8);
  vector selpar_A50b04_cPT_out(1,8);
  vector selpar_slopeb04_cPT_out(1,8);

// cTW        
  matrix sel_cTW(styr,endyr,1,nages);  //mirrors cGN pot sel
  vector sel_cTW_b01(1,nages);
  vector sel_cTW_b02(1,nages);
  vector sel_cTW_b03(1,nages);
  vector sel_cTW_b04(1,nages);
  vector sel_cTW_b05(1,nages); 
  
// rHB
  matrix sel_rHB(styr,endyr,1,nages);  
  vector sel_rHB_b01(1,nages);
  vector sel_rHB_b02(1,nages);
  vector sel_rHB_b03(1,nages);
  vector sel_rHB_b04(1,nages);
  vector sel_rHB_b05(1,nages);   
  
  init_bounded_number selpar_A50b01_rHB(selpar_A50b01_rHB_LO,selpar_A50b01_rHB_HI,selpar_A50b01_rHB_PH);
  init_bounded_number selpar_slopeb01_rHB(selpar_slopeb01_rHB_LO,selpar_slopeb01_rHB_HI,selpar_slopeb01_rHB_PH);
  init_bounded_number selpar_A50b02_rHB(selpar_A50b02_rHB_LO,selpar_A50b02_rHB_HI,selpar_A50b02_rHB_PH);
  init_bounded_number selpar_slopeb02_rHB(selpar_slopeb02_rHB_LO,selpar_slopeb02_rHB_HI,selpar_slopeb02_rHB_PH);    
  init_bounded_number selpar_A50b03_rHB(selpar_A50b03_rHB_LO,selpar_A50b03_rHB_HI,selpar_A50b03_rHB_PH);
  init_bounded_number selpar_slopeb03_rHB(selpar_slopeb03_rHB_LO,selpar_slopeb03_rHB_HI,selpar_slopeb03_rHB_PH);
  init_bounded_number selpar_A50b04_rHB(selpar_A50b04_rHB_LO,selpar_A50b04_rHB_HI,selpar_A50b04_rHB_PH);
  init_bounded_number selpar_slopeb04_rHB(selpar_slopeb04_rHB_LO,selpar_slopeb04_rHB_HI,selpar_slopeb04_rHB_PH);
  init_bounded_number selpar_A50b05_rHB(selpar_A50b05_rHB_LO,selpar_A50b05_rHB_HI,selpar_A50b05_rHB_PH);
  init_bounded_number selpar_slopeb05_rHB(selpar_slopeb05_rHB_LO,selpar_slopeb05_rHB_HI,selpar_slopeb05_rHB_PH);
  
  vector selpar_A50b01_rHB_out(1,8);
  vector selpar_slopeb01_rHB_out(1,8);
  vector selpar_A50b02_rHB_out(1,8);
  vector selpar_slopeb02_rHB_out(1,8);
  vector selpar_A50b03_rHB_out(1,8);
  vector selpar_slopeb03_rHB_out(1,8);
  vector selpar_A50b04_rHB_out(1,8);
  vector selpar_slopeb04_rHB_out(1,8);
  vector selpar_A50b05_rHB_out(1,8);
  vector selpar_slopeb05_rHB_out(1,8);
  
  number selpar_Age0_rHD;
  number selpar_Age1_rHD;
  number selpar_Age2_rHD;
  
  // rHD
  matrix sel_rHD(styr,endyr,1,nages); 
  vector sel_rHD_b01(1,nages);
  vector sel_rHD_b02(1,nages);
  vector sel_rHD_b03(1,nages);
  vector sel_rHD_b04(1,nages);
  vector sel_rHD_b05(1,nages); 
  
  init_bounded_number selpar_A50b04_rHD(selpar_A50b04_rHD_LO,selpar_A50b04_rHD_HI,selpar_A50b04_rHD_PH);
  init_bounded_number selpar_slopeb04_rHD(selpar_slopeb04_rHD_LO,selpar_slopeb04_rHD_HI,selpar_slopeb04_rHD_PH);
  init_bounded_number selpar_A502b04_rHD(selpar_A502b04_rHD_LO,selpar_A502b04_rHD_HI,selpar_A502b04_rHD_PH);
  init_bounded_number selpar_slope2b04_rHD(selpar_slope2b04_rHD_LO,selpar_slope2b04_rHD_HI,selpar_slope2b04_rHD_PH);
  
  init_bounded_number selpar_A50b05_rHD(selpar_A50b05_rHD_LO,selpar_A50b05_rHD_HI,selpar_A50b05_rHD_PH);
  init_bounded_number selpar_slopeb05_rHD(selpar_slopeb05_rHD_LO,selpar_slopeb05_rHD_HI,selpar_slopeb05_rHD_PH);
  init_bounded_number selpar_A502b05_rHD(selpar_A502b05_rHD_LO,selpar_A502b05_rHD_HI,selpar_A502b05_rHD_PH);
  init_bounded_number selpar_slope2b05_rHD(selpar_slope2b05_rHD_LO,selpar_slope2b05_rHD_HI,selpar_slope2b05_rHD_PH);
  
  vector selpar_A50b04_rHD_out(1,8);
  vector selpar_slopeb04_rHD_out(1,8);
  vector selpar_A502b04_rHD_out(1,8);
  vector selpar_slope2b04_rHD_out(1,8);
  vector selpar_A50b05_rHD_out(1,8);
  vector selpar_slopeb05_rHD_out(1,8);
  vector selpar_A502b05_rHD_out(1,8);
  vector selpar_slope2b05_rHD_out(1,8);

  init_bounded_number selpar_logitAge0_rHD(selpar_logitAge0_rHD_LO,selpar_logitAge0_rHD_HI,selpar_logitAge0_rHD_PH);
  init_bounded_number selpar_logitAge1_rHD(selpar_logitAge1_rHD_LO,selpar_logitAge1_rHD_HI,selpar_logitAge1_rHD_PH);
  init_bounded_number selpar_logitAge2_rHD(selpar_logitAge2_rHD_LO,selpar_logitAge2_rHD_HI,selpar_logitAge2_rHD_PH);
  
  vector selpar_logitAge0_rHD_out(1,8);                                  
  vector selpar_logitAge1_rHD_out(1,8);                                  
  vector selpar_logitAge2_rHD_out(1,8); 
  
  vector prob_belowsizelim_block1(1,nages);
  vector prob_belowsizelim_block2(1,nages);
  vector prob_belowsizelim_block3(1,nages);
  vector prob_belowsizelim_block4(1,nages);
  vector prob_belowsizelim_block5(1,nages);
  
  number zscore_lsizelim1;
  number zscore_lsizelim2;
  number zscore_lsizelim3;
  number zscore_lsizelim4;
  number zscore_lsizelim5;
  
  number cprob_lsizelim1;
  number cprob_lsizelim2;
  number cprob_lsizelim3;
  number cprob_lsizelim4;
  number cprob_lsizelim5;
  
  // rGN
  matrix sel_rGN(styr,endyr,1,nages);
  vector sel_rGN_b01(1,nages);
  vector sel_rGN_b02(1,nages);
  vector sel_rGN_b03(1,nages);
  vector sel_rGN_b04(1,nages);
  vector sel_rGN_b05(1,nages);
  
  // rIN
  matrix sel_rIN(styr,endyr,1,nages);
  vector sel_rIN_b01(1,nages);
  vector sel_rIN_b02(1,nages);
  vector sel_rIN_b03(1,nages);
  vector sel_rIN_b04(1,nages);
  vector sel_rIN_b05(1,nages); 
  
  // rGD
  matrix sel_rGD(styr,endyr,1,nages);
  vector sel_rGD_b01(1,nages);
  vector sel_rGD_b02(1,nages);
  vector sel_rGD_b03(1,nages);
  vector sel_rGD_b04(1,nages);
  vector sel_rGD_b05(1,nages);   
  
  // rID
  matrix sel_rID(styr,endyr,1,nages);
  vector sel_rID_b01(1,nages);
  vector sel_rID_b02(1,nages);
  vector sel_rID_b03(1,nages);
  vector sel_rID_b04(1,nages);
  vector sel_rID_b05(1,nages);   
  
//effort-weighted, recent selectivities
  vector sel_wgted_L(1,nages);  //toward landings
  vector sel_wgted_D(1,nages);  //toward discards  
  vector sel_wgted_tot(1,nages);//toward Z, landings plus deads discards

//-----------------------------------------------------------------------------------------------------------------------------------------------
//-------CPUE Predictions--------------------------------
  vector pred_cpue_sBT(styr_cpue_sBT,endyr_cpue_sBT);       //predicted sBT U (fish/trap-hour)
  matrix N_sBT(styr_cpue_sBT,endyr_cpue_sBT,1,nages);       //used to compute sBT index
  vector pred_cpue_sTV(styr_cpue_sTV,endyr_cpue_sTV);       //predicted sTV U (fish/trap-hour)
  vector obs_cpue_corrected_sTV(styr_cpue_sTV,endyr_cpue_sTV);       //predicted sTV U (fish/trap-hour)
  vector obs_cv_cpue_corrected_sTV(styr_cpue_sTV,endyr_cpue_sTV);       //predicted sTV U (fish/trap-hour)
  vector posobs_cpue_sTV(styr_cpue_sTV,endyr_cpue_sTV);     //whether observed sTV index is negative and should be included in objective function (if positive 1 if negative 0)
  matrix N_sTV(styr_cpue_sTV,endyr_cpue_sTV,1,nages);       //used to compute sTV index
  vector pred_cpue_cHL(styr_cpue_cHL,endyr_cpue_cHL);             //predicted cHL U (pounds/hook-hour)
  matrix N_cHL(styr_cpue_cHL,endyr_cpue_cHL,1,nages);             //used to compute cHL index
  vector pred_cpue_rHB(styr_cpue_rHB,endyr_cpue_rHB);             //predicted rHB U (pounds/hour)
  matrix N_rHB(styr_cpue_rHB,endyr_cpue_rHB,1,nages);             //used to compute rHB index
  //vector pred_cpue_rHD(styr_rHD_cpue,endyr_rHD_cpue);          //predicted rHD U (fish/angler-hour)
  //matrix N_rHD(styr_rHD_cpue,endyr_rHD_cpue,1,nages);          //used to compute rHD index

//---Catchability (CPUE q's)----------------------------------------------------------
  //init_bounded_number log_q_sBT(-20,-10,1);
  //init_bounded_number log_q_sTV(-20,-10,1);
  //init_bounded_number log_q_cHL(-20,-5,1);
  //init_bounded_number log_q_rHB(-20,-5,1);
  //init_bounded_number log_q_rHD(-20,-10,1);
  // init_bounded_number q_rate(0.001,0.1,set_q_rate_phase);
  init_bounded_number q_rate(q_rate_LO,q_rate_HI,q_rate_PH);
  vector q_rate_out(1,8);
  number q_rate_phase;
  
  init_bounded_number log_q_sBT(log_q_sBT_LO,log_q_sBT_HI,log_q_sBT_PH);
  init_bounded_number log_q_sTV(log_q_sTV_LO,log_q_sTV_HI,log_q_sTV_PH);
  init_bounded_number log_q_cHL(log_q_cHL_LO,log_q_cHL_HI,log_q_cHL_PH);
  init_bounded_number log_q_rHB(log_q_rHB_LO,log_q_rHB_HI,log_q_rHB_PH);
  //init_bounded_number log_q_rHD(log_q_rHD_LO,log_q_rHD_HI,log_q_rHD_PH);
  
  vector log_q_sBT_out(1,8);
  vector log_q_sTV_out(1,8);
  vector log_q_cHL_out(1,8);
  vector log_q_rHB_out(1,8);
  //vector log_q_rHD_out(1,8);
  

  vector q_rate_fcn_cHL(styr_cpue_cHL,endyr_cpue_cHL);       //increase due to technology creep (saturates in 2003)
  vector q_rate_fcn_rHB(styr_cpue_rHB,endyr_cpue_rHB);         //increase due to technology creep (saturates in 2003)
  //vector q_rate_fcn_rHD(styr_rHD_cpue,endyr_rHD_cpue);      //increase due to technology creep (saturates in 2003)
  
  //init_bounded_number q_DD_beta(0.1,0.9,set_q_DD_phase);  
  number q_DD_beta;
  vector q_DD_fcn(styr,endyr);    //density dependent function as a multiple of q (scaled a la Katsukawa and Matsuda. 2003)
  number B0_q_DD;                 //B0 of ages q_DD_age plus
  vector B_q_DD(styr,endyr);      //annual biomass of ages q_DD_age plus

  //init_bounded_vector q_RW_log_dev_cHL(styr_cpue_cHL,endyr_cpue_cHL-1,-3.0,3.0,set_q_RW_phase);
  //init_bounded_vector q_RW_log_dev_rHB(styr_cpue_rHB,endyr_cpue_rHB-1,-3.0,3.0,set_q_RW_phase);
  //init_bounded_vector q_RW_log_dev_rHD(styr_rHD_cpue,endyr_rHD_cpue-1,-3.0,3.0,set_q_RW_phase);
 
  //init_bounded_vector q_RW_log_dev_cHL(styr_cpue_cHL,endyr_cpue_cHL-1,log_RWq_LO,log_RWq_HI,log_RWq_PH);
  //init_bounded_vector q_RW_log_dev_rHB(styr_cpue_rHB,endyr_cpue_rHB-1,log_RWq_LO,log_RWq_HI,log_RWq_PH);
  //init_bounded_vector q_RW_log_dev_rHD(styr_rHD_cpue,endyr_rHD_cpue-1,log_RWq_LO,log_RWq_HI,log_RWq_PH);
  init_bounded_vector q_RW_log_dev_cHL(styr_cpue_cHL,endyr_cpue_cHL-1,log_RWq_LO,log_RWq_HI,log_RWq_PH);
  init_bounded_vector q_RW_log_dev_rHB(styr_cpue_rHB,endyr_cpue_rHB-1,log_RWq_LO,log_RWq_HI,log_RWq_PH);
  //init_bounded_vector q_RW_log_dev_rHD(styr_rHD_cpue,endyr_rHD_cpue-1,log_RWq_LO,log_RWq_HI,log_RWq_PH);
 
  vector q_cHL(styr_cpue_cHL,endyr_cpue_cHL);
  vector q_rHB(styr_cpue_rHB,endyr_cpue_rHB);
  //vector q_rHD(styr_rHD_cpue,endyr_rHD_cpue);

//---Landings in numbers (total or 1000 fish) and in wgt (klb)--------------------------------------------------
  matrix L_cHL_num(styr,endyr,1,nages);  //landings (numbers) at age
  matrix L_cHL_klb(styr,endyr,1,nages);  //landings (1000 lb whole weight) at age
  vector pred_L_cHL_knum(styr,endyr);    //yearly landings in 1000 fish summed over ages
  vector pred_L_cHL_klb(styr,endyr);     //yearly landings in 1000 lb summed over ages

  matrix L_cPT_num(styr,endyr,1,nages);   //landings (numbers) at age
  matrix L_cPT_klb(styr,endyr,1,nages);   //landings (1000 lb whole weight) at age    
  vector pred_L_cPT_knum(styr,endyr);     //yearly landings in 1000 fish summed over ages   
  vector pred_L_cPT_klb(styr,endyr);      //yearly landings in 1000 lb summed over ages 

  matrix L_cTW_num(styr,endyr,1,nages);   //landings (numbers) at age
  matrix L_cTW_klb(styr,endyr,1,nages);   //landings (1000 lb whole weight) at age    
  vector pred_L_cTW_knum(styr,endyr);     //yearly landings in 1000 fish summed over ages   
  vector pred_L_cTW_klb(styr,endyr);      //yearly landings in 1000 lb summed over ages 

  matrix L_rHB_num(styr,endyr,1,nages);   //landings (numbers) at age
  matrix L_rHB_klb(styr,endyr,1,nages);   //landings (1000 lb whole weight) at age    
  vector pred_L_rHB_knum(styr,endyr);     //yearly landings in 1000 fish summed over ages  
  vector pred_L_rHB_klb(styr,endyr);      //yearly landings in 1000 lb summed over ages

  matrix L_rGN_num(styr,endyr,1,nages);  //landings (numbers) at age
  matrix L_rGN_klb(styr,endyr,1,nages);  //landings (1000 lb whole weight) at age    
  vector pred_L_rGN_knum(styr,endyr);    //yearly landings in 1000 fish summed over ages  
  vector pred_L_rGN_klb(styr,endyr);     //yearly landings in 1000 lb summed over ages

  matrix L_rIN_num(styr,endyr,1,nages);  //landings (numbers) at age
  matrix L_rIN_klb(styr,endyr,1,nages);  //landings (1000 lb whole weight) at age    
  vector pred_L_rIN_knum(styr,endyr);    //yearly landings in 1000 fish summed over ages  
  vector pred_L_rIN_klb(styr,endyr);     //yearly landings in 1000 lb summed over ages  

  matrix L_total_num(styr,endyr,1,nages);//total landings in number at age
  matrix L_total_klb(styr,endyr,1,nages);//landings in klb at age 
  vector L_total_knum_yr(styr,endyr);    //total landings in 1000 fish by yr summed over ages  
  vector L_total_klb_yr(styr,endyr);     //total landings (klb) by yr summed over ages

//---Dead discards in numbers (total or 1000 fish) and in wgt (klb) --------------------------------------------------
  matrix D_cGD_num(styr,endyr,1,nages);     //discards (numbers) at age
  matrix D_cGD_klb(styr,endyr,1,nages);     //discards (1000 lb) at age
  vector pred_D_cGD_knum(styr,endyr);       //yearly discards summed over ages
  vector pred_D_cGD_klb(styr,endyr);        //yearly discards in klb summed over ages
  vector obs_D_cGD(styr_D_cGD,endyr_D_cGD); //observed releases multiplied by discard mortality
  // vector obs_D_cHD(styr_D_cHD,endyr_D_cHD); //observed releases multiplied by discard mortality
  // vector pred_D_cHD_knum(styr,endyr);       // Not used in modeling. empty vector passed to rdat
  // vector obs_D_cPD(styr_D_cHD,endyr_D_cHD); //observed releases multiplied by discard mortality
  // vector pred_D_cPD_knum(styr,endyr);       // Not used in modeling. empty vector passed to rdat
  // vector cGD_cv(styr_D_cGD,endyr_D_cGD);    //CVs for fitting combined discards
  
  matrix D_rHD_num(styr,endyr,1,nages);     //discards (numbers) at age
  matrix D_rHD_klb(styr,endyr,1,nages);     //discards (1000 lb) at age
  vector pred_D_rHD_knum(styr,endyr);       //yearly discards summed over ages
  vector pred_D_rHD_klb(styr,endyr);        //yearly discards in klb summed over ages
  vector obs_D_rHD(styr_D_rHD,endyr_D_rHD);   //observed releases multiplied by discard mortality

  matrix D_rGD_num(styr,endyr,1,nages);   //discards (numbers) at age
  matrix D_rGD_klb(styr,endyr,1,nages);   //discards (1000 lb) at age
  vector pred_D_rGD_knum(styr,endyr);     //yearly discards summed over ages
  vector pred_D_rGD_klb(styr,endyr);      //yearly discards in klb summed over ages
  vector obs_D_rGD(styr_D_rGD,endyr_D_rGD); //observed releases multiplied by discard mortality
  
  matrix D_rID_num(styr,endyr,1,nages);   //discards (numbers) at age
  matrix D_rID_klb(styr,endyr,1,nages);   //discards (1000 lb) at age
  vector pred_D_rID_knum(styr,endyr);     //yearly discards summed over ages
  vector pred_D_rID_klb(styr,endyr);      //yearly discards in klb summed over ages
  vector obs_D_rID(styr_D_rID,endyr_D_rID); //observed releases multiplied by discard mortality
  
  matrix D_total_num(styr,endyr,1,nages);  //total discards in number at age
  matrix D_total_klb(styr,endyr,1,nages);  //discards in klb at age 
  vector D_total_knum_yr(styr,endyr);      //total discards in 1000 fish by yr summed over ages  
  vector D_total_klb_yr(styr,endyr);       //total discards (klb) by yr summed over ages

////---MSY calcs----------------------------------------------------------------------------
  number prop_F_cHL;       //proportion of F_sum attributable to hal, last X=nyr_selpar_wgted yrs, used for avg body weights
  number prop_F_cPT;       //proportion of F_sum attributable to pots, last X yrs
  number prop_F_rHB;       //proportion of F_sum attributable to headboat, last X yrs
  number prop_F_rGN;     //proportion of F_sum attributable to rGN, last X yrs
  number prop_F_rIN;     //proportion of F_sum attributable to rIN, last X yrs
  number prop_F_cGD;   //proportion of F_sum attributable to cGD discards, last X yrs
  number prop_F_rHD;     //proportion of F_sum attributable to headboat discards, last X yrs
  number prop_F_rGD;   //proportion of F_sum attributable to rGD discards, last X yrs
  number prop_F_rID;   //proportion of F_sum attributable to rID discards, last X yrs
  number F_temp_sum;      //sum of geom mean Fsum's in last X yrs, used to compute F_fishery_prop

  vector F_end(1,nages);
  vector F_end_L(1,nages);
  vector F_end_D(1,nages);    
  number F_end_apex;
  
  number SSB_msy_n_out;           //SSB (popn fecudity) at msy
  number F_msy_n_out;             //F at msy
  number msy_n_klb_out;           //max sustainable yield (1000 lb)
  number msy_n_knum_out;          //max sustainable yield (1000 fish)  
  number B_msy_n_out;             //total biomass at MSY 
  number R_msy_n_out;             //equilibrium recruitment at F=Fmsy
  number D_msy_n_knum_out;        //equilibrium dead discards (1000 fish) at F=Fmsy
  number D_msy_n_klb_out;         //equilibrium dead discards (1000 lb) at F=Fmsy
  number T_msy_n_knum_out;        //equilibrium Total yield  (1000 fish; dead discards and landings) at F=Fmsy
  number T_msy_n_klb_out;         //equilibrium dead discards (1000 lb; dead discards and landings) at F=Fmsy  
  number spr_msy_n_out;           //spr at F=Fmsy
  
  number SSB_msy_lb_out;           //SSB (popn fecudity) at msy
  number F_msy_lb_out;             //F at msy
  number msy_lb_klb_out;           //max sustainable yield (1000 lb)
  number msy_lb_knum_out;          //max sustainable yield (1000 fish)  
  number B_msy_lb_out;             //total biomass at MSY 
  number R_msy_lb_out;             //equilibrium recruitment at F=Fmsy
  number D_msy_lb_knum_out;        //equilibrium dead discards (1000 fish) at F=Fmsy
  number D_msy_lb_klb_out;         //equilibrium dead discards (1000 lb) at F=Fmsy
  number T_msy_lb_knum_out;        //equilibrium Total yield  (1000 fish; dead discards and landings) at F=Fmsy
  number T_msy_lb_klb_out;         //equilibrium dead discards (1000 lb; dead discards and landings) at F=Fmsy  
  number spr_msy_lb_out;           //spr at F=Fmsy
  
  number F30_dum;				//intermediate calculation for F30
  number F40_dum;				//intermediate calculation for F40
  number F50_dum;				//intermediate calculation for F50
  number Fproxy_dum;		    //intermediate calculation for input spr_proxy level
  number Fproxy_out;            //FX%, defined by user  
  number SSB_Fproxy_out;
  number B_Fproxy_out;
  number R_Fproxy_out;
  number L_Fproxy_knum_out;
  number L_Fproxy_klb_out;
  number D_Fproxy_knum_out;
  number D_Fproxy_klb_out;
  number F50_out;            //F50%
  number SSB_F50_out;
  number B_F50_out;
  number R_F50_out;
  number L_F50_knum_out;
  number L_F50_klb_out;
  number D_F50_knum_out;
  number D_F50_klb_out;
  number F40_out;            //F40%
  number SSB_F40_out;
  number B_F40_out;
  number R_F40_out;
  number L_F40_knum_out;
  number L_F40_klb_out;
  number D_F40_knum_out;
  number D_F40_klb_out;
  number F30_out;            //F30%
  number SSB_F30_out;
  number B_F30_out;
  number R_F30_out;
  number L_F30_knum_out;
  number L_F30_klb_out;
  number D_F30_knum_out;
  number D_F30_klb_out;
  number rec_mean;  			//arithemetic average recruitment used in SPR-related quantities
  number rec_mu;  			    //geometric average recruitment deviate used to predict recruitment after endyr_rec_dev
  number rec_mu_temp;  			    //product of recruitment deviate to calculate recruitment after endyr_rec_dev
  
  vector N_age_msy(1,nages);         //numbers at age for MSY calculations: beginning of yr
  vector N_age_msy_spawn(1,nages);   //numbers at age for MSY calculations: time of peak spawning
  vector L_age_msy(1,nages);         //catch at age for MSY calculations
  vector Z_age_msy(1,nages);         //total mortality at age for MSY calculations
  vector D_age_msy(1,nages);         //discard mortality (dead discards) at age for MSY calculations
  vector T_age_msy(1,nages);         //total yeild (landings and dead discards) at age for MSY calculations
  vector F_L_age_msy(1,nages);       //fishing mortality landings (not discards) at age for MSY calculations
  vector F_D_age_msy(1,nages);       //fishing mortality of discards at age for MSY calculations
  vector F_msy(1,n_iter_msy);        //values of full F to be used in equilibrium calculations
  vector spr_msy(1,n_iter_msy);      //reproductive capacity-per-recruit values corresponding to F values in F_msy
  vector R_eq(1,n_iter_msy);         //equilibrium recruitment values corresponding to F values in F_msy
  vector L_eq_klb(1,n_iter_msy);     //equilibrium landings(klb) values corresponding to F values in F_msy
  vector L_eq_knum(1,n_iter_msy);    //equilibrium landings(1000 fish) values corresponding to F values in F_msy
  vector SSB_eq(1,n_iter_msy);       //equilibrium reproductive capacity values corresponding to F values in F_msy
  vector B_eq(1,n_iter_msy);         //equilibrium biomass values corresponding to F values in F_msy
  vector D_eq_klb(1,n_iter_msy);     //equilibrium discards (klb) corresponding to F values in F_msy
  vector D_eq_knum(1,n_iter_msy);    //equilibrium discards (1000s) corresponding to F values in F_msy
  vector T_eq_klb(1,n_iter_msy);     //equilibrium total yield (klb) landings and discards corresponding to F values in F_msy
  vector T_eq_knum(1,n_iter_msy);    //equilibrium total yield (1000s) landings and discards corresponding to F values in F_msy

  vector FdF_msy_n(styr,endyr);
  vector FdF_msy_lb(styr,endyr);
  vector FdFproxy(styr,endyr);
  vector FdF50(styr,endyr);
  vector FdF40(styr,endyr);
  vector FdF30(styr,endyr);
  vector SdSSB_msy_n(styr,endyr);
  number SdSSB_msy_n_end;
  number SSBend_MSST_n;
  vector SdSSB_msy_lb(styr,endyr);
  number SdSSB_msy_lb_end;
  number SSBend_MSST_lb;
  number FdF_msy_n_end;
  number FdF_msy_n_end_mean;          //geometric mean of last 3 yrs
  number FdF_msy_lb_end;
  number FdF_msy_lb_end_mean;          //geometric mean of last 3 yrs
  
  vector SdSSB_Fproxy(styr,endyr);
  vector Sdmsst_Fproxy(styr,endyr);
  number SdSSB_Fproxy_end;
  number Sdmsst_Fproxy_end;
  number FdFproxy_end_mean;          //geometric mean of last nyr_selpar_wgted yrs
  vector SdSSB_F50(styr,endyr);
  vector Sdmsst_F50(styr,endyr);
  number SdSSB_F50_end;
  number Sdmsst_F50_end;
  number FdF50_end_mean;          //geometric mean of last nyr_selpar_wgted yrs
  vector SdSSB_F40(styr,endyr);
  vector Sdmsst_F40(styr,endyr);
  number SdSSB_F40_end;
  number Sdmsst_F40_end;
  number FdF40_end_mean;          //geometric mean of last nyr_selpar_wgted yrs
  vector SdSSB_F30(styr,endyr);
  vector Sdmsst_F30(styr,endyr);
  number SdSSB_F30_end;
  number Sdmsst_F30_end;
  number FdF30_end_mean;          //geometric mean of last nyr_selpar_wgted yrs
  number Fend_mean_temp;			 //intermediate calc for geometric mean of last nyr_selpar_wgted yrs
  number Fend_mean;					 //geometric mean of last nyr_selpar_wgted yrs
  vector L_age_Fproxy(1,nages);      //landings at age for Fproxy calculations
  vector D_age_Fproxy(1,nages);      //discards at age for Fproxy calculations  

  vector wgt_wgted_L_klb(1,nages);  //fishery-weighted average weight at age of landings
  vector wgt_wgted_D_klb(1,nages);  //fishery-weighted average weight at age of discards  
  number wgt_wgted_L_denom;         //used in intermediate calculations
  number wgt_wgted_D_denom;         //used in intermediate calculations
  
  number iter_inc_msy;               //increments used to compute msy, equals 1/(n_iter_msy-1)
  
////--------Mortality------------------------------------------------------------------
  vector M(1,nages);                         //age-dependent natural mortality
  //number M_constant;                         //age-indpendent: used only for MSST
  init_bounded_number M_constant(M_constant_LO,M_constant_HI,M_constant_PH);   //age-indpendent: used only for MSST
  vector M_constant_out(1,8);
  // vector Mscale_len(1,nages);
  // vector M_lorenzen(1,nages);   
  // number cum_surv_1plus;
  // number nages_M_sum;

  number smsy2msstM;                         //scales Smsy to get msst using (1-M). Used only in output.
  
  matrix F(styr,endyr,1,nages);
  vector Fsum(styr,endyr);                   //Full fishing mortality rate by year
  vector Fapex(styr,endyr);                  //Max across ages, fishing mortality rate by year (may differ from Fsum bc of dome-shaped sel 
  matrix Z(styr,endyr,1,nages);
  
  init_bounded_number log_avg_F_L_cHL(log_avg_F_cHL_LO,log_avg_F_cHL_HI,log_avg_F_cHL_PH);
  vector out_log_avg_F_L_cHL(1,8);  
  init_bounded_dev_vector log_dev_F_L_cHL(styr_L_cHL,endyr_L_cHL,log_dev_F_cHL_LO,log_dev_F_cHL_HI,log_dev_F_cHL_PH);
  vector log_dev_F_L_out_cHL(styr_L_cHL,endyr_L_cHL);
  matrix F_cHL(styr,endyr,1,nages);
  vector F_cHL_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number log_dev_F_init_cHL;
  number log_dev_F_end_cHL;  

  init_bounded_number log_avg_F_L_cPT(log_avg_F_cPT_LO,log_avg_F_cPT_HI,log_avg_F_cPT_PH);
  vector out_log_avg_F_L_cPT(1,8);  
  init_bounded_dev_vector log_dev_F_L_cPT(styr_L_cPT,endyr_L_cPT,log_dev_F_cPT_LO,log_dev_F_cPT_HI,log_dev_F_cPT_PH);
  vector log_dev_F_L_out_cPT(styr_L_cPT,endyr_L_cPT);
  matrix F_cPT(styr,endyr,1,nages);
  vector F_cPT_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number log_dev_F_init_cPT;
  number log_dev_F_end_cPT;  

  init_bounded_number log_avg_F_L_cTW(log_avg_F_cTW_LO,log_avg_F_cTW_HI,log_avg_F_cTW_PH);
  vector out_log_avg_F_L_cTW(1,8);  
  init_bounded_dev_vector log_dev_F_L_cTW(styr_L_cTW,endyr_L_cTW,log_dev_F_cTW_LO,log_dev_F_cTW_HI,log_dev_F_cTW_PH);
  vector log_dev_F_L_out_cTW(styr_L_cTW,endyr_L_cTW);
  matrix F_cTW(styr,endyr,1,nages);
  vector F_cTW_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number log_dev_F_init_cTW;
  number log_dev_F_end_cTW;  
 
  init_bounded_number log_avg_F_L_rHB(log_avg_F_rHB_LO,log_avg_F_rHB_HI,log_avg_F_rHB_PH);
  vector out_log_avg_F_L_rHB(1,8); 
  init_bounded_dev_vector log_dev_F_L_rHB(styr_L_rHB,endyr_L_rHB,log_dev_F_rHB_LO,log_dev_F_rHB_HI,log_dev_F_rHB_PH);    
  vector log_dev_F_L_out_rHB(styr_L_rHB,endyr_L_rHB);
  matrix F_rHB(styr,endyr,1,nages);
  vector F_rHB_out(styr,endyr);        //used for intermediate calculations in fcn get_mortality
  number log_dev_F_init_rHB;
  number log_dev_F_end_rHB;    

  init_bounded_number log_avg_F_L_rGN(log_avg_F_rGN_LO,log_avg_F_rGN_HI,log_avg_F_rGN_PH);
  vector out_log_avg_F_L_rGN(1,8); 
  init_bounded_dev_vector log_dev_F_L_rGN(styr_L_rGN,endyr_L_rGN,log_dev_F_rGN_LO,log_dev_F_rGN_HI,log_dev_F_rGN_PH);    
  vector log_dev_F_L_out_rGN(styr_L_rGN,endyr_L_rGN);
  matrix F_rGN(styr,endyr,1,nages);
  vector F_rGN_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number log_dev_F_init_rGN;
  number log_dev_F_end_rGN;
  
  init_bounded_number log_avg_F_L_rIN(log_avg_F_rIN_LO,log_avg_F_rIN_HI,log_avg_F_rIN_PH);
  vector out_log_avg_F_L_rIN(1,8); 
  init_bounded_dev_vector log_dev_F_L_rIN(styr_L_rIN,endyr_L_rIN,log_dev_F_rIN_LO,log_dev_F_rIN_HI,log_dev_F_rIN_PH);    
  vector log_dev_F_L_out_rIN(styr_L_rIN,endyr_L_rIN);
  matrix F_rIN(styr,endyr,1,nages);
  vector F_rIN_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number log_dev_F_init_rIN;
  number log_dev_F_end_rIN;  
  
  init_bounded_number F_init_ratio(F_init_ratio_LO,F_init_ratio_HI,F_init_ratio_PH); //number F_init_ratio;
  vector F_init_ratio_out(1,8); 

  //--Discard mortality stuff------------------------------------------------------------------------------
  init_bounded_number log_avg_F_D_cGD(LO_log_avg_F_D_cGD,HI_log_avg_F_D_cGD,PH_log_avg_F_D_cGD);
  vector out_log_avg_F_D_cGD(1,8);  
  init_bounded_dev_vector log_dev_F_D_cGD(styr_D_cGD,endyr_D_cGD,log_dev_F_cGD_LO,log_dev_F_cGD_HI,log_dev_F_cGD_PH);
  vector log_dev_F_D_out_cGD(styr_D_cGD,endyr_D_cGD);
  matrix F_cGD(styr,endyr,1,nages);
  vector F_cGD_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number log_dev_F_cGD2;   //avg log deviations in reg period 2 (for estimation 1984-1992, prior to data) 
  number log_dev_F_end_cGD;  
  
  init_bounded_number log_avg_F_D_rHD(LO_log_avg_F_D_rHD,HI_log_avg_F_D_rHD,PH_log_avg_F_D_rHD);
  vector out_log_avg_F_D_rHD(1,8);  
  init_bounded_dev_vector log_dev_F_D_rHD(styr_D_rHD,endyr_D_rHD,log_dev_F_rHD_LO,log_dev_F_rHD_HI,log_dev_F_rHD_PH);
  vector log_dev_F_D_out_rHD(styr_D_rHD,endyr_D_rHD);
  matrix F_rHD(styr,endyr,1,nages);
  vector F_rHD_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number log_dev_F_end_rHD;  
       
  init_bounded_number log_avg_F_D_rGD(LO_log_avg_F_D_rGD,HI_log_avg_F_D_rGD,PH_log_avg_F_D_rGD);
  vector out_log_avg_F_D_rGD(1,8);  
  init_bounded_dev_vector log_dev_F_D_rGD(styr_D_rGD,endyr_D_rGD,log_dev_F_rGD_LO,log_dev_F_rGD_HI,log_dev_F_rGD_PH);
  vector log_dev_F_D_out_rGD(styr_D_rGD,endyr_D_rGD);
  matrix F_rGD(styr,endyr,1,nages);
  vector F_rGD_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number log_dev_F_init_rGD;
  number log_dev_F_end_rGD;
  
  init_bounded_number log_avg_F_D_rID(LO_log_avg_F_D_rID,HI_log_avg_F_D_rID,PH_log_avg_F_D_rID);
  vector out_log_avg_F_D_rID(1,8);  
  init_bounded_dev_vector log_dev_F_D_rID(styr_D_rID,endyr_D_rID,log_dev_F_rID_LO,log_dev_F_rID_HI,log_dev_F_rID_PH);
  vector log_dev_F_D_out_rID(styr_D_rID,endyr_D_rID);
  matrix F_rID(styr,endyr,1,nages);
  vector F_rID_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number log_dev_F_init_rID;
  number log_dev_F_end_rID;  

  number Dmort_cGD;
  // number Dmort_cPD1;
  // number Dmort_cPD2;   
  // number Dmort_cHD;
  number Dmort_rHD;
  number Dmort_rGD;
  number Dmort_rID;

//---Per-recruit stuff----------------------------------------------------------------------------------
  vector N_age_spr(1,nages);         //numbers at age for SPR calculations: beginning of year
  vector N_age_spr_spawn(1,nages);   //numbers at age for SPR calculations: at time of peak spawning  
  vector L_age_spr(1,nages);         //landings catch at age for SPR calculations
  vector T_age_spr(1,nages);         //total catch at age for SPR calculations
  vector Z_age_spr(1,nages);         //total mortality at age for SPR calculations
  vector spr_static(styr,endyr);     //vector of static SPR values by year
  vector F_L_age_spr(1,nages);       //fishing mortality of landings (not discards) at age for SPR calculations
  vector F_T_age_spr(1,nages);       //fishing mortality of landings and discards at age for SPR calculations
  vector F_D_age_spr(1,nages);       //fishing mortality of discards at age for SPR calculations
  vector F_spr(1,n_iter_spr);        //values of full F to be used in per-recruit calculations
  vector spr_spr(1,n_iter_spr);      //reproductive capacity-per-recruit values corresponding to F values in F_spr
  vector spr_ratio(1,n_iter_spr);    //reproductive capacity-per-recruit relative to spr_F0 values corresponding to F values in F_spr
  vector L_spr(1,n_iter_spr);        //landings(lb)-per-recruit (ypr) values corresponding to F values in F_spr
  vector T_spr(1,n_iter_spr);        //total catch(lb)-per-recruit (ypr) values corresponding to F values in F_spr (combined landings and discards)

  vector N_spr_F0(1,nages);          //Used to compute spr at F=0: at time of peak spawning
  vector N_bpr_F0(1,nages);          //Used to compute bpr at F=0: at start of year  
  vector N_spr_initial(1,nages);     //Initial spawners per recruit at age given initial F
  vector N_initial_eq(1,nages);      //Initial equilibrium abundance at age
  vector F_initial(1,nages);         //initial F at age
  vector Z_initial(1,nages);         //initial Z at age
  number spr_initial;                //initial spawners per recruit
  number spr_F0;                     //Spawning biomass per recruit at F=0
  number bpr_F0;                     //Biomass per recruit at F=0

  number iter_inc_spr;               //increments used to compute msy, equals max_F_spr_msy/(n_iter_spr-1)

////-------SDNR output-----------------------------------------------------------------------------
 
  number sdnr_lc_sBT;
  number sdnr_lc_sTV;
  number sdnr_lc_cHL;
  number sdnr_lc_cPT;
  number sdnr_lc_rHB;  
  number sdnr_lc_rHD;  
  number sdnr_lc_rGN;  
 
  number sdnr_ac_sBT;
  number sdnr_ac_sTV;
  number sdnr_ac_cHL; 
  number sdnr_ac_cPT;  
  number sdnr_ac_rHB;
  //number sdnr_ac_rGN;
  
  number sdnr_I_sBT;
  number sdnr_I_sTV;
  number sdnr_I_cHL;
  number sdnr_I_rHB;  
  
//-------Objective function components-----------------------------------------------------------------------------
  number w_L;
  number w_D;

  number w_lenc_sBT;
  number w_lenc_cHL;
  number w_lenc_cPT;
  number w_lenc_rHB;
  number w_lenc_rHD;
  number w_lenc_rGN;

  number w_agec_sBT; 
  number w_agec_sTV; 
  number w_agec_cHL; 
  number w_agec_cPT;  
  number w_agec_rHB;
  
  number w_cpue_sBT;  
  number w_cpue_sTV;
  number w_cpue_cHL;
  number w_cpue_rHB;
  
  number w_rec;
  number w_rec_early;
  number w_rec_end;
  
  number f_cpue_sBT;
  number f_cpue_sTV;
  number f_cpue_cHL;
  number f_cpue_rHB;

 
  number f_L_cHL;  
  number f_L_cPT;
  number f_L_cTW;
  number f_L_rHB;
  number f_L_rGN;
  number f_L_rIN;

  number f_D_cGD; 
  number f_D_rHD;
  number f_D_rGD;
  number f_D_rID;

  number f_lenc_sBT;   
  number f_lenc_cHL;
  number f_lenc_cPT;
  number f_lenc_rHB;
  number f_lenc_rHD;
  number f_lenc_rGN;

  number f_agec_sBT;
  number f_agec_sTV;  
  number f_agec_cHL;  
  number f_agec_cPT;    
  number f_agec_rHB; 
  number f_agec_rGN;
  
  number f_cpue_RW_cHL; //random walk component of indices
  number f_cpue_RW_rHB;
  number f_cpue_RW_rHD;   
  
//Penalties and constraints. Not all are used.
  number f_rec_dev;                //weight on recruitment deviations to fit S-R curve
  number f_rec_dev_early;          //extra weight on deviations in first recruitment stanza
  number f_rec_dev_end;            //extra weight on deviations in first recruitment stanza
  number f_priors;                  //prior information on parameters
  
  objective_function_value fval;
  number fval_data;
  number grad_max;
  
//--Dummy variables ----
  number denom;                   //denominator used in some calculations
  number numer;                   //numerator used in some calculations


  //---------- Projection quantities--------------------------------------------------------------
  number F_reg_proj;						   //value used to define the projections
  vector F_proj(styr_proj,endyr_proj);         //F by yr for projections (=F_reg_proj after regulations start, current F till then)
  vector L_knum_proj(styr_proj,endyr_proj);    //total landings in 1000 fish for projections
  vector L_klb_proj(styr_proj,endyr_proj);     //total landings in weight (1000 lb) for projections

  vector B_proj(styr_proj,endyr_proj);         //Biomass for projections
  vector SSB_proj(styr_proj,endyr_proj);       //SSB for projections
  vector R_proj(styr_proj,endyr_proj);     	   //recruits for projections
  vector FL_age_proj(1,nages);      		   //F (landings) by age for projections
  vector FD_age_proj(1,nages);      		   //F (landings) by age for projections
  
  matrix N_proj(styr_proj,endyr_proj,1,nages);           //Population numbers by year and age at start of yr
  matrix N_spawn_proj(styr_proj,endyr_proj,1,nages);     //Population numbers by year and age at peaking spawning: used for SSB in projections
  matrix Z_proj(styr_proj,endyr_proj,1,nages);           //Z by year and age for projections
  matrix L_age_proj(styr_proj,endyr_proj,1,nages);       //Projected landings at age in numbers
  matrix D_age_proj(styr_proj,endyr_proj,1,nages);       //Projected discards at age in numbers
  vector D_knum_proj(styr_proj,endyr_proj);    //total discards in 1000 fish for projections
  vector D_klb_proj(styr_proj,endyr_proj);     //total discards in weight (1000 lb) for projections




//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
// INITIALIZATION_SECTION


//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
GLOBALS_SECTION
  #include "admodel.h"          // Include AD class definitions
  #include "admb2r.cpp"    // Include S-compatible output functions (needs preceding)
	time_t start,finish;
	long hour,minute,second;	
	double elapsed_time;
	
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
RUNTIME_SECTION
 maximum_function_evaluations 10000, 20000, 30000, 50000, 100000, 100000, 100000;
 convergence_criteria 1e-2, 1e-2,1e-2, 1e-4;
 
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
PRELIMINARY_CALCS_SECTION

//// Set values of fixed parameters or set initial guess of estimated parameters
 sqrt2pi=sqrt(2.*3.14159265);
 g2mt=0.000001;         //conversion of grams to metric tons
 g2kg=0.001;            //conversion of grams to kg 
 mt2klb=2.20462;        //conversion of metric tons to 1000 lb 
 mt2lb=mt2klb*1000.0;   //conversion of metric tons to lb
 g2klb=g2mt*mt2klb;     //conversion of grams to 1000 lb 
 dzero=0.00001;         
 huge_number=1.0e+10;   
 onehalf=0.5;

  Dmort_cGD=set_Dmort_cGD;
  // Dmort_cHD=set_Dmort_cHD;
  // Dmort_cPD1=set_Dmort_cPD1;
  // Dmort_cPD2=set_Dmort_cPD2;
  Dmort_rHD=set_Dmort_rHD;
  Dmort_rGD=set_Dmort_rGD;
  Dmort_rID=set_Dmort_rID;


  // //Calculate values used for weighting selex and avg weight of cGN discards for when cHL is closed (2009-2013)
  // // cHL
  // Dopen_cHL4=Dmort_cHD*sum(obs_released_cHD(styr_closure_com,(endyr_sel_b04+1)))/((endyr_sel_b04+1)-styr_closure_com+1);
  // Dclosed_cHL4=Dmort_cHD*sum(obs_released_cHDclosed(styr_closure_com,(endyr_sel_b04+1)))/((endyr_sel_b04+1)-styr_closure_com+1);
  // Lopen_cHL4=Dmort_cHD*sum(obs_L_cHL(styr_closure_com,(endyr_sel_b04+1)))/((endyr_sel_b04+1)-styr_closure_com+1);
  // // cPT
  // Dopen_cPT4=Dmort_cPD2*sum(obs_released_cPD(styr_closure_com,(endyr_sel_b04+1)))/((endyr_sel_b04+1)-styr_closure_com+1);
  // Dclosed_cPT4=Dmort_cPD2*sum(obs_released_cPDclosed(styr_closure_com,(endyr_sel_b04+1)))/((endyr_sel_b04+1)-styr_closure_com+1);
  // Lopen_cPT4=Dmort_cPD2*sum(obs_L_cPT(styr_closure_com,(endyr_sel_b04+1)))/((endyr_sel_b04+1)-styr_closure_com+1);  //Average wt is close to 1 lb so didn't convert to numbers
  // 
  // D_sum_cHLcPT4=Dopen_cHL4+Dclosed_cHL4+Dopen_cPT4+Dclosed_cPT4;
  // 
  // Dprop_cGN_sel_D4=(Dopen_cHL4 + Dopen_cPT4 + Dclosed_cHL4*(Dopen_cHL4/(Dopen_cHL4+Lopen_cHL4)) + Dclosed_cPT4*(Dopen_cPT4/(Dopen_cPT4+Lopen_cPT4)))/D_sum_cHLcPT4; 
  // Dprop_cGN_sel_cHL4=Dclosed_cHL4*(Lopen_cHL4/(Dopen_cHL4+Lopen_cHL4))/D_sum_cHLcPT4; 
  // Dprop_cGN_sel_cPT4=Dclosed_cPT4*(Lopen_cPT4/(Dopen_cPT4+Lopen_cPT4))/D_sum_cHLcPT4; 
  // 
  // obs_D_cHD=Dmort_cHD*obs_released_cHD;
  // obs_D_cHD(styr_D_cHDclosed,endyr_D_cHDclosed)+=Dmort_cHD*obs_released_cHDclosed;
  // 
  // obs_D_cPD(styr_D_cPD,2006)=Dmort_cPD1*obs_released_cPD(styr_D_cPD,2006);
  // obs_D_cPD(2007,endyr_D_cPD)=Dmort_cPD2*obs_released_cPD(2007,endyr_D_cPD);
  // obs_D_cPD(styr_D_cPDclosed,endyr_D_cPDclosed)+=Dmort_cPD2*obs_released_cPDclosed;
  // 
  // obs_D_cGD=obs_D_cHD+obs_D_cPD;
  
  obs_D_cGD=Dmort_cGD*obs_released_cGD;
  obs_D_rHD=Dmort_rHD*obs_released_rHD;
  obs_D_rGD=Dmort_rGD*obs_released_rGD;
  obs_D_rID=Dmort_rID*obs_released_rID;
  // cGD_cv=obs_cv_D_cHD;
 
  Linf=set_Linf(1);
  K=set_K(1);
  t0=set_t0(1);

  M=set_M; 
  M_constant=set_M_constant(1);
  // Mscale_len=Linf*(1.0-mfexp(-K*(agebins-t0+0.5))); 
  // M_lorenzen=pow(Mscale_len,-1.0);
  // nages_M_sum=nages-min_M_age+1;
  
  
  steep=set_steep(1);
  R_autocorr=set_R_autocorr(1);
  rec_sigma=set_rec_sigma(1);
  
  log_dm_lenc_sBT=set_log_dm_lenc_sBT(1);
  // log_dm_sTV_lc=set_log_dm_sTV_lc(1);
  log_dm_lenc_cHL=set_log_dm_lenc_cHL(1);
  log_dm_lenc_cPT=set_log_dm_lenc_cPT(1);
  log_dm_lenc_rHB=set_log_dm_lenc_rHB(1);
  log_dm_lenc_rHD=set_log_dm_lenc_rHD(1);
  log_dm_lenc_rGN=set_log_dm_lenc_rGN(1);
  
  log_dm_agec_sBT=set_log_dm_agec_sBT(1);
  log_dm_agec_sTV=set_log_dm_agec_sTV(1);
  log_dm_agec_cHL=set_log_dm_agec_cHL(1);
  log_dm_agec_cPT=set_log_dm_agec_cPT(1);
  log_dm_agec_rHB=set_log_dm_agec_rHB(1);
    
  log_q_sBT=set_log_q_cpue_sBT(1);  
  log_q_sTV=set_log_q_cpue_sTV(1);
  log_q_cHL=set_log_q_cpue_cHL(1);
  log_q_rHB=set_log_q_cpue_rHB(1);
 
  q_rate=set_q_rate(1);
  q_rate_phase=set_q_rate(4);
  q_rate_fcn_cHL=1.0;
  q_rate_fcn_rHB=1.0;
 
  q_DD_beta=set_q_DD_beta;
  q_DD_fcn=1.0;
  q_RW_log_dev_cHL.initialize();
  q_RW_log_dev_rHB.initialize();

  if (q_rate_phase < 0.0 & q_rate!=0.0)
  {
      for (iyear=styr_cpue_cHL; iyear<=endyr_cpue_cHL; iyear++)
      {   if (iyear>styr_cpue_cHL & iyear <=2003) 
          {
             q_rate_fcn_cHL(iyear)=(1.0+(iyear-styr_cpue_cHL)*q_rate)*q_rate_fcn_cHL(styr_cpue_cHL);
          }
          if (iyear>2003) {q_rate_fcn_cHL(iyear)=q_rate_fcn_cHL(iyear-1);} 
      }   
      for (iyear=styr_cpue_rHB; iyear<=endyr_cpue_rHB; iyear++)
      {   if (iyear>styr_cpue_rHB & iyear <=2003) 
          {
            q_rate_fcn_rHB(iyear)=(1.0+(iyear-styr_cpue_rHB)*q_rate)*q_rate_fcn_rHB(styr_cpue_rHB); 
           }
          if (iyear>2003) {q_rate_fcn_rHB(iyear)=q_rate_fcn_rHB(iyear-1);} 
      }  
  } //end q_rate conditional
  
  // Fill in vector of whether the values are positive or not for calculating objective function
  // and replace negative values (e.g. -999) in observed time series with missing,
  // which shows up in the rdat file as NA
  for(int i=styr_cpue_sTV; i<=endyr_cpue_sTV; i++) {
    if (obs_cpue_sTV(i) <=0) {
      posobs_cpue_sTV(i)=0;
      obs_cpue_sTV(i)=missing;
    } else {posobs_cpue_sTV(i)=1;}
  }
  // make corrections to observed cpue and cv to remove negatives and set to 1
  // note than pow(missing,0) = 1
  obs_cv_cpue_corrected_sTV=pow(obs_cv_cpue_sTV,posobs_cpue_sTV);
  obs_cpue_corrected_sTV=pow(obs_cpue_sTV,posobs_cpue_sTV);

  w_L=set_w_L;
  w_D=set_w_D;
  
  w_lenc_sBT=set_w_lenc_sBT;
  w_lenc_cHL=set_w_lenc_cHL;
  w_lenc_cPT=set_w_lenc_cPT;
  w_lenc_rHB=set_w_lenc_rHB;
  w_lenc_rHD=set_w_lenc_rHD;
  w_lenc_rGN=set_w_lenc_rGN;
  
  w_agec_sBT=set_w_agec_sBT;
  w_agec_sTV=set_w_agec_sTV;
  w_agec_cHL=set_w_agec_cHL;
  w_agec_cPT=set_w_agec_cPT;
  w_agec_rHB=set_w_agec_rHB;

  w_cpue_sTV=set_w_cpue_sTV;
  w_cpue_sBT=set_w_cpue_sBT;
  w_cpue_cHL=set_w_cpue_cHL;
  w_cpue_rHB=set_w_cpue_rHB;

  w_rec=set_w_rec;
  w_rec_early=set_w_rec_early;
  w_rec_end=set_w_rec_end;

  log_avg_F_L_cHL=set_log_avg_F_L_cHL(1);  
  log_avg_F_L_cPT=set_log_avg_F_L_cPT(1);
  log_avg_F_L_cTW=set_log_avg_F_L_cTW(1);  
  log_avg_F_L_rHB=set_log_avg_F_L_rHB(1);
  log_avg_F_L_rGN=set_log_avg_F_L_rGN(1);
  log_avg_F_L_rIN=set_log_avg_F_L_rIN(1);
  F_init_ratio=set_F_init_ratio(1);

  
  log_avg_F_D_cGD=set_log_avg_F_D_cGD(1); 
  log_avg_F_D_rHD=set_log_avg_F_D_rHD(1);
  log_avg_F_D_rGD=set_log_avg_F_D_rGD(1);
  log_avg_F_D_rID=set_log_avg_F_D_rID(1);

  log_dev_F_L_cHL=set_log_dev_vals_F_L_cHL;
  log_dev_F_L_cPT=set_log_dev_vals_F_L_cPT;
  log_dev_F_L_cTW=set_log_dev_vals_F_L_cTW;
  log_dev_F_L_rHB=set_log_dev_vals_F_L_rHB;
  log_dev_F_L_rGN=set_log_dev_vals_F_L_rGN;
  log_dev_F_L_rIN=set_log_dev_vals_F_L_rIN;
  log_dev_F_D_cGD=set_log_dev_vals_F_D_cGD;
  log_dev_F_D_rHD=set_log_dev_vals_F_D_rHD;
  log_dev_F_D_rGD=set_log_dev_vals_F_D_rGD;
  log_dev_F_D_rID=set_log_dev_vals_F_D_rID;

  len_cv_val=set_len_cv(1);

  log_R0=set_log_R0(1);

  selpar_A50_sBT= set_selpar_A50_sBT(1);
  selpar_slope_sBT=set_selpar_slope_sBT(1); 
  selpar_A502_sBT=set_selpar_A502_sBT(1);
  selpar_slope2_sBT=set_selpar_slope2_sBT(1); 

  selpar_A50_sTV=set_selpar_A50_sTV(1);
  selpar_slope_sTV=set_selpar_slope_sTV(1); 
  selpar_A502_sTV=set_selpar_A502_sTV(1);
  selpar_slope2_sTV=set_selpar_slope2_sTV(1); 


  selpar_A50b02_cHL=set_selpar_A50b02_cHL(1);
  selpar_slopeb02_cHL=set_selpar_slopeb02_cHL(1); 
  selpar_A50b03_cHL=set_selpar_A50b03_cHL(1);
  selpar_slopeb03_cHL=set_selpar_slopeb03_cHL(1); 
  selpar_A50b04_cHL=set_selpar_A50b04_cHL(1);
  selpar_slopeb04_cHL=set_selpar_slopeb04_cHL(1); 

  selpar_A50b02_cPT=set_selpar_A50b02_cPT(1);
  selpar_slopeb02_cPT=set_selpar_slopeb02_cPT(1); 
  selpar_A50b03_cPT=set_selpar_A50b03_cPT(1);
  selpar_slopeb03_cPT=set_selpar_slopeb03_cPT(1); 
  selpar_A50b04_cPT=set_selpar_A50b04_cPT(1);
  selpar_slopeb04_cPT=set_selpar_slopeb04_cPT(1);

  selpar_A50b01_rHB=set_selpar_A50b01_rHB(1);
  selpar_slopeb01_rHB=set_selpar_slopeb01_rHB(1); 
  selpar_A50b02_rHB=set_selpar_A50b02_rHB(1);
  selpar_slopeb02_rHB=set_selpar_slopeb02_rHB(1); 
  selpar_A50b03_rHB=set_selpar_A50b03_rHB(1);
  selpar_slopeb03_rHB=set_selpar_slopeb03_rHB(1); 
  selpar_A50b04_rHB=set_selpar_A50b04_rHB(1);
  selpar_slopeb04_rHB=set_selpar_slopeb04_rHB(1); 
  selpar_A50b05_rHB=set_selpar_A50b05_rHB(1);
  selpar_slopeb05_rHB=set_selpar_slopeb05_rHB(1);

  selpar_logitAge0_rHD=set_selpar_logitAge0_rHD(1);  
  selpar_logitAge1_rHD=set_selpar_logitAge1_rHD(1); 
  selpar_logitAge2_rHD=set_selpar_logitAge2_rHD(1);

  selpar_A50b04_rHD=set_selpar_A50b04_rHD(1);
  selpar_slopeb04_rHD=set_selpar_slopeb04_rHD(1); 
  selpar_A502b04_rHD=set_selpar_A502b04_rHD(1);
  selpar_slope2b04_rHD=set_selpar_slope2b04_rHD(1); 
  selpar_A50b05_rHD=set_selpar_A50b05_rHD(1);
  selpar_slopeb05_rHD=set_selpar_slopeb05_rHD(1);
  selpar_A502b05_rHD=set_selpar_A502b05_rHD(1);
  selpar_slope2b05_rHD=set_selpar_slope2b05_rHD(1);
   
 SSB_msy_n_out=0.0;
 SSB_msy_lb_out=0.0;

 iter_inc_msy=max_F_spr_msy/(n_iter_msy-1);
 iter_inc_spr=max_F_spr_msy/(n_iter_spr-1); 

 maturity_f=obs_maturity_f;
 prop_f=obs_prop_f;

 
//Fill in sample sizes of comps, possibly sampled in nonconsec yrs 
//Used primarily for output in R object   
  nsamp_allyr_lenc_sBT=missing;//"missing" defined in admb2r.cpp
  nsamp_allyr_lenc_cHL=missing;  
  nsamp_allyr_lenc_cPT=missing;
  nsamp_allyr_lenc_rHB=missing;
  nsamp_allyr_lenc_rHD=missing;
  nsamp_allyr_lenc_rGN=missing;
  nsamp_allyr_agec_sBT=missing;
  nsamp_allyr_agec_sTV=missing;  
  nsamp_allyr_agec_cHL=missing;
  nsamp_allyr_agec_cPT=missing;      
  nsamp_allyr_agec_rHB=missing;

  nfish_allyr_lenc_sBT=missing;//"missing" defined in admb2r.cpp
  nfish_allyr_lenc_cHL=missing;  
  nfish_allyr_lenc_cPT=missing;
  nfish_allyr_lenc_rHB=missing;
  nfish_allyr_lenc_rHD=missing;
  nfish_allyr_lenc_rGN=missing;
  nfish_allyr_agec_sBT=missing;
  nfish_allyr_agec_sTV=missing;  
  nfish_allyr_agec_cHL=missing;
  nfish_allyr_agec_cPT=missing;      
  nfish_allyr_agec_rHB=missing;
          

      for (iyear=1; iyear<=nyr_lenc_sBT; iyear++)
         {  if (nsamp_lenc_sBT(iyear)>maxSS_lenc)
             {nsamp_lenc_sBT(iyear)=maxSS_lenc;}
            if (nsamp_lenc_sBT(iyear)>=minSS_lenc)
             {nsamp_allyr_lenc_sBT(yrs_lenc_sBT(iyear))=nsamp_lenc_sBT(iyear);
              nfish_allyr_lenc_sBT(yrs_lenc_sBT(iyear))=nfish_lenc_sBT(iyear);}}
            
      for (iyear=1; iyear<=nyr_lenc_cHL; iyear++)
         {  if (nsamp_lenc_cHL(iyear)>maxSS_lenc)
             {nsamp_lenc_cHL(iyear)=maxSS_lenc;}
            if (nsamp_lenc_cHL(iyear)>=minSS_lenc)
             {nsamp_allyr_lenc_cHL(yrs_lenc_cHL(iyear))=nsamp_lenc_cHL(iyear);
              nfish_allyr_lenc_cHL(yrs_lenc_cHL(iyear))=nfish_lenc_cHL(iyear);}}
         
      for (iyear=1; iyear<=nyr_lenc_cPT; iyear++)
         {  if (nsamp_lenc_cPT(iyear)>maxSS_lenc)
             {nsamp_lenc_cPT(iyear)=maxSS_lenc;}
            if (nsamp_lenc_cPT(iyear)>=minSS_lenc)
             {nsamp_allyr_lenc_cPT(yrs_lenc_cPT(iyear))=nsamp_lenc_cPT(iyear);
              nfish_allyr_lenc_cPT(yrs_lenc_cPT(iyear))=nfish_lenc_cPT(iyear);}}
      
      for (iyear=1; iyear<=nyr_lenc_rHB; iyear++)
         {  if (nsamp_lenc_rHB(iyear)>maxSS_lenc)
             {nsamp_lenc_rHB(iyear)=maxSS_lenc;}
            if (nsamp_lenc_rHB(iyear)>=minSS_lenc)
             {nsamp_allyr_lenc_rHB(yrs_lenc_rHB(iyear))=nsamp_lenc_rHB(iyear);
              nfish_allyr_lenc_rHB(yrs_lenc_rHB(iyear))=nfish_lenc_rHB(iyear);}}

      for (iyear=1; iyear<=nyr_lenc_rHD; iyear++)
         {  if (nsamp_lenc_rHD(iyear)>maxSS_lenc)
             {nsamp_lenc_rHD(iyear)=maxSS_lenc;}
            if (nsamp_lenc_rHD(iyear)>=minSS_lenc)
             {nsamp_allyr_lenc_rHD(yrs_lenc_rHD(iyear))=nsamp_lenc_rHD(iyear);
              nfish_allyr_lenc_rHD(yrs_lenc_rHD(iyear))=nfish_lenc_rHD(iyear);}}
         
      for (iyear=1; iyear<=nyr_lenc_rGN; iyear++)
         {  if (nsamp_lenc_rGN(iyear)>maxSS_lenc)
             {nsamp_lenc_rGN(iyear)=maxSS_lenc;}
            if (nsamp_lenc_rGN(iyear)>=minSS_lenc)
             {nsamp_allyr_lenc_rGN(yrs_lenc_rGN(iyear))=nsamp_lenc_rGN(iyear);
              nfish_allyr_lenc_rGN(yrs_lenc_rGN(iyear))=nfish_lenc_rGN(iyear);}}

                  
      for (iyear=1; iyear<=nyr_agec_sBT; iyear++)
         {  if (nsamp_agec_sBT(iyear)>maxSS_agec)
             {nsamp_agec_sBT(iyear)=maxSS_agec;}
            if (nsamp_agec_sBT(iyear)>=minSS_agec)
             {nsamp_allyr_agec_sBT(yrs_agec_sBT(iyear))=nsamp_agec_sBT(iyear);
              nfish_allyr_agec_sBT(yrs_agec_sBT(iyear))=nfish_agec_sBT(iyear);}} 

      for (iyear=1; iyear<=nyr_agec_sTV; iyear++)
         {  if (nsamp_agec_sTV(iyear)>maxSS_agec)
             {nsamp_agec_sTV(iyear)=maxSS_agec;}
            if (nsamp_agec_sTV(iyear)>=minSS_agec)
             {nsamp_allyr_agec_sTV(yrs_agec_sTV(iyear))=nsamp_agec_sTV(iyear);
              nfish_allyr_agec_sTV(yrs_agec_sTV(iyear))=nfish_agec_sTV(iyear);}} 

      for (iyear=1; iyear<=nyr_agec_cHL; iyear++)
         {  if (nsamp_agec_cHL(iyear)>maxSS_agec)
             {nsamp_agec_cHL(iyear)=maxSS_agec;}
            if (nsamp_agec_cHL(iyear)>=minSS_agec)
             {nsamp_allyr_agec_cHL(yrs_agec_cHL(iyear))=nsamp_agec_cHL(iyear);
              nfish_allyr_agec_cHL(yrs_agec_cHL(iyear))=nfish_agec_cHL(iyear);}} 
            
      for (iyear=1; iyear<=nyr_agec_cPT; iyear++)
         {  if (nsamp_agec_cPT(iyear)>maxSS_agec)
             {nsamp_agec_cPT(iyear)=maxSS_agec;}
            if (nsamp_agec_cPT(iyear)>=minSS_agec)
             {nsamp_allyr_agec_cPT(yrs_agec_cPT(iyear))=nsamp_agec_cPT(iyear);
              nfish_allyr_agec_cPT(yrs_agec_cPT(iyear))=nfish_agec_cPT(iyear);}}
             
      for (iyear=1; iyear<=nyr_agec_rHB; iyear++)
         {  if (nsamp_agec_rHB(iyear)>maxSS_agec)
             {nsamp_agec_rHB(iyear)=maxSS_agec;}
            if (nsamp_agec_rHB(iyear)>=minSS_agec)
             {nsamp_allyr_agec_rHB(yrs_agec_rHB(iyear))=nsamp_agec_rHB(iyear);
              nfish_allyr_agec_rHB(yrs_agec_rHB(iyear))=nfish_agec_rHB(iyear);}}


//fill in Fs for msy and per-recruit analyses
  F_msy(1)=0.0;  
    for (ff=2;ff<=n_iter_msy;ff++) {F_msy(ff)=F_msy(ff-1)+iter_inc_msy;}
  F_spr(1)=0.0;  
    for (ff=2;ff<=n_iter_spr;ff++){F_spr(ff)=F_spr(ff-1)+iter_inc_spr;}


//fill in F's, Catch matrices, and log rec dev with zero's
  F_cHL.initialize(); 
  F_cPT.initialize(); 
  F_cTW.initialize(); 
  F_rHB.initialize(); 
  F_rGN.initialize();
  F_rIN.initialize();
  
  F_cGD.initialize(); 
  F_rHD.initialize(); 
  F_rGD.initialize();  
  F_rID.initialize();  
      
  L_cHL_num.initialize(); 
  L_cPT_num.initialize(); 
  L_cTW_num.initialize(); 
  L_rHB_num.initialize(); 
  L_rGN_num.initialize();
  L_rIN_num.initialize();
  
  D_cGD_num.initialize(); 
  D_rHD_num.initialize(); 
  D_rGD_num.initialize();
  D_rID_num.initialize();
      
  F_cHL_out.initialize(); 
  F_cPT_out.initialize(); 
  F_cTW_out.initialize(); 
  F_rHB_out.initialize(); 
  F_rGN_out.initialize();
  F_rIN_out.initialize();

  F_cGD_out.initialize(); 
  F_rHD_out.initialize(); 
  F_rGD_out.initialize();
  F_rID_out.initialize();
  
  pred_L_cHL_klb.initialize();
  pred_L_cPT_klb.initialize();
  pred_L_cTW_klb.initialize(); 
  pred_L_rHB_klb.initialize(); 
  pred_L_rGN_klb.initialize();
  pred_L_rIN_klb.initialize();
  pred_L_cHL_knum.initialize();
  pred_L_cPT_knum.initialize();
  pred_L_cTW_knum.initialize(); 
  pred_L_rHB_knum.initialize(); 
  pred_L_rGN_knum.initialize();
  pred_L_rIN_knum.initialize();
  
  pred_D_cGD_klb.initialize(); 
  pred_D_rHD_klb.initialize();
  pred_D_rGD_klb.initialize();
  pred_D_rID_klb.initialize();
  pred_D_cGD_knum.initialize();
  // pred_D_cHD_knum.initialize(); // Not used in modeling. empty vector passed to rdat
  // pred_D_cPD_knum.initialize(); // Not used in modeling. empty vector passed to rdat
  pred_D_rHD_knum.initialize(); 
  pred_D_rGD_knum.initialize();
  pred_D_rID_knum.initialize();
 
  sel_sTV.initialize(); 
  sel_sBT.initialize();  
  sel_cHL.initialize(); 
  sel_cPT.initialize(); 
  sel_cTW.initialize();  
  sel_rHB.initialize(); 
  sel_rGN.initialize();
  sel_rIN.initialize();
  sel_cGD.initialize(); 
  sel_rHD.initialize(); 
  sel_rGD.initialize();  
  sel_rID.initialize();  
  
  prob_belowsizelim_block1.initialize();
  prob_belowsizelim_block2.initialize();
  prob_belowsizelim_block3.initialize();
  prob_belowsizelim_block4.initialize();
  prob_belowsizelim_block5.initialize();
  
  log_rec_dev_output.initialize();
  log_rec_dev_output=missing;
  log_Nage_dev_output.initialize();
  log_dev_rec.initialize();
  log_dev_Nage.initialize();
  
 
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
TOP_OF_MAIN_SECTION
  time(&start);
  arrmblsize=20000000;
  gradient_structure::set_MAX_NVAR_OFFSET(1600);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(2000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(2000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(500);
  

//>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
PROCEDURE_SECTION 
 //cout<<"start"<<endl;

 // get_M_at_age(); 
 
 get_length_weight_at_age(); 
 //cout << "got length, weight, fecundity transitions" <<endl;
 get_reprod();
// cout << "got reprod" << endl;
 get_length_at_age_dist(); 
// cout<< "got predicted length at age distribution"<<endl;
 get_weight_at_age_landings();
// cout<< "got weight at age of landings"<<endl; 
 get_spr_F0();
// cout << "got F0 spr" << endl;
 get_selectivity(); 
// cout << "got selectivity" << endl;
 get_mortality(); 
// cout << "got mortalities" << endl;
 get_bias_corr(); 
// cout<< "got recruitment bias correction" << endl;
 get_numbers_at_age(); 
// cout << "got numbers at age" << endl;
 get_landings_numbers();
// cout << "got catch at age" << endl;
 get_landings_wgt();
// cout << "got landings" << endl;
 get_dead_discards();
 //cout << "got discards" << endl;
 get_catchability_fcns();
 //cout << "got catchability_fcns" << endl;
 get_indices();
 //cout << "got indices" << endl;
 get_length_comps();
 //cout<< "got length comps"<< endl;
 get_age_comps();
// cout<< "got age comps"<< endl;
 evaluate_objective_function();
 //cout << "objective function calculations complete" << endl;


// FUNCTION get_M_at_age
   // cum_surv_1plus=nages_M_sum*M_constant;  
   // M=M_lorenzen*(cum_surv_1plus/sum(M_lorenzen(min_M_age,nages)));

FUNCTION get_length_weight_at_age
  //compute mean length (mm) and weight (whole) at age
    meanlen_TL=Linf*(1.0-mfexp(-K*(agebins-t0+0.5)));    //total length in mm
    wgt_g=wgtpar_a*pow(meanlen_TL,wgtpar_b);             //wgt in grams 
    wgt_kg=g2kg*wgt_g;                                   //wgt in kilograms 
    wgt_mt=g2mt*wgt_g;                                   //mt of whole wgt: g2mt converts g to mt
    wgt_klb=mt2klb*wgt_mt;                               //1000 lb of whole wgt
    wgt_lb=mt2lb*wgt_mt;                                 //1000 lb of whole wgt
    fecundity=fecpar_batches*mfexp(fecpar_a+wgt_g*fecpar_b)/fecpar_scale;    //fecundity at age, scaled

FUNCTION get_reprod 
   //reprod is product of stuff going into reproductive capacity calcs
   reprod=elem_prod(wgt_klb,maturity_f);  
   reprod2=elem_prod(elem_prod(prop_f,maturity_f),wgt_klb);  
 
FUNCTION get_length_at_age_dist
  
//compute matrix of length at age, based on the normal distribution
  //len_cv=len_cv_val; 
  //len_sd=elem_prod(len_cv, meanlen_TL);
  for (iage=1;iage<=nages;iage++)
   {len_cv(iage)=len_cv_val;
    len_sd(iage)=meanlen_TL(iage)*len_cv(iage);
    zscore_lzero=(0.0-meanlen_TL(iage))/len_sd(iage); 
	cprob_lzero=cumd_norm(zscore_lzero);	 
	//population
    zscore_len=((lenbins(1)+0.5*lenbins_width)-meanlen_TL(iage)) / len_sd(iage);
    cprob_lenvec(1)=cumd_norm(zscore_len);          //includes any probability mass below zero
    lenprob(iage,1)=cprob_lenvec(1)-cprob_lzero;    //removes any probability mass below zero	

	//First size limit 8" Period 1 for both commercial and recreational discards (1984 through 1998)
	zscore_lsizelim1=(sizelim1-meanlen_TL(iage)) / len_sd(iage);
	cprob_lsizelim1=cumd_norm(zscore_lsizelim1);   //includes any probability mass below zero
	prob_belowsizelim_block1(iage)=	cprob_lsizelim1-cprob_lzero; //removes any probability mass below zero

	//Second size limit 10" Period 2 for both commercial and recreational discards and Period 3 for recreational landings (1999 through 2012 for commercial) (1999 through 2006 for recreational)
	zscore_lsizelim2=(sizelim2-meanlen_TL(iage)) / len_sd(iage);
	cprob_lsizelim2=cumd_norm(zscore_lsizelim2);   //includes any probability mass below zero
	prob_belowsizelim_block2(iage)=	cprob_lsizelim2-cprob_lzero; //removes any probability mass below zero

	//Third size limit 11" Period 3 for commercial (2013 through terminal year)
	zscore_lsizelim3=(sizelim3-meanlen_TL(iage)) / len_sd(iage);
	cprob_lsizelim3=cumd_norm(zscore_lsizelim3);                 //includes any probability mass below zero
	prob_belowsizelim_block3(iage)=	cprob_lsizelim3-cprob_lzero; //removes any probability mass below zero

	//Fourth size limit 12" Period 3 for recreational discards and Period 4 for recreational landings (2007 through 2012)
	zscore_lsizelim4=(sizelim4-meanlen_TL(iage)) / len_sd(iage);
	cprob_lsizelim4=cumd_norm(zscore_lsizelim4);                 //includes any probability mass below zero
	prob_belowsizelim_block4(iage)=	cprob_lsizelim4-cprob_lzero; //removes any probability mass below zero

	//Fifth size limit 13" Period 4 for recreational discards and Period 5 for recreational landings (2013 through the terminal year)
	zscore_lsizelim5=(sizelim5-meanlen_TL(iage)) / len_sd(iage);
	cprob_lsizelim5=cumd_norm(zscore_lsizelim5);                 //includes any probability mass below zero
	prob_belowsizelim_block5(iage)=	cprob_lsizelim5-cprob_lzero; //removes any probability mass below zero


	//most other length bins  
    //population
    for (ilen=2;ilen<nlenbins;ilen++)
      {
        zscore_len=((lenbins(ilen)+0.5*lenbins_width)-meanlen_TL(iage)) / len_sd(iage); 
		cprob_lenvec(ilen)=cumd_norm(zscore_len);
        lenprob(iage,ilen)=cprob_lenvec(ilen)-cprob_lenvec(ilen-1);
      }
	//last length bin is a plus group
	//population
    zscore_len=((lenbins(nlenbins)-0.5*lenbins_width)-meanlen_TL(iage)) / len_sd(iage); 
	lenprob(iage,nlenbins)=1.0-cumd_norm(zscore_len);
      lenprob(iage)=lenprob(iage)/(1.0-cprob_lzero);  //renormalize to account for any prob mass below size=0
 }
 
FUNCTION get_weight_at_age_landings
  // initialize temporary variable to speed code (NPK 2026-10-09)
  dvariable mean_len_mm = 0.0;

  //fleets under identical size limits are set equal at end of fcn
  for (iyear=styr; iyear<=endyr_sel_b01; iyear++)
  {
    len_cHL_mm(iyear)=meanlen_TL;
    wgt_cHL_klb(iyear)=wgt_klb;       

    for (iage=1;iage<=nages; iage++)
    {
    len_rGD_mm(iyear,iage)=sum(elem_prod(lenprob(iage),lenbins)); //assumes same size distn in period 1 as in period 2       
    }
    wgt_rGD_klb(iyear)=g2klb*wgtpar_a*pow(len_rGD_mm(iyear),wgtpar_b);
  } // end iyear loop
  
  
  for (iyear=(endyr_sel_b01+1); iyear<=endyr; iyear++)
  {
    for (iage=1;iage<=nages; iage++)
    {
    mean_len_mm = sum(elem_prod(lenprob(iage),lenbins));
      
    len_cHL_mm(iyear,iage)= mean_len_mm; 
    len_cGD_mm(iyear,iage)= mean_len_mm;
    len_rGD_mm(iyear,iage)= mean_len_mm;
    len_rID_mm(iyear,iage)= mean_len_mm;
    }
    wgt_cHL_klb(iyear)=g2klb*wgtpar_a*pow(len_cHL_mm(iyear),wgtpar_b);
    wgt_cGD_klb(iyear)=g2klb*wgtpar_a*pow(len_cGD_mm(iyear),wgtpar_b);  
    wgt_rGD_klb(iyear)=g2klb*wgtpar_a*pow(len_rGD_mm(iyear),wgtpar_b);
    wgt_rID_klb(iyear)=g2klb*wgtpar_a*pow(len_rID_mm(iyear),wgtpar_b);
  }

    //identical fleets set equal here (for speed)      
    len_cPT_mm=len_cHL_mm;  wgt_cPT_klb=wgt_cHL_klb;  
    len_cTW_mm=len_cHL_mm;  wgt_cTW_klb=wgt_cHL_klb;  
    len_rGN_mm=len_cHL_mm;  wgt_rGN_klb=wgt_cHL_klb;  
    len_rHB_mm=len_cHL_mm;  wgt_rHB_klb=wgt_cHL_klb; 
    len_rHD_mm=len_rGD_mm;  wgt_rHD_klb=wgt_rGD_klb; 
    
FUNCTION get_spr_F0
  //at mdyr, apply half this yr's mortality, half next yr's
  N_spr_F0(1)=1.0*mfexp(-1.0*M(1)*spawn_time_frac); //at peak spawning time
  N_bpr_F0(1)=1.0;      //at start of year
  for (iage=2; iage<=nages; iage++)
  {
    //N_spr_F0(iage)=N_spr_F0(iage-1)*mfexp(-1.0*(M(iage-1));
    N_spr_F0(iage)=N_spr_F0(iage-1)*
                   mfexp(-1.0*(M(iage-1)*(1.0-spawn_time_frac) + M(iage)*spawn_time_frac)); 
    N_bpr_F0(iage)=N_bpr_F0(iage-1)*mfexp(-1.0*(M(iage-1)));    
  }
  N_spr_F0(nages)=N_spr_F0(nages)/(1.0-mfexp(-1.0*M(nages))); //plus group (sum of geometric series)
  N_bpr_F0(nages)=N_bpr_F0(nages)/(1.0-mfexp(-1.0*M(nages)));
   
  spr_F0=sum(elem_prod(N_spr_F0,reprod));
  bpr_F0=sum(elem_prod(N_bpr_F0,wgt_klb));    
 

FUNCTION get_selectivity

// Logistic transformation of selectivity parameters
    selpar_Age0_rHD=1.0/(1.0+mfexp(-selpar_logitAge0_rHD));
    selpar_Age1_rHD=1.0/(1.0+mfexp(-selpar_logitAge1_rHD));
    selpar_Age2_rHD=1.0/(1.0+mfexp(-selpar_logitAge2_rHD));
    
// define selectivity functions (vectors by agebins) for each block    
  // sel b01
    sel_sTV_b01 = logistic_double(agebins, selpar_A50_sTV, selpar_slope_sTV,selpar_A502_sTV, selpar_slope2_sTV); 
    sel_sBT_b01 = logistic_double(agebins, selpar_A50_sBT, selpar_slope_sBT,selpar_A502_sBT, selpar_slope2_sBT);  
    sel_cHL_b01 = logistic(agebins, selpar_A50b02_cHL, selpar_slopeb02_cHL);
    sel_cPT_b01 = logistic(agebins, selpar_A50b02_cPT, selpar_slopeb02_cPT);
    sel_cTW_b01 = sel_cPT_b01;
    sel_rHB_b01 = logistic(agebins, selpar_A50b01_rHB, selpar_slopeb01_rHB); 
    sel_rGN_b01 = sel_rHB_b01;
    sel_rIN_b01 = sel_rHB_b01;
    // sel_cGD_b01; // Initialized above as vector of zeros
    sel_rHD_b01 = prob_belowsizelim_block1; // Default for all ages
      sel_rHD_b01(1) = selpar_Age0_rHD;     // Overwrite Age 0
      sel_rHD_b01(2) = selpar_Age1_rHD;     // Overwrite Age 1
      sel_rHD_b01(3) = selpar_Age2_rHD;     // Overwrite Age 2
      sel_rHD_b01(4) = 1.0;                 // Overwrite Age 3 
    sel_rGD_b01 = sel_rHD_b01;
    sel_rID_b01 = sel_rHD_b01;
    
  // sel b02
    sel_sTV_b02 = sel_sTV_b01; 
    sel_sBT_b02 = sel_sBT_b01;  
    sel_cHL_b02 = sel_cHL_b01;
    sel_cPT_b02 = sel_cPT_b01;
    sel_cTW_b02 = sel_cPT_b02;
    sel_rHB_b02 = logistic(agebins, selpar_A50b02_rHB, selpar_slopeb02_rHB);
    sel_rGN_b02 = sel_rHB_b02;
    sel_rIN_b02 = sel_rHB_b02;
    sel_cGD_b02 = sel_rHD_b01;
    sel_rHD_b02 = sel_rHD_b01;
    sel_rGD_b02 = sel_rHD_b02;
    sel_rID_b02 = sel_rHD_b02;

    // sel b03
    sel_sTV_b03 = sel_sTV_b01; 
    sel_sBT_b03 = sel_sBT_b01;  
    sel_cHL_b03 = logistic(agebins, selpar_A50b03_cHL, selpar_slopeb03_cHL);
    sel_cPT_b03 = logistic(agebins, selpar_A50b03_cPT, selpar_slopeb03_cPT);
    sel_cTW_b03 = sel_cPT_b03;
    sel_rHB_b03 = logistic(agebins, selpar_A50b03_rHB, selpar_slopeb03_rHB);
    sel_rGN_b03 = sel_rHB_b03;
    sel_rIN_b03 = sel_rHB_b03;
    sel_cGD_b03 = prob_belowsizelim_block3;
      sel_cGD_b03(1) = selpar_Age0_rHD;
      sel_cGD_b03(2) = selpar_Age1_rHD;
      sel_cGD_b03(3) = selpar_Age2_rHD;
      sel_cGD_b03(4) = 1.0;
    sel_rHD_b03 = prob_belowsizelim_block2;
      sel_rHD_b03(1) = selpar_Age0_rHD;
      sel_rHD_b03(2) = selpar_Age1_rHD;
      sel_rHD_b03(3) = selpar_Age2_rHD;
      sel_rHD_b03(4) = 1.0;
      sel_rGD_b03 = sel_rHD_b03;
      sel_rID_b03 = sel_rHD_b03;
    
    // sel b04
    sel_sTV_b04 = sel_sTV_b01; 
    sel_sBT_b04 = sel_sBT_b01;  
    sel_cHL_b04 = sel_cHL_b03;
    sel_cPT_b04 = sel_cPT_b03;
    sel_cTW_b04 = sel_cPT_b04;
    sel_rHB_b04 = logistic(agebins, selpar_A50b04_rHB, selpar_slopeb04_rHB);
    sel_rGN_b04 = sel_rHB_b04;
    sel_rIN_b04 = sel_rHB_b04;
    sel_cGD_b04 = sel_rHD_b03;
    sel_rHD_b04 = logistic_exponential(agebins, selpar_A50b04_rHD, selpar_slopeb04_rHD, selpar_A502b04_rHD, selpar_slope2b04_rHD);
    sel_rGD_b04 = sel_rHD_b04;
    sel_rID_b04 = sel_rHD_b04;
    
    // sel b05
    sel_sTV_b05 = sel_sTV_b01; 
    sel_sBT_b05 = sel_sBT_b01;  
    sel_cHL_b05 = logistic(agebins, selpar_A50b04_cHL, selpar_slopeb04_cHL);
    sel_cPT_b05 = logistic(agebins, selpar_A50b04_cPT, selpar_slopeb04_cPT);
    sel_cTW_b05 = sel_cPT_b05;
    sel_rHB_b05 = logistic(agebins, selpar_A50b05_rHB, selpar_slopeb05_rHB);
    sel_rGN_b05 = sel_rHB_b05;
    sel_rIN_b05 = sel_rHB_b05;
    sel_cGD_b05 = sel_rHD_b03;
    sel_rHD_b05 = logistic_exponential(agebins, selpar_A50b05_rHD, selpar_slopeb05_rHD, selpar_A502b05_rHD, selpar_slope2b05_rHD);
    sel_rGD_b05 = sel_rHD_b05; 
    sel_rID_b05 = sel_rHD_b05; 

  
  // sel b01
  for (iyear=styr; iyear<=endyr_sel_b01; iyear++)
  {
    sel_sBT(iyear) = sel_sBT_b01; 
    sel_sTV(iyear) = sel_sTV_b01;
    sel_cHL(iyear) = sel_cHL_b01;
    sel_cPT(iyear) = sel_cPT_b01;
    sel_cTW(iyear) = sel_cTW_b01;  
    sel_rHB(iyear) = sel_rHB_b01; 
    sel_rGN(iyear) = sel_rGN_b01;
    sel_rIN(iyear) = sel_rIN_b01;
    
    sel_rHD(iyear) = sel_rHD_b01; 
    sel_rGD(iyear) = sel_rGD_b01;
    sel_rID(iyear) = sel_rID_b01;
    sel_cGD(iyear) = sel_cGD_b01; 
  }

  // sel b02
  for (iyear=endyr_sel_b01+1; iyear<=endyr_sel_b02; iyear++)
  {
     sel_sBT(iyear) = sel_sBT_b02; 
     sel_sTV(iyear) = sel_sTV_b02;    
     sel_cHL(iyear) = sel_cHL_b02;
     sel_cPT(iyear) = sel_cPT_b02;
     sel_cTW(iyear) = sel_cTW_b02;
     sel_rHB(iyear) = sel_rHB_b02;
     sel_rGN(iyear) = sel_rGN_b02;
     sel_rIN(iyear) = sel_rIN_b02;

     sel_rHD(iyear) = sel_rHD_b02;
     sel_rGD(iyear) = sel_rGD_b02;
     sel_rID(iyear) = sel_rID_b02;
     sel_cGD(iyear) = sel_cGD_b02;
     
  }

  // sel b03
  for (iyear=endyr_sel_b02+1; iyear<=endyr_sel_b03; iyear++)
  {
     sel_sBT(iyear) = sel_sBT_b03; 
     sel_sTV(iyear) = sel_sTV_b03;    
     sel_cHL(iyear) = sel_cHL_b03; 
     sel_cPT(iyear) = sel_cPT_b03;
     sel_cTW(iyear) = sel_cTW_b03;     
     sel_rHB(iyear) = sel_rHB_b03;
     sel_rGN(iyear) = sel_rGN_b03;
     sel_rIN(iyear) = sel_rIN_b03;

     sel_rHD(iyear) = sel_rHD_b03;
     sel_rGD(iyear) = sel_rGD_b03;
     sel_rID(iyear) = sel_rID_b03;
     sel_cGD(iyear) = sel_rHD_b03;     
  }  
  
  // sel b04
  for (iyear=endyr_sel_b03+1; iyear<=endyr_sel_b04; iyear++)
  {
     sel_sBT(iyear) = sel_sBT_b04; 
     sel_sTV(iyear) = sel_sTV_b04;    
     sel_cHL(iyear) = sel_cHL_b04; 
     sel_cPT(iyear) = sel_cPT_b04;
     sel_cTW(iyear) = sel_cTW_b04;     
     sel_rHB(iyear) = sel_rHB_b04;
     sel_rGN(iyear) = sel_rGN_b04;
     sel_rIN(iyear) = sel_rIN_b04;
     
     sel_rHD(iyear) = sel_rHD_b04;
     sel_rGD(iyear) = sel_rGD_b04;
     sel_rID(iyear) = sel_rID_b04;
     sel_cGD(iyear) = sel_cGD_b04;       
  }
  // sel b05
  for (iyear=endyr_sel_b04+1; iyear<=endyr; iyear++)
  {
     sel_sBT(iyear) = sel_sBT_b05; 
     sel_sTV(iyear) = sel_sTV_b05;    
     sel_cHL(iyear) = sel_cHL_b05; 
     sel_cPT(iyear) = sel_cPT_b05;
     sel_cTW(iyear) = sel_cTW_b05;     
     sel_rHB(iyear) = sel_rHB_b05;
     sel_rGN(iyear) = sel_rGN_b05;
     sel_rIN(iyear) = sel_rIN_b05;
     
     sel_rHD(iyear) = sel_rHD_b05;
     sel_rGD(iyear) = sel_rGD_b05; 
     sel_rID(iyear) = sel_rID_b05; 
     sel_cGD(iyear) = sel_cGD_b05;           
  }
  
  //// Selectivity blocks that only apply to cGD.
  //// 2026-05-12 NPK It would be more consistent with BAM practices to create two more little selectivity blocks for all fleets, and most fleets just wouldn't really use them. But I'm hoping to simplify the selectivity structure of this model and just eliminate these blocks since commercial discards are so minimal and selectivity is mostly based on rHD, so I'm not going to do that at this point.
  // //Discard quota: overwrite with wghted average for period 3 styr cGD closed=2009 through 2013
  // for (iyear=styr_closure_com; iyear<= (endyr_sel_b04+1); iyear++)
  // {sel_cGD(iyear)=Dprop_cGN_sel_D4*sel_rHD_b03 + Dprop_cGN_sel_cHL4*sel_cHL_b03 +
  //                    Dprop_cGN_sel_cPT4*sel_rHD_b03; //Treat the discards from the cPT closed season as if they were discards during open season 
  //  sel_cGD(iyear)=sel_cGD(iyear)/max(sel_cGD(iyear));}
  //  
  // // For period after close of cHL starts in 2014
  // for (iyear=(endyr_sel_b04+2); iyear<=endyr; iyear++)
  // {
  //    sel_cGD(iyear)= sel_cGD_b03;
  //  } 

  //set selectivities that mirror others
  //sel_cTW=sel_cPT;  
  //sel_rGD=sel_rHD;
  //sel_cGD=sel_rHD; // Mirroring sel_rHD has almost no effect on the TotalLikelihood. It changes by 0.3 units.


FUNCTION get_mortality
  Fsum.initialize();
  Fapex.initialize();
  F.initialize();
  //initialization F is avg from first 3 yrs of observed landings 
  log_dev_F_init_cHL=sum(log_dev_F_L_cHL(styr_L_cHL,(styr_L_cHL+2)))/3.0;
  log_dev_F_init_cPT=sum(log_dev_F_L_cPT(styr_L_cPT,(styr_L_cPT+2)))/3.0;
  log_dev_F_init_cTW=sum(log_dev_F_L_cTW(styr_L_cTW,(styr_L_cTW+2)))/3.0;
  log_dev_F_init_rHB=sum(log_dev_F_L_rHB(styr_L_rHB,(styr_L_rHB+2)))/3.0;
  log_dev_F_init_rGN=sum(log_dev_F_L_rGN(styr_L_rGN,(styr_L_rGN+2)))/3.0;
  log_dev_F_init_rIN=sum(log_dev_F_L_rIN(styr_L_rIN,(styr_L_rIN+2)))/3.0;
  log_dev_F_init_rGD=sum(log_dev_F_D_rGD(styr_D_rGD,(styr_D_rGD+2)))/3.0;  
  log_dev_F_init_rID=sum(log_dev_F_D_rID(styr_D_rID,(styr_D_rID+2)))/3.0;
  
  log_dev_F_cGD2=sum(log_dev_F_D_cGD(styr_D_cGD,(styr_D_cGD+5)))/6.0; //for cGN D 1984-1992
  
  //cout<<styr<<endl;  
  for (iyear=styr; iyear<=endyr; iyear++) 
  {
    //-------------
    if(iyear>=styr_L_cHL & iyear<=endyr_L_cHL)
    {  F_cHL_out(iyear)=mfexp(log_avg_F_L_cHL+log_dev_F_L_cHL(iyear)); //}    
       //if (iyear<styr_L_cHL){F_cHL_out(iyear)=mfexp(log_avg_F_L_cHL+log_dev_F_init_cHL);}        
       F_cHL(iyear)=sel_cHL(iyear)*F_cHL_out(iyear);
       Fsum(iyear)+=F_cHL_out(iyear);
    }
    
    //-------------
    if(iyear>=styr_L_cPT & iyear<=endyr_L_cPT)
    {  F_cPT_out(iyear)=mfexp(log_avg_F_L_cPT+log_dev_F_L_cPT(iyear)); //}
       //if (iyear<styr_L_cPT) {F_cPT_out(iyear)=0.0;}
       F_cPT(iyear)=sel_cPT(iyear)*F_cPT_out(iyear);
       Fsum(iyear)+=F_cPT_out(iyear);
    }

    //-------------
    if(iyear>=styr_L_cTW & iyear<=endyr_L_cTW)
    {  F_cTW_out(iyear)=mfexp(log_avg_F_L_cTW+log_dev_F_L_cTW(iyear)); //}
      // if (iyear<styr_L_cTW) {F_cTW_out(iyear)=0.0;}
       F_cTW(iyear)=sel_cTW(iyear)*F_cTW_out(iyear);
       Fsum(iyear)+=F_cTW_out(iyear);
    } 
        
    //-------------
    if(iyear>=styr_L_rHB & iyear<=endyr_L_rHB)
    {  F_rHB_out(iyear)=mfexp(log_avg_F_L_rHB+log_dev_F_L_rHB(iyear));//}
    // if (iyear<styr_L_rHB){F_rHB_out(iyear)=mfexp(log_avg_F_L_rHB+log_F_init_rHB);}
       F_rHB(iyear)=sel_rHB(iyear)*F_rHB_out(iyear);    
       Fsum(iyear)+=F_rHB_out(iyear);
    }
    
    //-------------
    if(iyear>=styr_L_rGN & iyear<=endyr_L_rGN)
       {F_rGN_out(iyear)=mfexp(log_avg_F_L_rGN+log_dev_F_L_rGN(iyear));}
    if (iyear<styr_L_rGN){F_rGN_out(iyear)=mfexp(log_avg_F_L_rGN+log_dev_F_init_rGN);}
    F_rGN(iyear)=sel_rGN(iyear)*F_rGN_out(iyear);    
    Fsum(iyear)+=F_rGN_out(iyear);
    
    //-------------
    if(iyear>=styr_L_rIN & iyear<=endyr_L_rIN)
    {F_rIN_out(iyear)=mfexp(log_avg_F_L_rIN+log_dev_F_L_rIN(iyear));}
    if (iyear<styr_L_rIN){F_rIN_out(iyear)=mfexp(log_avg_F_L_rIN+log_dev_F_init_rIN);}
    F_rIN(iyear)=sel_rIN(iyear)*F_rIN_out(iyear);    
    Fsum(iyear)+=F_rIN_out(iyear);    
    
    //discards-------------

    if(iyear>=styr_D_cGD & iyear<=endyr_D_cGD)
      {F_cGD_out(iyear)=mfexp(log_avg_F_D_cGD+log_dev_F_D_cGD(iyear));}
    if(iyear > endyr_sel_b01 & iyear < styr_D_cGD)
      {F_cGD_out(iyear)=mfexp(log_avg_F_D_cGD+log_dev_F_cGD2);}
    F_cGD(iyear)=sel_cGD(iyear)*F_cGD_out(iyear);  
    Fsum(iyear)+=F_cGD_out(iyear);    


    if(iyear>=styr_D_rHD & iyear<=endyr_D_rHD)
    {  F_rHD_out(iyear)=mfexp(log_avg_F_D_rHD+log_dev_F_D_rHD(iyear)); 
       F_rHD(iyear)=sel_rHD(iyear)*F_rHD_out(iyear); 
       Fsum(iyear)+=F_rHD_out(iyear); 
    }
    
    if(iyear<styr_D_rGD)
       {F_rGD_out(iyear)=mfexp(log_avg_F_D_rGD+log_dev_F_init_rGD);}
    if(iyear>=styr_D_rGD& iyear<=endyr_D_rGD)
       { F_rGD_out(iyear)=mfexp(log_avg_F_D_rGD+log_dev_F_D_rGD(iyear));}
    F_rGD(iyear)=sel_rGD(iyear)*F_rGD_out(iyear);  
    Fsum(iyear)+=F_rGD_out(iyear);
    
    if(iyear<styr_D_rID)
    {F_rID_out(iyear)=mfexp(log_avg_F_D_rID+log_dev_F_init_rID);}
    if(iyear>=styr_D_rID& iyear<=endyr_D_rID)
    { F_rID_out(iyear)=mfexp(log_avg_F_D_rID+log_dev_F_D_rID(iyear));}
    F_rID(iyear)=sel_rID(iyear)*F_rID_out(iyear);  
    Fsum(iyear)+=F_rID_out(iyear);    

    //Total F at age
    F(iyear)=F_cHL(iyear); //first in additive series (NO +=)
    F(iyear)+=F_cPT(iyear);
    F(iyear)+=F_cTW(iyear);    
    F(iyear)+=F_rHB(iyear);
    F(iyear)+=F_rGN(iyear);
    F(iyear)+=F_rIN(iyear);
   
    F(iyear)+=F_cGD(iyear);
    F(iyear)+=F_rHD(iyear);
    F(iyear)+=F_rGD(iyear);
    F(iyear)+=F_rID(iyear);
    
    Fapex(iyear)=max(F(iyear));
   // cout <<"Fapex" << Fapex <<endl;
    Z(iyear)=M+F(iyear);
  }  //end iyear
 
 
FUNCTION get_bias_corr
  var_rec_dev=norm2(log_dev_rec(styr_rec_dev,endyr_rec_dev)-
              sum(log_dev_rec(styr_rec_dev,endyr_rec_dev))
              /(nyr_rec_dev))/(nyr_rec_dev); 
			  //cout<<"var rec dev"<<var_rec_dev<<endl;
  //if (set_BiasCor <= 0.0) {BiasCor=mfexp(var_rec_dev/2.0);}   //bias correction
  rec_sigma_sq=square(rec_sigma);
  rec_sigma_sqd2=rec_sigma_sq/2.0;
  if (set_BiasCor <= 0.0) {BiasCor=mfexp(rec_sigma_sqd2);}   //bias correction               
  else {BiasCor=set_BiasCor;}

//  cout<< "get bias corr" << var_rec_dev <<endl;

FUNCTION get_numbers_at_age
//Initialization
  R0=mfexp(log_R0);
  S0=spr_F0*R0;
  R_virgin=SR_eq_func(R0, steep, spr_F0, spr_F0, BiasCor, SR_switch);
  B0=bpr_F0*R_virgin;   
  B0_q_DD=R_virgin*sum(elem_prod(N_bpr_F0(set_q_DD_stage,nages),wgt_klb(set_q_DD_stage,nages))); 

  F_initial=sel_cHL(styr)*mfexp(log_avg_F_L_cHL+log_dev_F_init_cHL)+
            sel_cPT(styr)*mfexp(log_avg_F_L_cPT+log_dev_F_init_cPT)+
            sel_cTW(styr)*mfexp(log_avg_F_L_cTW+log_dev_F_init_cTW)+
            sel_rHB(styr)*mfexp(log_avg_F_L_rHB+log_dev_F_init_rHB)+
            sel_rGN(styr)*mfexp(log_avg_F_L_rGN+log_dev_F_init_rGN)+
            sel_rIN(styr)*mfexp(log_avg_F_L_rIN+log_dev_F_init_rIN)+            
            sel_rGD(styr)*mfexp(log_avg_F_D_rGD+log_dev_F_init_rGD)+
            sel_rID(styr)*mfexp(log_avg_F_D_rID+log_dev_F_init_rID);
  Z_initial=M+F_init_ratio*F_initial;


//Initial equilibrium age structure
  N_spr_initial(1)=1.0*mfexp(-1.0*Z_initial(1)*spawn_time_frac); //at peak spawning time;
  for (iage=2; iage<=nages; iage++)
    {
      N_spr_initial(iage)=N_spr_initial(iage-1)*
                   mfexp(-1.0*(Z_initial(iage-1)*(1.0-spawn_time_frac) + Z_initial(iage)*spawn_time_frac)); 
    }
  N_spr_initial(nages)=N_spr_initial(nages)/(1.0-mfexp(-1.0*Z_initial(nages))); //plus group
//    N_spr_F_init_mdyr(1,(nages-1))=elem_prod(N_spr_initial(1,(nages-1)),
//                                   mfexp((-1.*(M(nages-1)+ F_initial))/2.0));   

  spr_initial=sum(elem_prod(N_spr_initial,reprod));

  if (styr==styr_rec_dev) {R1=SR_eq_func(R0, steep, spr_F0, spr_initial, 1.0, SR_switch);} //without bias correction (deviation added later)
  else {R1=SR_eq_func(R0, steep, spr_F0, spr_initial, BiasCor, SR_switch);} //with bias correction                 

  if(R1<10.0) {R1=10.0;} //Avoid negative (or unreasonably low) popn sizes during search algorithm

  
   
//Compute equilibrium age structure for first year
  N_initial_eq(1)=R1;
  for (iage=2; iage<=nages; iage++)
  {
    N_initial_eq(iage)=N_initial_eq(iage-1)*
        mfexp(-1.0*(Z_initial(iage-1)));    
  }
  //plus group calculation
  N_initial_eq(nages)=N_initial_eq(nages)/(1.0-mfexp(-1.0*Z_initial(nages))); //plus group
  
//Add deviations to initial equilibrium N
  N(styr)(2,nages)=elem_prod(N_initial_eq(2,nages),mfexp(log_dev_Nage));
   
  if (styr==styr_rec_dev) {N(styr,1)=N_initial_eq(1)*mfexp(log_dev_rec(styr_rec_dev));}
  else {N(styr,1)=N_initial_eq(1);}
  
  N_mdyr(styr)(1,nages)=elem_prod(N(styr)(1,nages),(mfexp(-1.*(Z_initial(1,nages))*0.5))); //mid year 
  N_spawn(styr)(2,nages)=elem_prod(N(styr)(2,nages),(mfexp(-1.*(Z_initial(2,nages))*spawn_time_frac))); //peak spawning time 

  SSB(styr)=sum(elem_prod(N_spawn(styr),reprod));
  MatFemB(styr)=sum(elem_prod(N_spawn(styr),reprod2));
  B_q_DD(styr)=sum(elem_prod(N(styr)(set_q_DD_stage,nages),wgt_klb(set_q_DD_stage,nages)));
    
//Rest of years 
  for (iyear=styr; iyear<endyr; iyear++)
  {
    if(iyear<(styr_rec_dev-1)) //recruitment follows S-R curve exactly
    {
        N(iyear+1,1)=0.0;   //no age 0's mature in SSB calculations, value replaced below for abundance calcs
        N(iyear+1)(2,nages)=++elem_prod(N(iyear)(1,nages-1),(mfexp(-1.*Z(iyear)(1,nages-1))));
        N(iyear+1,nages)+=N(iyear,nages)*mfexp(-1.*Z(iyear,nages));//plus group
        N_spawn(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*spawn_time_frac))); //peak spawning time 
        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod)); 
        MatFemB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod2));   
        B_q_DD(iyear+1)=sum(elem_prod(N(iyear+1)(set_q_DD_stage,nages),wgt_klb(set_q_DD_stage,nages)));   
        N(iyear+1,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB(iyear),SR_switch);
        N_mdyr(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.5))); //mid year                    
    }
    else if(iyear>(endyr_rec_dev-1)) //recruitment follows S-R curve (with bias correction) exactly or else recent mean recruitment
    {   
        if (iyear==endyr_rec_dev){// rec_mu=sum(trans(N)(1)(styr_rec_mu,endyr_rec_dev))/nyr_rec_mu;  //Arithmetic mean of recruitment values
	    rec_mu=sum(log_dev_rec(styr_rec_mu,endyr_rec_dev))/nyr_rec_mu; //arithmetic mean of deviates
	    // //geometric mean of deviates Doesn't work!!
	    // rec_mu_temp=1.0;
	    // for (iyear=styr_rec_mu; iyear<=endyr_rec_dev; iyear++) {rec_mu_temp*=log_dev_rec(iyear);}
	    // if (rec_mu_temp ==0.0){
	    //     rec_mu=0.0;
	    // } else {
 	    //     rec_mu=pow(rec_mu_temp,(1.0/nyr_rec_mu));
	    // }
	}		
		N(iyear+1,1)=0.0;   //no age 0's mature in SSB calculations, value replaced below for abundance calcs
		N(iyear+1)(2,nages)=++elem_prod(N(iyear)(1,nages-1),(mfexp(-1.*Z(iyear)(1,nages-1))));
        N(iyear+1,nages)+=N(iyear,nages)*mfexp(-1.*Z(iyear,nages)); //plus group
        N_spawn(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*spawn_time_frac))); //peak spawning time
        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod));
        MatFemB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod2));  
	B_q_DD(iyear+1)=sum(elem_prod(N(iyear+1)(set_q_DD_stage,nages),wgt_klb(set_q_DD_stage,nages)));
	N(iyear+1,1)=SR_func(R0, steep, spr_F0, SSB(iyear),SR_switch)*mfexp(rec_mu);
        N_mdyr(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.5))); //mid year
    } 
    else   //recruitment follows S-R curve with lognormal deviation
    {       
        N(iyear+1,1)=0.0;   //no age 0's mature in SSB calculations, value replaced below for abundance calcs
        N(iyear+1)(2,nages)=++elem_prod(N(iyear)(1,nages-1),(mfexp(-1.*Z(iyear)(1,nages-1))));
        N(iyear+1,nages)+=N(iyear,nages)*mfexp(-1.*Z(iyear,nages));//plus group
        N_spawn(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*spawn_time_frac))); //peak spawning time 
        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod));
        MatFemB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod2));
        B_q_DD(iyear+1)=sum(elem_prod(N(iyear+1)(set_q_DD_stage,nages),wgt_klb(set_q_DD_stage,nages)));
        N(iyear+1,1)=SR_func(R0, steep, spr_F0, SSB(iyear),SR_switch)*mfexp(log_dev_rec(iyear+1));
        N_mdyr(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.5))); //mid year 
        
    }
  } //end iyear

  rec_mean=sum(rec(styr_rec_spr, endyr_rec_spr))/nyr_rec_spr;

//Time series of interest
  rec=column(N,1);
  SdS0=SSB/S0;
  
FUNCTION get_landings_numbers //Baranov catch eqn
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    for (iage=1; iage<=nages; iage++)
    {
      L_cHL_num(iyear,iage)=N(iyear,iage)*F_cHL(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      L_cPT_num(iyear,iage)=N(iyear,iage)*F_cPT(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      L_cTW_num(iyear,iage)=N(iyear,iage)*F_cTW(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      L_rHB_num(iyear,iage)=N(iyear,iage)*F_rHB(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      L_rGN_num(iyear,iage)=N(iyear,iage)*F_rGN(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      L_rIN_num(iyear,iage)=N(iyear,iage)*F_rIN(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);      
    }
                
    pred_L_cHL_knum(iyear)=sum(L_cHL_num(iyear))/1000.0;
    pred_L_cPT_knum(iyear)=sum(L_cPT_num(iyear))/1000.0;
    pred_L_cTW_knum(iyear)=sum(L_cTW_num(iyear))/1000.0;    
    pred_L_rHB_knum(iyear)=sum(L_rHB_num(iyear))/1000.0;
    pred_L_rGN_knum(iyear)=sum(L_rGN_num(iyear))/1000.0;
    pred_L_rIN_knum(iyear)=sum(L_rIN_num(iyear))/1000.0;
  }
 
 
FUNCTION get_landings_wgt

////---Predicted landings------------------------
  for (iyear=styr; iyear<=endyr; iyear++)
  {    
    L_cHL_klb(iyear)=elem_prod(L_cHL_num(iyear),wgt_cHL_klb(iyear));     //in 1000 lb
    L_cPT_klb(iyear)=elem_prod(L_cPT_num(iyear),wgt_cPT_klb(iyear));     //in 1000 lb
    L_cTW_klb(iyear)=elem_prod(L_cTW_num(iyear),wgt_cTW_klb(iyear));     //in 1000 lb    
    L_rHB_klb(iyear)=elem_prod(L_rHB_num(iyear),wgt_rHB_klb(iyear));     //in 1000 lb
    L_rGN_klb(iyear)=elem_prod(L_rGN_num(iyear),wgt_rGN_klb(iyear));     //in 1000 lb  
    L_rIN_klb(iyear)=elem_prod(L_rIN_num(iyear),wgt_rIN_klb(iyear));     //in 1000 lb  

    pred_L_cHL_klb(iyear)=sum(L_cHL_klb(iyear));
    pred_L_cPT_klb(iyear)=sum(L_cPT_klb(iyear));
    pred_L_cTW_klb(iyear)=sum(L_cTW_klb(iyear));
    pred_L_rHB_klb(iyear)=sum(L_rHB_klb(iyear));
    pred_L_rGN_klb(iyear)=sum(L_rGN_klb(iyear));
    pred_L_rIN_klb(iyear)=sum(L_rIN_klb(iyear));    
  }
 
   
    
FUNCTION get_dead_discards //Baranov catch eqn
  //dead discards at age (number fish) 

  for (iyear=styr; iyear<=endyr; iyear++)
  {
    for (iage=1; iage<=nages; iage++)
    {
      D_cGD_num(iyear,iage)=N(iyear,iage)*F_cGD(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      D_rHD_num(iyear,iage)=N(iyear,iage)*F_rHD(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      D_rGD_num(iyear,iage)=N(iyear,iage)*F_rGD(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      D_rID_num(iyear,iage)=N(iyear,iage)*F_rID(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);        
    }
    pred_D_cGD_knum(iyear)=sum(D_cGD_num(iyear))/1000.0;         //pred annual dead discards in 1000s (for matching data)
    pred_D_cGD_klb(iyear)=sum(elem_prod(D_cGD_num(iyear),wgt_cGD_klb(iyear)));  //annual dead discards in 1000 lb (for output only)
    
    // pred_D_cHD_knum(iyear)=missing; // Not used in modeling. missing fills empty vector with NA in rdat
    // pred_D_cPD_knum(iyear)=missing; // Not used in modeling. missing empty vector with NA in rdat

    pred_D_rHD_knum(iyear)=sum(D_rHD_num(iyear))/1000.0;             //pred annual dead discards in 1000s (for matching data)
    pred_D_rHD_klb(iyear)=sum(elem_prod(D_rHD_num(iyear),wgt_rHD_klb(iyear)));        //annual dead discards in 1000 lb (for output only)

    pred_D_rGD_knum(iyear)=sum(D_rGD_num(iyear))/1000.0;         //pred annual dead discards in 1000s (for matching data)
    pred_D_rGD_klb(iyear)=sum(elem_prod(D_rGD_num(iyear),wgt_rGD_klb(iyear)));  //annual dead discards in 1000 lb (for output only)
    
    pred_D_rID_knum(iyear)=sum(D_rID_num(iyear))/1000.0;         //pred annual dead discards in 1000s (for matching data)
    pred_D_rID_klb(iyear)=sum(elem_prod(D_rID_num(iyear),wgt_rID_klb(iyear)));  //annual dead discards in 1000 lb (for output only)
  }

FUNCTION get_catchability_fcns    
 //Get rate increase if estimated, otherwise fixed above
  if (q_rate_phase > 0.0)
  {
      for (iyear=styr_cpue_cHL; iyear<=endyr_cpue_cHL; iyear++)
      {   if (iyear>styr_cpue_cHL & iyear <=2003) 
          {
             q_rate_fcn_cHL(iyear)=(1.0+(iyear-styr_cpue_cHL)*q_rate)*q_rate_fcn_cHL(styr_cpue_cHL); 
          }
          if (iyear>2003) {q_rate_fcn_cHL(iyear)=q_rate_fcn_cHL(iyear-1);} 
      }   
      for (iyear=styr_cpue_rHB; iyear<=endyr_cpue_rHB; iyear++)
      {   if (iyear>styr_cpue_rHB & iyear <=2003) 
          {
             q_rate_fcn_rHB(iyear)=(1.0+(iyear-styr_cpue_rHB)*q_rate)*q_rate_fcn_rHB(styr_cpue_rHB);  
          }
          if (iyear>2003) {q_rate_fcn_rHB(iyear)=q_rate_fcn_rHB(iyear-1);} 
      }
  } //end q_rate conditional      

 //Get density dependence scalar (=1.0 if density independent model is used)   
  if (q_DD_beta>0.0) 
  {
    B_q_DD+=dzero;
    for (iyear=styr;iyear<=endyr;iyear++)
        {q_DD_fcn(iyear)=pow(B0_q_DD,q_DD_beta)*pow(B_q_DD(iyear),-q_DD_beta);}
          //{q_DD_fcn(iyear)=1.0+4.0/(1.0+mfexp(0.75*(B_q_DD(iyear)-0.1*B0_q_DD))); }
  }  
  
      
FUNCTION get_indices
//---Predicted CPUEs------------------------
 //Survey 1: sBT
  for (iyear=styr_cpue_sBT; iyear<=endyr_cpue_sBT; iyear++)
  {   //index in number units
      N_sBT(iyear)=elem_prod(N_mdyr(iyear),sel_sBT(iyear)); 
      pred_cpue_sBT(iyear)=mfexp(log_q_sBT)*sum(N_sBT(iyear));
  }

 //Survey 2: MCVT
  for (iyear=styr_cpue_sTV; iyear<=endyr_cpue_sTV; iyear++)
  {   //index in number units
      N_sTV(iyear)=elem_prod(N_mdyr(iyear),sel_sTV(iyear)); 
      pred_cpue_sTV(iyear)=mfexp(log_q_sTV)*sum(N_sTV(iyear));
  }
  
  
 //Commercial handline cpue
  q_cHL(styr_cpue_cHL)=mfexp(log_q_cHL); 
  for (iyear=styr_cpue_cHL; iyear<=endyr_cpue_cHL; iyear++)
  {   //index in weight units. original index in lb and re-scaled. predicted in klb, but difference is absorbed by q
      N_cHL(iyear)=elem_prod(elem_prod(N_mdyr(iyear),sel_cHL(iyear)),wgt_cHL_klb(iyear)); 
      pred_cpue_cHL(iyear)=q_cHL(iyear)*q_rate_fcn_cHL(iyear)*q_DD_fcn(iyear)*sum(N_cHL(iyear));
      if (iyear<endyr_cpue_cHL){q_cHL(iyear+1)=q_cHL(iyear)*mfexp(q_RW_log_dev_cHL(iyear));}
  }

 //Headboat cpue
  q_rHB(styr_cpue_rHB)=mfexp(log_q_rHB);
  for (iyear=styr_cpue_rHB; iyear<=endyr_cpue_rHB; iyear++)
  {   //index in weight units. original index in lb and re-scaled. predicted in klb, but difference is absorbed by q
      N_rHB(iyear)=elem_prod(elem_prod(N_mdyr(iyear),sel_rHB(iyear)),wgt_rHB_klb(iyear)); 
      pred_cpue_rHB(iyear)=q_rHB(iyear)*q_rate_fcn_rHB(iyear)*q_DD_fcn(iyear)*sum(N_rHB(iyear));
      if (iyear<endyr_cpue_rHB){q_rHB(iyear+1)=q_rHB(iyear)*mfexp(q_RW_log_dev_rHB(iyear));}
  }
  //rHD cpue
  //q_rHD(styr_rHD_cpue)=mfexp(log_q_rHD);
  //for (iyear=styr_rHD_cpue; iyear<=endyr_rHD_cpue; iyear++)
  //{   //index in number units
  //    N_rHD(iyear)=elem_prod(N_mdyr(iyear),sel_rHD(iyear)); 
  //    pred_cpue_rHD(iyear)=q_rHD(iyear)*q_rate_fcn_rHD(iyear)*q_DD_fcn(iyear)*sum(N_rHD(iyear));
  //    if (iyear<endyr_rHD_cpue){q_rHD(iyear+1)=q_rHD(iyear)*mfexp(q_RW_log_dev_rHD(iyear));}
  //}
  
FUNCTION get_length_comps
  
  for (iyear=1;iyear<=nyr_lenc_sBT;iyear++)
    {pred_lenc_sBT(iyear)=(N_sBT(yrs_lenc_sBT(iyear))*lenprob)/sum(N_sBT(yrs_lenc_sBT(iyear)));} 

  for (iyear=1;iyear<=nyr_lenc_cHL;iyear++)
    {pred_lenc_cHL(iyear)=(L_cHL_num(yrs_lenc_cHL(iyear))*lenprob)/sum(L_cHL_num(yrs_lenc_cHL(iyear)));}
   
  for (iyear=1;iyear<=nyr_lenc_cPT;iyear++)
    {pred_lenc_cPT(iyear)=(L_cPT_num(yrs_lenc_cPT(iyear))*lenprob)/sum(L_cPT_num(yrs_lenc_cPT(iyear)));}  
 
  for (iyear=1;iyear<=nyr_lenc_rHB;iyear++)
    {pred_lenc_rHB(iyear)=(L_rHB_num(yrs_lenc_rHB(iyear))*lenprob)/sum(L_rHB_num(yrs_lenc_rHB(iyear)));}
  
  for (iyear=1;iyear<=nyr_lenc_rHD;iyear++)
    {pred_lenc_rHD(iyear)=(D_rHD_num(yrs_lenc_rHD(iyear))*lenprob)/sum(D_rHD_num(yrs_lenc_rHD(iyear)));}

  for (iyear=1;iyear<=nyr_lenc_rGN;iyear++)
    {pred_lenc_rGN(iyear)=(L_rGN_num(yrs_lenc_rGN(iyear))*lenprob)/sum(L_rGN_num(yrs_lenc_rGN(iyear)));} 
 
  
FUNCTION get_age_comps

  for (iyear=1;iyear<=nyr_agec_sBT;iyear++)
  {
    predNoErr_agec_sBT(iyear)=N_sBT(yrs_agec_sBT(iyear))/sum(N_sBT(yrs_agec_sBT(iyear)));
    pred_agec_sBT(iyear)=age_error*predNoErr_agec_sBT(iyear);                      
  }

  for (iyear=1;iyear<=nyr_agec_sTV;iyear++)
  {
    predNoErr_agec_sTV(iyear)=N_sTV(yrs_agec_sTV(iyear))/sum(N_sTV(yrs_agec_sTV(iyear)));
    pred_agec_sTV(iyear)=age_error*predNoErr_agec_sTV(iyear);                      
  }
    
  for (iyear=1;iyear<=nyr_agec_cHL;iyear++)
  {
    predNoErr_agec_cHL(iyear)=L_cHL_num(yrs_agec_cHL(iyear))/sum(L_cHL_num(yrs_agec_cHL(iyear)));
    pred_agec_cHL(iyear)=age_error*predNoErr_agec_cHL(iyear);                      
  }

  for (iyear=1;iyear<=nyr_agec_cPT;iyear++)
  {
    predNoErr_agec_cPT(iyear)=L_cPT_num(yrs_agec_cPT(iyear))/sum(L_cPT_num(yrs_agec_cPT(iyear)));
    pred_agec_cPT(iyear)=age_error*predNoErr_agec_cPT(iyear);                      
  }
  
  for (iyear=1;iyear<=nyr_agec_rHB;iyear++)
  {
    predNoErr_agec_rHB(iyear)=L_rHB_num(yrs_agec_rHB(iyear))/sum(L_rHB_num(yrs_agec_rHB(iyear)));
    pred_agec_rHB(iyear)=age_error*predNoErr_agec_rHB(iyear);                    
  }
    
////--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
FUNCTION get_weighted_current 
  F_temp_sum=0.0;
  F_temp_sum+=mfexp((nyr_selpar_wgted*log_avg_F_L_cHL+sum(log_dev_F_L_cHL((endyr-nyr_selpar_wgted+1),endyr)))/nyr_selpar_wgted);
  F_temp_sum+=mfexp((nyr_selpar_wgted*log_avg_F_L_cPT+sum(log_dev_F_L_cPT((endyr-nyr_selpar_wgted+1),endyr)))/nyr_selpar_wgted);
  F_temp_sum+=mfexp((nyr_selpar_wgted*log_avg_F_L_rHB+sum(log_dev_F_L_rHB((endyr-nyr_selpar_wgted+1),endyr)))/nyr_selpar_wgted);
  F_temp_sum+=mfexp((nyr_selpar_wgted*log_avg_F_L_rGN+sum(log_dev_F_L_rGN((endyr-nyr_selpar_wgted+1),endyr)))/nyr_selpar_wgted);
  F_temp_sum+=mfexp((nyr_selpar_wgted*log_avg_F_L_rIN+sum(log_dev_F_L_rIN((endyr-nyr_selpar_wgted+1),endyr)))/nyr_selpar_wgted);
  F_temp_sum+=mfexp((nyr_selpar_wgted*log_avg_F_D_cGD+sum(log_dev_F_D_cGD((endyr-nyr_selpar_wgted+1),endyr)))/nyr_selpar_wgted); 
  F_temp_sum+=mfexp((nyr_selpar_wgted*log_avg_F_D_rHD+sum(log_dev_F_D_rHD((endyr-nyr_selpar_wgted+1),endyr)))/nyr_selpar_wgted);
  F_temp_sum+=mfexp((nyr_selpar_wgted*log_avg_F_D_rGD+sum(log_dev_F_D_rGD((endyr-nyr_selpar_wgted+1),endyr)))/nyr_selpar_wgted);
  F_temp_sum+=mfexp((nyr_selpar_wgted*log_avg_F_D_rID+sum(log_dev_F_D_rID((endyr-nyr_selpar_wgted+1),endyr)))/nyr_selpar_wgted);

  prop_F_cHL=mfexp((nyr_selpar_wgted*log_avg_F_L_cHL+sum(log_dev_F_L_cHL((endyr-nyr_selpar_wgted+1),endyr)))/nyr_selpar_wgted)/F_temp_sum;
  prop_F_cPT=mfexp((nyr_selpar_wgted*log_avg_F_L_cPT+sum(log_dev_F_L_cPT((endyr-nyr_selpar_wgted+1),endyr)))/nyr_selpar_wgted)/F_temp_sum;
  prop_F_rHB=mfexp((nyr_selpar_wgted*log_avg_F_L_rHB+sum(log_dev_F_L_rHB((endyr-nyr_selpar_wgted+1),endyr)))/nyr_selpar_wgted)/F_temp_sum;
  prop_F_rGN=mfexp((nyr_selpar_wgted*log_avg_F_L_rGN+sum(log_dev_F_L_rGN((endyr-nyr_selpar_wgted+1),endyr)))/nyr_selpar_wgted)/F_temp_sum;
  prop_F_rIN=mfexp((nyr_selpar_wgted*log_avg_F_L_rIN+sum(log_dev_F_L_rIN((endyr-nyr_selpar_wgted+1),endyr)))/nyr_selpar_wgted)/F_temp_sum;
  prop_F_cGD=mfexp((nyr_selpar_wgted*log_avg_F_D_cGD+sum(log_dev_F_D_cGD((endyr-nyr_selpar_wgted+1),endyr)))/nyr_selpar_wgted)/F_temp_sum;
  prop_F_rHD=mfexp((nyr_selpar_wgted*log_avg_F_D_rHD+sum(log_dev_F_D_rHD((endyr-nyr_selpar_wgted+1),endyr)))/nyr_selpar_wgted)/F_temp_sum;
  prop_F_rGD=mfexp((nyr_selpar_wgted*log_avg_F_D_rGD+sum(log_dev_F_D_rGD((endyr-nyr_selpar_wgted+1),endyr)))/nyr_selpar_wgted)/F_temp_sum;
  prop_F_rID=mfexp((nyr_selpar_wgted*log_avg_F_D_rID+sum(log_dev_F_D_rID((endyr-nyr_selpar_wgted+1),endyr)))/nyr_selpar_wgted)/F_temp_sum;


  log_dev_F_end_cHL=sum(log_dev_F_L_cHL((endyr-nyr_selpar_wgted+1),endyr))/nyr_selpar_wgted;
  log_dev_F_end_cPT=sum(log_dev_F_L_cPT((endyr-nyr_selpar_wgted+1),endyr))/nyr_selpar_wgted;  
  log_dev_F_end_rHB=sum(log_dev_F_L_rHB((endyr-nyr_selpar_wgted+1),endyr))/nyr_selpar_wgted;
  log_dev_F_end_rGN=sum(log_dev_F_L_rGN((endyr-nyr_selpar_wgted+1),endyr))/nyr_selpar_wgted;
  log_dev_F_end_rIN=sum(log_dev_F_L_rIN((endyr-nyr_selpar_wgted+1),endyr))/nyr_selpar_wgted;

  log_dev_F_end_cGD=sum(log_dev_F_D_cGD((endyr-nyr_selpar_wgted+1),endyr))/nyr_selpar_wgted;
  log_dev_F_end_rHD=sum(log_dev_F_D_rHD((endyr-nyr_selpar_wgted+1),endyr))/nyr_selpar_wgted;
  log_dev_F_end_rGD=sum(log_dev_F_D_rGD((endyr-nyr_selpar_wgted+1),endyr))/nyr_selpar_wgted;
  log_dev_F_end_rID=sum(log_dev_F_D_rID((endyr-nyr_selpar_wgted+1),endyr))/nyr_selpar_wgted;

  F_end_L=sel_cHL(endyr)*mfexp(log_avg_F_L_cHL+log_dev_F_end_cHL)+
          sel_cPT(endyr)*mfexp(log_avg_F_L_cPT+log_dev_F_end_cPT)+
          sel_rHB(endyr)*mfexp(log_avg_F_L_rHB+log_dev_F_end_rHB)+
          sel_rGN(endyr)*mfexp(log_avg_F_L_rGN+log_dev_F_end_rGN)+
          sel_rIN(endyr)*mfexp(log_avg_F_L_rIN+log_dev_F_end_rIN);
                 
  F_end_D=sel_cGD(endyr)*mfexp(log_avg_F_D_cGD+log_dev_F_end_cGD)+
          sel_rHD(endyr)*mfexp(log_avg_F_D_rHD+log_dev_F_end_rHD)+
          sel_rGD(endyr)*mfexp(log_avg_F_D_rGD+log_dev_F_end_rGD)+
          sel_rID(endyr)*mfexp(log_avg_F_D_rID+log_dev_F_end_rID);    

  F_end=F_end_L+F_end_D;
  F_end_apex=max(F_end);

  sel_wgted_tot=F_end/F_end_apex;
  sel_wgted_L=elem_prod(sel_wgted_tot, elem_div(F_end_L,F_end));
  sel_wgted_D=elem_prod(sel_wgted_tot, elem_div(F_end_D,F_end));
  
  wgt_wgted_L_denom=prop_F_cHL+prop_F_cPT+prop_F_rHB+prop_F_rGN;
  wgt_wgted_L_klb=prop_F_cHL/wgt_wgted_L_denom*wgt_cHL_klb(endyr)+
              prop_F_cPT/wgt_wgted_L_denom*wgt_cPT_klb(endyr)+
              prop_F_rHB/wgt_wgted_L_denom*wgt_rHB_klb(endyr)+
              prop_F_rGN/wgt_wgted_L_denom*wgt_rGN_klb(endyr)+
              prop_F_rIN/wgt_wgted_L_denom*wgt_rIN_klb(endyr);                

  wgt_wgted_D_denom=prop_F_cGD+prop_F_rHD+prop_F_rGD;
  wgt_wgted_D_klb=prop_F_cGD/wgt_wgted_D_denom*wgt_cGD_klb(endyr)+
              prop_F_rHD/wgt_wgted_D_denom*wgt_rHD_klb(endyr)+
              prop_F_rGD/wgt_wgted_D_denom*wgt_rGD_klb(endyr)+
              prop_F_rID/wgt_wgted_D_denom*wgt_rID_klb(endyr);
  
FUNCTION get_msy
  
  //compute values as functions of F
  for(ff=1; ff<=n_iter_msy; ff++)
  {
    //uses fishery-weighted F's
    Z_age_msy=0.0;
    F_L_age_msy=0.0;
    F_D_age_msy=0.0;
      
    F_L_age_msy=F_msy(ff)*sel_wgted_L;
    F_D_age_msy=F_msy(ff)*sel_wgted_D;
    Z_age_msy=M+F_L_age_msy+F_D_age_msy;         
    
    N_age_msy(1)=1.0;
    for (iage=2; iage<=nages; iage++)
    {
      N_age_msy(iage)=N_age_msy(iage-1)*mfexp(-1.*Z_age_msy(iage-1));
    }
    N_age_msy(nages)=N_age_msy(nages)/(1.0-mfexp(-1.*Z_age_msy(nages)));
    N_age_msy_spawn(1,(nages-1))=elem_prod(N_age_msy(1,(nages-1)),
                                   mfexp((-1.*Z_age_msy(1,(nages-1)))*spawn_time_frac));                 
    N_age_msy_spawn(nages)=(N_age_msy_spawn(nages-1)*
                          (mfexp(-1.*(Z_age_msy(nages-1)*(1.0-spawn_time_frac) + 
                                 Z_age_msy(nages)*spawn_time_frac) )))
                          /(1.0-mfexp(-1.*Z_age_msy(nages)));
                     
    spr_msy(ff)=sum(elem_prod(N_age_msy_spawn,reprod));
        

    //Compute equilibrium values of R (including bias correction), SSB and Yield at each F
    R_eq(ff)=SR_eq_func(R0, steep, spr_msy(1), spr_msy(ff), BiasCor, SR_switch);
    if (R_eq(ff)<dzero) {R_eq(ff)=dzero;}    
    N_age_msy*=R_eq(ff);
    N_age_msy_spawn*=R_eq(ff);
    
    for (iage=1; iage<=nages; iage++)
    {
      L_age_msy(iage)=N_age_msy(iage)*(F_L_age_msy(iage)/Z_age_msy(iage))*
                      (1.-mfexp(-1.*Z_age_msy(iage)));
      D_age_msy(iage)=N_age_msy(iage)*(F_D_age_msy(iage)/Z_age_msy(iage))*
                      (1.-mfexp(-1.0*Z_age_msy(iage)));
      T_age_msy(iage)=N_age_msy(iage)*((F_D_age_msy(iage)+F_L_age_msy(iage))/Z_age_msy(iage))*(1.-mfexp(-1.0*Z_age_msy(iage)));
    }
    
    
    SSB_eq(ff)=sum(elem_prod(N_age_msy_spawn,reprod));
    B_eq(ff)=sum(elem_prod(N_age_msy,wgt_klb));
    L_eq_klb(ff)=sum(elem_prod(L_age_msy,wgt_wgted_L_klb));
    L_eq_knum(ff)=sum(L_age_msy)/1000.0;  
    D_eq_klb(ff)=sum(elem_prod(D_age_msy,wgt_wgted_D_klb));    
    D_eq_knum(ff)=sum(D_age_msy)/1000.0;  
    T_eq_klb(ff)=sum(elem_prod(T_age_msy,wgt_wgted_L_klb));    
    T_eq_knum(ff)=sum(T_age_msy)/1000.0;
  }  //end ff loop
  
  msy_lb_klb_out=max(T_eq_klb);
  msy_n_knum_out=max(T_eq_knum);
  
  for(ff=1; ff<=n_iter_msy; ff++)
  {
   if(T_eq_klb(ff) == msy_lb_klb_out) 
      {    
        SSB_msy_lb_out=SSB_eq(ff);
        B_msy_lb_out=B_eq(ff);
        R_msy_lb_out=R_eq(ff);
        msy_lb_knum_out=L_eq_knum(ff);
        D_msy_lb_knum_out=D_eq_knum(ff);
        D_msy_lb_klb_out=D_eq_klb(ff);
        F_msy_lb_out=F_msy(ff);  
        spr_msy_lb_out=spr_msy(ff);  
        T_msy_lb_knum_out=D_eq_knum(ff);
        T_msy_lb_klb_out=D_eq_klb(ff);    
      }

   if(T_eq_knum(ff) == msy_n_knum_out) 
      {    
        SSB_msy_n_out=SSB_eq(ff);
        B_msy_n_out=B_eq(ff);
        R_msy_n_out=R_eq(ff);
        msy_n_klb_out=L_eq_klb(ff);
        D_msy_n_knum_out=D_eq_knum(ff);
        D_msy_n_klb_out=D_eq_klb(ff);
        F_msy_n_out=F_msy(ff);  
        spr_msy_n_out=spr_msy(ff);  
        T_msy_n_knum_out=D_eq_knum(ff);
        T_msy_n_klb_out=D_eq_klb(ff);     
      }
  }


//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
FUNCTION get_per_recruit_stuff

  //static per-recruit stuff
 
  for(iyear=styr; iyear<=endyr; iyear++)
  {
    N_age_spr(1)=1.0;
    for(iage=2; iage<=nages; iage++)
    {
      N_age_spr(iage)=N_age_spr(iage-1)*mfexp(-1.*Z(iyear,iage-1));
    }
    N_age_spr(nages)=N_age_spr(nages)/(1.0-mfexp(-1.*Z(iyear,nages)));    
    N_age_spr_spawn(1,(nages-1))=elem_prod(N_age_spr(1,(nages-1)),
                                mfexp(-1.*Z(iyear)(1,(nages-1))*spawn_time_frac));
    N_age_spr_spawn(nages)=(N_age_spr_spawn(nages-1)*
                          (mfexp(-1.*(Z(iyear)(nages-1)*(1.0-spawn_time_frac) + Z(iyear)(nages)*spawn_time_frac) )))
                          /(1.0-mfexp(-1.*Z(iyear)(nages)));           
    spr_static(iyear)=sum(elem_prod(N_age_spr_spawn,reprod))/spr_F0;
  }
  

  //compute SSB/R and YPR as functions of F
  for(ff=1; ff<=n_iter_spr; ff++)
  {
    //uses fishery-weighted F's, same as in MSY calculations
    Z_age_spr=0.0;
    F_L_age_spr=0.0;

    F_L_age_spr=F_spr(ff)*sel_wgted_L;
    F_T_age_spr=F_L_age_spr+F_spr(ff)*sel_wgted_D;
    
    Z_age_spr=M+F_L_age_spr+F_spr(ff)*sel_wgted_D;

    N_age_spr(1)=1.0;
    for (iage=2; iage<=nages; iage++)
    {
      N_age_spr(iage)=N_age_spr(iage-1)*mfexp(-1.*Z_age_spr(iage-1));
    }
    N_age_spr(nages)=N_age_spr(nages)/(1-mfexp(-1.*Z_age_spr(nages)));
    N_age_spr_spawn(1,(nages-1))=elem_prod(N_age_spr(1,(nages-1)),
                                   mfexp((-1.*Z_age_spr(1,(nages-1)))*spawn_time_frac));                 
    N_age_spr_spawn(nages)=(N_age_spr_spawn(nages-1)*
                          (mfexp(-1.*(Z_age_spr(nages-1)*(1.0-spawn_time_frac) + Z_age_spr(nages)*spawn_time_frac) )))
                          /(1.0-mfexp(-1.*Z_age_spr(nages)));
    
    spr_spr(ff)=sum(elem_prod(N_age_spr_spawn,reprod));
    L_spr(ff)=0.0;
    T_spr(ff)=0.0;
    for (iage=1; iage<=nages; iage++)
    {
      L_age_spr(iage)=N_age_spr(iage)*(F_L_age_spr(iage)/Z_age_spr(iage))*
                      (1.-mfexp(-1.*Z_age_spr(iage)));
      L_spr(ff)+=L_age_spr(iage)*wgt_wgted_L_klb(iage)*1000.0; //in lb
      // SPR for total catch
      T_age_spr(iage)=N_age_spr(iage)*(F_T_age_spr(iage)/Z_age_spr(iage))*
                      (1.-mfexp(-1.*Z_age_spr(iage)));
      T_spr(ff)+=T_age_spr(iage)*wgt_wgted_L_klb(iage)*1000.0; //in lb
    }   
  }
  spr_ratio=spr_spr/spr_F0;
  F30_dum=min(fabs(spr_ratio-0.3));
  F40_dum=min(fabs(spr_ratio-0.4));
  F50_dum=min(fabs(spr_ratio-0.5));
  Fproxy_dum=min(fabs(spr_ratio-spr_proxy));
  for(ff=1; ff<=n_iter_spr; ff++)
  {   
      if (fabs(spr_ratio(ff)-0.3)==F30_dum) {F30_out=F_spr(ff);}
      if (fabs(spr_ratio(ff)-0.4)==F40_dum) {F40_out=F_spr(ff);}
      if (fabs(spr_ratio(ff)-0.5)==F50_dum) {F50_out=F_spr(ff);}
      if (fabs(spr_ratio(ff)-spr_proxy)==Fproxy_dum) {Fproxy_out=F_spr(ff);}
  }

    
  rec=column(N,1);	  
  switch(SPR_rec_switch){
       case 1: //Observed mean over years specified
          rec_mean=sum(rec(styr_rec_spr, endyr_rec_spr))/nyr_rec_spr;
          break;
       case 2: //Expected recruitment is R0 with bias correction (if that feature is turned on)
          rec_mean=BiasCor*R0;
          break;
       default: // no such switch available
          cout << "Error in input: SPR_rec_switch must be set to 1 or 2." << endl;
          cout << "Presently it is set to " << SPR_rec_switch <<"."<< endl;
          exit(0);          
	}	  
    R_Fproxy_out=rec_mean; 	  
    R_F50_out=rec_mean; 	  
    R_F40_out=rec_mean; 	  
    R_F30_out=rec_mean; 

// proxy
  F_L_age_spr=Fproxy_out*sel_wgted_L;
  F_D_age_spr=Fproxy_out*sel_wgted_D;
  Z_age_spr=M+F_L_age_spr+F_D_age_spr;

  N_age_spr(1)=R_Fproxy_out;
  for (iage=2; iage<=nages; iage++)
     {N_age_spr(iage)=N_age_spr(iage-1)*mfexp(-1.*Z_age_spr(iage-1));}
  N_age_spr(nages)=N_age_spr(nages)/(1-mfexp(-1.*Z_age_spr(nages)));
  N_age_spr_spawn(1,(nages-1))=elem_prod(N_age_spr(1,(nages-1)),
                                   mfexp((-1.*Z_age_spr(1,(nages-1)))*spawn_time_frac));                 
  N_age_spr_spawn(nages)=(N_age_spr_spawn(nages-1)*
                          (mfexp(-1.*(Z_age_spr(nages-1)*(1.0-spawn_time_frac) + Z_age_spr(nages)*spawn_time_frac) )))
                          /(1.0-mfexp(-1.*Z_age_spr(nages)));

  for (iage=1; iage<=nages; iage++)
    {
      L_age_Fproxy(iage)=N_age_spr(iage)*(F_L_age_spr(iage)/Z_age_spr(iage))*
                      (1.-mfexp(-1.*Z_age_spr(iage)));
      D_age_Fproxy(iage)=N_age_spr(iage)*(F_D_age_spr(iage)/Z_age_spr(iage))*
                      (1.-mfexp(-1.0*Z_age_spr(iage)));                      
    }
		
  SSB_Fproxy_out=sum(elem_prod(N_age_spr_spawn,reprod));
  B_Fproxy_out=sum(elem_prod(N_age_spr,wgt_klb));
  L_Fproxy_klb_out=sum(elem_prod(L_age_Fproxy,wgt_wgted_L_klb)); //in whole weight
  L_Fproxy_knum_out=sum(L_age_Fproxy)/1000.0;
  D_Fproxy_klb_out=sum(elem_prod(D_age_Fproxy,wgt_wgted_D_klb)); //in whole weight   
  D_Fproxy_knum_out=sum(D_age_Fproxy)/1000.0;
  
// 50
  F_L_age_spr=F50_out*sel_wgted_L;
  F_D_age_spr=F50_out*sel_wgted_D;
  Z_age_spr=M+F_L_age_spr+F_D_age_spr;

  N_age_spr(1)=R_F50_out;
  for (iage=2; iage<=nages; iage++)
     {N_age_spr(iage)=N_age_spr(iage-1)*mfexp(-1.*Z_age_spr(iage-1));}
  N_age_spr(nages)=N_age_spr(nages)/(1-mfexp(-1.*Z_age_spr(nages)));
  N_age_spr_spawn(1,(nages-1))=elem_prod(N_age_spr(1,(nages-1)),
                                   mfexp((-1.*Z_age_spr(1,(nages-1)))*spawn_time_frac));                 
  N_age_spr_spawn(nages)=(N_age_spr_spawn(nages-1)*
                          (mfexp(-1.*(Z_age_spr(nages-1)*(1.0-spawn_time_frac) + Z_age_spr(nages)*spawn_time_frac) )))
                          /(1.0-mfexp(-1.*Z_age_spr(nages)));

  for (iage=1; iage<=nages; iage++)
    {
      L_age_Fproxy(iage)=N_age_spr(iage)*(F_L_age_spr(iage)/Z_age_spr(iage))*
                      (1.-mfexp(-1.*Z_age_spr(iage)));
      D_age_Fproxy(iage)=N_age_spr(iage)*(F_D_age_spr(iage)/Z_age_spr(iage))*
                      (1.-mfexp(-1.0*Z_age_spr(iage)));                      
    }
		
  SSB_F50_out=sum(elem_prod(N_age_spr_spawn,reprod));
  B_F50_out=sum(elem_prod(N_age_spr,wgt_klb));
  L_F50_klb_out=sum(elem_prod(L_age_Fproxy,wgt_wgted_L_klb)); //in whole weight
  L_F50_knum_out=sum(L_age_Fproxy)/1000.0;
  D_F50_klb_out=sum(elem_prod(D_age_Fproxy,wgt_wgted_D_klb)); //in whole weight   
  D_F50_knum_out=sum(D_age_Fproxy)/1000.0;

// 40
  F_L_age_spr=F40_out*sel_wgted_L;
  F_D_age_spr=F40_out*sel_wgted_D;
  Z_age_spr=M+F_L_age_spr+F_D_age_spr;
  N_age_spr(1)=R_F40_out;
  for (iage=2; iage<=nages; iage++)
     {N_age_spr(iage)=N_age_spr(iage-1)*mfexp(-1.*Z_age_spr(iage-1));}
  N_age_spr(nages)=N_age_spr(nages)/(1-mfexp(-1.*Z_age_spr(nages)));
  N_age_spr_spawn(1,(nages-1))=elem_prod(N_age_spr(1,(nages-1)),
                                   mfexp((-1.*Z_age_spr(1,(nages-1)))*spawn_time_frac));
  N_age_spr_spawn(nages)=(N_age_spr_spawn(nages-1)*
                          (mfexp(-1.*(Z_age_spr(nages-1)*(1.0-spawn_time_frac) + Z_age_spr(nages)*spawn_time_frac) )))
                          /(1.0-mfexp(-1.*Z_age_spr(nages)));
  for (iage=1; iage<=nages; iage++)
    {
      L_age_Fproxy(iage)=N_age_spr(iage)*(F_L_age_spr(iage)/Z_age_spr(iage))*
                      (1.-mfexp(-1.*Z_age_spr(iage)));
      D_age_Fproxy(iage)=N_age_spr(iage)*(F_D_age_spr(iage)/Z_age_spr(iage))*
                      (1.-mfexp(-1.0*Z_age_spr(iage)));                      
    }
  SSB_F40_out=sum(elem_prod(N_age_spr_spawn,reprod));
  B_F40_out=sum(elem_prod(N_age_spr,wgt_klb));
  L_F40_klb_out=sum(elem_prod(L_age_Fproxy,wgt_wgted_L_klb)); //in whole weight
  L_F40_knum_out=sum(L_age_Fproxy)/1000.0;
  D_F40_klb_out=sum(elem_prod(D_age_Fproxy,wgt_wgted_D_klb)); //in whole weight   
  D_F40_knum_out=sum(D_age_Fproxy)/1000.0;

// 30
  F_L_age_spr=F30_out*sel_wgted_L;
  F_D_age_spr=F30_out*sel_wgted_D;
  Z_age_spr=M+F_L_age_spr+F_D_age_spr;
  N_age_spr(1)=R_F30_out;
  for (iage=2; iage<=nages; iage++)
     {N_age_spr(iage)=N_age_spr(iage-1)*mfexp(-1.*Z_age_spr(iage-1));}
  N_age_spr(nages)=N_age_spr(nages)/(1-mfexp(-1.*Z_age_spr(nages)));
  N_age_spr_spawn(1,(nages-1))=elem_prod(N_age_spr(1,(nages-1)),
                                   mfexp((-1.*Z_age_spr(1,(nages-1)))*spawn_time_frac));
  N_age_spr_spawn(nages)=(N_age_spr_spawn(nages-1)*
                          (mfexp(-1.*(Z_age_spr(nages-1)*(1.0-spawn_time_frac) + Z_age_spr(nages)*spawn_time_frac) )))
                          /(1.0-mfexp(-1.*Z_age_spr(nages)));
  for (iage=1; iage<=nages; iage++)
    {
      L_age_Fproxy(iage)=N_age_spr(iage)*(F_L_age_spr(iage)/Z_age_spr(iage))*
                      (1.-mfexp(-1.*Z_age_spr(iage)));
      D_age_Fproxy(iage)=N_age_spr(iage)*(F_D_age_spr(iage)/Z_age_spr(iage))*
                      (1.-mfexp(-1.0*Z_age_spr(iage)));                      
    }
  SSB_F30_out=sum(elem_prod(N_age_spr_spawn,reprod));
  B_F30_out=sum(elem_prod(N_age_spr,wgt_klb));
  L_F30_klb_out=sum(elem_prod(L_age_Fproxy,wgt_wgted_L_klb)); //in whole weight
  L_F30_knum_out=sum(L_age_Fproxy)/1000.0;
  D_F30_klb_out=sum(elem_prod(D_age_Fproxy,wgt_wgted_D_klb)); //in whole weight   
  D_F30_knum_out=sum(D_age_Fproxy)/1000.0;   

  //--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
FUNCTION get_miscellaneous_stuff

//switch here if var_rec_dev <=dzero
  if(var_rec_dev>0.0)
  {
    sigma_rec_dev=sqrt(var_rec_dev); //sample SD of predicted residuals (may not equal rec_sigma)
  } else{sigma_rec_dev=0.0;}
  
  //compute total landings- and discards-at-age in 1000 fish and klb
  L_total_num.initialize();
  L_total_klb.initialize();
  D_total_num.initialize();
  D_total_klb.initialize();
  L_total_knum_yr.initialize();
  L_total_klb_yr.initialize();  
  D_total_knum_yr.initialize();
  D_total_klb_yr.initialize();
  
  for(iyear=styr; iyear<=endyr; iyear++)
  {
        L_total_klb_yr(iyear)=
          pred_L_cHL_klb(iyear)+
          pred_L_cPT_klb(iyear)+
          pred_L_cTW_klb(iyear)+
          pred_L_rHB_klb(iyear)+
          pred_L_rGN_klb(iyear)+
          pred_L_rIN_klb(iyear);
    
        L_total_knum_yr(iyear)=
          pred_L_cHL_knum(iyear)+
          pred_L_cPT_knum(iyear)+
          pred_L_cTW_knum(iyear)+
          pred_L_rHB_knum(iyear)+
          pred_L_rGN_knum(iyear)+
          pred_L_rIN_knum(iyear);
                
        D_total_knum_yr(iyear)+=
          pred_D_cGD_knum(iyear)+
          pred_D_rHD_knum(iyear)+
          pred_D_rGD_knum(iyear)+
          pred_D_rID_knum(iyear);
        
        D_total_klb_yr(iyear)+=
          pred_D_cGD_klb(iyear)+
          pred_D_rHD_klb(iyear)+
          pred_D_rGD_klb(iyear)+
          pred_D_rID_klb(iyear);

        D_cGD_klb(iyear)=elem_prod(D_cGD_num(iyear),wgt_cGD_klb(iyear));    //in 1000 lb
        D_rHD_klb(iyear)=elem_prod(D_rHD_num(iyear),wgt_rHD_klb(iyear));    //in 1000 lb
        D_rGD_klb(iyear)=elem_prod(D_rGD_num(iyear),wgt_rGD_klb(iyear));    //in 1000 lb 
        D_rID_klb(iyear)=elem_prod(D_rID_num(iyear),wgt_rID_klb(iyear));    //in 1000 lb 
        
        B(iyear)=elem_prod(N(iyear),wgt_klb);
        totN(iyear)=sum(N(iyear));
        totB(iyear)=sum(B(iyear));              
  }
  
  //landings at age in number fish
  L_total_num=
    L_cHL_num+
    L_cPT_num+
    L_cTW_num+
    L_rHB_num+
    L_rGN_num+
    L_rIN_num; 
  
  //landings at age in klb whole weight
  L_total_klb=
    L_cHL_klb+
    L_cPT_klb+
    L_cTW_klb+
    L_rHB_klb+
    L_rGN_klb+
    L_rIN_klb;

  //discards at age in number fish
  D_total_num=
    D_cGD_num+
    D_rHD_num+
    D_rGD_num+
    D_rID_num;
  
  //discards at age in klb whole weight
  D_total_klb=
    D_cGD_klb+
    D_rHD_klb+
    D_rGD_klb+            
    D_rID_klb;            
  
  Fend_mean_temp=1.0;
  for (iyear=1; iyear<=nyr_selpar_wgted; iyear++) {Fend_mean_temp*=Fapex(endyr-iyear+1);}
  Fend_mean=pow(Fend_mean_temp,(1.0/nyr_selpar_wgted));
  if(F_msy_lb_out>0)
    {
      FdF_msy_lb=Fapex/F_msy_lb_out;
      FdF_msy_lb_end=FdF_msy_lb(endyr);
      FdF_msy_lb_end_mean=Fend_mean/F_msy_lb_out;
    }
  if(F_msy_n_out>0)
    {
      FdF_msy_n=Fapex/F_msy_n_out;
      FdF_msy_n_end=FdF_msy_n(endyr);
      FdF_msy_n_end_mean=Fend_mean/F_msy_n_out;
    }
  if(SSB_msy_lb_out>0)
    {
      SdSSB_msy_lb=SSB/SSB_msy_lb_out;
      SdSSB_msy_lb_end=SdSSB_msy_lb(endyr);
    }  
  if(SSB_msy_n_out>0)
    {
      SdSSB_msy_n=SSB/SSB_msy_n_out;
      SdSSB_msy_n_end=SdSSB_msy_n(endyr);
    }  
  if(Fproxy_out>0)
    {
      FdFproxy=Fapex/Fproxy_out;
      FdFproxy_end_mean=Fend_mean/Fproxy_out;
    }
  if(SSB_Fproxy_out>0)
    {
      SdSSB_Fproxy=SSB/SSB_Fproxy_out;
      Sdmsst_Fproxy=SSB/((1-M_constant)*SSB_Fproxy_out);
      SdSSB_Fproxy_end=SdSSB_Fproxy(endyr);
      Sdmsst_Fproxy_end=Sdmsst_Fproxy(endyr);
    }
      if(F50_out>0)
    {
      FdF50=Fapex/F50_out;
      FdF50_end_mean=Fend_mean/F50_out;
    }
  if(SSB_F50_out>0)
    {
      SdSSB_F50=SSB/SSB_F50_out;
      Sdmsst_F50=SSB/((1-M_constant)*SSB_F50_out);
      SdSSB_F50_end=SdSSB_F50(endyr);
      Sdmsst_F50_end=Sdmsst_F50(endyr);
    }
      if(F40_out>0)
    {
      FdF40=Fapex/F40_out;
      FdF40_end_mean=Fend_mean/F40_out;
    }
  if(SSB_F40_out>0)
    {
      SdSSB_F40=SSB/SSB_F40_out;
      Sdmsst_F40=SSB/((1-M_constant)*SSB_F40_out);
      SdSSB_F40_end=SdSSB_F40(endyr);
      Sdmsst_F40_end=Sdmsst_F40(endyr);
    }
  if(F30_out>0)
    {
      FdF30=Fapex/F30_out;
      FdF30_end_mean=Fend_mean/F30_out;
    }
  if(SSB_F30_out>0)
    {
      SdSSB_F30=SSB/SSB_F30_out;
      Sdmsst_F30=SSB/((1-M_constant)*SSB_F30_out);
      SdSSB_F30_end=SdSSB_F30(endyr);
      Sdmsst_F30_end=Sdmsst_F30(endyr);
    }
    
   //fill in log recruitment deviations for yrs they are nonzero
   for(iyear=styr_rec_dev; iyear<=endyr_rec_dev; iyear++)
     {log_rec_dev_output(iyear)=log_dev_rec(iyear);}
   //fill in log Nage deviations for ages they are nonzero (ages2+)
   for(iage=2; iage<=nages; iage++)
   { log_Nage_dev_output(iage)=log_dev_Nage(iage);}


FUNCTION get_projection

    switch(Fproj_switch){
       case 1: //F=Fcurrent
          F_reg_proj=Fend_mean;
          break;
       case 2: //F=Fmsy
          F_reg_proj=F_msy_lb_out;
          break;
       case 3: //F=F30
          F_reg_proj=F30_out;
          break;
       case 4: //F=F40
          F_reg_proj=F40_out;
          break;
       case 5: //F=F50
          F_reg_proj=F50_out;
          break;
       case 6: //F=Fproxy
          F_reg_proj=Fproxy_out;
          break;  	  
       default: // no such switch available
          cout << "Error in input: Projection switch Fproj_switch must be set to 1, 2, 3, 4, 5, or 6." << endl;
          cout << "Presently it is set to " << Fproj_switch <<"."<< endl;
          exit(0);
   }

  //N_proj initial conditions
      switch(Rproj_switch){
       case 1: //SR curve
          N_proj(styr_proj,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB(endyr),SR_switch);
          break;
       case 2: //mean from SPR calcs
          N_proj(styr_proj,1)=rec_mean;
          break;
       case 3: //mean from end of assessment years
          N_proj(styr_proj,1)=rec_mu;
          break;
       default: // no such switch available
          cout << "Error in input: Projection switch Rproj_switch must be set to 1, 2, or 3." << endl;
          cout << "Presently it is set to " << Rproj_switch <<"."<< endl;
          exit(0);
   }
   
  N_proj(styr_proj)(2,nages)=++elem_prod(N(endyr)(1,nages-1),(mfexp(-1.*Z(endyr)(1,nages-1))));
  N_proj(styr_proj,nages)+=N(endyr,nages)*mfexp(-1.*Z(endyr,nages)); //plus group

 for (iyear=styr_proj; iyear<=endyr_proj; iyear++)
 {     
       if (iyear<styr_regs) {F_proj(iyear)=Fend_mean;}
		else {F_proj(iyear)=Fproj_mult*F_reg_proj;}
		
		FL_age_proj=sel_wgted_L*F_proj(iyear);
		FD_age_proj=sel_wgted_D*F_proj(iyear);
		
       Z_proj(iyear)=M+FL_age_proj+FD_age_proj;
       N_spawn_proj(iyear)(1,nages)=elem_prod(N_proj(iyear)(1,nages),(mfexp(-1.*(Z_proj(iyear)(1,nages))*spawn_time_frac))); //peak spawning time
		SSB_proj(iyear)= sum(elem_prod(N_spawn_proj(iyear),reprod));
       B_proj(iyear)=sum(elem_prod(N_proj(iyear),wgt_klb)); //uses spawning weight
	     
		for (iage=1; iage<=nages; iage++)
			{L_age_proj(iyear,iage)=N_proj(iyear,iage)*FL_age_proj(iage)*(1.-mfexp(-1.*Z_proj(iyear,iage)))/Z_proj(iyear,iage);
		     D_age_proj(iyear,iage)=N_proj(iyear,iage)*FD_age_proj(iage)*(1.-mfexp(-1.*Z_proj(iyear,iage)))/Z_proj(iyear,iage);
			}          
       L_knum_proj(iyear)=sum(L_age_proj(iyear))/1000.0;
	    D_knum_proj(iyear)=sum(D_age_proj(iyear))/1000.0;
	    L_klb_proj(iyear)=sum(elem_prod(L_age_proj(iyear),wgt_wgted_L_klb));     //in 1000 lb
       D_klb_proj(iyear)=sum(elem_prod(D_age_proj(iyear),wgt_wgted_D_klb));     //in 1000 lb
		
		if (iyear<endyr_proj) {
			N_proj(iyear+1,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB_proj(iyear),SR_switch);
			N_proj(iyear+1)(2,nages)=++elem_prod(N_proj(iyear)(1,nages-1),(mfexp(-1.*Z_proj(iyear)(1,nages-1))));
			N_proj(iyear+1,nages)+=N_proj(iyear,nages)*mfexp(-1.*Z_proj(iyear,nages)); //plus group		
		}
 }
  R_proj=column(N_proj,1);                          
// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  



FUNCTION evaluate_objective_function
  fval=0.0;
  fval_data=0.0;

//---likelihoods---------------------------

//---Indices-------------------------------

  f_cpue_sBT=0.0;
  f_cpue_sBT=lk_lognormal(pred_cpue_sBT, obs_cpue_sBT, obs_cv_cpue_sBT, w_cpue_sBT);
  fval+=f_cpue_sBT;
  fval_data+=f_cpue_sBT;

  f_cpue_sTV=0.0;
  f_cpue_sTV=lk_lognormal(pow(pred_cpue_sTV,posobs_cpue_sTV), obs_cpue_corrected_sTV, obs_cv_cpue_corrected_sTV, w_cpue_sTV);
  fval+=f_cpue_sTV;
  fval_data+=f_cpue_sTV;
  
  f_cpue_cHL=0.0;
  f_cpue_cHL=lk_lognormal(pred_cpue_cHL, obs_cpue_cHL, obs_cv_cpue_cHL, w_cpue_cHL);
  fval+=f_cpue_cHL;
  fval_data+=f_cpue_cHL;
  
  f_cpue_rHB=0.0;
  f_cpue_rHB=lk_lognormal(pred_cpue_rHB, obs_cpue_rHB, obs_cv_cpue_rHB, w_cpue_rHB);
  fval+=f_cpue_rHB;
  fval_data+=f_cpue_rHB;
  
  //f_cpue_rHD=0.0;
  //f_cpue_rHD=lk_lognormal(pred_cpue_rHD, obs_rHD_cpue, rHD_cpue_cv, w_I_rHD);
  //fval+=f_cpue_rHD;
  //fval_data+=f_cpue_rHD;

////---Landings-------------------------------
 
  //f_L_cHL in 1000 lb ww
  f_L_cHL=lk_lognormal(pred_L_cHL_klb(styr_L_cHL,endyr_L_cHL), obs_L_cHL(styr_L_cHL,endyr_L_cHL), 
                      obs_cv_L_cHL(styr_L_cHL,endyr_L_cHL), w_L);
  fval+=f_L_cHL;
  fval_data+=f_L_cHL;
  
  //f_L_cPT in 1000 lb ww
  f_L_cPT=lk_lognormal(pred_L_cPT_klb(styr_L_cPT,endyr_L_cPT), obs_L_cPT(styr_L_cPT,endyr_L_cPT), 
                      obs_cv_L_cPT(styr_L_cPT,endyr_L_cPT), w_L);
  fval+=f_L_cPT;
  fval_data+=f_L_cPT;
  
  //f_L_cTW in 1000 lb ww
  f_L_cTW=lk_lognormal(pred_L_cTW_klb(styr_L_cTW,endyr_L_cTW), obs_L_cTW(styr_L_cTW,endyr_L_cTW), 
                      obs_cv_L_cTW(styr_L_cTW,endyr_L_cTW), w_L);
  fval+=f_L_cTW;
  fval_data+=f_L_cTW;
  
  //f_L_rHB in 1000 fish
  f_L_rHB=lk_lognormal(pred_L_rHB_knum(styr_L_rHB,endyr_L_rHB), obs_L_rHB(styr_L_rHB,endyr_L_rHB),
                      obs_cv_L_rHB(styr_L_rHB,endyr_L_rHB), w_L);
  fval+=f_L_rHB;
  fval_data+=f_L_rHB;

  //f_L_rGN in 1000 fish
  f_L_rGN=lk_lognormal(pred_L_rGN_knum(styr_L_rGN,endyr_L_rGN), obs_L_rGN(styr_L_rGN,endyr_L_rGN), 
                       obs_cv_L_rGN(styr_L_rGN,endyr_L_rGN), w_L);
  fval+=f_L_rGN;
  fval_data+=f_L_rGN;
  
  //f_L_rIN in 1000 fish
  f_L_rIN=lk_lognormal(pred_L_rIN_knum(styr_L_rIN,endyr_L_rIN), obs_L_rIN(styr_L_rIN,endyr_L_rIN), 
                       obs_cv_L_rIN(styr_L_rIN,endyr_L_rIN), w_L);
  fval+=f_L_rIN;
  fval_data+=f_L_rIN;  


//---Discards-------------------------------

  //f_D_cGD in 1000 fish
  f_D_cGD=lk_lognormal(pred_D_cGD_knum(styr_D_cGD,endyr_D_cGD), obs_D_cGD(styr_D_cGD,endyr_D_cGD), 
                       obs_cv_D_cGD(styr_D_cGD,endyr_D_cGD), w_D);
  fval+=f_D_cGD;
  fval_data+=f_D_cGD;
  
  
  //f_D_rHD in 1000 fish
  f_D_rHD=lk_lognormal(pred_D_rHD_knum(styr_D_rHD,endyr_D_rHD), obs_D_rHD(styr_D_rHD,endyr_D_rHD), 
                       obs_cv_D_rHD(styr_D_rHD,endyr_D_rHD), w_D);
  fval+=f_D_rHD;
  fval_data+=f_D_rHD;

  //f_D_rGD in 1000 fish
  f_D_rGD=lk_lognormal(pred_D_rGD_knum(styr_D_rGD,endyr_D_rGD), obs_D_rGD(styr_D_rGD,endyr_D_rGD), 
                       obs_cv_D_rGD(styr_D_rGD,endyr_D_rGD), w_D);
  fval+=f_D_rGD;
  fval_data+=f_D_rGD;
  
  //f_D_rID in 1000 fish
  f_D_rID=lk_lognormal(pred_D_rID_knum(styr_D_rID,endyr_D_rID), obs_D_rID(styr_D_rID,endyr_D_rID), 
                       obs_cv_D_rID(styr_D_rID,endyr_D_rID), w_D);
  fval+=f_D_rID;
  fval_data+=f_D_rID;
  


//---Length comps-------------------------------

  //f_lenc_sBT
  f_lenc_sBT=lk_dirichlet_multinomial(nsamp_lenc_sBT, pred_lenc_sBT, obs_lenc_sBT, nyr_lenc_sBT, double(nlenbins), minSS_lenc, log_dm_lenc_sBT); //w_lenc_sBT);
  fval+=f_lenc_sBT;
  fval_data+=f_lenc_sBT;

  //f_lenc_cHL
  f_lenc_cHL=lk_dirichlet_multinomial(nsamp_lenc_cHL, pred_lenc_cHL, obs_lenc_cHL, nyr_lenc_cHL, double(nlenbins), minSS_lenc, log_dm_lenc_cHL); //w_lenc_cHL);
  fval+=f_lenc_cHL;
  fval_data+=f_lenc_cHL;
  
  //f_lenc_cPT
  f_lenc_cPT=lk_dirichlet_multinomial(nsamp_lenc_cPT, pred_lenc_cPT, obs_lenc_cPT, nyr_lenc_cPT, double(nlenbins), minSS_lenc, log_dm_lenc_cPT); //w_lenc_cPT);
  fval+=f_lenc_cPT;
  fval_data+=f_lenc_cPT;
  
  //f_lenc_rHB
  f_lenc_rHB=lk_dirichlet_multinomial(nsamp_lenc_rHB, pred_lenc_rHB, obs_lenc_rHB, nyr_lenc_rHB, double(nlenbins), minSS_lenc, log_dm_lenc_rHB); //w_lenc_rHB);
  fval+=f_lenc_rHB;
  fval_data+=f_lenc_rHB;  
  
  //f_lenc_rHD
  f_lenc_rHD=lk_dirichlet_multinomial(nsamp_lenc_rHD, pred_lenc_rHD, obs_lenc_rHD, nyr_lenc_rHD, double(nlenbins), minSS_lenc, log_dm_lenc_rHD); //w_lenc_rHD);
  fval+=f_lenc_rHD;
  fval_data+=f_lenc_rHD;
    
  //f_lenc_rGN
  f_lenc_rGN=lk_dirichlet_multinomial(nsamp_lenc_rGN, pred_lenc_rGN, obs_lenc_rGN, nyr_lenc_rGN, double(nlenbins), minSS_lenc, log_dm_lenc_rGN); //w_lenc_rGN);
  fval+=f_lenc_rGN;
  fval_data+=f_lenc_rGN;


//---Age comps-------------------------------

  //f_agec_sBT
  f_agec_sBT=lk_dirichlet_multinomial(nsamp_agec_sBT, pred_agec_sBT, obs_agec_sBT, nyr_agec_sBT, double(nages), minSS_agec, log_dm_agec_sBT); //w_agec_sBT);
  fval+=f_agec_sBT;
  fval_data+=f_agec_sBT;

  //f_agec_sTV
  f_agec_sTV=lk_dirichlet_multinomial(nsamp_agec_sTV, pred_agec_sTV, obs_agec_sTV, nyr_agec_sTV, double(nages), minSS_agec, log_dm_agec_sTV); //w_agec_sTV);
  fval+=f_agec_sTV;
  fval_data+=f_agec_sTV;

  //f_agec_cHL
  f_agec_cHL=lk_dirichlet_multinomial(nsamp_agec_cHL, pred_agec_cHL, obs_agec_cHL, nyr_agec_cHL, double(nages), minSS_agec, log_dm_agec_cHL); //w_agec_cHL);
  fval+=f_agec_cHL;
  fval_data+=f_agec_cHL;
  
  //f_agec_cPT
  f_agec_cPT=lk_dirichlet_multinomial(nsamp_agec_cPT, pred_agec_cPT, obs_agec_cPT, nyr_agec_cPT, double(nages), minSS_agec, log_dm_agec_cPT); //w_agec_cPT);
  fval+=f_agec_cPT;
  fval_data+=f_agec_cPT;
    
  //f_agec_rHB
  f_agec_rHB=lk_dirichlet_multinomial(nsamp_agec_rHB, pred_agec_rHB, obs_agec_rHB, nyr_agec_rHB, double(nages), minSS_agec, log_dm_agec_rHB); //w_agec_rHB);
  fval+=f_agec_rHB;
  fval_data+=f_agec_rHB;  
  

//-----------Constraints and penalties--------------------------------

  f_rec_dev=0.0;
  rec_logL_add=nyr_rec_dev*log(rec_sigma);
  f_rec_dev=(square(log_dev_rec(styr_rec_dev) + rec_sigma_sqd2)/(2.0*rec_sigma_sq));
  for(iyear=(styr_rec_dev+1); iyear<=endyr_rec_dev; iyear++)
  {f_rec_dev+=(square(log_dev_rec(iyear)-R_autocorr*log_dev_rec(iyear-1) + rec_sigma_sqd2)/
               (2.0*rec_sigma_sq));}
  f_rec_dev+=rec_logL_add;            
  fval+=w_rec*f_rec_dev;

  f_rec_dev_early=0.0; //possible extra constraint on early rec deviations
  if (w_rec_early>0.0)
    { if (styr_rec_dev<endyr_rec_phase1)
        {  
          f_rec_dev_early=(square(log_dev_rec(styr_rec_dev) + rec_sigma_sq/2.0)/(2.0*rec_sigma_sq)) + rec_logL_add;
          for(iyear=(styr_rec_dev+1); iyear<=endyr_rec_phase1; iyear++)
          {f_rec_dev_early+=(square(log_dev_rec(iyear)-R_autocorr*log_dev_rec(iyear-1) + rec_sigma_sqd2)/
                            (2.0*rec_sigma_sq)) + rec_logL_add;}
        }
  fval+=w_rec_early*f_rec_dev_early;
  }
  
  f_rec_dev_end=0.0; //possible extra constraint on ending rec deviations
  if (w_rec_end>0.0)
  { if (endyr_rec_phase2<endyr_rec_dev)
        {  
          for(iyear=(endyr_rec_phase2+1); iyear<=endyr_rec_dev; iyear++)
          {f_rec_dev_end+=(square(log_dev_rec(iyear)-R_autocorr*log_dev_rec(iyear-1) + rec_sigma_sqd2)/
                            (2.0*rec_sigma_sq)) + rec_logL_add;}
        }
      fval+=w_rec_end*f_rec_dev_end;
   }
//	cout<< "log rec devs" << log_dev_rec << endl;
  
  //applies if initial age structure is estimated
  fval+=norm2(log_dev_Nage); 

  //Random walk components of fishery dependent indices: these components equal zero if RW turned off
  f_cpue_RW_cHL=0.0;
  for (iyear=styr_cpue_cHL; iyear<endyr_cpue_cHL; iyear++)
      {f_cpue_RW_cHL+=square(q_RW_log_dev_cHL(iyear))/(2.0*set_q_RW_var_cHL);}
  fval+=f_cpue_RW_cHL;  
  
  f_cpue_RW_rHB=0.0;
  for (iyear=styr_cpue_rHB; iyear<endyr_cpue_rHB; iyear++)
      {f_cpue_RW_rHB+=square(q_RW_log_dev_rHB(iyear))/(2.0*set_q_RW_var_rHB);}      
  fval+=f_cpue_RW_rHB;   
  
//---Priors---------------------------------------------------
//neg_log_prior arguments: estimate, prior, variance, pdf type
//Variance input as a negative value is considered to be CV in arithmetic space (CV=-1 implies loose prior) 
//pdf type(1=none, 2=lognormal, 3=normal, 4=beta)

  f_priors=0.0; 
  f_priors+=neg_log_prior(steep, set_steep(5), set_steep(6), set_steep(7)); 
//	cout<<f_priors<<endl;
  f_priors+=neg_log_prior(rec_sigma,set_rec_sigma(5),set_rec_sigma(6),set_rec_sigma(7));
  f_priors+=neg_log_prior(R_autocorr,set_R_autocorr(5),set_R_autocorr(6),set_R_autocorr(7));
  //f_priors+=neg_log_prior(q_DD_beta, set_q_DD_beta, square(set_q_DD_beta_se), 3);	
  f_priors+=neg_log_prior(q_rate, set_q_rate(5), set_q_rate(6), set_q_rate(7));  
  f_priors+=neg_log_prior(F_init_ratio, set_F_init_ratio(5),set_F_init_ratio(6),set_F_init_ratio(7));
  f_priors+=neg_log_prior(M_constant, set_M_constant(5), set_M_constant(6),set_M_constant(7));	

  //f_priors+=neg_log_prior(Linf,set_Linf,square(set_Linf_se),3);
  //f_priors+=neg_log_prior(K,set_K,square(set_K_se),3);
  //f_priors+=neg_log_prior(t0,set_K,square(set_t0_se),3);
  f_priors+=neg_log_prior(len_cv_val,set_len_cv(5),set_len_cv(6),set_len_cv(7));

  f_priors+=neg_log_prior(log_dm_lenc_sBT,set_log_dm_lenc_sBT(5),set_log_dm_lenc_sBT(6),set_log_dm_lenc_sBT(7));
  // f_priors+=neg_log_prior(log_dm_sTV_lc,set_log_dm_sTV_lc(5),set_log_dm_sTV_lc(6),set_log_dm_sTV_lc(7));
  f_priors+=neg_log_prior(log_dm_lenc_cHL,set_log_dm_lenc_cHL(5),set_log_dm_lenc_cHL(6),set_log_dm_lenc_cHL(7));
  f_priors+=neg_log_prior(log_dm_lenc_cPT,set_log_dm_lenc_cPT(5),set_log_dm_lenc_cPT(6),set_log_dm_lenc_cPT(7));
  f_priors+=neg_log_prior(log_dm_lenc_rHB,set_log_dm_lenc_rHB(5),set_log_dm_lenc_rHB(6),set_log_dm_lenc_rHB(7));
  f_priors+=neg_log_prior(log_dm_lenc_rHD,set_log_dm_lenc_rHD(5),set_log_dm_lenc_rHD(6),set_log_dm_lenc_rHD(7));
  f_priors+=neg_log_prior(log_dm_lenc_rGN,set_log_dm_lenc_rGN(5),set_log_dm_lenc_rGN(6),set_log_dm_lenc_rGN(7));
  f_priors+=neg_log_prior(log_dm_agec_sBT,set_log_dm_agec_sBT(5),set_log_dm_agec_sBT(6),set_log_dm_agec_sBT(7));  
  f_priors+=neg_log_prior(log_dm_agec_sTV,set_log_dm_agec_sTV(5),set_log_dm_agec_sTV(6),set_log_dm_agec_sTV(7));
  f_priors+=neg_log_prior(log_dm_agec_cHL,set_log_dm_agec_cHL(5),set_log_dm_agec_cHL(6),set_log_dm_agec_cHL(7));
  f_priors+=neg_log_prior(log_dm_agec_cPT,set_log_dm_agec_cPT(5),set_log_dm_agec_cPT(6),set_log_dm_agec_cPT(7));
  f_priors+=neg_log_prior(log_dm_agec_rHB,set_log_dm_agec_rHB(5),set_log_dm_agec_rHB(6),set_log_dm_agec_rHB(7));
  // f_priors+=neg_log_prior(log_dm_rGN_ac,set_log_dm_rGN_ac(5),set_log_dm_rGN_ac(6),set_log_dm_rGN_ac(7));


  f_priors+=neg_log_prior(selpar_A50_sBT, set_selpar_A50_sBT(5),set_selpar_A50_sBT(6),set_selpar_A50_sBT(7));  
  f_priors+=neg_log_prior(selpar_slope_sBT, set_selpar_slope_sBT(5),set_selpar_slope_sBT(6),set_selpar_slope_sBT(7));                  
  f_priors+=neg_log_prior(selpar_A502_sBT, set_selpar_A502_sBT(5),set_selpar_A502_sBT(6),set_selpar_A502_sBT(7));  
  f_priors+=neg_log_prior(selpar_slope2_sBT, set_selpar_slope2_sBT(5),set_selpar_slope2_sBT(6),set_selpar_slope2_sBT(7));                  

  f_priors+=neg_log_prior(selpar_A50_sTV, set_selpar_A50_sTV(5),set_selpar_A50_sTV(6),set_selpar_A50_sTV(7));  
  f_priors+=neg_log_prior(selpar_slope_sTV, set_selpar_slope_sTV(5),set_selpar_slope_sTV(6),set_selpar_slope_sTV(7));                  
  f_priors+=neg_log_prior(selpar_A502_sTV, set_selpar_A502_sTV(5),set_selpar_A502_sTV(6),set_selpar_A502_sTV(7));  
  f_priors+=neg_log_prior(selpar_slope2_sTV, set_selpar_slope2_sTV(5),set_selpar_slope2_sTV(6),set_selpar_slope2_sTV(7));                  

  f_priors+=neg_log_prior(selpar_A50b02_cHL, set_selpar_A50b02_cHL(5),set_selpar_A50b02_cHL(6),set_selpar_A50b02_cHL(7));  
  f_priors+=neg_log_prior(selpar_slopeb02_cHL, set_selpar_slopeb02_cHL(5),set_selpar_slopeb02_cHL(6),set_selpar_slopeb02_cHL(7));                  
  f_priors+=neg_log_prior(selpar_A50b03_cHL, set_selpar_A50b03_cHL(5),set_selpar_A50b03_cHL(6),set_selpar_A50b03_cHL(7));  
  f_priors+=neg_log_prior(selpar_slopeb03_cHL, set_selpar_slopeb03_cHL(5),set_selpar_slopeb03_cHL(6),set_selpar_slopeb03_cHL(7)); 
  f_priors+=neg_log_prior(selpar_A50b04_cHL, set_selpar_A50b04_cHL(5),set_selpar_A50b04_cHL(6),set_selpar_A50b04_cHL(7));  
  f_priors+=neg_log_prior(selpar_slopeb04_cHL, set_selpar_slopeb04_cHL(5),set_selpar_slopeb04_cHL(6),set_selpar_slopeb04_cHL(7));   

  f_priors+=neg_log_prior(selpar_A50b02_cPT, set_selpar_A50b02_cPT(5),set_selpar_A50b02_cPT(6),set_selpar_A50b02_cPT(7));  
  f_priors+=neg_log_prior(selpar_slopeb02_cPT, set_selpar_slopeb02_cPT(5),set_selpar_slopeb02_cPT(6),set_selpar_slopeb02_cPT(7));                  
  f_priors+=neg_log_prior(selpar_A50b03_cPT, set_selpar_A50b03_cPT(5),set_selpar_A50b03_cPT(6),set_selpar_A50b03_cPT(7));  
  f_priors+=neg_log_prior(selpar_slopeb03_cPT, set_selpar_slopeb03_cPT(5),set_selpar_slopeb03_cPT(6),set_selpar_slopeb03_cPT(7)); 
  f_priors+=neg_log_prior(selpar_A50b04_cPT, set_selpar_A50b04_cPT(5),set_selpar_A50b04_cPT(6),set_selpar_A50b04_cPT(7));  
  f_priors+=neg_log_prior(selpar_slopeb04_cPT, set_selpar_slopeb04_cPT(5),set_selpar_slopeb04_cPT(6),set_selpar_slopeb04_cPT(7));   

  f_priors+=neg_log_prior(selpar_A50b01_rHB, set_selpar_A50b01_rHB(5),set_selpar_A50b01_rHB(6),set_selpar_A50b01_rHB(7));  
  f_priors+=neg_log_prior(selpar_slopeb01_rHB, set_selpar_slopeb01_rHB(5),set_selpar_slopeb01_rHB(6),set_selpar_slopeb01_rHB(7));                  
  f_priors+=neg_log_prior(selpar_A50b02_rHB, set_selpar_A50b02_rHB(5),set_selpar_A50b02_rHB(6),set_selpar_A50b02_rHB(7));  
  f_priors+=neg_log_prior(selpar_slopeb02_rHB, set_selpar_slopeb02_rHB(5),set_selpar_slopeb02_rHB(6),set_selpar_slopeb02_rHB(7));                  
  f_priors+=neg_log_prior(selpar_A50b03_rHB, set_selpar_A50b03_rHB(5),set_selpar_A50b03_rHB(6),set_selpar_A50b03_rHB(7));  
  f_priors+=neg_log_prior(selpar_slopeb03_rHB, set_selpar_slopeb03_rHB(5),set_selpar_slopeb03_rHB(6),set_selpar_slopeb03_rHB(7));      
  f_priors+=neg_log_prior(selpar_A50b04_rHB, set_selpar_A50b04_rHB(5),set_selpar_A50b04_rHB(6),set_selpar_A50b04_rHB(7)); 
  f_priors+=neg_log_prior(selpar_slopeb04_rHB, set_selpar_slopeb04_rHB(5),set_selpar_slopeb04_rHB(6),set_selpar_slopeb04_rHB(7)); 
  f_priors+=neg_log_prior(selpar_A50b05_rHB, set_selpar_A50b05_rHB(5),set_selpar_A50b05_rHB(6),set_selpar_A50b05_rHB(7));   
  f_priors+=neg_log_prior(selpar_slopeb05_rHB, set_selpar_slopeb05_rHB(5),set_selpar_slopeb05_rHB(6),set_selpar_slopeb05_rHB(7));   

  f_priors+=neg_log_prior(selpar_logitAge0_rHD, set_selpar_logitAge0_rHD(5),set_selpar_logitAge0_rHD(6),set_selpar_logitAge0_rHD(7));
  f_priors+=neg_log_prior(selpar_logitAge1_rHD, set_selpar_logitAge1_rHD(5),set_selpar_logitAge1_rHD(6),set_selpar_logitAge1_rHD(7));
  f_priors+=neg_log_prior(selpar_logitAge2_rHD, set_selpar_logitAge2_rHD(5),set_selpar_logitAge2_rHD(6), set_selpar_logitAge2_rHD(7));
    
  f_priors+=neg_log_prior(selpar_A50b04_rHD, set_selpar_A50b04_rHD(5),set_selpar_A50b04_rHD(6),set_selpar_A50b04_rHD(7)); 
  f_priors+=neg_log_prior(selpar_slopeb04_rHD, set_selpar_slopeb04_rHD(5),set_selpar_slopeb04_rHD(6),set_selpar_slopeb04_rHD(7)); 
  f_priors+=neg_log_prior(selpar_A502b04_rHD, set_selpar_A502b04_rHD(5),set_selpar_A502b04_rHD(6),set_selpar_A502b04_rHD(7)); f_priors+=neg_log_prior(selpar_A50b05_rHD, set_selpar_A50b05_rHD(5),set_selpar_A50b05_rHD(6),set_selpar_A50b05_rHD(7));   
  f_priors+=neg_log_prior(selpar_slopeb05_rHD, set_selpar_slopeb05_rHD(5),set_selpar_slopeb05_rHD(6),set_selpar_slopeb05_rHD(7)); 
  f_priors+=neg_log_prior(selpar_A502b05_rHD, set_selpar_A502b05_rHD(5),set_selpar_A502b05_rHD(6),set_selpar_A502b05_rHD(7));

  f_priors+=neg_log_prior(log_q_sBT,set_log_q_cpue_sBT(5),set_log_q_cpue_sBT(6),set_log_q_cpue_sBT(7));
  f_priors+=neg_log_prior(log_q_sTV,set_log_q_cpue_sTV(5),set_log_q_cpue_sTV(6),set_log_q_cpue_sTV(7));
  f_priors+=neg_log_prior(log_q_cHL,set_log_q_cpue_cHL(5),set_log_q_cpue_cHL(6),set_log_q_cpue_cHL(7));
  f_priors+=neg_log_prior(log_q_rHB,set_log_q_cpue_rHB(5),set_log_q_cpue_rHB(6),set_log_q_cpue_rHB(7));

  f_priors+=neg_log_prior(log_avg_F_L_cHL,set_log_avg_F_L_cHL(5),set_log_avg_F_L_cHL(6),set_log_avg_F_L_cHL(7));
  f_priors+=neg_log_prior(log_avg_F_L_cPT,set_log_avg_F_L_cPT(5),set_log_avg_F_L_cPT(6),set_log_avg_F_L_cPT(7));
  f_priors+=neg_log_prior(log_avg_F_L_cTW,set_log_avg_F_L_cTW(5),set_log_avg_F_L_cTW(6),set_log_avg_F_L_cTW(7));
  f_priors+=neg_log_prior(log_avg_F_L_rHB,set_log_avg_F_L_rHB(5),set_log_avg_F_L_rHB(6),set_log_avg_F_L_rHB(7));
  f_priors+=neg_log_prior(log_avg_F_L_rGN,set_log_avg_F_L_rGN(5),set_log_avg_F_L_rGN(6),set_log_avg_F_L_rGN(7));
  f_priors+=neg_log_prior(log_avg_F_L_rIN,set_log_avg_F_L_rIN(5),set_log_avg_F_L_rIN(6),set_log_avg_F_L_rIN(7));
  f_priors+=neg_log_prior(log_avg_F_D_cGD,set_log_avg_F_D_cGD(5),set_log_avg_F_D_cGD(6),set_log_avg_F_D_cGD(7));
  f_priors+=neg_log_prior(log_avg_F_D_rHD,set_log_avg_F_D_rHD(5),set_log_avg_F_D_rHD(6),set_log_avg_F_D_rHD(7));
  f_priors+=neg_log_prior(log_avg_F_D_rGD,set_log_avg_F_D_rGD(5),set_log_avg_F_D_rGD(6),set_log_avg_F_D_rGD(7));
  f_priors+=neg_log_prior(log_avg_F_D_rID,set_log_avg_F_D_rID(5),set_log_avg_F_D_rID(6),set_log_avg_F_D_rID(7));


  fval+=f_priors;

  //cout << "fval = " << fval << "  fval_data = " << fval_data << endl;


//----------------------------------------------------------------------------------
//Logistic function: 2 parameters
FUNCTION dvar_vector logistic(const dvar_vector& ages, const dvariable& A50, const dvariable& slope)
  //ages=vector of ages, A50=age at 50% selectivity, slope=rate of increase
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=1./(1.+mfexp(-1.*slope*(ages-A50))); //logistic;  
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;

//-----------------------------------------------------------------------------------
//Logistic-exponential: 4 parameters (but 1 is fixed)
FUNCTION dvar_vector logistic_exponential(const dvar_vector& ages, const dvariable& A50, const dvariable& slope, const dvariable& sigma, const dvariable& joint)
  //ages=vector of ages, A50=age at 50% sel (ascending limb), slope=rate of increase, sigma=controls rate of descent (descending)
  //joint=age to join curves
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=1.0;
  for (iage=ages.indexmin(); iage<=ages.indexmax(); iage++)
  {
   if (ages(iage)<joint) {Sel_Tmp(iage)=1./(1.+mfexp(-1.*slope*(ages(iage)-A50)));}
   if (ages(iage)>joint){Sel_Tmp(iage)=mfexp(-1.*square((ages(iage)-joint)/sigma));}
  }
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;
 
//-----------------------------------------------------------------------------------
//Logistic function: 4 parameters
FUNCTION dvar_vector logistic_double(const dvar_vector& ages, const dvariable& A501, const dvariable& slope1, const dvariable& A502, const dvariable& slope2)
  //ages=vector of ages, A50=age at 50% selectivity, slope=rate of increase, A502=age at 50% decrease additive to A501, slope2=slope of decrease
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=elem_prod( (1./(1.+mfexp(-1.*slope1*(ages-A501)))),(1.-(1./(1.+mfexp(-1.*slope2*(ages-(A501+A502)))))) );     
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;

//-----------------------------------------------------------------------------------
//Jointed logistic function: 6 parameters (increasing and decreasing logistics joined at peak selectivity)
FUNCTION dvar_vector logistic_joint(const dvar_vector& ages, const dvariable& A501, const dvariable& slope1, const dvariable& A502, const dvariable& slope2, const dvariable& satval, const dvariable& joint)
  //ages=vector of ages, A501=age at 50% sel (ascending limb), slope1=rate of increase,A502=age at 50% sel (descending), slope1=rate of increase (ascending), 
  //satval=saturation value of descending limb, joint=location in age vector to join curves (may equal age or age + 1 if age-0 is included)
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=1.0; 
  for (iage=ages.indexmin(); iage<=ages.indexmax(); iage++)
  {
   if (double(iage)<joint) {Sel_Tmp(iage)=1./(1.+mfexp(-1.*slope1*(ages(iage)-A501)));}  
   if (double(iage)>joint){Sel_Tmp(iage)=1.0-(1.0-satval)/(1.+mfexp(-1.*slope2*(ages(iage)-A502)));}  
  }  
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;

//-----------------------------------------------------------------------------------  
//Double Gaussian function: 6 parameters (as in SS3)
FUNCTION dvar_vector gaussian_double(const dvar_vector& ages, const dvariable& peak, const dvariable& top, const dvariable& ascwid, const dvariable& deswid, const dvariable& init, const dvariable& final)
  //ages=vector of ages, peak=ascending inflection location (as logistic), top=width of plateau, ascwid=ascent width (as log(width))
  //deswid=descent width (as log(width))
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step1(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step2(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step3(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step4(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step5(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step6(ages.indexmin(),ages.indexmax());
  dvar_vector pars_tmp(1,6); dvar_vector sel_tmp_iq(1,2);
  
  pars_tmp(1)=peak;
  pars_tmp(2)=peak+1.0+(0.99*ages(ages.indexmin())-peak-1.0)/(1.0+mfexp(-top));
  pars_tmp(3)=mfexp(ascwid);
  pars_tmp(4)=mfexp(deswid);
  pars_tmp(5)=1.0/(1.0+mfexp(-init));
  pars_tmp(6)=1.0/(1.0+mfexp(-final));
       
  sel_tmp_iq(1)=mfexp(-(square(ages(ages.indexmin())-pars_tmp(1))/pars_tmp(3)));
  sel_tmp_iq(2)=mfexp(-(square(ages(ages.indexmax())-pars_tmp(2))/pars_tmp(4)));
  
  sel_step1=mfexp(-(square(ages-pars_tmp(1))/pars_tmp(3)));
  sel_step2=pars_tmp(5)+(1.0-pars_tmp(5))*(sel_step1-sel_tmp_iq(1))/(1.0-sel_tmp_iq(1));  
  sel_step3=mfexp(-(square(ages-pars_tmp(2))/pars_tmp(4)));
  sel_step4=1.0+(pars_tmp(6)-1.0)*(sel_step3-1.0)/(sel_tmp_iq(2)-1.0);
  sel_step5=1.0/ (1.0+mfexp(-(20.0* elem_div((ages-pars_tmp(1)), (1.0+sfabs(ages-pars_tmp(1)))) )));
  sel_step6=1.0/(1.0+mfexp(-(20.0*elem_div((ages-pars_tmp(2)),(1.0+sfabs(ages-pars_tmp(2)))) )));  

  Sel_Tmp=elem_prod(sel_step2,(1.0-sel_step5))+ 
          elem_prod(sel_step5,((1.0-sel_step6)+ elem_prod(sel_step4,sel_step6)) ); 
 
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;
    
//-----------------------------------------------------------------------------------    
//Spawner-recruit function (Beverton-Holt, Ricker, or mean)
FUNCTION dvariable SR_func(const dvariable& R0, const dvariable& h, const dvariable& spr_F0, const dvariable& SSB, int func)
  //R0=virgin recruitment, h=steepness, spr_F0=spawners per recruit @ F=0, SSB=spawning biomass
  //func=1 for Beverton-Holt, 2 for Ricker
  RETURN_ARRAYS_INCREMENT();
  dvariable Recruits_Tmp;
  switch(func) {
    case 1: //Beverton-Holt
      Recruits_Tmp=((0.8*R0*h*SSB)/(0.2*R0*spr_F0*(1.0-h)+(h-0.2)*SSB));       
    break;
    case 2: //Ricker
      Recruits_Tmp=((SSB/spr_F0)*mfexp(h*(1-SSB/(R0*spr_F0))));       
    break;
	case 3: //Mean
      Recruits_Tmp=R0;       
    break;
  }
  RETURN_ARRAYS_DECREMENT();
  return Recruits_Tmp;

//-----------------------------------------------------------------------------------    
//Spawner-recruit equilibrium function (Beverton-Holt, Ricker, or mean)
FUNCTION dvariable SR_eq_func(const dvariable& R0, const dvariable& h, const dvariable& spr_F0, const dvariable& spr_F, const dvariable& BC, int func)
  //R0=virgin recruitment, h=steepness, spr_F0=spawners per recruit @ F=0, spr_F=spawners per recruit @ F, BC=bias correction
  //func=1 for Beverton-Holt, 2 for Ricker
  RETURN_ARRAYS_INCREMENT();
  dvariable Recruits_Tmp;
  switch(func) {
    case 1: //Beverton-Holt
      Recruits_Tmp=(R0/((5.0*h-1.0)*spr_F))*(BC*4.0*h*spr_F-spr_F0*(1.0-h));    
    break;
    case 2: //Ricker
      Recruits_Tmp=R0/(spr_F/spr_F0)*(1.0+log(BC*spr_F/spr_F0)/h);      
    break;
    case 3: //Mean
      Recruits_Tmp=BC*R0;      
    break;	
  }
  RETURN_ARRAYS_DECREMENT();
  return Recruits_Tmp;
  
//-----------------------------------------------------------------------------------
//Likelihood contribution: lognormal
FUNCTION dvariable lk_lognormal(const dvar_vector& pred, const dvar_vector& obs, const dvar_vector& cv, const dvariable& wgt_dat)
  //pred=vector of predicted vals, obs=vector of observed vals, cv=vector of CVs in arithmetic space, wgt_dat=constant scaling of CVs
  //dzero is small value to avoid log(0) during search
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvar_vector var(cv.indexmin(),cv.indexmax()); //variance in log space
  var=log(1.0+square(cv/wgt_dat));   // convert cv in arithmetic space to variance in log space
  LkvalTmp=sum(0.5*elem_div(square(log(elem_div((pred+dzero),(obs+dzero)))),var) );
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;

//-----------------------------------------------------------------------------------
//Likelihood contribution: multinomial
FUNCTION dvariable lk_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const double& minSS, const dvariable& wgt_dat)
  //nsamp=vector of N's, pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, ncomp = number of yrs in matrix, minSS=min N threshold, wgt_dat=scaling of N's
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  LkvalTmp=0.0;
  for (int ii=1; ii<=ncomp; ii++)
  {if (nsamp(ii)>=minSS)
    {LkvalTmp-=wgt_dat*nsamp(ii)*sum(elem_prod((obs_comp(ii)+dzero),
               log(elem_div((pred_comp(ii)+dzero), (obs_comp(ii)+dzero)))));
    }
  }  
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;

  //-----------------------------------------------------------------------------------
//Likelihood contribution: Dirichlet-multinomial
FUNCTION dvariable lk_dirichlet_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const dvariable& mbin, const double& minSS, const dvariable& log_dir_par)
  //nsamp=vector of N's, pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, ncomp = number of yrs in matrix, mbin=number of bins, minSS=min N threshold, wgt_dat=scaling of N's
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvariable small_number=0.00001;
  LkvalTmp=0.0; 
  dvar_vector nsamp_adjust=nsamp*mfexp(log_dir_par);
   //dvar_vector nsamp_adjust=mfexp(log_dir_par);
  for (int ii=1; ii<=ncomp; ii++)
  {
	if (nsamp(ii)>=minSS)
    {
		LkvalTmp-=gammln(nsamp_adjust(ii))-gammln(nsamp(ii)+nsamp_adjust(ii));
		LkvalTmp-=sum(gammln(nsamp(ii)*obs_comp(ii)+nsamp_adjust(ii)*pred_comp(ii)+small_number));
		LkvalTmp+=sum(gammln(nsamp_adjust(ii)*pred_comp(ii)+small_number));
    }
  }  
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;

//-----------------------------------------------------------------------------------
//Likelihood contribution: priors
FUNCTION  dvariable neg_log_prior(dvariable pred, const double& prior, dvariable var, int pdf)
  //prior=prior point estimate, var=variance (if negative, treated as CV in arithmetic space), pred=predicted value, pdf=prior type (1=none, 2=lognormal, 3=normal, 4=beta)
    dvariable LkvalTmp;
    dvariable alpha, beta, ab_iq;
    LkvalTmp=0.0;
    // compute generic pdf's
    switch(pdf) {
        case 1: //option to turn off prior
          LkvalTmp=0.0;
          break;
        case 2: // lognormal 
          if(prior<=0.0) cout << "YIKES: Don't use a lognormal distn for a negative prior" << endl;
          else if(pred<=0) LkvalTmp=huge_number;
          else {
            if(var<0.0) var=log(1.0+var*var) ;      // convert cv to variance on log scale
            LkvalTmp= 0.5*( square(log(pred/prior))/var + log(var) );
          }
	    break;
        case 3: // normal
          if(var<0.0 && prior!=0.0) var=square(var*prior);       // convert cv to variance on observation scale
          else if(var<0.0 && prior==0.0) var=-var;               // cv not really appropriate if prior value equals zero
          LkvalTmp= 0.5*( square(pred-prior)/var + log(var) );
          break;
        case 4: // beta
          if(var<0.0) var=square(var*prior);          // convert cv to variance on observation scale
          if(prior<=0.0 || prior>=1.0) cout << "YIKES: Don't use a beta distn for a prior outside (0,1)" << endl;
          ab_iq=prior*(1.0-prior)/var - 1.0; alpha=prior*ab_iq; beta=(1.0-prior)*ab_iq;
          if(pred>=0 && pred<=1) LkvalTmp= (1.0-alpha)*log(pred)+(1.0-beta)*log(1.0-pred)-gammln(alpha+beta)+gammln(alpha)+gammln(beta);
          else LkvalTmp=huge_number;
          break;
        default: // no such prior pdf currently available
          cout << "The prior must be either 1(lognormal), 2(normal), or 3(beta)." << endl;
          cout << "Presently it is " << pdf << endl;
          exit(0);
    }
    return LkvalTmp;
 
//-----------------------------------------------------------------------------------
//SDNR: age comp likelihood (assumes fits are done with the robust multinomial function)
FUNCTION dvariable sdnr_multinomial(const double& ncomp, const dvar_vector& ages, const dvar_vector& nsamp, 
                                    const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const dvariable& wgt_dat)
  //ncomp=number of years of data, ages=vector of ages, nsamp=vector of N's, 
  //pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, wgt_dat=likelihood weight for data source
  RETURN_ARRAYS_INCREMENT();
  dvariable SdnrTmp;
  dvar_vector o(1,ncomp);  
  dvar_vector p(1,ncomp);  
  dvar_vector ose(1,ncomp);  
  dvar_vector res(1,ncomp);
  SdnrTmp=0.0;
  for (int ii=1; ii<=ncomp; ii++)
  {
    o(ii)=sum(elem_prod(ages,obs_comp(ii)));
    p(ii)=sum(elem_prod(ages,pred_comp(ii)));
    ose(ii)=sqrt((sum(elem_prod(square(ages),pred_comp(ii)))-square(p(ii)))/(nsamp(ii)*wgt_dat));
  }
  res=elem_div((o-p),ose); 
  SdnrTmp=sqrt(sum(square(res-(sum(res)/ncomp))/(ncomp-1.0))); 
  RETURN_ARRAYS_DECREMENT();
  return SdnrTmp;

//-----------------------------------------------------------------------------------
//SDNR: lognormal likelihood
FUNCTION dvariable sdnr_lognormal(const dvar_vector& pred, const dvar_vector& obs, const dvar_vector& cv, const dvariable& wgt_dat)
  //nyr=number of years of data, pred=vector of predicted data, obs=vector of observed data, cv=vector of cv's, wgt_dat=likelihood weight for data source
  RETURN_ARRAYS_INCREMENT();
  dvariable SdnrTmp;
  dvariable small_number=0.00001;
  dvariable n;
  dvar_vector res(cv.indexmin(),cv.indexmax());
  SdnrTmp=0.0;
  res=elem_div(log(elem_div(obs+small_number,pred+small_number)),sqrt(log(1+square(cv/wgt_dat))));
  n=cv.indexmax()-cv.indexmin()+1;
  SdnrTmp=sqrt(sum(square(res-(sum(res)/n))/(n-1.0))); 
  RETURN_ARRAYS_DECREMENT();
  return SdnrTmp; 

//-----------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------
REPORT_SECTION

  if (last_phase())  
  {

      cout<<"start report"<<endl;
      get_weighted_current();
      cout<<"got weighted"<<endl;
      get_msy();
      cout<<"got msy"<<endl;
      get_per_recruit_stuff();
      cout<<"got per recruit"<<endl;  
      get_miscellaneous_stuff();
      cout<<"got misc stuff"<<endl;
      get_projection();
      cout<<"got projection"<<endl;
	  
      grad_max=objective_function_value::pobjfun->gmax;
      time(&finish);
	  elapsed_time=difftime(finish,start);
	  hour=long(elapsed_time)/3600;
	  minute=long(elapsed_time)%3600/60;
	  second=(long(elapsed_time)%3600)%60;
	  cout<<endl<<endl<<"*******************************************"<<endl;
	  cout<<"--Start time: "<<ctime(&start)<<endl;
	  cout<<"--Finish time: "<<ctime(&finish)<<endl;
	  cout<<"--Runtime: ";
	  cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
	  cout << "--TotalLikelihood: " << fval << endl;
      cout<<"--Final gradient: "<<objective_function_value::pobjfun->gmax << endl;
	  cout<<"*******************************************"<<endl;
      cout <<endl;     
      cout << "><>--><>--><>--><>--><>--><>--><>--><>--><>--><>"  <<endl;
      cout << "BC Fmsy=" << F_msy_lb_out<< "   BC SSBmsy=" << SSB_msy_lb_out <<endl;
      cout <<"F status="<<FdF_msy_lb_end<<endl;
      cout <<"Pop status="<<SdSSB_msy_lb_end<<endl;
      cout << "h="<<steep<<"   R0="<<R0<<endl;
      cout << "><>--><>--><>--><>--><>--><>--><>--><>--><>--><>"  <<endl;  
       
      report << "TotalLikelihood " << fval << endl;
      report << "N" << endl;
      report << N<<endl;      
      report << "SSB" << endl;       
      report << SSB << endl;
	    
	    sdnr_I_sBT=sdnr_lognormal(pred_cpue_sBT, obs_cpue_sBT, obs_cv_cpue_sBT, w_cpue_sBT);  
	    sdnr_I_sTV=sdnr_lognormal(pow(pred_cpue_sTV,posobs_cpue_sTV), obs_cpue_corrected_sTV, obs_cv_cpue_corrected_sTV, w_cpue_sTV);  
	    sdnr_I_cHL=sdnr_lognormal(pred_cpue_cHL, obs_cpue_cHL, obs_cv_cpue_cHL, w_cpue_cHL);
	    sdnr_I_rHB=sdnr_lognormal(pred_cpue_rHB, obs_cpue_rHB, obs_cv_cpue_rHB, w_cpue_rHB);
	    
	    
	    Linf_out(8)=Linf; Linf_out(1,7)=set_Linf; 
	    K_out(8)=K; K_out(1,7)=set_K;
	    t0_out(8)=t0; t0_out(1,7)=set_t0;
	    len_cv_val_out(8)=len_cv_val; len_cv_val_out(1,7)=set_len_cv;
	    
	    log_R0_out(8)=log_R0; log_R0_out(1,7)=set_log_R0;
	    M_constant_out(8)=M_constant; M_constant_out(1,7)=set_M_constant;
	    steep_out(8)=steep; steep_out(1,7)=set_steep;
	    rec_sigma_out(8)=rec_sigma; rec_sigma_out(1,7)=set_rec_sigma;
	    R_autocorr_out(8)=R_autocorr; R_autocorr_out(1,7)=set_R_autocorr;
	    
	    log_dm_sBT_lc_out(8)=log_dm_lenc_sBT; log_dm_sBT_lc_out(1,7)=set_log_dm_lenc_sBT;
	    // log_dm_sTV_lc_out(8)=log_dm_sTV_lc; log_dm_sTV_lc_out(1,7)=set_log_dm_sTV_lc;
	    log_dm_cHL_lc_out(8)=log_dm_lenc_cHL; log_dm_cHL_lc_out(1,7)=set_log_dm_lenc_cHL;
	    log_dm_cPT_lc_out(8)=log_dm_lenc_cPT; log_dm_cPT_lc_out(1,7)=set_log_dm_lenc_cPT;
	    log_dm_rHB_lc_out(8)=log_dm_lenc_rHB; log_dm_rHB_lc_out(1,7)=set_log_dm_lenc_rHB;
	    log_dm_rHD_lc_out(8)=log_dm_lenc_rHD; log_dm_rHD_lc_out(1,7)=set_log_dm_lenc_rHD;
	    log_dm_rGN_lc_out(8)=log_dm_lenc_rGN; log_dm_rGN_lc_out(1,7)=set_log_dm_lenc_rGN;
	    
	    log_dm_sBT_ac_out(8)=log_dm_agec_sBT; log_dm_sBT_ac_out(1,7)=set_log_dm_agec_sBT;
	    log_dm_sTV_ac_out(8)=log_dm_agec_sTV; log_dm_sTV_ac_out(1,7)=set_log_dm_agec_sTV;
	    log_dm_cHL_ac_out(8)=log_dm_agec_cHL; log_dm_cHL_ac_out(1,7)=set_log_dm_agec_cHL;
	    log_dm_cPT_ac_out(8)=log_dm_agec_cPT; log_dm_cPT_ac_out(1,7)=set_log_dm_agec_cPT;
	    log_dm_rHB_ac_out(8)=log_dm_agec_rHB; log_dm_rHB_ac_out(1,7)=set_log_dm_agec_rHB;
	    //log_dm_rGN_ac_out(8)=log_dm_rGN_ac; log_dm_rGN_ac_out(1,7)=set_log_dm_rGN_ac;
	   
	   
       selpar_A50_sBT_out(8)=selpar_A50_sBT; selpar_A50_sBT_out(1,7)=set_selpar_A50_sBT;
       selpar_slope_sBT_out(8)=selpar_slope_sBT; selpar_slope_sBT_out(1,7)=set_selpar_slope_sBT;
	     selpar_A502_sBT_out(8)=selpar_A502_sBT; selpar_A502_sBT_out(1,7)=set_selpar_A502_sBT;
       selpar_slope2_sBT_out(8)=selpar_slope2_sBT; selpar_slope2_sBT_out(1,7)=set_selpar_slope2_sBT;
       
	     selpar_A50_sTV_out(8)=selpar_A50_sTV; selpar_A50_sTV_out(1,7)=set_selpar_A50_sTV;
       selpar_slope_sTV_out(8)=selpar_slope_sTV; selpar_slope_sTV_out(1,7)=set_selpar_slope_sTV;
	     selpar_A502_sTV_out(8)=selpar_A502_sTV; selpar_A502_sTV_out(1,7)=set_selpar_A502_sTV;
       selpar_slope2_sTV_out(8)=selpar_slope2_sTV; selpar_slope2_sTV_out(1,7)=set_selpar_slope2_sTV;

       
       //selpar_A50b01_cHL_out(8)=selpar_A50b01_cHL; selpar_A50b01_cHL_out(1,7)=set_selpar_A50b01_cHL;
       //selpar_slopeb01_cHL_out(8)=selpar_slopeb01_cHL; selpar_slopeb01_cHL_out(1,7)=set_selpar_slopeb01_cHL;
       selpar_A50b02_cHL_out(8)=selpar_A50b02_cHL; selpar_A50b02_cHL_out(1,7)=set_selpar_A50b02_cHL;
       selpar_slopeb02_cHL_out(8)=selpar_slopeb02_cHL; selpar_slopeb02_cHL_out(1,7)=set_selpar_slopeb02_cHL;
       selpar_A50b03_cHL_out(8)=selpar_A50b03_cHL; selpar_A50b03_cHL_out(1,7)=set_selpar_A50b03_cHL;
       selpar_slopeb03_cHL_out(8)=selpar_slopeb03_cHL; selpar_slopeb03_cHL_out(1,7)=set_selpar_slopeb03_cHL;
       selpar_A50b04_cHL_out(8)=selpar_A50b04_cHL; selpar_A50b04_cHL_out(1,7)=set_selpar_A50b04_cHL;
       selpar_slopeb04_cHL_out(8)=selpar_slopeb04_cHL; selpar_slopeb04_cHL_out(1,7)=set_selpar_slopeb04_cHL;
	   
	     selpar_A50b02_cPT_out(8)=selpar_A50b02_cPT; selpar_A50b02_cPT_out(1,7)=set_selpar_A50b02_cPT;
       selpar_slopeb02_cPT_out(8)=selpar_slopeb02_cPT; selpar_slopeb02_cPT_out(1,7)=set_selpar_slopeb02_cPT;
       selpar_A50b03_cPT_out(8)=selpar_A50b03_cPT; selpar_A50b03_cPT_out(1,7)=set_selpar_A50b03_cPT;
       selpar_slopeb03_cPT_out(8)=selpar_slopeb03_cPT; selpar_slopeb03_cPT_out(1,7)=set_selpar_slopeb03_cPT;
	     selpar_A50b04_cPT_out(8)=selpar_A50b04_cPT; selpar_A50b04_cPT_out(1,7)=set_selpar_A50b04_cPT;
       selpar_slopeb04_cPT_out(8)=selpar_slopeb04_cPT; selpar_slopeb04_cPT_out(1,7)=set_selpar_slopeb04_cPT;
          
       selpar_A50b01_rHB_out(8)=selpar_A50b01_rHB; selpar_A50b01_rHB_out(1,7)=set_selpar_A50b01_rHB;
       selpar_slopeb01_rHB_out(8)=selpar_slopeb01_rHB; selpar_slopeb01_rHB_out(1,7)=set_selpar_slopeb01_rHB;
       selpar_A50b02_rHB_out(8)=selpar_A50b02_rHB; selpar_A50b02_rHB_out(1,7)=set_selpar_A50b02_rHB;
       selpar_slopeb02_rHB_out(8)=selpar_slopeb02_rHB; selpar_slopeb02_rHB_out(1,7)=set_selpar_slopeb02_rHB;
       selpar_A50b03_rHB_out(8)=selpar_A50b03_rHB; selpar_A50b03_rHB_out(1,7)=set_selpar_A50b03_rHB;
       selpar_slopeb03_rHB_out(8)=selpar_slopeb03_rHB; selpar_slopeb03_rHB_out(1,7)=set_selpar_slopeb03_rHB;
	     selpar_A50b04_rHB_out(8)=selpar_A50b04_rHB; selpar_A50b04_rHB_out(1,7)=set_selpar_A50b04_rHB;
       selpar_slopeb04_rHB_out(8)=selpar_slopeb04_rHB; selpar_slopeb04_rHB_out(1,7)=set_selpar_slopeb04_rHB;
	     selpar_A50b05_rHB_out(8)=selpar_A50b05_rHB; selpar_A50b05_rHB_out(1,7)=set_selpar_A50b05_rHB;
       selpar_slopeb05_rHB_out(8)=selpar_slopeb05_rHB; selpar_slopeb05_rHB_out(1,7)=set_selpar_slopeb05_rHB;

	     selpar_logitAge0_rHD_out(8)=selpar_logitAge0_rHD; selpar_logitAge0_rHD_out(1,7)=set_selpar_logitAge0_rHD;
       selpar_logitAge1_rHD_out(8)=selpar_logitAge1_rHD; selpar_logitAge1_rHD_out(1,7)=set_selpar_logitAge1_rHD;
	     selpar_logitAge2_rHD_out(8)=selpar_logitAge2_rHD; selpar_logitAge2_rHD_out(1,7)=set_selpar_logitAge2_rHD;

	     selpar_A50b04_rHD_out(8)=selpar_A50b04_rHD; selpar_A50b04_rHD_out(1,7)=set_selpar_A50b04_rHD;
       selpar_slopeb04_rHD_out(8)=selpar_slopeb04_rHD; selpar_slopeb04_rHD_out(1,7)=set_selpar_slopeb04_rHD;
	     selpar_A502b04_rHD_out(8)=selpar_A502b04_rHD; selpar_A502b04_rHD_out(1,7)=set_selpar_A502b04_rHD;
       selpar_slope2b04_rHD_out(8)=selpar_slope2b04_rHD; selpar_slope2b04_rHD_out(1,7)=set_selpar_slope2b04_rHD;
	     selpar_A50b05_rHD_out(8)=selpar_A50b05_rHD; selpar_A50b05_rHD_out(1,7)=set_selpar_A50b05_rHD;
       selpar_slopeb05_rHD_out(8)=selpar_slopeb05_rHD; selpar_slopeb05_rHD_out(1,7)=set_selpar_slopeb05_rHD;
	     selpar_A502b05_rHD_out(8)=selpar_A502b05_rHD; selpar_A502b05_rHD_out(1,7)=set_selpar_A502b05_rHD;
       selpar_slope2b05_rHD_out(8)=selpar_slope2b05_rHD; selpar_slope2b05_rHD_out(1,7)=set_selpar_slope2b05_rHD;
	   
       log_q_sBT_out(8)=log_q_sBT; log_q_sBT_out(1,7)=set_log_q_cpue_sBT;
	     log_q_sTV_out(8)=log_q_sTV; log_q_sTV_out(1,7)=set_log_q_cpue_sTV;
       log_q_cHL_out(8)=log_q_cHL; log_q_cHL_out(1,7)=set_log_q_cpue_cHL;
	     log_q_rHB_out(8)=log_q_rHB; log_q_rHB_out(1,7)=set_log_q_cpue_rHB;
                            
       out_log_avg_F_L_cHL(8)=log_avg_F_L_cHL; out_log_avg_F_L_cHL(1,7)=set_log_avg_F_L_cHL;
	     out_log_avg_F_L_cPT(8)=log_avg_F_L_cPT; out_log_avg_F_L_cPT(1,7)=set_log_avg_F_L_cPT;
       out_log_avg_F_L_rHB(8)=log_avg_F_L_rHB; out_log_avg_F_L_rHB(1,7)=set_log_avg_F_L_rHB;
       out_log_avg_F_L_cTW(8)=log_avg_F_L_cTW; out_log_avg_F_L_cTW(1,7)=set_log_avg_F_L_cTW;       
       out_log_avg_F_L_rGN(8)=log_avg_F_L_rGN; out_log_avg_F_L_rGN(1,7)=set_log_avg_F_L_rGN;
       out_log_avg_F_L_rIN(8)=log_avg_F_L_rIN; out_log_avg_F_L_rIN(1,7)=set_log_avg_F_L_rIN;
       out_log_avg_F_D_cGD(8)=log_avg_F_D_cGD; out_log_avg_F_D_cGD(1,7)=set_log_avg_F_D_cGD;
       out_log_avg_F_D_rHD(8)=log_avg_F_D_rHD; out_log_avg_F_D_rHD(1,7)=set_log_avg_F_D_rHD;
       out_log_avg_F_D_rGD(8)=log_avg_F_D_rGD; out_log_avg_F_D_rGD(1,7)=set_log_avg_F_D_rGD;
       out_log_avg_F_D_rID(8)=log_avg_F_D_rID; out_log_avg_F_D_rID(1,7)=set_log_avg_F_D_rID;
              
       F_init_ratio_out(8)=F_init_ratio; F_init_ratio_out(1,7)=set_F_init_ratio;
       q_rate_out(8)=q_rate; q_rate_out(1,7)=set_q_rate;       
       
       log_rec_dev_out(styr_rec_dev, endyr_rec_dev)=log_dev_rec;
       log_dev_F_L_out_cHL(styr_L_cHL,endyr_L_cHL)=log_dev_F_L_cHL;
       log_dev_F_L_out_cPT(styr_L_cPT,endyr_L_cPT)=log_dev_F_L_cPT;
       log_dev_F_L_out_cTW(styr_L_cTW,endyr_L_cTW)=log_dev_F_L_cTW;
       log_dev_F_L_out_rHB(styr_L_rHB,endyr_L_rHB)=log_dev_F_L_rHB;
       log_dev_F_L_out_rGN(styr_L_rGN,endyr_L_rGN)=log_dev_F_L_rGN;
       log_dev_F_L_out_rIN(styr_L_rIN,endyr_L_rIN)=log_dev_F_L_rIN;
       log_dev_F_D_out_cGD(styr_D_cGD,endyr_D_cGD)=log_dev_F_D_cGD;
       log_dev_F_D_out_rHD(styr_D_rHD,endyr_D_rHD)=log_dev_F_D_rHD;
       log_dev_F_D_out_rGD(styr_D_rGD,endyr_D_rGD)=log_dev_F_D_rGD;
       log_dev_F_D_out_rID(styr_D_rID,endyr_D_rID)=log_dev_F_D_rID;
       
       out_neff_allyr_lenc_sBT=(1.0+nsamp_allyr_lenc_sBT*mfexp(log_dm_sBT_lc_out(8))) / (1.0 + mfexp(log_dm_sBT_lc_out(8)));
       out_neff_allyr_lenc_cHL=(1.0+nsamp_allyr_lenc_cHL*mfexp(log_dm_cHL_lc_out(8))) / (1.0 + mfexp(log_dm_cHL_lc_out(8)));
       out_neff_allyr_lenc_cPT=(1.0+nsamp_allyr_lenc_cPT*mfexp(log_dm_cPT_lc_out(8))) / (1.0 + mfexp(log_dm_cPT_lc_out(8)));
       out_neff_allyr_lenc_rHB=(1.0+nsamp_allyr_lenc_rHB*mfexp(log_dm_rHB_lc_out(8))) / (1.0 + mfexp(log_dm_rHB_lc_out(8)));
       out_neff_allyr_lenc_rGN=(1.0+nsamp_allyr_lenc_rGN*mfexp(log_dm_rGN_lc_out(8))) / (1.0 + mfexp(log_dm_rGN_lc_out(8)));
       out_neff_allyr_lenc_rHD=(1.0+nsamp_allyr_lenc_rHD*mfexp(log_dm_rHD_lc_out(8))) / (1.0 + mfexp(log_dm_rHD_lc_out(8)));
       
       out_neff_allyr_agec_sBT=(1.0+nsamp_allyr_agec_sBT*mfexp(log_dm_sBT_ac_out(8))) / (1.0 + mfexp(log_dm_sBT_ac_out(8)));
       out_neff_allyr_agec_sTV=(1.0+nsamp_allyr_agec_sTV*mfexp(log_dm_sTV_ac_out(8))) / (1.0 + mfexp(log_dm_sTV_ac_out(8)));
       out_neff_allyr_agec_cHL=(1.0+nsamp_allyr_agec_cHL*mfexp(log_dm_cHL_ac_out(8))) / (1.0 + mfexp(log_dm_cHL_ac_out(8)));
       out_neff_allyr_agec_cPT=(1.0+nsamp_allyr_agec_cPT*mfexp(log_dm_cPT_ac_out(8))) / (1.0 + mfexp(log_dm_cPT_ac_out(8)));
       out_neff_allyr_agec_rHB=(1.0+nsamp_allyr_agec_rHB*mfexp(log_dm_rHB_ac_out(8))) / (1.0 + mfexp(log_dm_rHB_ac_out(8)));

// Fill in years with samples less than minimum cutoff for Neff calculations for MCBE resampling

      for (iyear=1; iyear<=nyr_lenc_sBT; iyear++)
         { if(nsamp_lenc_sBT(iyear)>=0 & nsamp_lenc_sBT(iyear)<minSS_lenc) {
	      nsamp_allyr_lenc_sBT(yrs_lenc_sBT(iyear))=0;
              nfish_allyr_lenc_sBT(yrs_lenc_sBT(iyear))=nfish_lenc_sBT(iyear);}}
            
      for (iyear=1; iyear<=nyr_lenc_cHL; iyear++)
         {  if(nsamp_lenc_cHL(iyear)>=0 & nsamp_lenc_cHL(iyear)<minSS_lenc) {
	      nsamp_allyr_lenc_cHL(yrs_lenc_cHL(iyear))=0;
              nfish_allyr_lenc_cHL(yrs_lenc_cHL(iyear))=nfish_lenc_cHL(iyear);}}
         
      for (iyear=1; iyear<=nyr_lenc_cPT; iyear++)
         { if(nsamp_lenc_cPT(iyear)>=0 & nsamp_lenc_cPT(iyear)<minSS_lenc) {
	      nsamp_allyr_lenc_cPT(yrs_lenc_cPT(iyear))=0;
              nfish_allyr_lenc_cPT(yrs_lenc_cPT(iyear))=nfish_lenc_cPT(iyear);}}
      
      for (iyear=1; iyear<=nyr_lenc_rHB; iyear++)
         { if(nsamp_lenc_rHB(iyear)>=0 & nsamp_lenc_rHB(iyear)<minSS_lenc) {
	      nsamp_allyr_lenc_rHB(yrs_lenc_rHB(iyear))=0;
              nfish_allyr_lenc_rHB(yrs_lenc_rHB(iyear))=nfish_lenc_rHB(iyear);}}

      for (iyear=1; iyear<=nyr_lenc_rHD; iyear++)
         {  if(nsamp_lenc_rHD(iyear)>=0 & nsamp_lenc_rHD(iyear)<minSS_lenc) {
	      nsamp_allyr_lenc_rHD(yrs_lenc_rHD(iyear))=0;
              nfish_allyr_lenc_rHD(yrs_lenc_rHD(iyear))=nfish_lenc_rHD(iyear);}}
         
      for (iyear=1; iyear<=nyr_lenc_rGN; iyear++)
         {  if(nsamp_lenc_rGN(iyear)>=0 & nsamp_lenc_rGN(iyear)<minSS_lenc) {
	      nsamp_allyr_lenc_rGN(yrs_lenc_rGN(iyear))=0;
              nfish_allyr_lenc_rGN(yrs_lenc_rGN(iyear))=nfish_lenc_rGN(iyear);}}

                  
      for (iyear=1; iyear<=nyr_agec_sBT; iyear++)
         { if (nsamp_agec_sBT(iyear)>=0 & nsamp_agec_sBT(iyear)<minSS_agec)
             {nsamp_allyr_agec_sBT(yrs_agec_sBT(iyear))=0;
              nfish_allyr_agec_sBT(yrs_agec_sBT(iyear))=nfish_agec_sBT(iyear);}} 

      for (iyear=1; iyear<=nyr_agec_sTV; iyear++)
         {  if (nsamp_agec_sTV(iyear)>=0 & nsamp_agec_sTV(iyear)<minSS_agec)
             {nsamp_allyr_agec_sTV(yrs_agec_sTV(iyear))=0;
              nfish_allyr_agec_sTV(yrs_agec_sTV(iyear))=nfish_agec_sTV(iyear);}} 

      for (iyear=1; iyear<=nyr_agec_cHL; iyear++)
         { if (nsamp_agec_cHL(iyear)>=0 & nsamp_agec_cHL(iyear)<minSS_agec)
             {nsamp_allyr_agec_cHL(yrs_agec_cHL(iyear))=0;
              nfish_allyr_agec_cHL(yrs_agec_cHL(iyear))=nfish_agec_cHL(iyear);}} 
            
      for (iyear=1; iyear<=nyr_agec_cPT; iyear++)
         {  if (nsamp_agec_cPT(iyear)>=0 & nsamp_agec_cPT(iyear)<minSS_agec)
             {nsamp_allyr_agec_cPT(yrs_agec_cPT(iyear))=0;
              nfish_allyr_agec_cPT(yrs_agec_cPT(iyear))=nfish_agec_cPT(iyear);}}
             
      for (iyear=1; iyear<=nyr_agec_rHB; iyear++)
         { if (nsamp_agec_rHB(iyear)>=0 & nsamp_agec_rHB(iyear)<minSS_agec)
             {nsamp_allyr_agec_rHB(yrs_agec_rHB(iyear))=0;
              nfish_allyr_agec_rHB(yrs_agec_rHB(iyear))=nfish_agec_rHB(iyear);}}

    SSBend_MSST_lb=SdSSB_msy_lb_end/(1.0-M_constant);
    SSBend_MSST_n=SdSSB_msy_n_end/(1.0-M_constant);
     #include "bam.cxx"   // write the S-compatible report

  } //endl last phase loop     
  
  
