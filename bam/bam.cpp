#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
  #include "admodel.h"          // Include AD class definitions
  #include "admb2r.cpp"    // Include S-compatible output functions (needs preceding)
	time_t start,finish;
	long hour,minute,second;	
	double elapsed_time;
#ifdef DEBUG
  #include <chrono>
#endif
#include <admodel.h>
#ifdef USE_ADMB_CONTRIBS
#include <contrib.h>

#endif
  extern "C"  {
    void ad_boundf(int i);
  }
#include <bam.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  adstring tmpstring;
  tmpstring=adprogram_name + adstring(".dat");
  if (argc > 1)
  {
    int on=0;
    if ( (on=option_match(argc,argv,"-ind"))>-1)
    {
      if (on>argc-2 || argv[on+1][0] == '-')
      {
        cerr << "Invalid input data command line option"
                " -- ignored" << endl;
      }
      else
      {
        tmpstring = adstring(argv[on+1]);
      }
    }
  }
  global_datafile = new cifstream(tmpstring);
  if (!global_datafile)
  {
    cerr << "Error: Unable to allocate global_datafile in model_data constructor.";
    ad_exit(1);
  }
  if (!(*global_datafile))
  {
    delete global_datafile;
    global_datafile=NULL;
  }
cout << "Starting Beaufort Assessment Model" << endl;
cout << endl;
cout << "                BAM!" << endl;
cout << endl;
  styr.allocate("styr");
  endyr.allocate("endyr");
  styr_rec_dev.allocate("styr_rec_dev");
  endyr_rec_dev.allocate("endyr_rec_dev");
  styr_rec_mu.allocate("styr_rec_mu");
  endyr_rec_phase1.allocate("endyr_rec_phase1");
  endyr_rec_phase2.allocate("endyr_rec_phase2");
  nyr.allocate("nyr");
  nyr_rec_dev.allocate("nyr_rec_dev");
  nyr_rec_mu.allocate("nyr_rec_mu");
  endyr_sel_b01.allocate("endyr_sel_b01");
  endyr_sel_b02.allocate("endyr_sel_b02");
  endyr_sel_b03.allocate("endyr_sel_b03");
  endyr_sel_b04.allocate("endyr_sel_b04");
  styr_closure_com.allocate("styr_closure_com");
  sizelim1.allocate("sizelim1");
  sizelim2.allocate("sizelim2");
  sizelim3.allocate("sizelim3");
  sizelim4.allocate("sizelim4");
  sizelim5.allocate("sizelim5");
  nages.allocate("nages");
  agebins.allocate(1,nages,"agebins");
  nlenbins.allocate("nlenbins");
  lenbins_width.allocate("lenbins_width");
  lenbins.allocate(1,nlenbins,"lenbins");
  max_F_spr_msy.allocate("max_F_spr_msy");
  n_iter_msy.allocate("n_iter_msy");
  n_iter_spr.allocate("n_iter_spr");
  styr_rec_spr.allocate("styr_rec_spr");
  endyr_rec_spr.allocate("endyr_rec_spr");
  nyr_rec_spr.allocate("nyr_rec_spr");
  SPR_rec_switch.allocate("SPR_rec_switch");
  spr_proxy.allocate("spr_proxy");
  SR_switch.allocate("SR_switch");
  nyr_selpar_wgted.allocate("nyr_selpar_wgted");
  set_BiasCor.allocate("set_BiasCor");
  styr_cpue_sBT.allocate("styr_cpue_sBT");
  endyr_cpue_sBT.allocate("endyr_cpue_sBT");
  obs_cpue_sBT.allocate(styr_cpue_sBT,endyr_cpue_sBT,"obs_cpue_sBT");
  obs_cv_cpue_sBT.allocate(styr_cpue_sBT,endyr_cpue_sBT,"obs_cv_cpue_sBT");
  nyr_lenc_sBT.allocate("nyr_lenc_sBT");
  yrs_lenc_sBT.allocate(1,nyr_lenc_sBT,"yrs_lenc_sBT");
  nsamp_lenc_sBT.allocate(1,nyr_lenc_sBT,"nsamp_lenc_sBT");
  nfish_lenc_sBT.allocate(1,nyr_lenc_sBT,"nfish_lenc_sBT");
  obs_lenc_sBT.allocate(1,nyr_lenc_sBT,1,nlenbins,"obs_lenc_sBT");
  nyr_agec_sBT.allocate("nyr_agec_sBT");
  yrs_agec_sBT.allocate(1,nyr_agec_sBT,"yrs_agec_sBT");
  nsamp_agec_sBT.allocate(1,nyr_agec_sBT,"nsamp_agec_sBT");
  nfish_agec_sBT.allocate(1,nyr_agec_sBT,"nfish_agec_sBT");
  obs_agec_sBT.allocate(1,nyr_agec_sBT,1,nages,"obs_agec_sBT");
  styr_cpue_sTV.allocate("styr_cpue_sTV");
  endyr_cpue_sTV.allocate("endyr_cpue_sTV");
  obs_cpue_sTV.allocate(styr_cpue_sTV,endyr_cpue_sTV,"obs_cpue_sTV");
  obs_cv_cpue_sTV.allocate(styr_cpue_sTV,endyr_cpue_sTV,"obs_cv_cpue_sTV");
  nyr_agec_sTV.allocate("nyr_agec_sTV");
  yrs_agec_sTV.allocate(1,nyr_agec_sTV,"yrs_agec_sTV");
  nsamp_agec_sTV.allocate(1,nyr_agec_sTV,"nsamp_agec_sTV");
  nfish_agec_sTV.allocate(1,nyr_agec_sTV,"nfish_agec_sTV");
  obs_agec_sTV.allocate(1,nyr_agec_sTV,1,nages,"obs_agec_sTV");
  styr_cpue_cHL.allocate("styr_cpue_cHL");
  endyr_cpue_cHL.allocate("endyr_cpue_cHL");
  obs_cpue_cHL.allocate(styr_cpue_cHL,endyr_cpue_cHL,"obs_cpue_cHL");
  obs_cv_cpue_cHL.allocate(styr_cpue_cHL,endyr_cpue_cHL,"obs_cv_cpue_cHL");
  styr_L_cHL.allocate("styr_L_cHL");
  endyr_L_cHL.allocate("endyr_L_cHL");
  obs_L_cHL.allocate(styr_L_cHL,endyr_L_cHL,"obs_L_cHL");
  obs_cv_L_cHL.allocate(styr_L_cHL,endyr_L_cHL,"obs_cv_L_cHL");
  nyr_lenc_cHL.allocate("nyr_lenc_cHL");
  yrs_lenc_cHL.allocate(1,nyr_lenc_cHL,"yrs_lenc_cHL");
  nsamp_lenc_cHL.allocate(1,nyr_lenc_cHL,"nsamp_lenc_cHL");
  nfish_lenc_cHL.allocate(1,nyr_lenc_cHL,"nfish_lenc_cHL");
  obs_lenc_cHL.allocate(1,nyr_lenc_cHL,1,nlenbins,"obs_lenc_cHL");
  nyr_agec_cHL.allocate("nyr_agec_cHL");
  yrs_agec_cHL.allocate(1,nyr_agec_cHL,"yrs_agec_cHL");
  nsamp_agec_cHL.allocate(1,nyr_agec_cHL,"nsamp_agec_cHL");
  nfish_agec_cHL.allocate(1,nyr_agec_cHL,"nfish_agec_cHL");
  obs_agec_cHL.allocate(1,nyr_agec_cHL,1,nages,"obs_agec_cHL");
  styr_L_cPT.allocate("styr_L_cPT");
  endyr_L_cPT.allocate("endyr_L_cPT");
  obs_L_cPT.allocate(styr_L_cPT,endyr_L_cPT,"obs_L_cPT");
  obs_cv_L_cPT.allocate(styr_L_cPT,endyr_L_cPT,"obs_cv_L_cPT");
  nyr_lenc_cPT.allocate("nyr_lenc_cPT");
  yrs_lenc_cPT.allocate(1,nyr_lenc_cPT,"yrs_lenc_cPT");
  nsamp_lenc_cPT.allocate(1,nyr_lenc_cPT,"nsamp_lenc_cPT");
  nfish_lenc_cPT.allocate(1,nyr_lenc_cPT,"nfish_lenc_cPT");
  obs_lenc_cPT.allocate(1,nyr_lenc_cPT,1,nlenbins,"obs_lenc_cPT");
  nyr_agec_cPT.allocate("nyr_agec_cPT");
  yrs_agec_cPT.allocate(1,nyr_agec_cPT,"yrs_agec_cPT");
  nsamp_agec_cPT.allocate(1,nyr_agec_cPT,"nsamp_agec_cPT");
  nfish_agec_cPT.allocate(1,nyr_agec_cPT,"nfish_agec_cPT");
  obs_agec_cPT.allocate(1,nyr_agec_cPT,1,nages,"obs_agec_cPT");
  styr_D_cGD.allocate("styr_D_cGD");
  endyr_D_cGD.allocate("endyr_D_cGD");
  obs_released_cGD.allocate(styr_D_cGD,endyr_D_cGD,"obs_released_cGD");
  obs_cv_D_cGD.allocate(styr_D_cGD,endyr_D_cGD,"obs_cv_D_cGD");
  styr_L_cTW.allocate("styr_L_cTW");
  endyr_L_cTW.allocate("endyr_L_cTW");
  obs_L_cTW.allocate(styr_L_cTW,endyr_L_cTW,"obs_L_cTW");
  obs_cv_L_cTW.allocate(styr_L_cTW,endyr_L_cTW,"obs_cv_L_cTW");
  styr_cpue_rHB.allocate("styr_cpue_rHB");
  endyr_cpue_rHB.allocate("endyr_cpue_rHB");
  obs_cpue_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,"obs_cpue_rHB");
  obs_cv_cpue_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,"obs_cv_cpue_rHB");
  styr_L_rHB.allocate("styr_L_rHB");
  endyr_L_rHB.allocate("endyr_L_rHB");
  obs_L_rHB.allocate(styr_L_rHB,endyr_L_rHB,"obs_L_rHB");
  obs_cv_L_rHB.allocate(styr_L_rHB,endyr_L_rHB,"obs_cv_L_rHB");
  nyr_lenc_rHB.allocate("nyr_lenc_rHB");
  yrs_lenc_rHB.allocate(1,nyr_lenc_rHB,"yrs_lenc_rHB");
  nsamp_lenc_rHB.allocate(1,nyr_lenc_rHB,"nsamp_lenc_rHB");
  nfish_lenc_rHB.allocate(1,nyr_lenc_rHB,"nfish_lenc_rHB");
  obs_lenc_rHB.allocate(1,nyr_lenc_rHB,1,nlenbins,"obs_lenc_rHB");
  nyr_agec_rHB.allocate("nyr_agec_rHB");
  yrs_agec_rHB.allocate(1,nyr_agec_rHB,"yrs_agec_rHB");
  nsamp_agec_rHB.allocate(1,nyr_agec_rHB,"nsamp_agec_rHB");
  nfish_agec_rHB.allocate(1,nyr_agec_rHB,"nfish_agec_rHB");
  obs_agec_rHB.allocate(1,nyr_agec_rHB,1,nages,"obs_agec_rHB");
  styr_D_rHD.allocate("styr_D_rHD");
  endyr_D_rHD.allocate("endyr_D_rHD");
  obs_released_rHD.allocate(styr_D_rHD,endyr_D_rHD,"obs_released_rHD");
  obs_cv_D_rHD.allocate(styr_D_rHD,endyr_D_rHD,"obs_cv_D_rHD");
  nyr_lenc_rHD.allocate("nyr_lenc_rHD");
  yrs_lenc_rHD.allocate(1,nyr_lenc_rHD,"yrs_lenc_rHD");
  nsamp_lenc_rHD.allocate(1,nyr_lenc_rHD,"nsamp_lenc_rHD");
  nfish_lenc_rHD.allocate(1,nyr_lenc_rHD,"nfish_lenc_rHD");
  obs_lenc_rHD.allocate(1,nyr_lenc_rHD,1,nlenbins,"obs_lenc_rHD");
  styr_L_rGN.allocate("styr_L_rGN");
  endyr_L_rGN.allocate("endyr_L_rGN");
  obs_L_rGN.allocate(styr_L_rGN,endyr_L_rGN,"obs_L_rGN");
  obs_cv_L_rGN.allocate(styr_L_rGN,endyr_L_rGN,"obs_cv_L_rGN");
  nyr_lenc_rGN.allocate("nyr_lenc_rGN");
  yrs_lenc_rGN.allocate(1,nyr_lenc_rGN,"yrs_lenc_rGN");
  nsamp_lenc_rGN.allocate(1,nyr_lenc_rGN,"nsamp_lenc_rGN");
  nfish_lenc_rGN.allocate(1,nyr_lenc_rGN,"nfish_lenc_rGN");
  obs_lenc_rGN.allocate(1,nyr_lenc_rGN,1,nlenbins,"obs_lenc_rGN");
  styr_L_rIN.allocate("styr_L_rIN");
  endyr_L_rIN.allocate("endyr_L_rIN");
  obs_L_rIN.allocate(styr_L_rIN,endyr_L_rIN,"obs_L_rIN");
  obs_cv_L_rIN.allocate(styr_L_rIN,endyr_L_rIN,"obs_cv_L_rIN");
  styr_D_rGD.allocate("styr_D_rGD");
  endyr_D_rGD.allocate("endyr_D_rGD");
  obs_released_rGD.allocate(styr_D_rGD,endyr_D_rGD,"obs_released_rGD");
  obs_cv_D_rGD.allocate(styr_D_rGD,endyr_D_rGD,"obs_cv_D_rGD");
  styr_D_rID.allocate("styr_D_rID");
  endyr_D_rID.allocate("endyr_D_rID");
  obs_released_rID.allocate(styr_D_rID,endyr_D_rID,"obs_released_rID");
  obs_cv_D_rID.allocate(styr_D_rID,endyr_D_rID,"obs_cv_D_rID");
  set_Dmort_cGD.allocate("set_Dmort_cGD");
  set_Dmort_rHD.allocate("set_Dmort_rHD");
  set_Dmort_rGD.allocate("set_Dmort_rGD");
  set_Dmort_rID.allocate("set_Dmort_rID");
  set_Linf.allocate(1,7,"set_Linf");
  set_K.allocate(1,7,"set_K");
  set_t0.allocate(1,7,"set_t0");
  set_len_cv.allocate(1,7,"set_len_cv");
  set_M_constant.allocate(1,7,"set_M_constant");
  set_steep.allocate(1,7,"set_steep");
  set_log_R0.allocate(1,7,"set_log_R0");
  set_R_autocorr.allocate(1,7,"set_R_autocorr");
  set_rec_sigma.allocate(1,7,"set_rec_sigma");
  set_log_dm_lenc_sBT.allocate(1,7,"set_log_dm_lenc_sBT");
  set_log_dm_lenc_cHL.allocate(1,7,"set_log_dm_lenc_cHL");
  set_log_dm_lenc_cPT.allocate(1,7,"set_log_dm_lenc_cPT");
  set_log_dm_lenc_rHB.allocate(1,7,"set_log_dm_lenc_rHB");
  set_log_dm_lenc_rHD.allocate(1,7,"set_log_dm_lenc_rHD");
  set_log_dm_lenc_rGN.allocate(1,7,"set_log_dm_lenc_rGN");
  set_log_dm_agec_sBT.allocate(1,7,"set_log_dm_agec_sBT");
  set_log_dm_agec_sTV.allocate(1,7,"set_log_dm_agec_sTV");
  set_log_dm_agec_cHL.allocate(1,7,"set_log_dm_agec_cHL");
  set_log_dm_agec_cPT.allocate(1,7,"set_log_dm_agec_cPT");
  set_log_dm_agec_rHB.allocate(1,7,"set_log_dm_agec_rHB");
  set_selpar_A50_sBT.allocate(1,7,"set_selpar_A50_sBT");
  set_selpar_slope_sBT.allocate(1,7,"set_selpar_slope_sBT");
  set_selpar_A502_sBT.allocate(1,7,"set_selpar_A502_sBT");
  set_selpar_slope2_sBT.allocate(1,7,"set_selpar_slope2_sBT");
  set_selpar_A50_sTV.allocate(1,7,"set_selpar_A50_sTV");
  set_selpar_slope_sTV.allocate(1,7,"set_selpar_slope_sTV");
  set_selpar_A502_sTV.allocate(1,7,"set_selpar_A502_sTV");
  set_selpar_slope2_sTV.allocate(1,7,"set_selpar_slope2_sTV");
  set_selpar_A50b02_cHL.allocate(1,7,"set_selpar_A50b02_cHL");
  set_selpar_slopeb02_cHL.allocate(1,7,"set_selpar_slopeb02_cHL");
  set_selpar_A50b03_cHL.allocate(1,7,"set_selpar_A50b03_cHL");
  set_selpar_slopeb03_cHL.allocate(1,7,"set_selpar_slopeb03_cHL");
  set_selpar_A50b04_cHL.allocate(1,7,"set_selpar_A50b04_cHL");
  set_selpar_slopeb04_cHL.allocate(1,7,"set_selpar_slopeb04_cHL");
  set_selpar_A50b02_cPT.allocate(1,7,"set_selpar_A50b02_cPT");
  set_selpar_slopeb02_cPT.allocate(1,7,"set_selpar_slopeb02_cPT");
  set_selpar_A50b03_cPT.allocate(1,7,"set_selpar_A50b03_cPT");
  set_selpar_slopeb03_cPT.allocate(1,7,"set_selpar_slopeb03_cPT");
  set_selpar_A50b04_cPT.allocate(1,7,"set_selpar_A50b04_cPT");
  set_selpar_slopeb04_cPT.allocate(1,7,"set_selpar_slopeb04_cPT");
  set_selpar_A50b01_rHB.allocate(1,7,"set_selpar_A50b01_rHB");
  set_selpar_slopeb01_rHB.allocate(1,7,"set_selpar_slopeb01_rHB");
  set_selpar_A50b02_rHB.allocate(1,7,"set_selpar_A50b02_rHB");
  set_selpar_slopeb02_rHB.allocate(1,7,"set_selpar_slopeb02_rHB");
  set_selpar_A50b03_rHB.allocate(1,7,"set_selpar_A50b03_rHB");
  set_selpar_slopeb03_rHB.allocate(1,7,"set_selpar_slopeb03_rHB");
  set_selpar_A50b04_rHB.allocate(1,7,"set_selpar_A50b04_rHB");
  set_selpar_slopeb04_rHB.allocate(1,7,"set_selpar_slopeb04_rHB");
  set_selpar_A50b05_rHB.allocate(1,7,"set_selpar_A50b05_rHB");
  set_selpar_slopeb05_rHB.allocate(1,7,"set_selpar_slopeb05_rHB");
  set_selpar_logitAge0_rHD.allocate(1,7,"set_selpar_logitAge0_rHD");
  set_selpar_logitAge1_rHD.allocate(1,7,"set_selpar_logitAge1_rHD");
  set_selpar_logitAge2_rHD.allocate(1,7,"set_selpar_logitAge2_rHD");
  set_selpar_A50b04_rHD.allocate(1,7,"set_selpar_A50b04_rHD");
  set_selpar_slopeb04_rHD.allocate(1,7,"set_selpar_slopeb04_rHD");
  set_selpar_A502b04_rHD.allocate(1,7,"set_selpar_A502b04_rHD");
  set_selpar_slope2b04_rHD.allocate(1,7,"set_selpar_slope2b04_rHD");
  set_selpar_A50b05_rHD.allocate(1,7,"set_selpar_A50b05_rHD");
  set_selpar_slopeb05_rHD.allocate(1,7,"set_selpar_slopeb05_rHD");
  set_selpar_A502b05_rHD.allocate(1,7,"set_selpar_A502b05_rHD");
  set_selpar_slope2b05_rHD.allocate(1,7,"set_selpar_slope2b05_rHD");
  set_F_init_ratio.allocate(1,7,"set_F_init_ratio");
  set_log_q_cpue_sBT.allocate(1,7,"set_log_q_cpue_sBT");
  set_log_q_cpue_sTV.allocate(1,7,"set_log_q_cpue_sTV");
  set_log_q_cpue_cHL.allocate(1,7,"set_log_q_cpue_cHL");
  set_log_q_cpue_rHB.allocate(1,7,"set_log_q_cpue_rHB");
  set_q_rate.allocate(1,7,"set_q_rate");
  set_log_avg_F_L_cHL.allocate(1,7,"set_log_avg_F_L_cHL");
  set_log_avg_F_L_cPT.allocate(1,7,"set_log_avg_F_L_cPT");
  set_log_avg_F_L_cTW.allocate(1,7,"set_log_avg_F_L_cTW");
  set_log_avg_F_L_rHB.allocate(1,7,"set_log_avg_F_L_rHB");
  set_log_avg_F_L_rGN.allocate(1,7,"set_log_avg_F_L_rGN");
  set_log_avg_F_L_rIN.allocate(1,7,"set_log_avg_F_L_rIN");
  set_log_avg_F_D_cGD.allocate(1,7,"set_log_avg_F_D_cGD");
  set_log_avg_F_D_rHD.allocate(1,7,"set_log_avg_F_D_rHD");
  set_log_avg_F_D_rGD.allocate(1,7,"set_log_avg_F_D_rGD");
  set_log_avg_F_D_rID.allocate(1,7,"set_log_avg_F_D_rID");
  set_log_dev_F_L_cHL.allocate(1,3,"set_log_dev_F_L_cHL");
  set_log_dev_F_L_cPT.allocate(1,3,"set_log_dev_F_L_cPT");
  set_log_dev_F_L_cTW.allocate(1,3,"set_log_dev_F_L_cTW");
  set_log_dev_F_L_rHB.allocate(1,3,"set_log_dev_F_L_rHB");
  set_log_dev_F_L_rGN.allocate(1,3,"set_log_dev_F_L_rGN");
  set_log_dev_F_L_rIN.allocate(1,3,"set_log_dev_F_L_rIN");
  set_log_dev_F_D_cGD.allocate(1,3,"set_log_dev_F_D_cGD");
  set_log_dev_F_D_rHD.allocate(1,3,"set_log_dev_F_D_rHD");
  set_log_dev_F_D_rGD.allocate(1,3,"set_log_dev_F_D_rGD");
  set_log_dev_F_D_rID.allocate(1,3,"set_log_dev_F_D_rID");
  set_log_dev_RWq.allocate(1,3,"set_log_dev_RWq");
  set_log_dev_rec.allocate(1,3,"set_log_dev_rec");
  set_log_dev_Nage.allocate(1,3,"set_log_dev_Nage");
  set_log_dev_vals_F_L_cHL.allocate(styr_L_cHL,endyr_L_cHL,"set_log_dev_vals_F_L_cHL");
  set_log_dev_vals_F_L_cPT.allocate(styr_L_cPT,endyr_L_cPT,"set_log_dev_vals_F_L_cPT");
  set_log_dev_vals_F_L_cTW.allocate(styr_L_cTW,endyr_L_cTW,"set_log_dev_vals_F_L_cTW");
  set_log_dev_vals_F_L_rHB.allocate(styr_L_rHB,endyr_L_rHB,"set_log_dev_vals_F_L_rHB");
  set_log_dev_vals_F_L_rGN.allocate(styr_L_rGN,endyr_L_rGN,"set_log_dev_vals_F_L_rGN");
  set_log_dev_vals_F_L_rIN.allocate(styr_L_rGN,endyr_L_rIN,"set_log_dev_vals_F_L_rIN");
  set_log_dev_vals_F_D_cGD.allocate(styr_D_cGD,endyr_D_cGD,"set_log_dev_vals_F_D_cGD");
  set_log_dev_vals_F_D_rHD.allocate(styr_D_rHD,endyr_D_rHD,"set_log_dev_vals_F_D_rHD");
  set_log_dev_vals_F_D_rGD.allocate(styr_D_rGD,endyr_D_rGD,"set_log_dev_vals_F_D_rGD");
  set_log_dev_vals_F_D_rID.allocate(styr_D_rGD,endyr_D_rID,"set_log_dev_vals_F_D_rID");
  set_log_dev_vals_rec.allocate(styr_rec_dev,endyr_rec_dev,"set_log_dev_vals_rec");
  set_log_dev_vals_Nage.allocate(2,nages,"set_log_dev_vals_Nage");
  set_w_L.allocate("set_w_L");
  set_w_D.allocate("set_w_D");
  set_w_lenc_sBT.allocate("set_w_lenc_sBT");
  set_w_lenc_cHL.allocate("set_w_lenc_cHL");
  set_w_lenc_cPT.allocate("set_w_lenc_cPT");
  set_w_lenc_rHB.allocate("set_w_lenc_rHB");
  set_w_lenc_rHD.allocate("set_w_lenc_rHD");
  set_w_lenc_rGN.allocate("set_w_lenc_rGN");
  set_w_agec_sBT.allocate("set_w_agec_sBT");
  set_w_agec_sTV.allocate("set_w_agec_sTV");
  set_w_agec_cHL.allocate("set_w_agec_cHL");
  set_w_agec_cPT.allocate("set_w_agec_cPT");
  set_w_agec_rHB.allocate("set_w_agec_rHB");
  set_w_cpue_sBT.allocate("set_w_cpue_sBT");
  set_w_cpue_sTV.allocate("set_w_cpue_sTV");
  set_w_cpue_cHL.allocate("set_w_cpue_cHL");
  set_w_cpue_rHB.allocate("set_w_cpue_rHB");
  set_w_rec.allocate("set_w_rec");
  set_w_rec_early.allocate("set_w_rec_early");
  set_w_rec_end.allocate("set_w_rec_end");
  wgtpar_a.allocate("wgtpar_a");
  wgtpar_b.allocate("wgtpar_b");
  fecpar_a.allocate("fecpar_a");
  fecpar_b.allocate("fecpar_b");
  fecpar_batches.allocate("fecpar_batches");
  fecpar_scale.allocate("fecpar_scale");
  obs_maturity_f.allocate(1,nages,"obs_maturity_f");
  obs_prop_f.allocate(1,nages,"obs_prop_f");
  set_M.allocate(1,nages,"set_M");
  spawn_time_frac.allocate("spawn_time_frac");
  set_q_DD_phase.allocate("set_q_DD_phase");
  set_q_DD_beta.allocate("set_q_DD_beta");
  set_q_DD_beta_se.allocate("set_q_DD_beta_se");
  set_q_DD_stage.allocate("set_q_DD_stage");
  set_q_RW_phase.allocate("set_q_RW_phase");
  set_q_RW_var_cHL.allocate("set_q_RW_var_cHL");
  set_q_RW_var_rHB.allocate("set_q_RW_var_rHB");
  minSS_lenc.allocate("minSS_lenc");
  minSS_agec.allocate("minSS_agec");
  maxSS_lenc.allocate("maxSS_lenc");
  maxSS_agec.allocate("maxSS_agec");
  endyr_proj.allocate("endyr_proj");
  styr_regs.allocate("styr_regs");
  Fproj_switch.allocate("Fproj_switch");
  Fproj_mult.allocate("Fproj_mult");
  Rproj_switch.allocate("Rproj_switch");
   styr_proj=endyr+1;
  age_error.allocate(1,nages,1,nages,"age_error");
  end_of_data_file.allocate("end_of_data_file");
    if(end_of_data_file!=999)
    {
        cout << "*** WARNING: Data File NOT READ CORRECTLY ****" << endl;
        exit(0);  
    }
    else
    {cout << "Data File read correctly" << endl;} 
  if (global_datafile)
  {
    delete global_datafile;
    global_datafile = NULL;
  }
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
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
  
  Linf.allocate(Linf_LO,Linf_HI,Linf_PH,"Linf");
  K.allocate(K_LO,K_HI,K_PH,"K");
  t0.allocate(t0_LO,t0_HI,t0_PH,"t0");
  len_cv_val.allocate(len_cv_LO,len_cv_HI,len_cv_PH,"len_cv_val");
  Linf_out.allocate(1,8,"Linf_out");
  #ifndef NO_AD_INITIALIZE
    Linf_out.initialize();
  #endif
  K_out.allocate(1,8,"K_out");
  #ifndef NO_AD_INITIALIZE
    K_out.initialize();
  #endif
  t0_out.allocate(1,8,"t0_out");
  #ifndef NO_AD_INITIALIZE
    t0_out.initialize();
  #endif
  len_cv_val_out.allocate(1,8,"len_cv_val_out");
  #ifndef NO_AD_INITIALIZE
    len_cv_val_out.initialize();
  #endif
  meanlen_TL.allocate(1,nages,"meanlen_TL");
  #ifndef NO_AD_INITIALIZE
    meanlen_TL.initialize();
  #endif
  wgt_g.allocate(1,nages,"wgt_g");
  #ifndef NO_AD_INITIALIZE
    wgt_g.initialize();
  #endif
  wgt_kg.allocate(1,nages,"wgt_kg");
  #ifndef NO_AD_INITIALIZE
    wgt_kg.initialize();
  #endif
  wgt_mt.allocate(1,nages,"wgt_mt");
  #ifndef NO_AD_INITIALIZE
    wgt_mt.initialize();
  #endif
  wgt_klb.allocate(1,nages,"wgt_klb");
  #ifndef NO_AD_INITIALIZE
    wgt_klb.initialize();
  #endif
  wgt_lb.allocate(1,nages,"wgt_lb");
  #ifndef NO_AD_INITIALIZE
    wgt_lb.initialize();
  #endif
  fecundity.allocate(1,nages,"fecundity");
  #ifndef NO_AD_INITIALIZE
    fecundity.initialize();
  #endif
  len_cHL_mm.allocate(styr,endyr,1,nages,"len_cHL_mm");
  #ifndef NO_AD_INITIALIZE
    len_cHL_mm.initialize();
  #endif
  wgt_cHL_klb.allocate(styr,endyr,1,nages,"wgt_cHL_klb");
  #ifndef NO_AD_INITIALIZE
    wgt_cHL_klb.initialize();
  #endif
  len_cPT_mm.allocate(styr,endyr,1,nages,"len_cPT_mm");
  #ifndef NO_AD_INITIALIZE
    len_cPT_mm.initialize();
  #endif
  wgt_cPT_klb.allocate(styr,endyr,1,nages,"wgt_cPT_klb");
  #ifndef NO_AD_INITIALIZE
    wgt_cPT_klb.initialize();
  #endif
  len_cTW_mm.allocate(styr,endyr,1,nages,"len_cTW_mm");
  #ifndef NO_AD_INITIALIZE
    len_cTW_mm.initialize();
  #endif
  wgt_cTW_klb.allocate(styr,endyr,1,nages,"wgt_cTW_klb");
  #ifndef NO_AD_INITIALIZE
    wgt_cTW_klb.initialize();
  #endif
  len_rHB_mm.allocate(styr,endyr,1,nages,"len_rHB_mm");
  #ifndef NO_AD_INITIALIZE
    len_rHB_mm.initialize();
  #endif
  wgt_rHB_klb.allocate(styr,endyr,1,nages,"wgt_rHB_klb");
  #ifndef NO_AD_INITIALIZE
    wgt_rHB_klb.initialize();
  #endif
  len_rGN_mm.allocate(styr,endyr,1,nages,"len_rGN_mm");
  #ifndef NO_AD_INITIALIZE
    len_rGN_mm.initialize();
  #endif
  wgt_rGN_klb.allocate(styr,endyr,1,nages,"wgt_rGN_klb");
  #ifndef NO_AD_INITIALIZE
    wgt_rGN_klb.initialize();
  #endif
  len_rIN_mm.allocate(styr,endyr,1,nages,"len_rIN_mm");
  #ifndef NO_AD_INITIALIZE
    len_rIN_mm.initialize();
  #endif
  wgt_rIN_klb.allocate(styr,endyr,1,nages,"wgt_rIN_klb");
  #ifndef NO_AD_INITIALIZE
    wgt_rIN_klb.initialize();
  #endif
  len_cGD_mm.allocate(styr,endyr,1,nages,"len_cGD_mm");
  #ifndef NO_AD_INITIALIZE
    len_cGD_mm.initialize();
  #endif
  wgt_cGD_klb.allocate(styr,endyr,1,nages,"wgt_cGD_klb");
  #ifndef NO_AD_INITIALIZE
    wgt_cGD_klb.initialize();
  #endif
  len_rHD_mm.allocate(styr,endyr,1,nages,"len_rHD_mm");
  #ifndef NO_AD_INITIALIZE
    len_rHD_mm.initialize();
  #endif
  wgt_rHD_klb.allocate(styr,endyr,1,nages,"wgt_rHD_klb");
  #ifndef NO_AD_INITIALIZE
    wgt_rHD_klb.initialize();
  #endif
  len_rGD_mm.allocate(styr,endyr,1,nages,"len_rGD_mm");
  #ifndef NO_AD_INITIALIZE
    len_rGD_mm.initialize();
  #endif
  wgt_rGD_klb.allocate(styr,endyr,1,nages,"wgt_rGD_klb");
  #ifndef NO_AD_INITIALIZE
    wgt_rGD_klb.initialize();
  #endif
  len_rID_mm.allocate(styr,endyr,1,nages,"len_rID_mm");
  #ifndef NO_AD_INITIALIZE
    len_rID_mm.initialize();
  #endif
  wgt_rID_klb.allocate(styr,endyr,1,nages,"wgt_rID_klb");
  #ifndef NO_AD_INITIALIZE
    wgt_rID_klb.initialize();
  #endif
  lenprob.allocate(1,nages,1,nlenbins,"lenprob");
  #ifndef NO_AD_INITIALIZE
    lenprob.initialize();
  #endif
  zscore_len.allocate("zscore_len");
  #ifndef NO_AD_INITIALIZE
  zscore_len.initialize();
  #endif
  cprob_lenvec.allocate(1,nlenbins,"cprob_lenvec");
  #ifndef NO_AD_INITIALIZE
    cprob_lenvec.initialize();
  #endif
  zscore_lzero.allocate("zscore_lzero");
  #ifndef NO_AD_INITIALIZE
  zscore_lzero.initialize();
  #endif
  cprob_lzero.allocate("cprob_lzero");
  #ifndef NO_AD_INITIALIZE
  cprob_lzero.initialize();
  #endif
  len_sd.allocate(1,nages,"len_sd");
  #ifndef NO_AD_INITIALIZE
    len_sd.initialize();
  #endif
  len_cv.allocate(1,nages,"len_cv");
  #ifndef NO_AD_INITIALIZE
    len_cv.initialize();
  #endif
  pred_lenc_sBT.allocate(1,nyr_lenc_sBT,1,nlenbins,"pred_lenc_sBT");
  #ifndef NO_AD_INITIALIZE
    pred_lenc_sBT.initialize();
  #endif
  pred_lenc_cHL.allocate(1,nyr_lenc_cHL,1,nlenbins,"pred_lenc_cHL");
  #ifndef NO_AD_INITIALIZE
    pred_lenc_cHL.initialize();
  #endif
  pred_lenc_cPT.allocate(1,nyr_lenc_cPT,1,nlenbins,"pred_lenc_cPT");
  #ifndef NO_AD_INITIALIZE
    pred_lenc_cPT.initialize();
  #endif
  pred_lenc_rHB.allocate(1,nyr_lenc_rHB,1,nlenbins,"pred_lenc_rHB");
  #ifndef NO_AD_INITIALIZE
    pred_lenc_rHB.initialize();
  #endif
  pred_lenc_rHD.allocate(1,nyr_lenc_rHD,1,nlenbins,"pred_lenc_rHD");
  #ifndef NO_AD_INITIALIZE
    pred_lenc_rHD.initialize();
  #endif
  pred_lenc_rGN.allocate(1,nyr_lenc_rGN,1,nlenbins,"pred_lenc_rGN");
  #ifndef NO_AD_INITIALIZE
    pred_lenc_rGN.initialize();
  #endif
  predNoErr_agec_sBT.allocate(1,nyr_agec_sBT,1,nages,"predNoErr_agec_sBT");
  #ifndef NO_AD_INITIALIZE
    predNoErr_agec_sBT.initialize();
  #endif
  predNoErr_agec_sTV.allocate(1,nyr_agec_sTV,1,nages,"predNoErr_agec_sTV");
  #ifndef NO_AD_INITIALIZE
    predNoErr_agec_sTV.initialize();
  #endif
  predNoErr_agec_cHL.allocate(1,nyr_agec_cHL,1,nages,"predNoErr_agec_cHL");
  #ifndef NO_AD_INITIALIZE
    predNoErr_agec_cHL.initialize();
  #endif
  predNoErr_agec_cPT.allocate(1,nyr_agec_cPT,1,nages,"predNoErr_agec_cPT");
  #ifndef NO_AD_INITIALIZE
    predNoErr_agec_cPT.initialize();
  #endif
  predNoErr_agec_rHB.allocate(1,nyr_agec_rHB,1,nages,"predNoErr_agec_rHB");
  #ifndef NO_AD_INITIALIZE
    predNoErr_agec_rHB.initialize();
  #endif
  pred_agec_sBT.allocate(1,nyr_agec_sBT,1,nages,"pred_agec_sBT");
  #ifndef NO_AD_INITIALIZE
    pred_agec_sBT.initialize();
  #endif
  pred_agec_sTV.allocate(1,nyr_agec_sTV,1,nages,"pred_agec_sTV");
  #ifndef NO_AD_INITIALIZE
    pred_agec_sTV.initialize();
  #endif
  pred_agec_cHL.allocate(1,nyr_agec_cHL,1,nages,"pred_agec_cHL");
  #ifndef NO_AD_INITIALIZE
    pred_agec_cHL.initialize();
  #endif
  pred_agec_cPT.allocate(1,nyr_agec_cPT,1,nages,"pred_agec_cPT");
  #ifndef NO_AD_INITIALIZE
    pred_agec_cPT.initialize();
  #endif
  pred_agec_rHB.allocate(1,nyr_agec_rHB,1,nages,"pred_agec_rHB");
  #ifndef NO_AD_INITIALIZE
    pred_agec_rHB.initialize();
  #endif
  nsamp_allyr_lenc_sBT.allocate(styr,endyr,"nsamp_allyr_lenc_sBT");
  #ifndef NO_AD_INITIALIZE
    nsamp_allyr_lenc_sBT.initialize();
  #endif
  nsamp_allyr_lenc_cHL.allocate(styr,endyr,"nsamp_allyr_lenc_cHL");
  #ifndef NO_AD_INITIALIZE
    nsamp_allyr_lenc_cHL.initialize();
  #endif
  nsamp_allyr_lenc_cPT.allocate(styr,endyr,"nsamp_allyr_lenc_cPT");
  #ifndef NO_AD_INITIALIZE
    nsamp_allyr_lenc_cPT.initialize();
  #endif
  nsamp_allyr_lenc_rHB.allocate(styr,endyr,"nsamp_allyr_lenc_rHB");
  #ifndef NO_AD_INITIALIZE
    nsamp_allyr_lenc_rHB.initialize();
  #endif
  nsamp_allyr_lenc_rHD.allocate(styr,endyr,"nsamp_allyr_lenc_rHD");
  #ifndef NO_AD_INITIALIZE
    nsamp_allyr_lenc_rHD.initialize();
  #endif
  nsamp_allyr_lenc_rGN.allocate(styr,endyr,"nsamp_allyr_lenc_rGN");
  #ifndef NO_AD_INITIALIZE
    nsamp_allyr_lenc_rGN.initialize();
  #endif
  nsamp_allyr_agec_sBT.allocate(styr,endyr,"nsamp_allyr_agec_sBT");
  #ifndef NO_AD_INITIALIZE
    nsamp_allyr_agec_sBT.initialize();
  #endif
  nsamp_allyr_agec_sTV.allocate(styr,endyr,"nsamp_allyr_agec_sTV");
  #ifndef NO_AD_INITIALIZE
    nsamp_allyr_agec_sTV.initialize();
  #endif
  nsamp_allyr_agec_cHL.allocate(styr,endyr,"nsamp_allyr_agec_cHL");
  #ifndef NO_AD_INITIALIZE
    nsamp_allyr_agec_cHL.initialize();
  #endif
  nsamp_allyr_agec_cPT.allocate(styr,endyr,"nsamp_allyr_agec_cPT");
  #ifndef NO_AD_INITIALIZE
    nsamp_allyr_agec_cPT.initialize();
  #endif
  nsamp_allyr_agec_rHB.allocate(styr,endyr,"nsamp_allyr_agec_rHB");
  #ifndef NO_AD_INITIALIZE
    nsamp_allyr_agec_rHB.initialize();
  #endif
  nfish_allyr_lenc_sBT.allocate(styr,endyr,"nfish_allyr_lenc_sBT");
  #ifndef NO_AD_INITIALIZE
    nfish_allyr_lenc_sBT.initialize();
  #endif
  nfish_allyr_lenc_cHL.allocate(styr,endyr,"nfish_allyr_lenc_cHL");
  #ifndef NO_AD_INITIALIZE
    nfish_allyr_lenc_cHL.initialize();
  #endif
  nfish_allyr_lenc_cPT.allocate(styr,endyr,"nfish_allyr_lenc_cPT");
  #ifndef NO_AD_INITIALIZE
    nfish_allyr_lenc_cPT.initialize();
  #endif
  nfish_allyr_lenc_rHB.allocate(styr,endyr,"nfish_allyr_lenc_rHB");
  #ifndef NO_AD_INITIALIZE
    nfish_allyr_lenc_rHB.initialize();
  #endif
  nfish_allyr_lenc_rHD.allocate(styr,endyr,"nfish_allyr_lenc_rHD");
  #ifndef NO_AD_INITIALIZE
    nfish_allyr_lenc_rHD.initialize();
  #endif
  nfish_allyr_lenc_rGN.allocate(styr,endyr,"nfish_allyr_lenc_rGN");
  #ifndef NO_AD_INITIALIZE
    nfish_allyr_lenc_rGN.initialize();
  #endif
  nfish_allyr_agec_sBT.allocate(styr,endyr,"nfish_allyr_agec_sBT");
  #ifndef NO_AD_INITIALIZE
    nfish_allyr_agec_sBT.initialize();
  #endif
  nfish_allyr_agec_sTV.allocate(styr,endyr,"nfish_allyr_agec_sTV");
  #ifndef NO_AD_INITIALIZE
    nfish_allyr_agec_sTV.initialize();
  #endif
  nfish_allyr_agec_cHL.allocate(styr,endyr,"nfish_allyr_agec_cHL");
  #ifndef NO_AD_INITIALIZE
    nfish_allyr_agec_cHL.initialize();
  #endif
  nfish_allyr_agec_cPT.allocate(styr,endyr,"nfish_allyr_agec_cPT");
  #ifndef NO_AD_INITIALIZE
    nfish_allyr_agec_cPT.initialize();
  #endif
  nfish_allyr_agec_rHB.allocate(styr,endyr,"nfish_allyr_agec_rHB");
  #ifndef NO_AD_INITIALIZE
    nfish_allyr_agec_rHB.initialize();
  #endif
  out_neff_allyr_lenc_sBT.allocate(styr,endyr,"out_neff_allyr_lenc_sBT");
  #ifndef NO_AD_INITIALIZE
    out_neff_allyr_lenc_sBT.initialize();
  #endif
  out_neff_allyr_lenc_cHL.allocate(styr,endyr,"out_neff_allyr_lenc_cHL");
  #ifndef NO_AD_INITIALIZE
    out_neff_allyr_lenc_cHL.initialize();
  #endif
  out_neff_allyr_lenc_cPT.allocate(styr,endyr,"out_neff_allyr_lenc_cPT");
  #ifndef NO_AD_INITIALIZE
    out_neff_allyr_lenc_cPT.initialize();
  #endif
  out_neff_allyr_lenc_rHB.allocate(styr,endyr,"out_neff_allyr_lenc_rHB");
  #ifndef NO_AD_INITIALIZE
    out_neff_allyr_lenc_rHB.initialize();
  #endif
  out_neff_allyr_lenc_rHD.allocate(styr,endyr,"out_neff_allyr_lenc_rHD");
  #ifndef NO_AD_INITIALIZE
    out_neff_allyr_lenc_rHD.initialize();
  #endif
  out_neff_allyr_lenc_rGN.allocate(styr,endyr,"out_neff_allyr_lenc_rGN");
  #ifndef NO_AD_INITIALIZE
    out_neff_allyr_lenc_rGN.initialize();
  #endif
  out_neff_allyr_agec_sBT.allocate(styr,endyr,"out_neff_allyr_agec_sBT");
  #ifndef NO_AD_INITIALIZE
    out_neff_allyr_agec_sBT.initialize();
  #endif
  out_neff_allyr_agec_sTV.allocate(styr,endyr,"out_neff_allyr_agec_sTV");
  #ifndef NO_AD_INITIALIZE
    out_neff_allyr_agec_sTV.initialize();
  #endif
  out_neff_allyr_agec_cHL.allocate(styr,endyr,"out_neff_allyr_agec_cHL");
  #ifndef NO_AD_INITIALIZE
    out_neff_allyr_agec_cHL.initialize();
  #endif
  out_neff_allyr_agec_cPT.allocate(styr,endyr,"out_neff_allyr_agec_cPT");
  #ifndef NO_AD_INITIALIZE
    out_neff_allyr_agec_cPT.initialize();
  #endif
  out_neff_allyr_agec_rHB.allocate(styr,endyr,"out_neff_allyr_agec_rHB");
  #ifndef NO_AD_INITIALIZE
    out_neff_allyr_agec_rHB.initialize();
  #endif
  N.allocate(styr,endyr,1,nages,"N");
  #ifndef NO_AD_INITIALIZE
    N.initialize();
  #endif
  N_mdyr.allocate(styr,endyr,1,nages,"N_mdyr");
  #ifndef NO_AD_INITIALIZE
    N_mdyr.initialize();
  #endif
  N_spawn.allocate(styr,endyr,1,nages,"N_spawn");
  #ifndef NO_AD_INITIALIZE
    N_spawn.initialize();
  #endif
  log_dev_Nage.allocate(2,nages,log_Nage_dev_LO,log_Nage_dev_HI,log_Nage_dev_PH,"log_dev_Nage");
  log_Nage_dev_output.allocate(1,nages,"log_Nage_dev_output");
  #ifndef NO_AD_INITIALIZE
    log_Nage_dev_output.initialize();
  #endif
  B.allocate(styr,endyr,1,nages,"B");
  #ifndef NO_AD_INITIALIZE
    B.initialize();
  #endif
  totB.allocate(styr,endyr,"totB");
  #ifndef NO_AD_INITIALIZE
    totB.initialize();
  #endif
  totN.allocate(styr,endyr,"totN");
  #ifndef NO_AD_INITIALIZE
    totN.initialize();
  #endif
  SSB.allocate(styr,endyr,"SSB");
  #ifndef NO_AD_INITIALIZE
    SSB.initialize();
  #endif
  MatFemB.allocate(styr,endyr,"MatFemB");
  #ifndef NO_AD_INITIALIZE
    MatFemB.initialize();
  #endif
  rec.allocate(styr,endyr,"rec");
  #ifndef NO_AD_INITIALIZE
    rec.initialize();
  #endif
  prop_f.allocate(1,nages,"prop_f");
  #ifndef NO_AD_INITIALIZE
    prop_f.initialize();
  #endif
  maturity_f.allocate(1,nages,"maturity_f");
  #ifndef NO_AD_INITIALIZE
    maturity_f.initialize();
  #endif
  reprod.allocate(1,nages,"reprod");
  #ifndef NO_AD_INITIALIZE
    reprod.initialize();
  #endif
  reprod2.allocate(1,nages,"reprod2");
  #ifndef NO_AD_INITIALIZE
    reprod2.initialize();
  #endif
  log_R0.allocate(log_R0_LO,log_R0_HI,log_R0_PH,"log_R0");
  log_R0_out.allocate(1,8,"log_R0_out");
  #ifndef NO_AD_INITIALIZE
    log_R0_out.initialize();
  #endif
  R0.allocate("R0");
  #ifndef NO_AD_INITIALIZE
  R0.initialize();
  #endif
  steep.allocate(steep_LO,steep_HI,steep_PH,"steep");
  steep_out.allocate(1,8,"steep_out");
  #ifndef NO_AD_INITIALIZE
    steep_out.initialize();
  #endif
  rec_sigma.allocate(rec_sigma_LO,rec_sigma_HI,rec_sigma_PH,"rec_sigma");
  rec_sigma_out.allocate(1,8,"rec_sigma_out");
  #ifndef NO_AD_INITIALIZE
    rec_sigma_out.initialize();
  #endif
  R_autocorr.allocate(R_autocorr_LO,R_autocorr_HI,R_autocorr_PH,"R_autocorr");
  R_autocorr_out.allocate(1,8,"R_autocorr_out");
  #ifndef NO_AD_INITIALIZE
    R_autocorr_out.initialize();
  #endif
  rec_sigma_sq.allocate("rec_sigma_sq");
  #ifndef NO_AD_INITIALIZE
  rec_sigma_sq.initialize();
  #endif
  rec_sigma_sqd2.allocate("rec_sigma_sqd2");
  #ifndef NO_AD_INITIALIZE
  rec_sigma_sqd2.initialize();
  #endif
  rec_logL_add.allocate("rec_logL_add");
  #ifndef NO_AD_INITIALIZE
  rec_logL_add.initialize();
  #endif
  log_dev_rec.allocate(styr_rec_dev,endyr_rec_dev,log_rec_dev_LO,log_rec_dev_HI,log_rec_dev_PH,"log_dev_rec");
  log_rec_dev_output.allocate(styr_rec_dev,endyr,"log_rec_dev_output");
  #ifndef NO_AD_INITIALIZE
    log_rec_dev_output.initialize();
  #endif
  log_rec_dev_out.allocate(styr_rec_dev,endyr_rec_dev,"log_rec_dev_out");
  #ifndef NO_AD_INITIALIZE
    log_rec_dev_out.initialize();
  #endif
  var_rec_dev.allocate("var_rec_dev");
  #ifndef NO_AD_INITIALIZE
  var_rec_dev.initialize();
  #endif
  sigma_rec_dev.allocate("sigma_rec_dev");
  #ifndef NO_AD_INITIALIZE
  sigma_rec_dev.initialize();
  #endif
  BiasCor.allocate("BiasCor");
  #ifndef NO_AD_INITIALIZE
  BiasCor.initialize();
  #endif
  S0.allocate("S0");
  #ifndef NO_AD_INITIALIZE
  S0.initialize();
  #endif
  B0.allocate("B0");
  #ifndef NO_AD_INITIALIZE
  B0.initialize();
  #endif
  R1.allocate("R1");
  #ifndef NO_AD_INITIALIZE
  R1.initialize();
  #endif
  R_virgin.allocate("R_virgin");
  #ifndef NO_AD_INITIALIZE
  R_virgin.initialize();
  #endif
  SdS0.allocate(styr,endyr,"SdS0");
  #ifndef NO_AD_INITIALIZE
    SdS0.initialize();
  #endif
  log_dm_lenc_sBT.allocate(log_dm_sBT_lc_LO,log_dm_sBT_lc_HI,log_dm_sBT_lc_PH,"log_dm_lenc_sBT");
  log_dm_lenc_cHL.allocate(log_dm_cHL_lc_LO,log_dm_cHL_lc_HI,log_dm_cHL_lc_PH,"log_dm_lenc_cHL");
  log_dm_lenc_cPT.allocate(log_dm_cPT_lc_LO,log_dm_cPT_lc_HI,log_dm_cPT_lc_PH,"log_dm_lenc_cPT");
  log_dm_lenc_rHB.allocate(log_dm_rHB_lc_LO,log_dm_rHB_lc_HI,log_dm_rHB_lc_PH,"log_dm_lenc_rHB");
  log_dm_lenc_rHD.allocate(log_dm_rHD_lc_LO,log_dm_rHD_lc_HI,log_dm_rHD_lc_PH,"log_dm_lenc_rHD");
  log_dm_lenc_rGN.allocate(log_dm_rGN_lc_LO,log_dm_rGN_lc_HI,log_dm_rGN_lc_PH,"log_dm_lenc_rGN");
  log_dm_agec_sBT.allocate(log_dm_sBT_ac_LO,log_dm_sBT_ac_HI,log_dm_sBT_ac_PH,"log_dm_agec_sBT");
  log_dm_agec_sTV.allocate(log_dm_sTV_ac_LO,log_dm_sTV_ac_HI,log_dm_sTV_ac_PH,"log_dm_agec_sTV");
  log_dm_agec_cHL.allocate(log_dm_cHL_ac_LO,log_dm_cHL_ac_HI,log_dm_cHL_ac_PH,"log_dm_agec_cHL");
  log_dm_agec_cPT.allocate(log_dm_cPT_ac_LO,log_dm_cPT_ac_HI,log_dm_cPT_ac_PH,"log_dm_agec_cPT");
  log_dm_agec_rHB.allocate(log_dm_rHB_ac_LO,log_dm_rHB_ac_HI,log_dm_rHB_ac_PH,"log_dm_agec_rHB");
  log_dm_sBT_lc_out.allocate(1,8,"log_dm_sBT_lc_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_sBT_lc_out.initialize();
  #endif
  log_dm_cHL_lc_out.allocate(1,8,"log_dm_cHL_lc_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_cHL_lc_out.initialize();
  #endif
  log_dm_cPT_lc_out.allocate(1,8,"log_dm_cPT_lc_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_cPT_lc_out.initialize();
  #endif
  log_dm_rHB_lc_out.allocate(1,8,"log_dm_rHB_lc_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_rHB_lc_out.initialize();
  #endif
  log_dm_rGN_lc_out.allocate(1,8,"log_dm_rGN_lc_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_rGN_lc_out.initialize();
  #endif
  log_dm_rHD_lc_out.allocate(1,8,"log_dm_rHD_lc_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_rHD_lc_out.initialize();
  #endif
  log_dm_sBT_ac_out.allocate(1,8,"log_dm_sBT_ac_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_sBT_ac_out.initialize();
  #endif
  log_dm_sTV_ac_out.allocate(1,8,"log_dm_sTV_ac_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_sTV_ac_out.initialize();
  #endif
  log_dm_cHL_ac_out.allocate(1,8,"log_dm_cHL_ac_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_cHL_ac_out.initialize();
  #endif
  log_dm_cPT_ac_out.allocate(1,8,"log_dm_cPT_ac_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_cPT_ac_out.initialize();
  #endif
  log_dm_rHB_ac_out.allocate(1,8,"log_dm_rHB_ac_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_rHB_ac_out.initialize();
  #endif
  sel_sBT.allocate(styr,endyr,1,nages,"sel_sBT");
  #ifndef NO_AD_INITIALIZE
    sel_sBT.initialize();
  #endif
  sel_sBT_b01.allocate(1,nages,"sel_sBT_b01");
  #ifndef NO_AD_INITIALIZE
    sel_sBT_b01.initialize();
  #endif
  sel_sBT_b02.allocate(1,nages,"sel_sBT_b02");
  #ifndef NO_AD_INITIALIZE
    sel_sBT_b02.initialize();
  #endif
  sel_sBT_b03.allocate(1,nages,"sel_sBT_b03");
  #ifndef NO_AD_INITIALIZE
    sel_sBT_b03.initialize();
  #endif
  sel_sBT_b04.allocate(1,nages,"sel_sBT_b04");
  #ifndef NO_AD_INITIALIZE
    sel_sBT_b04.initialize();
  #endif
  sel_sBT_b05.allocate(1,nages,"sel_sBT_b05");
  #ifndef NO_AD_INITIALIZE
    sel_sBT_b05.initialize();
  #endif
  selpar_A50_sBT.allocate(selpar_A50_sBT_LO,selpar_A50_sBT_HI,selpar_A50_sBT_PH,"selpar_A50_sBT");
  selpar_slope_sBT.allocate(selpar_slope_sBT_LO,selpar_slope_sBT_HI,selpar_slope_sBT_PH,"selpar_slope_sBT");
  selpar_A502_sBT.allocate(selpar_A502_sBT_LO,selpar_A502_sBT_HI,selpar_A502_sBT_PH,"selpar_A502_sBT");
  selpar_slope2_sBT.allocate(selpar_slope2_sBT_LO,selpar_slope2_sBT_HI,selpar_slope2_sBT_PH,"selpar_slope2_sBT");
  selpar_A50_sBT_out.allocate(1,8,"selpar_A50_sBT_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_sBT_out.initialize();
  #endif
  selpar_slope_sBT_out.allocate(1,8,"selpar_slope_sBT_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_sBT_out.initialize();
  #endif
  selpar_A502_sBT_out.allocate(1,8,"selpar_A502_sBT_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502_sBT_out.initialize();
  #endif
  selpar_slope2_sBT_out.allocate(1,8,"selpar_slope2_sBT_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_sBT_out.initialize();
  #endif
  sel_sTV.allocate(styr,endyr,1,nages,"sel_sTV");
  #ifndef NO_AD_INITIALIZE
    sel_sTV.initialize();
  #endif
  sel_sTV_b01.allocate(1,nages,"sel_sTV_b01");
  #ifndef NO_AD_INITIALIZE
    sel_sTV_b01.initialize();
  #endif
  sel_sTV_b02.allocate(1,nages,"sel_sTV_b02");
  #ifndef NO_AD_INITIALIZE
    sel_sTV_b02.initialize();
  #endif
  sel_sTV_b03.allocate(1,nages,"sel_sTV_b03");
  #ifndef NO_AD_INITIALIZE
    sel_sTV_b03.initialize();
  #endif
  sel_sTV_b04.allocate(1,nages,"sel_sTV_b04");
  #ifndef NO_AD_INITIALIZE
    sel_sTV_b04.initialize();
  #endif
  sel_sTV_b05.allocate(1,nages,"sel_sTV_b05");
  #ifndef NO_AD_INITIALIZE
    sel_sTV_b05.initialize();
  #endif
  selpar_A50_sTV.allocate(selpar_A50_sTV_LO,selpar_A50_sTV_HI,selpar_A50_sTV_PH,"selpar_A50_sTV");
  selpar_slope_sTV.allocate(selpar_slope_sTV_LO,selpar_slope_sTV_HI,selpar_slope_sTV_PH,"selpar_slope_sTV");
  selpar_A502_sTV.allocate(selpar_A502_sTV_LO,selpar_A502_sTV_HI,selpar_A502_sTV_PH,"selpar_A502_sTV");
  selpar_slope2_sTV.allocate(selpar_slope2_sTV_LO,selpar_slope2_sTV_HI,selpar_slope2_sTV_PH,"selpar_slope2_sTV");
  selpar_A50_sTV_out.allocate(1,8,"selpar_A50_sTV_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_sTV_out.initialize();
  #endif
  selpar_slope_sTV_out.allocate(1,8,"selpar_slope_sTV_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_sTV_out.initialize();
  #endif
  selpar_A502_sTV_out.allocate(1,8,"selpar_A502_sTV_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502_sTV_out.initialize();
  #endif
  selpar_slope2_sTV_out.allocate(1,8,"selpar_slope2_sTV_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_sTV_out.initialize();
  #endif
  sel_cHL.allocate(styr,endyr,1,nages,"sel_cHL");
  #ifndef NO_AD_INITIALIZE
    sel_cHL.initialize();
  #endif
  sel_cHL_b01.allocate(1,nages,"sel_cHL_b01");
  #ifndef NO_AD_INITIALIZE
    sel_cHL_b01.initialize();
  #endif
  sel_cHL_b02.allocate(1,nages,"sel_cHL_b02");
  #ifndef NO_AD_INITIALIZE
    sel_cHL_b02.initialize();
  #endif
  sel_cHL_b03.allocate(1,nages,"sel_cHL_b03");
  #ifndef NO_AD_INITIALIZE
    sel_cHL_b03.initialize();
  #endif
  sel_cHL_b04.allocate(1,nages,"sel_cHL_b04");
  #ifndef NO_AD_INITIALIZE
    sel_cHL_b04.initialize();
  #endif
  sel_cHL_b05.allocate(1,nages,"sel_cHL_b05");
  #ifndef NO_AD_INITIALIZE
    sel_cHL_b05.initialize();
  #endif
  selpar_A50b02_cHL.allocate(selpar_A50b02_cHL_LO,selpar_A50b02_cHL_HI,selpar_A50b02_cHL_PH,"selpar_A50b02_cHL");
  selpar_slopeb02_cHL.allocate(selpar_slopeb02_cHL_LO,selpar_slopeb02_cHL_HI,selpar_slopeb02_cHL_PH,"selpar_slopeb02_cHL");
  selpar_A50b03_cHL.allocate(selpar_A50b03_cHL_LO,selpar_A50b03_cHL_HI,selpar_A50b03_cHL_PH,"selpar_A50b03_cHL");
  selpar_slopeb03_cHL.allocate(selpar_slopeb03_cHL_LO,selpar_slopeb03_cHL_HI,selpar_slopeb03_cHL_PH,"selpar_slopeb03_cHL");
  selpar_A50b04_cHL.allocate(selpar_A50b04_cHL_LO,selpar_A50b04_cHL_HI,selpar_A50b04_cHL_PH,"selpar_A50b04_cHL");
  selpar_slopeb04_cHL.allocate(selpar_slopeb04_cHL_LO,selpar_slopeb04_cHL_HI,selpar_slopeb04_cHL_PH,"selpar_slopeb04_cHL");
  selpar_A50b02_cHL_out.allocate(1,8,"selpar_A50b02_cHL_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50b02_cHL_out.initialize();
  #endif
  selpar_slopeb02_cHL_out.allocate(1,8,"selpar_slopeb02_cHL_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slopeb02_cHL_out.initialize();
  #endif
  selpar_A50b03_cHL_out.allocate(1,8,"selpar_A50b03_cHL_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50b03_cHL_out.initialize();
  #endif
  selpar_slopeb03_cHL_out.allocate(1,8,"selpar_slopeb03_cHL_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slopeb03_cHL_out.initialize();
  #endif
  selpar_A50b04_cHL_out.allocate(1,8,"selpar_A50b04_cHL_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50b04_cHL_out.initialize();
  #endif
  selpar_slopeb04_cHL_out.allocate(1,8,"selpar_slopeb04_cHL_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slopeb04_cHL_out.initialize();
  #endif
  sel_cGD.allocate(styr,endyr,1,nages,"sel_cGD");
  #ifndef NO_AD_INITIALIZE
    sel_cGD.initialize();
  #endif
  sel_cGD_b01.allocate(1,nages,"sel_cGD_b01");
  #ifndef NO_AD_INITIALIZE
    sel_cGD_b01.initialize();
  #endif
  sel_cGD_b02.allocate(1,nages,"sel_cGD_b02");
  #ifndef NO_AD_INITIALIZE
    sel_cGD_b02.initialize();
  #endif
  sel_cGD_b03.allocate(1,nages,"sel_cGD_b03");
  #ifndef NO_AD_INITIALIZE
    sel_cGD_b03.initialize();
  #endif
  sel_cGD_b04.allocate(1,nages,"sel_cGD_b04");
  #ifndef NO_AD_INITIALIZE
    sel_cGD_b04.initialize();
  #endif
  sel_cGD_b05.allocate(1,nages,"sel_cGD_b05");
  #ifndef NO_AD_INITIALIZE
    sel_cGD_b05.initialize();
  #endif
  Dopen_cHL4.allocate("Dopen_cHL4");
  #ifndef NO_AD_INITIALIZE
  Dopen_cHL4.initialize();
  #endif
  Dclosed_cHL4.allocate("Dclosed_cHL4");
  #ifndef NO_AD_INITIALIZE
  Dclosed_cHL4.initialize();
  #endif
  Lopen_cHL4.allocate("Lopen_cHL4");
  #ifndef NO_AD_INITIALIZE
  Lopen_cHL4.initialize();
  #endif
  Dopen_cPT4.allocate("Dopen_cPT4");
  #ifndef NO_AD_INITIALIZE
  Dopen_cPT4.initialize();
  #endif
  Dclosed_cPT4.allocate("Dclosed_cPT4");
  #ifndef NO_AD_INITIALIZE
  Dclosed_cPT4.initialize();
  #endif
  Lopen_cPT4.allocate("Lopen_cPT4");
  #ifndef NO_AD_INITIALIZE
  Lopen_cPT4.initialize();
  #endif
  D_sum_cHLcPT4.allocate("D_sum_cHLcPT4");
  #ifndef NO_AD_INITIALIZE
  D_sum_cHLcPT4.initialize();
  #endif
  Dprop_cGN_sel_D4.allocate("Dprop_cGN_sel_D4");
  #ifndef NO_AD_INITIALIZE
  Dprop_cGN_sel_D4.initialize();
  #endif
  Dprop_cGN_sel_cHL4.allocate("Dprop_cGN_sel_cHL4");
  #ifndef NO_AD_INITIALIZE
  Dprop_cGN_sel_cHL4.initialize();
  #endif
  Dprop_cGN_sel_cPT4.allocate("Dprop_cGN_sel_cPT4");
  #ifndef NO_AD_INITIALIZE
  Dprop_cGN_sel_cPT4.initialize();
  #endif
  sel_cPT.allocate(styr,endyr,1,nages,"sel_cPT");
  #ifndef NO_AD_INITIALIZE
    sel_cPT.initialize();
  #endif
  sel_cPT_b01.allocate(1,nages,"sel_cPT_b01");
  #ifndef NO_AD_INITIALIZE
    sel_cPT_b01.initialize();
  #endif
  sel_cPT_b02.allocate(1,nages,"sel_cPT_b02");
  #ifndef NO_AD_INITIALIZE
    sel_cPT_b02.initialize();
  #endif
  sel_cPT_b03.allocate(1,nages,"sel_cPT_b03");
  #ifndef NO_AD_INITIALIZE
    sel_cPT_b03.initialize();
  #endif
  sel_cPT_b04.allocate(1,nages,"sel_cPT_b04");
  #ifndef NO_AD_INITIALIZE
    sel_cPT_b04.initialize();
  #endif
  sel_cPT_b05.allocate(1,nages,"sel_cPT_b05");
  #ifndef NO_AD_INITIALIZE
    sel_cPT_b05.initialize();
  #endif
  selpar_A50b02_cPT.allocate(selpar_A50b02_cPT_LO,selpar_A50b02_cPT_HI,selpar_A50b02_cPT_PH,"selpar_A50b02_cPT");
  selpar_slopeb02_cPT.allocate(selpar_slopeb02_cPT_LO,selpar_slopeb02_cPT_HI,selpar_slopeb02_cPT_PH,"selpar_slopeb02_cPT");
  selpar_A50b03_cPT.allocate(selpar_A50b03_cPT_LO,selpar_A50b03_cPT_HI,selpar_A50b03_cPT_PH,"selpar_A50b03_cPT");
  selpar_slopeb03_cPT.allocate(selpar_slopeb03_cPT_LO,selpar_slopeb03_cPT_HI,selpar_slopeb03_cPT_PH,"selpar_slopeb03_cPT");
  selpar_A50b04_cPT.allocate(selpar_A50b04_cPT_LO,selpar_A50b04_cPT_HI,selpar_A50b04_cPT_PH,"selpar_A50b04_cPT");
  selpar_slopeb04_cPT.allocate(selpar_slopeb04_cPT_LO,selpar_slopeb04_cPT_HI,selpar_slopeb04_cPT_PH,"selpar_slopeb04_cPT");
  selpar_A50b02_cPT_out.allocate(1,8,"selpar_A50b02_cPT_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50b02_cPT_out.initialize();
  #endif
  selpar_slopeb02_cPT_out.allocate(1,8,"selpar_slopeb02_cPT_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slopeb02_cPT_out.initialize();
  #endif
  selpar_A50b03_cPT_out.allocate(1,8,"selpar_A50b03_cPT_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50b03_cPT_out.initialize();
  #endif
  selpar_slopeb03_cPT_out.allocate(1,8,"selpar_slopeb03_cPT_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slopeb03_cPT_out.initialize();
  #endif
  selpar_A50b04_cPT_out.allocate(1,8,"selpar_A50b04_cPT_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50b04_cPT_out.initialize();
  #endif
  selpar_slopeb04_cPT_out.allocate(1,8,"selpar_slopeb04_cPT_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slopeb04_cPT_out.initialize();
  #endif
  sel_cTW.allocate(styr,endyr,1,nages,"sel_cTW");
  #ifndef NO_AD_INITIALIZE
    sel_cTW.initialize();
  #endif
  sel_cTW_b01.allocate(1,nages,"sel_cTW_b01");
  #ifndef NO_AD_INITIALIZE
    sel_cTW_b01.initialize();
  #endif
  sel_cTW_b02.allocate(1,nages,"sel_cTW_b02");
  #ifndef NO_AD_INITIALIZE
    sel_cTW_b02.initialize();
  #endif
  sel_cTW_b03.allocate(1,nages,"sel_cTW_b03");
  #ifndef NO_AD_INITIALIZE
    sel_cTW_b03.initialize();
  #endif
  sel_cTW_b04.allocate(1,nages,"sel_cTW_b04");
  #ifndef NO_AD_INITIALIZE
    sel_cTW_b04.initialize();
  #endif
  sel_cTW_b05.allocate(1,nages,"sel_cTW_b05");
  #ifndef NO_AD_INITIALIZE
    sel_cTW_b05.initialize();
  #endif
  sel_rHB.allocate(styr,endyr,1,nages,"sel_rHB");
  #ifndef NO_AD_INITIALIZE
    sel_rHB.initialize();
  #endif
  sel_rHB_b01.allocate(1,nages,"sel_rHB_b01");
  #ifndef NO_AD_INITIALIZE
    sel_rHB_b01.initialize();
  #endif
  sel_rHB_b02.allocate(1,nages,"sel_rHB_b02");
  #ifndef NO_AD_INITIALIZE
    sel_rHB_b02.initialize();
  #endif
  sel_rHB_b03.allocate(1,nages,"sel_rHB_b03");
  #ifndef NO_AD_INITIALIZE
    sel_rHB_b03.initialize();
  #endif
  sel_rHB_b04.allocate(1,nages,"sel_rHB_b04");
  #ifndef NO_AD_INITIALIZE
    sel_rHB_b04.initialize();
  #endif
  sel_rHB_b05.allocate(1,nages,"sel_rHB_b05");
  #ifndef NO_AD_INITIALIZE
    sel_rHB_b05.initialize();
  #endif
  selpar_A50b01_rHB.allocate(selpar_A50b01_rHB_LO,selpar_A50b01_rHB_HI,selpar_A50b01_rHB_PH,"selpar_A50b01_rHB");
  selpar_slopeb01_rHB.allocate(selpar_slopeb01_rHB_LO,selpar_slopeb01_rHB_HI,selpar_slopeb01_rHB_PH,"selpar_slopeb01_rHB");
  selpar_A50b02_rHB.allocate(selpar_A50b02_rHB_LO,selpar_A50b02_rHB_HI,selpar_A50b02_rHB_PH,"selpar_A50b02_rHB");
  selpar_slopeb02_rHB.allocate(selpar_slopeb02_rHB_LO,selpar_slopeb02_rHB_HI,selpar_slopeb02_rHB_PH,"selpar_slopeb02_rHB");
  selpar_A50b03_rHB.allocate(selpar_A50b03_rHB_LO,selpar_A50b03_rHB_HI,selpar_A50b03_rHB_PH,"selpar_A50b03_rHB");
  selpar_slopeb03_rHB.allocate(selpar_slopeb03_rHB_LO,selpar_slopeb03_rHB_HI,selpar_slopeb03_rHB_PH,"selpar_slopeb03_rHB");
  selpar_A50b04_rHB.allocate(selpar_A50b04_rHB_LO,selpar_A50b04_rHB_HI,selpar_A50b04_rHB_PH,"selpar_A50b04_rHB");
  selpar_slopeb04_rHB.allocate(selpar_slopeb04_rHB_LO,selpar_slopeb04_rHB_HI,selpar_slopeb04_rHB_PH,"selpar_slopeb04_rHB");
  selpar_A50b05_rHB.allocate(selpar_A50b05_rHB_LO,selpar_A50b05_rHB_HI,selpar_A50b05_rHB_PH,"selpar_A50b05_rHB");
  selpar_slopeb05_rHB.allocate(selpar_slopeb05_rHB_LO,selpar_slopeb05_rHB_HI,selpar_slopeb05_rHB_PH,"selpar_slopeb05_rHB");
  selpar_A50b01_rHB_out.allocate(1,8,"selpar_A50b01_rHB_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50b01_rHB_out.initialize();
  #endif
  selpar_slopeb01_rHB_out.allocate(1,8,"selpar_slopeb01_rHB_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slopeb01_rHB_out.initialize();
  #endif
  selpar_A50b02_rHB_out.allocate(1,8,"selpar_A50b02_rHB_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50b02_rHB_out.initialize();
  #endif
  selpar_slopeb02_rHB_out.allocate(1,8,"selpar_slopeb02_rHB_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slopeb02_rHB_out.initialize();
  #endif
  selpar_A50b03_rHB_out.allocate(1,8,"selpar_A50b03_rHB_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50b03_rHB_out.initialize();
  #endif
  selpar_slopeb03_rHB_out.allocate(1,8,"selpar_slopeb03_rHB_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slopeb03_rHB_out.initialize();
  #endif
  selpar_A50b04_rHB_out.allocate(1,8,"selpar_A50b04_rHB_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50b04_rHB_out.initialize();
  #endif
  selpar_slopeb04_rHB_out.allocate(1,8,"selpar_slopeb04_rHB_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slopeb04_rHB_out.initialize();
  #endif
  selpar_A50b05_rHB_out.allocate(1,8,"selpar_A50b05_rHB_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50b05_rHB_out.initialize();
  #endif
  selpar_slopeb05_rHB_out.allocate(1,8,"selpar_slopeb05_rHB_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slopeb05_rHB_out.initialize();
  #endif
  selpar_Age0_rHD.allocate("selpar_Age0_rHD");
  #ifndef NO_AD_INITIALIZE
  selpar_Age0_rHD.initialize();
  #endif
  selpar_Age1_rHD.allocate("selpar_Age1_rHD");
  #ifndef NO_AD_INITIALIZE
  selpar_Age1_rHD.initialize();
  #endif
  selpar_Age2_rHD.allocate("selpar_Age2_rHD");
  #ifndef NO_AD_INITIALIZE
  selpar_Age2_rHD.initialize();
  #endif
  sel_rHD.allocate(styr,endyr,1,nages,"sel_rHD");
  #ifndef NO_AD_INITIALIZE
    sel_rHD.initialize();
  #endif
  sel_rHD_b01.allocate(1,nages,"sel_rHD_b01");
  #ifndef NO_AD_INITIALIZE
    sel_rHD_b01.initialize();
  #endif
  sel_rHD_b02.allocate(1,nages,"sel_rHD_b02");
  #ifndef NO_AD_INITIALIZE
    sel_rHD_b02.initialize();
  #endif
  sel_rHD_b03.allocate(1,nages,"sel_rHD_b03");
  #ifndef NO_AD_INITIALIZE
    sel_rHD_b03.initialize();
  #endif
  sel_rHD_b04.allocate(1,nages,"sel_rHD_b04");
  #ifndef NO_AD_INITIALIZE
    sel_rHD_b04.initialize();
  #endif
  sel_rHD_b05.allocate(1,nages,"sel_rHD_b05");
  #ifndef NO_AD_INITIALIZE
    sel_rHD_b05.initialize();
  #endif
  selpar_A50b04_rHD.allocate(selpar_A50b04_rHD_LO,selpar_A50b04_rHD_HI,selpar_A50b04_rHD_PH,"selpar_A50b04_rHD");
  selpar_slopeb04_rHD.allocate(selpar_slopeb04_rHD_LO,selpar_slopeb04_rHD_HI,selpar_slopeb04_rHD_PH,"selpar_slopeb04_rHD");
  selpar_A502b04_rHD.allocate(selpar_A502b04_rHD_LO,selpar_A502b04_rHD_HI,selpar_A502b04_rHD_PH,"selpar_A502b04_rHD");
  selpar_slope2b04_rHD.allocate(selpar_slope2b04_rHD_LO,selpar_slope2b04_rHD_HI,selpar_slope2b04_rHD_PH,"selpar_slope2b04_rHD");
  selpar_A50b05_rHD.allocate(selpar_A50b05_rHD_LO,selpar_A50b05_rHD_HI,selpar_A50b05_rHD_PH,"selpar_A50b05_rHD");
  selpar_slopeb05_rHD.allocate(selpar_slopeb05_rHD_LO,selpar_slopeb05_rHD_HI,selpar_slopeb05_rHD_PH,"selpar_slopeb05_rHD");
  selpar_A502b05_rHD.allocate(selpar_A502b05_rHD_LO,selpar_A502b05_rHD_HI,selpar_A502b05_rHD_PH,"selpar_A502b05_rHD");
  selpar_slope2b05_rHD.allocate(selpar_slope2b05_rHD_LO,selpar_slope2b05_rHD_HI,selpar_slope2b05_rHD_PH,"selpar_slope2b05_rHD");
  selpar_A50b04_rHD_out.allocate(1,8,"selpar_A50b04_rHD_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50b04_rHD_out.initialize();
  #endif
  selpar_slopeb04_rHD_out.allocate(1,8,"selpar_slopeb04_rHD_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slopeb04_rHD_out.initialize();
  #endif
  selpar_A502b04_rHD_out.allocate(1,8,"selpar_A502b04_rHD_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502b04_rHD_out.initialize();
  #endif
  selpar_slope2b04_rHD_out.allocate(1,8,"selpar_slope2b04_rHD_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2b04_rHD_out.initialize();
  #endif
  selpar_A50b05_rHD_out.allocate(1,8,"selpar_A50b05_rHD_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50b05_rHD_out.initialize();
  #endif
  selpar_slopeb05_rHD_out.allocate(1,8,"selpar_slopeb05_rHD_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slopeb05_rHD_out.initialize();
  #endif
  selpar_A502b05_rHD_out.allocate(1,8,"selpar_A502b05_rHD_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502b05_rHD_out.initialize();
  #endif
  selpar_slope2b05_rHD_out.allocate(1,8,"selpar_slope2b05_rHD_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2b05_rHD_out.initialize();
  #endif
  selpar_logitAge0_rHD.allocate(selpar_logitAge0_rHD_LO,selpar_logitAge0_rHD_HI,selpar_logitAge0_rHD_PH,"selpar_logitAge0_rHD");
  selpar_logitAge1_rHD.allocate(selpar_logitAge1_rHD_LO,selpar_logitAge1_rHD_HI,selpar_logitAge1_rHD_PH,"selpar_logitAge1_rHD");
  selpar_logitAge2_rHD.allocate(selpar_logitAge2_rHD_LO,selpar_logitAge2_rHD_HI,selpar_logitAge2_rHD_PH,"selpar_logitAge2_rHD");
  selpar_logitAge0_rHD_out.allocate(1,8,"selpar_logitAge0_rHD_out");
  #ifndef NO_AD_INITIALIZE
    selpar_logitAge0_rHD_out.initialize();
  #endif
  selpar_logitAge1_rHD_out.allocate(1,8,"selpar_logitAge1_rHD_out");
  #ifndef NO_AD_INITIALIZE
    selpar_logitAge1_rHD_out.initialize();
  #endif
  selpar_logitAge2_rHD_out.allocate(1,8,"selpar_logitAge2_rHD_out");
  #ifndef NO_AD_INITIALIZE
    selpar_logitAge2_rHD_out.initialize();
  #endif
  prob_belowsizelim_block1.allocate(1,nages,"prob_belowsizelim_block1");
  #ifndef NO_AD_INITIALIZE
    prob_belowsizelim_block1.initialize();
  #endif
  prob_belowsizelim_block2.allocate(1,nages,"prob_belowsizelim_block2");
  #ifndef NO_AD_INITIALIZE
    prob_belowsizelim_block2.initialize();
  #endif
  prob_belowsizelim_block3.allocate(1,nages,"prob_belowsizelim_block3");
  #ifndef NO_AD_INITIALIZE
    prob_belowsizelim_block3.initialize();
  #endif
  prob_belowsizelim_block4.allocate(1,nages,"prob_belowsizelim_block4");
  #ifndef NO_AD_INITIALIZE
    prob_belowsizelim_block4.initialize();
  #endif
  prob_belowsizelim_block5.allocate(1,nages,"prob_belowsizelim_block5");
  #ifndef NO_AD_INITIALIZE
    prob_belowsizelim_block5.initialize();
  #endif
  zscore_lsizelim1.allocate("zscore_lsizelim1");
  #ifndef NO_AD_INITIALIZE
  zscore_lsizelim1.initialize();
  #endif
  zscore_lsizelim2.allocate("zscore_lsizelim2");
  #ifndef NO_AD_INITIALIZE
  zscore_lsizelim2.initialize();
  #endif
  zscore_lsizelim3.allocate("zscore_lsizelim3");
  #ifndef NO_AD_INITIALIZE
  zscore_lsizelim3.initialize();
  #endif
  zscore_lsizelim4.allocate("zscore_lsizelim4");
  #ifndef NO_AD_INITIALIZE
  zscore_lsizelim4.initialize();
  #endif
  zscore_lsizelim5.allocate("zscore_lsizelim5");
  #ifndef NO_AD_INITIALIZE
  zscore_lsizelim5.initialize();
  #endif
  cprob_lsizelim1.allocate("cprob_lsizelim1");
  #ifndef NO_AD_INITIALIZE
  cprob_lsizelim1.initialize();
  #endif
  cprob_lsizelim2.allocate("cprob_lsizelim2");
  #ifndef NO_AD_INITIALIZE
  cprob_lsizelim2.initialize();
  #endif
  cprob_lsizelim3.allocate("cprob_lsizelim3");
  #ifndef NO_AD_INITIALIZE
  cprob_lsizelim3.initialize();
  #endif
  cprob_lsizelim4.allocate("cprob_lsizelim4");
  #ifndef NO_AD_INITIALIZE
  cprob_lsizelim4.initialize();
  #endif
  cprob_lsizelim5.allocate("cprob_lsizelim5");
  #ifndef NO_AD_INITIALIZE
  cprob_lsizelim5.initialize();
  #endif
  sel_rGN.allocate(styr,endyr,1,nages,"sel_rGN");
  #ifndef NO_AD_INITIALIZE
    sel_rGN.initialize();
  #endif
  sel_rGN_b01.allocate(1,nages,"sel_rGN_b01");
  #ifndef NO_AD_INITIALIZE
    sel_rGN_b01.initialize();
  #endif
  sel_rGN_b02.allocate(1,nages,"sel_rGN_b02");
  #ifndef NO_AD_INITIALIZE
    sel_rGN_b02.initialize();
  #endif
  sel_rGN_b03.allocate(1,nages,"sel_rGN_b03");
  #ifndef NO_AD_INITIALIZE
    sel_rGN_b03.initialize();
  #endif
  sel_rGN_b04.allocate(1,nages,"sel_rGN_b04");
  #ifndef NO_AD_INITIALIZE
    sel_rGN_b04.initialize();
  #endif
  sel_rGN_b05.allocate(1,nages,"sel_rGN_b05");
  #ifndef NO_AD_INITIALIZE
    sel_rGN_b05.initialize();
  #endif
  sel_rIN.allocate(styr,endyr,1,nages,"sel_rIN");
  #ifndef NO_AD_INITIALIZE
    sel_rIN.initialize();
  #endif
  sel_rIN_b01.allocate(1,nages,"sel_rIN_b01");
  #ifndef NO_AD_INITIALIZE
    sel_rIN_b01.initialize();
  #endif
  sel_rIN_b02.allocate(1,nages,"sel_rIN_b02");
  #ifndef NO_AD_INITIALIZE
    sel_rIN_b02.initialize();
  #endif
  sel_rIN_b03.allocate(1,nages,"sel_rIN_b03");
  #ifndef NO_AD_INITIALIZE
    sel_rIN_b03.initialize();
  #endif
  sel_rIN_b04.allocate(1,nages,"sel_rIN_b04");
  #ifndef NO_AD_INITIALIZE
    sel_rIN_b04.initialize();
  #endif
  sel_rIN_b05.allocate(1,nages,"sel_rIN_b05");
  #ifndef NO_AD_INITIALIZE
    sel_rIN_b05.initialize();
  #endif
  sel_rGD.allocate(styr,endyr,1,nages,"sel_rGD");
  #ifndef NO_AD_INITIALIZE
    sel_rGD.initialize();
  #endif
  sel_rGD_b01.allocate(1,nages,"sel_rGD_b01");
  #ifndef NO_AD_INITIALIZE
    sel_rGD_b01.initialize();
  #endif
  sel_rGD_b02.allocate(1,nages,"sel_rGD_b02");
  #ifndef NO_AD_INITIALIZE
    sel_rGD_b02.initialize();
  #endif
  sel_rGD_b03.allocate(1,nages,"sel_rGD_b03");
  #ifndef NO_AD_INITIALIZE
    sel_rGD_b03.initialize();
  #endif
  sel_rGD_b04.allocate(1,nages,"sel_rGD_b04");
  #ifndef NO_AD_INITIALIZE
    sel_rGD_b04.initialize();
  #endif
  sel_rGD_b05.allocate(1,nages,"sel_rGD_b05");
  #ifndef NO_AD_INITIALIZE
    sel_rGD_b05.initialize();
  #endif
  sel_rID.allocate(styr,endyr,1,nages,"sel_rID");
  #ifndef NO_AD_INITIALIZE
    sel_rID.initialize();
  #endif
  sel_rID_b01.allocate(1,nages,"sel_rID_b01");
  #ifndef NO_AD_INITIALIZE
    sel_rID_b01.initialize();
  #endif
  sel_rID_b02.allocate(1,nages,"sel_rID_b02");
  #ifndef NO_AD_INITIALIZE
    sel_rID_b02.initialize();
  #endif
  sel_rID_b03.allocate(1,nages,"sel_rID_b03");
  #ifndef NO_AD_INITIALIZE
    sel_rID_b03.initialize();
  #endif
  sel_rID_b04.allocate(1,nages,"sel_rID_b04");
  #ifndef NO_AD_INITIALIZE
    sel_rID_b04.initialize();
  #endif
  sel_rID_b05.allocate(1,nages,"sel_rID_b05");
  #ifndef NO_AD_INITIALIZE
    sel_rID_b05.initialize();
  #endif
  sel_wgted_L.allocate(1,nages,"sel_wgted_L");
  #ifndef NO_AD_INITIALIZE
    sel_wgted_L.initialize();
  #endif
  sel_wgted_D.allocate(1,nages,"sel_wgted_D");
  #ifndef NO_AD_INITIALIZE
    sel_wgted_D.initialize();
  #endif
  sel_wgted_tot.allocate(1,nages,"sel_wgted_tot");
  #ifndef NO_AD_INITIALIZE
    sel_wgted_tot.initialize();
  #endif
  pred_cpue_sBT.allocate(styr_cpue_sBT,endyr_cpue_sBT,"pred_cpue_sBT");
  #ifndef NO_AD_INITIALIZE
    pred_cpue_sBT.initialize();
  #endif
  N_sBT.allocate(styr_cpue_sBT,endyr_cpue_sBT,1,nages,"N_sBT");
  #ifndef NO_AD_INITIALIZE
    N_sBT.initialize();
  #endif
  pred_cpue_sTV.allocate(styr_cpue_sTV,endyr_cpue_sTV,"pred_cpue_sTV");
  #ifndef NO_AD_INITIALIZE
    pred_cpue_sTV.initialize();
  #endif
  obs_cpue_corrected_sTV.allocate(styr_cpue_sTV,endyr_cpue_sTV,"obs_cpue_corrected_sTV");
  #ifndef NO_AD_INITIALIZE
    obs_cpue_corrected_sTV.initialize();
  #endif
  obs_cv_cpue_corrected_sTV.allocate(styr_cpue_sTV,endyr_cpue_sTV,"obs_cv_cpue_corrected_sTV");
  #ifndef NO_AD_INITIALIZE
    obs_cv_cpue_corrected_sTV.initialize();
  #endif
  posobs_cpue_sTV.allocate(styr_cpue_sTV,endyr_cpue_sTV,"posobs_cpue_sTV");
  #ifndef NO_AD_INITIALIZE
    posobs_cpue_sTV.initialize();
  #endif
  N_sTV.allocate(styr_cpue_sTV,endyr_cpue_sTV,1,nages,"N_sTV");
  #ifndef NO_AD_INITIALIZE
    N_sTV.initialize();
  #endif
  pred_cpue_cHL.allocate(styr_cpue_cHL,endyr_cpue_cHL,"pred_cpue_cHL");
  #ifndef NO_AD_INITIALIZE
    pred_cpue_cHL.initialize();
  #endif
  N_cHL.allocate(styr_cpue_cHL,endyr_cpue_cHL,1,nages,"N_cHL");
  #ifndef NO_AD_INITIALIZE
    N_cHL.initialize();
  #endif
  pred_cpue_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,"pred_cpue_rHB");
  #ifndef NO_AD_INITIALIZE
    pred_cpue_rHB.initialize();
  #endif
  N_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,1,nages,"N_rHB");
  #ifndef NO_AD_INITIALIZE
    N_rHB.initialize();
  #endif
  q_rate.allocate(q_rate_LO,q_rate_HI,q_rate_PH,"q_rate");
  q_rate_out.allocate(1,8,"q_rate_out");
  #ifndef NO_AD_INITIALIZE
    q_rate_out.initialize();
  #endif
  q_rate_phase.allocate("q_rate_phase");
  #ifndef NO_AD_INITIALIZE
  q_rate_phase.initialize();
  #endif
  log_q_sBT.allocate(log_q_sBT_LO,log_q_sBT_HI,log_q_sBT_PH,"log_q_sBT");
  log_q_sTV.allocate(log_q_sTV_LO,log_q_sTV_HI,log_q_sTV_PH,"log_q_sTV");
  log_q_cHL.allocate(log_q_cHL_LO,log_q_cHL_HI,log_q_cHL_PH,"log_q_cHL");
  log_q_rHB.allocate(log_q_rHB_LO,log_q_rHB_HI,log_q_rHB_PH,"log_q_rHB");
  log_q_sBT_out.allocate(1,8,"log_q_sBT_out");
  #ifndef NO_AD_INITIALIZE
    log_q_sBT_out.initialize();
  #endif
  log_q_sTV_out.allocate(1,8,"log_q_sTV_out");
  #ifndef NO_AD_INITIALIZE
    log_q_sTV_out.initialize();
  #endif
  log_q_cHL_out.allocate(1,8,"log_q_cHL_out");
  #ifndef NO_AD_INITIALIZE
    log_q_cHL_out.initialize();
  #endif
  log_q_rHB_out.allocate(1,8,"log_q_rHB_out");
  #ifndef NO_AD_INITIALIZE
    log_q_rHB_out.initialize();
  #endif
  q_rate_fcn_cHL.allocate(styr_cpue_cHL,endyr_cpue_cHL,"q_rate_fcn_cHL");
  #ifndef NO_AD_INITIALIZE
    q_rate_fcn_cHL.initialize();
  #endif
  q_rate_fcn_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,"q_rate_fcn_rHB");
  #ifndef NO_AD_INITIALIZE
    q_rate_fcn_rHB.initialize();
  #endif
  q_DD_beta.allocate("q_DD_beta");
  #ifndef NO_AD_INITIALIZE
  q_DD_beta.initialize();
  #endif
  q_DD_fcn.allocate(styr,endyr,"q_DD_fcn");
  #ifndef NO_AD_INITIALIZE
    q_DD_fcn.initialize();
  #endif
  B0_q_DD.allocate("B0_q_DD");
  #ifndef NO_AD_INITIALIZE
  B0_q_DD.initialize();
  #endif
  B_q_DD.allocate(styr,endyr,"B_q_DD");
  #ifndef NO_AD_INITIALIZE
    B_q_DD.initialize();
  #endif
  q_RW_log_dev_cHL.allocate(styr_cpue_cHL,endyr_cpue_cHL-1,log_RWq_LO,log_RWq_HI,log_RWq_PH,"q_RW_log_dev_cHL");
  q_RW_log_dev_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB-1,log_RWq_LO,log_RWq_HI,log_RWq_PH,"q_RW_log_dev_rHB");
  q_cHL.allocate(styr_cpue_cHL,endyr_cpue_cHL,"q_cHL");
  #ifndef NO_AD_INITIALIZE
    q_cHL.initialize();
  #endif
  q_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,"q_rHB");
  #ifndef NO_AD_INITIALIZE
    q_rHB.initialize();
  #endif
  L_cHL_num.allocate(styr,endyr,1,nages,"L_cHL_num");
  #ifndef NO_AD_INITIALIZE
    L_cHL_num.initialize();
  #endif
  L_cHL_klb.allocate(styr,endyr,1,nages,"L_cHL_klb");
  #ifndef NO_AD_INITIALIZE
    L_cHL_klb.initialize();
  #endif
  pred_L_cHL_knum.allocate(styr,endyr,"pred_L_cHL_knum");
  #ifndef NO_AD_INITIALIZE
    pred_L_cHL_knum.initialize();
  #endif
  pred_L_cHL_klb.allocate(styr,endyr,"pred_L_cHL_klb");
  #ifndef NO_AD_INITIALIZE
    pred_L_cHL_klb.initialize();
  #endif
  L_cPT_num.allocate(styr,endyr,1,nages,"L_cPT_num");
  #ifndef NO_AD_INITIALIZE
    L_cPT_num.initialize();
  #endif
  L_cPT_klb.allocate(styr,endyr,1,nages,"L_cPT_klb");
  #ifndef NO_AD_INITIALIZE
    L_cPT_klb.initialize();
  #endif
  pred_L_cPT_knum.allocate(styr,endyr,"pred_L_cPT_knum");
  #ifndef NO_AD_INITIALIZE
    pred_L_cPT_knum.initialize();
  #endif
  pred_L_cPT_klb.allocate(styr,endyr,"pred_L_cPT_klb");
  #ifndef NO_AD_INITIALIZE
    pred_L_cPT_klb.initialize();
  #endif
  L_cTW_num.allocate(styr,endyr,1,nages,"L_cTW_num");
  #ifndef NO_AD_INITIALIZE
    L_cTW_num.initialize();
  #endif
  L_cTW_klb.allocate(styr,endyr,1,nages,"L_cTW_klb");
  #ifndef NO_AD_INITIALIZE
    L_cTW_klb.initialize();
  #endif
  pred_L_cTW_knum.allocate(styr,endyr,"pred_L_cTW_knum");
  #ifndef NO_AD_INITIALIZE
    pred_L_cTW_knum.initialize();
  #endif
  pred_L_cTW_klb.allocate(styr,endyr,"pred_L_cTW_klb");
  #ifndef NO_AD_INITIALIZE
    pred_L_cTW_klb.initialize();
  #endif
  L_rHB_num.allocate(styr,endyr,1,nages,"L_rHB_num");
  #ifndef NO_AD_INITIALIZE
    L_rHB_num.initialize();
  #endif
  L_rHB_klb.allocate(styr,endyr,1,nages,"L_rHB_klb");
  #ifndef NO_AD_INITIALIZE
    L_rHB_klb.initialize();
  #endif
  pred_L_rHB_knum.allocate(styr,endyr,"pred_L_rHB_knum");
  #ifndef NO_AD_INITIALIZE
    pred_L_rHB_knum.initialize();
  #endif
  pred_L_rHB_klb.allocate(styr,endyr,"pred_L_rHB_klb");
  #ifndef NO_AD_INITIALIZE
    pred_L_rHB_klb.initialize();
  #endif
  L_rGN_num.allocate(styr,endyr,1,nages,"L_rGN_num");
  #ifndef NO_AD_INITIALIZE
    L_rGN_num.initialize();
  #endif
  L_rGN_klb.allocate(styr,endyr,1,nages,"L_rGN_klb");
  #ifndef NO_AD_INITIALIZE
    L_rGN_klb.initialize();
  #endif
  pred_L_rGN_knum.allocate(styr,endyr,"pred_L_rGN_knum");
  #ifndef NO_AD_INITIALIZE
    pred_L_rGN_knum.initialize();
  #endif
  pred_L_rGN_klb.allocate(styr,endyr,"pred_L_rGN_klb");
  #ifndef NO_AD_INITIALIZE
    pred_L_rGN_klb.initialize();
  #endif
  L_rIN_num.allocate(styr,endyr,1,nages,"L_rIN_num");
  #ifndef NO_AD_INITIALIZE
    L_rIN_num.initialize();
  #endif
  L_rIN_klb.allocate(styr,endyr,1,nages,"L_rIN_klb");
  #ifndef NO_AD_INITIALIZE
    L_rIN_klb.initialize();
  #endif
  pred_L_rIN_knum.allocate(styr,endyr,"pred_L_rIN_knum");
  #ifndef NO_AD_INITIALIZE
    pred_L_rIN_knum.initialize();
  #endif
  pred_L_rIN_klb.allocate(styr,endyr,"pred_L_rIN_klb");
  #ifndef NO_AD_INITIALIZE
    pred_L_rIN_klb.initialize();
  #endif
  L_total_num.allocate(styr,endyr,1,nages,"L_total_num");
  #ifndef NO_AD_INITIALIZE
    L_total_num.initialize();
  #endif
  L_total_klb.allocate(styr,endyr,1,nages,"L_total_klb");
  #ifndef NO_AD_INITIALIZE
    L_total_klb.initialize();
  #endif
  L_total_knum_yr.allocate(styr,endyr,"L_total_knum_yr");
  #ifndef NO_AD_INITIALIZE
    L_total_knum_yr.initialize();
  #endif
  L_total_klb_yr.allocate(styr,endyr,"L_total_klb_yr");
  #ifndef NO_AD_INITIALIZE
    L_total_klb_yr.initialize();
  #endif
  D_cGD_num.allocate(styr,endyr,1,nages,"D_cGD_num");
  #ifndef NO_AD_INITIALIZE
    D_cGD_num.initialize();
  #endif
  D_cGD_klb.allocate(styr,endyr,1,nages,"D_cGD_klb");
  #ifndef NO_AD_INITIALIZE
    D_cGD_klb.initialize();
  #endif
  pred_D_cGD_knum.allocate(styr,endyr,"pred_D_cGD_knum");
  #ifndef NO_AD_INITIALIZE
    pred_D_cGD_knum.initialize();
  #endif
  pred_D_cGD_klb.allocate(styr,endyr,"pred_D_cGD_klb");
  #ifndef NO_AD_INITIALIZE
    pred_D_cGD_klb.initialize();
  #endif
  obs_D_cGD.allocate(styr_D_cGD,endyr_D_cGD,"obs_D_cGD");
  #ifndef NO_AD_INITIALIZE
    obs_D_cGD.initialize();
  #endif
  D_rHD_num.allocate(styr,endyr,1,nages,"D_rHD_num");
  #ifndef NO_AD_INITIALIZE
    D_rHD_num.initialize();
  #endif
  D_rHD_klb.allocate(styr,endyr,1,nages,"D_rHD_klb");
  #ifndef NO_AD_INITIALIZE
    D_rHD_klb.initialize();
  #endif
  pred_D_rHD_knum.allocate(styr,endyr,"pred_D_rHD_knum");
  #ifndef NO_AD_INITIALIZE
    pred_D_rHD_knum.initialize();
  #endif
  pred_D_rHD_klb.allocate(styr,endyr,"pred_D_rHD_klb");
  #ifndef NO_AD_INITIALIZE
    pred_D_rHD_klb.initialize();
  #endif
  obs_D_rHD.allocate(styr_D_rHD,endyr_D_rHD,"obs_D_rHD");
  #ifndef NO_AD_INITIALIZE
    obs_D_rHD.initialize();
  #endif
  D_rGD_num.allocate(styr,endyr,1,nages,"D_rGD_num");
  #ifndef NO_AD_INITIALIZE
    D_rGD_num.initialize();
  #endif
  D_rGD_klb.allocate(styr,endyr,1,nages,"D_rGD_klb");
  #ifndef NO_AD_INITIALIZE
    D_rGD_klb.initialize();
  #endif
  pred_D_rGD_knum.allocate(styr,endyr,"pred_D_rGD_knum");
  #ifndef NO_AD_INITIALIZE
    pred_D_rGD_knum.initialize();
  #endif
  pred_D_rGD_klb.allocate(styr,endyr,"pred_D_rGD_klb");
  #ifndef NO_AD_INITIALIZE
    pred_D_rGD_klb.initialize();
  #endif
  obs_D_rGD.allocate(styr_D_rGD,endyr_D_rGD,"obs_D_rGD");
  #ifndef NO_AD_INITIALIZE
    obs_D_rGD.initialize();
  #endif
  D_rID_num.allocate(styr,endyr,1,nages,"D_rID_num");
  #ifndef NO_AD_INITIALIZE
    D_rID_num.initialize();
  #endif
  D_rID_klb.allocate(styr,endyr,1,nages,"D_rID_klb");
  #ifndef NO_AD_INITIALIZE
    D_rID_klb.initialize();
  #endif
  pred_D_rID_knum.allocate(styr,endyr,"pred_D_rID_knum");
  #ifndef NO_AD_INITIALIZE
    pred_D_rID_knum.initialize();
  #endif
  pred_D_rID_klb.allocate(styr,endyr,"pred_D_rID_klb");
  #ifndef NO_AD_INITIALIZE
    pred_D_rID_klb.initialize();
  #endif
  obs_D_rID.allocate(styr_D_rID,endyr_D_rID,"obs_D_rID");
  #ifndef NO_AD_INITIALIZE
    obs_D_rID.initialize();
  #endif
  D_total_num.allocate(styr,endyr,1,nages,"D_total_num");
  #ifndef NO_AD_INITIALIZE
    D_total_num.initialize();
  #endif
  D_total_klb.allocate(styr,endyr,1,nages,"D_total_klb");
  #ifndef NO_AD_INITIALIZE
    D_total_klb.initialize();
  #endif
  D_total_knum_yr.allocate(styr,endyr,"D_total_knum_yr");
  #ifndef NO_AD_INITIALIZE
    D_total_knum_yr.initialize();
  #endif
  D_total_klb_yr.allocate(styr,endyr,"D_total_klb_yr");
  #ifndef NO_AD_INITIALIZE
    D_total_klb_yr.initialize();
  #endif
  prop_F_cHL.allocate("prop_F_cHL");
  #ifndef NO_AD_INITIALIZE
  prop_F_cHL.initialize();
  #endif
  prop_F_cPT.allocate("prop_F_cPT");
  #ifndef NO_AD_INITIALIZE
  prop_F_cPT.initialize();
  #endif
  prop_F_rHB.allocate("prop_F_rHB");
  #ifndef NO_AD_INITIALIZE
  prop_F_rHB.initialize();
  #endif
  prop_F_rGN.allocate("prop_F_rGN");
  #ifndef NO_AD_INITIALIZE
  prop_F_rGN.initialize();
  #endif
  prop_F_rIN.allocate("prop_F_rIN");
  #ifndef NO_AD_INITIALIZE
  prop_F_rIN.initialize();
  #endif
  prop_F_cGD.allocate("prop_F_cGD");
  #ifndef NO_AD_INITIALIZE
  prop_F_cGD.initialize();
  #endif
  prop_F_rHD.allocate("prop_F_rHD");
  #ifndef NO_AD_INITIALIZE
  prop_F_rHD.initialize();
  #endif
  prop_F_rGD.allocate("prop_F_rGD");
  #ifndef NO_AD_INITIALIZE
  prop_F_rGD.initialize();
  #endif
  prop_F_rID.allocate("prop_F_rID");
  #ifndef NO_AD_INITIALIZE
  prop_F_rID.initialize();
  #endif
  F_temp_sum.allocate("F_temp_sum");
  #ifndef NO_AD_INITIALIZE
  F_temp_sum.initialize();
  #endif
  F_end.allocate(1,nages,"F_end");
  #ifndef NO_AD_INITIALIZE
    F_end.initialize();
  #endif
  F_end_L.allocate(1,nages,"F_end_L");
  #ifndef NO_AD_INITIALIZE
    F_end_L.initialize();
  #endif
  F_end_D.allocate(1,nages,"F_end_D");
  #ifndef NO_AD_INITIALIZE
    F_end_D.initialize();
  #endif
  F_end_apex.allocate("F_end_apex");
  #ifndef NO_AD_INITIALIZE
  F_end_apex.initialize();
  #endif
  SSB_msy_n_out.allocate("SSB_msy_n_out");
  #ifndef NO_AD_INITIALIZE
  SSB_msy_n_out.initialize();
  #endif
  F_msy_n_out.allocate("F_msy_n_out");
  #ifndef NO_AD_INITIALIZE
  F_msy_n_out.initialize();
  #endif
  msy_n_klb_out.allocate("msy_n_klb_out");
  #ifndef NO_AD_INITIALIZE
  msy_n_klb_out.initialize();
  #endif
  msy_n_knum_out.allocate("msy_n_knum_out");
  #ifndef NO_AD_INITIALIZE
  msy_n_knum_out.initialize();
  #endif
  B_msy_n_out.allocate("B_msy_n_out");
  #ifndef NO_AD_INITIALIZE
  B_msy_n_out.initialize();
  #endif
  R_msy_n_out.allocate("R_msy_n_out");
  #ifndef NO_AD_INITIALIZE
  R_msy_n_out.initialize();
  #endif
  D_msy_n_knum_out.allocate("D_msy_n_knum_out");
  #ifndef NO_AD_INITIALIZE
  D_msy_n_knum_out.initialize();
  #endif
  D_msy_n_klb_out.allocate("D_msy_n_klb_out");
  #ifndef NO_AD_INITIALIZE
  D_msy_n_klb_out.initialize();
  #endif
  T_msy_n_knum_out.allocate("T_msy_n_knum_out");
  #ifndef NO_AD_INITIALIZE
  T_msy_n_knum_out.initialize();
  #endif
  T_msy_n_klb_out.allocate("T_msy_n_klb_out");
  #ifndef NO_AD_INITIALIZE
  T_msy_n_klb_out.initialize();
  #endif
  spr_msy_n_out.allocate("spr_msy_n_out");
  #ifndef NO_AD_INITIALIZE
  spr_msy_n_out.initialize();
  #endif
  SSB_msy_lb_out.allocate("SSB_msy_lb_out");
  #ifndef NO_AD_INITIALIZE
  SSB_msy_lb_out.initialize();
  #endif
  F_msy_lb_out.allocate("F_msy_lb_out");
  #ifndef NO_AD_INITIALIZE
  F_msy_lb_out.initialize();
  #endif
  msy_lb_klb_out.allocate("msy_lb_klb_out");
  #ifndef NO_AD_INITIALIZE
  msy_lb_klb_out.initialize();
  #endif
  msy_lb_knum_out.allocate("msy_lb_knum_out");
  #ifndef NO_AD_INITIALIZE
  msy_lb_knum_out.initialize();
  #endif
  B_msy_lb_out.allocate("B_msy_lb_out");
  #ifndef NO_AD_INITIALIZE
  B_msy_lb_out.initialize();
  #endif
  R_msy_lb_out.allocate("R_msy_lb_out");
  #ifndef NO_AD_INITIALIZE
  R_msy_lb_out.initialize();
  #endif
  D_msy_lb_knum_out.allocate("D_msy_lb_knum_out");
  #ifndef NO_AD_INITIALIZE
  D_msy_lb_knum_out.initialize();
  #endif
  D_msy_lb_klb_out.allocate("D_msy_lb_klb_out");
  #ifndef NO_AD_INITIALIZE
  D_msy_lb_klb_out.initialize();
  #endif
  T_msy_lb_knum_out.allocate("T_msy_lb_knum_out");
  #ifndef NO_AD_INITIALIZE
  T_msy_lb_knum_out.initialize();
  #endif
  T_msy_lb_klb_out.allocate("T_msy_lb_klb_out");
  #ifndef NO_AD_INITIALIZE
  T_msy_lb_klb_out.initialize();
  #endif
  spr_msy_lb_out.allocate("spr_msy_lb_out");
  #ifndef NO_AD_INITIALIZE
  spr_msy_lb_out.initialize();
  #endif
  F30_dum.allocate("F30_dum");
  #ifndef NO_AD_INITIALIZE
  F30_dum.initialize();
  #endif
  F40_dum.allocate("F40_dum");
  #ifndef NO_AD_INITIALIZE
  F40_dum.initialize();
  #endif
  F50_dum.allocate("F50_dum");
  #ifndef NO_AD_INITIALIZE
  F50_dum.initialize();
  #endif
  Fproxy_dum.allocate("Fproxy_dum");
  #ifndef NO_AD_INITIALIZE
  Fproxy_dum.initialize();
  #endif
  Fproxy_out.allocate("Fproxy_out");
  #ifndef NO_AD_INITIALIZE
  Fproxy_out.initialize();
  #endif
  SSB_Fproxy_out.allocate("SSB_Fproxy_out");
  #ifndef NO_AD_INITIALIZE
  SSB_Fproxy_out.initialize();
  #endif
  B_Fproxy_out.allocate("B_Fproxy_out");
  #ifndef NO_AD_INITIALIZE
  B_Fproxy_out.initialize();
  #endif
  R_Fproxy_out.allocate("R_Fproxy_out");
  #ifndef NO_AD_INITIALIZE
  R_Fproxy_out.initialize();
  #endif
  L_Fproxy_knum_out.allocate("L_Fproxy_knum_out");
  #ifndef NO_AD_INITIALIZE
  L_Fproxy_knum_out.initialize();
  #endif
  L_Fproxy_klb_out.allocate("L_Fproxy_klb_out");
  #ifndef NO_AD_INITIALIZE
  L_Fproxy_klb_out.initialize();
  #endif
  D_Fproxy_knum_out.allocate("D_Fproxy_knum_out");
  #ifndef NO_AD_INITIALIZE
  D_Fproxy_knum_out.initialize();
  #endif
  D_Fproxy_klb_out.allocate("D_Fproxy_klb_out");
  #ifndef NO_AD_INITIALIZE
  D_Fproxy_klb_out.initialize();
  #endif
  F50_out.allocate("F50_out");
  #ifndef NO_AD_INITIALIZE
  F50_out.initialize();
  #endif
  SSB_F50_out.allocate("SSB_F50_out");
  #ifndef NO_AD_INITIALIZE
  SSB_F50_out.initialize();
  #endif
  B_F50_out.allocate("B_F50_out");
  #ifndef NO_AD_INITIALIZE
  B_F50_out.initialize();
  #endif
  R_F50_out.allocate("R_F50_out");
  #ifndef NO_AD_INITIALIZE
  R_F50_out.initialize();
  #endif
  L_F50_knum_out.allocate("L_F50_knum_out");
  #ifndef NO_AD_INITIALIZE
  L_F50_knum_out.initialize();
  #endif
  L_F50_klb_out.allocate("L_F50_klb_out");
  #ifndef NO_AD_INITIALIZE
  L_F50_klb_out.initialize();
  #endif
  D_F50_knum_out.allocate("D_F50_knum_out");
  #ifndef NO_AD_INITIALIZE
  D_F50_knum_out.initialize();
  #endif
  D_F50_klb_out.allocate("D_F50_klb_out");
  #ifndef NO_AD_INITIALIZE
  D_F50_klb_out.initialize();
  #endif
  F40_out.allocate("F40_out");
  #ifndef NO_AD_INITIALIZE
  F40_out.initialize();
  #endif
  SSB_F40_out.allocate("SSB_F40_out");
  #ifndef NO_AD_INITIALIZE
  SSB_F40_out.initialize();
  #endif
  B_F40_out.allocate("B_F40_out");
  #ifndef NO_AD_INITIALIZE
  B_F40_out.initialize();
  #endif
  R_F40_out.allocate("R_F40_out");
  #ifndef NO_AD_INITIALIZE
  R_F40_out.initialize();
  #endif
  L_F40_knum_out.allocate("L_F40_knum_out");
  #ifndef NO_AD_INITIALIZE
  L_F40_knum_out.initialize();
  #endif
  L_F40_klb_out.allocate("L_F40_klb_out");
  #ifndef NO_AD_INITIALIZE
  L_F40_klb_out.initialize();
  #endif
  D_F40_knum_out.allocate("D_F40_knum_out");
  #ifndef NO_AD_INITIALIZE
  D_F40_knum_out.initialize();
  #endif
  D_F40_klb_out.allocate("D_F40_klb_out");
  #ifndef NO_AD_INITIALIZE
  D_F40_klb_out.initialize();
  #endif
  F30_out.allocate("F30_out");
  #ifndef NO_AD_INITIALIZE
  F30_out.initialize();
  #endif
  SSB_F30_out.allocate("SSB_F30_out");
  #ifndef NO_AD_INITIALIZE
  SSB_F30_out.initialize();
  #endif
  B_F30_out.allocate("B_F30_out");
  #ifndef NO_AD_INITIALIZE
  B_F30_out.initialize();
  #endif
  R_F30_out.allocate("R_F30_out");
  #ifndef NO_AD_INITIALIZE
  R_F30_out.initialize();
  #endif
  L_F30_knum_out.allocate("L_F30_knum_out");
  #ifndef NO_AD_INITIALIZE
  L_F30_knum_out.initialize();
  #endif
  L_F30_klb_out.allocate("L_F30_klb_out");
  #ifndef NO_AD_INITIALIZE
  L_F30_klb_out.initialize();
  #endif
  D_F30_knum_out.allocate("D_F30_knum_out");
  #ifndef NO_AD_INITIALIZE
  D_F30_knum_out.initialize();
  #endif
  D_F30_klb_out.allocate("D_F30_klb_out");
  #ifndef NO_AD_INITIALIZE
  D_F30_klb_out.initialize();
  #endif
  rec_mean.allocate("rec_mean");
  #ifndef NO_AD_INITIALIZE
  rec_mean.initialize();
  #endif
  rec_mu.allocate("rec_mu");
  #ifndef NO_AD_INITIALIZE
  rec_mu.initialize();
  #endif
  rec_mu_temp.allocate("rec_mu_temp");
  #ifndef NO_AD_INITIALIZE
  rec_mu_temp.initialize();
  #endif
  N_age_msy.allocate(1,nages,"N_age_msy");
  #ifndef NO_AD_INITIALIZE
    N_age_msy.initialize();
  #endif
  N_age_msy_spawn.allocate(1,nages,"N_age_msy_spawn");
  #ifndef NO_AD_INITIALIZE
    N_age_msy_spawn.initialize();
  #endif
  L_age_msy.allocate(1,nages,"L_age_msy");
  #ifndef NO_AD_INITIALIZE
    L_age_msy.initialize();
  #endif
  Z_age_msy.allocate(1,nages,"Z_age_msy");
  #ifndef NO_AD_INITIALIZE
    Z_age_msy.initialize();
  #endif
  D_age_msy.allocate(1,nages,"D_age_msy");
  #ifndef NO_AD_INITIALIZE
    D_age_msy.initialize();
  #endif
  T_age_msy.allocate(1,nages,"T_age_msy");
  #ifndef NO_AD_INITIALIZE
    T_age_msy.initialize();
  #endif
  F_L_age_msy.allocate(1,nages,"F_L_age_msy");
  #ifndef NO_AD_INITIALIZE
    F_L_age_msy.initialize();
  #endif
  F_D_age_msy.allocate(1,nages,"F_D_age_msy");
  #ifndef NO_AD_INITIALIZE
    F_D_age_msy.initialize();
  #endif
  F_msy.allocate(1,n_iter_msy,"F_msy");
  #ifndef NO_AD_INITIALIZE
    F_msy.initialize();
  #endif
  spr_msy.allocate(1,n_iter_msy,"spr_msy");
  #ifndef NO_AD_INITIALIZE
    spr_msy.initialize();
  #endif
  R_eq.allocate(1,n_iter_msy,"R_eq");
  #ifndef NO_AD_INITIALIZE
    R_eq.initialize();
  #endif
  L_eq_klb.allocate(1,n_iter_msy,"L_eq_klb");
  #ifndef NO_AD_INITIALIZE
    L_eq_klb.initialize();
  #endif
  L_eq_knum.allocate(1,n_iter_msy,"L_eq_knum");
  #ifndef NO_AD_INITIALIZE
    L_eq_knum.initialize();
  #endif
  SSB_eq.allocate(1,n_iter_msy,"SSB_eq");
  #ifndef NO_AD_INITIALIZE
    SSB_eq.initialize();
  #endif
  B_eq.allocate(1,n_iter_msy,"B_eq");
  #ifndef NO_AD_INITIALIZE
    B_eq.initialize();
  #endif
  D_eq_klb.allocate(1,n_iter_msy,"D_eq_klb");
  #ifndef NO_AD_INITIALIZE
    D_eq_klb.initialize();
  #endif
  D_eq_knum.allocate(1,n_iter_msy,"D_eq_knum");
  #ifndef NO_AD_INITIALIZE
    D_eq_knum.initialize();
  #endif
  T_eq_klb.allocate(1,n_iter_msy,"T_eq_klb");
  #ifndef NO_AD_INITIALIZE
    T_eq_klb.initialize();
  #endif
  T_eq_knum.allocate(1,n_iter_msy,"T_eq_knum");
  #ifndef NO_AD_INITIALIZE
    T_eq_knum.initialize();
  #endif
  FdF_msy_n.allocate(styr,endyr,"FdF_msy_n");
  #ifndef NO_AD_INITIALIZE
    FdF_msy_n.initialize();
  #endif
  FdF_msy_lb.allocate(styr,endyr,"FdF_msy_lb");
  #ifndef NO_AD_INITIALIZE
    FdF_msy_lb.initialize();
  #endif
  FdFproxy.allocate(styr,endyr,"FdFproxy");
  #ifndef NO_AD_INITIALIZE
    FdFproxy.initialize();
  #endif
  FdF50.allocate(styr,endyr,"FdF50");
  #ifndef NO_AD_INITIALIZE
    FdF50.initialize();
  #endif
  FdF40.allocate(styr,endyr,"FdF40");
  #ifndef NO_AD_INITIALIZE
    FdF40.initialize();
  #endif
  FdF30.allocate(styr,endyr,"FdF30");
  #ifndef NO_AD_INITIALIZE
    FdF30.initialize();
  #endif
  SdSSB_msy_n.allocate(styr,endyr,"SdSSB_msy_n");
  #ifndef NO_AD_INITIALIZE
    SdSSB_msy_n.initialize();
  #endif
  SdSSB_msy_n_end.allocate("SdSSB_msy_n_end");
  #ifndef NO_AD_INITIALIZE
  SdSSB_msy_n_end.initialize();
  #endif
  SSBend_MSST_n.allocate("SSBend_MSST_n");
  #ifndef NO_AD_INITIALIZE
  SSBend_MSST_n.initialize();
  #endif
  SdSSB_msy_lb.allocate(styr,endyr,"SdSSB_msy_lb");
  #ifndef NO_AD_INITIALIZE
    SdSSB_msy_lb.initialize();
  #endif
  SdSSB_msy_lb_end.allocate("SdSSB_msy_lb_end");
  #ifndef NO_AD_INITIALIZE
  SdSSB_msy_lb_end.initialize();
  #endif
  SSBend_MSST_lb.allocate("SSBend_MSST_lb");
  #ifndef NO_AD_INITIALIZE
  SSBend_MSST_lb.initialize();
  #endif
  FdF_msy_n_end.allocate("FdF_msy_n_end");
  #ifndef NO_AD_INITIALIZE
  FdF_msy_n_end.initialize();
  #endif
  FdF_msy_n_end_mean.allocate("FdF_msy_n_end_mean");
  #ifndef NO_AD_INITIALIZE
  FdF_msy_n_end_mean.initialize();
  #endif
  FdF_msy_lb_end.allocate("FdF_msy_lb_end");
  #ifndef NO_AD_INITIALIZE
  FdF_msy_lb_end.initialize();
  #endif
  FdF_msy_lb_end_mean.allocate("FdF_msy_lb_end_mean");
  #ifndef NO_AD_INITIALIZE
  FdF_msy_lb_end_mean.initialize();
  #endif
  SdSSB_Fproxy.allocate(styr,endyr,"SdSSB_Fproxy");
  #ifndef NO_AD_INITIALIZE
    SdSSB_Fproxy.initialize();
  #endif
  Sdmsst_Fproxy.allocate(styr,endyr,"Sdmsst_Fproxy");
  #ifndef NO_AD_INITIALIZE
    Sdmsst_Fproxy.initialize();
  #endif
  SdSSB_Fproxy_end.allocate("SdSSB_Fproxy_end");
  #ifndef NO_AD_INITIALIZE
  SdSSB_Fproxy_end.initialize();
  #endif
  Sdmsst_Fproxy_end.allocate("Sdmsst_Fproxy_end");
  #ifndef NO_AD_INITIALIZE
  Sdmsst_Fproxy_end.initialize();
  #endif
  FdFproxy_end_mean.allocate("FdFproxy_end_mean");
  #ifndef NO_AD_INITIALIZE
  FdFproxy_end_mean.initialize();
  #endif
  SdSSB_F50.allocate(styr,endyr,"SdSSB_F50");
  #ifndef NO_AD_INITIALIZE
    SdSSB_F50.initialize();
  #endif
  Sdmsst_F50.allocate(styr,endyr,"Sdmsst_F50");
  #ifndef NO_AD_INITIALIZE
    Sdmsst_F50.initialize();
  #endif
  SdSSB_F50_end.allocate("SdSSB_F50_end");
  #ifndef NO_AD_INITIALIZE
  SdSSB_F50_end.initialize();
  #endif
  Sdmsst_F50_end.allocate("Sdmsst_F50_end");
  #ifndef NO_AD_INITIALIZE
  Sdmsst_F50_end.initialize();
  #endif
  FdF50_end_mean.allocate("FdF50_end_mean");
  #ifndef NO_AD_INITIALIZE
  FdF50_end_mean.initialize();
  #endif
  SdSSB_F40.allocate(styr,endyr,"SdSSB_F40");
  #ifndef NO_AD_INITIALIZE
    SdSSB_F40.initialize();
  #endif
  Sdmsst_F40.allocate(styr,endyr,"Sdmsst_F40");
  #ifndef NO_AD_INITIALIZE
    Sdmsst_F40.initialize();
  #endif
  SdSSB_F40_end.allocate("SdSSB_F40_end");
  #ifndef NO_AD_INITIALIZE
  SdSSB_F40_end.initialize();
  #endif
  Sdmsst_F40_end.allocate("Sdmsst_F40_end");
  #ifndef NO_AD_INITIALIZE
  Sdmsst_F40_end.initialize();
  #endif
  FdF40_end_mean.allocate("FdF40_end_mean");
  #ifndef NO_AD_INITIALIZE
  FdF40_end_mean.initialize();
  #endif
  SdSSB_F30.allocate(styr,endyr,"SdSSB_F30");
  #ifndef NO_AD_INITIALIZE
    SdSSB_F30.initialize();
  #endif
  Sdmsst_F30.allocate(styr,endyr,"Sdmsst_F30");
  #ifndef NO_AD_INITIALIZE
    Sdmsst_F30.initialize();
  #endif
  SdSSB_F30_end.allocate("SdSSB_F30_end");
  #ifndef NO_AD_INITIALIZE
  SdSSB_F30_end.initialize();
  #endif
  Sdmsst_F30_end.allocate("Sdmsst_F30_end");
  #ifndef NO_AD_INITIALIZE
  Sdmsst_F30_end.initialize();
  #endif
  FdF30_end_mean.allocate("FdF30_end_mean");
  #ifndef NO_AD_INITIALIZE
  FdF30_end_mean.initialize();
  #endif
  Fend_mean_temp.allocate("Fend_mean_temp");
  #ifndef NO_AD_INITIALIZE
  Fend_mean_temp.initialize();
  #endif
  Fend_mean.allocate("Fend_mean");
  #ifndef NO_AD_INITIALIZE
  Fend_mean.initialize();
  #endif
  L_age_Fproxy.allocate(1,nages,"L_age_Fproxy");
  #ifndef NO_AD_INITIALIZE
    L_age_Fproxy.initialize();
  #endif
  D_age_Fproxy.allocate(1,nages,"D_age_Fproxy");
  #ifndef NO_AD_INITIALIZE
    D_age_Fproxy.initialize();
  #endif
  wgt_wgted_L_klb.allocate(1,nages,"wgt_wgted_L_klb");
  #ifndef NO_AD_INITIALIZE
    wgt_wgted_L_klb.initialize();
  #endif
  wgt_wgted_D_klb.allocate(1,nages,"wgt_wgted_D_klb");
  #ifndef NO_AD_INITIALIZE
    wgt_wgted_D_klb.initialize();
  #endif
  wgt_wgted_L_denom.allocate("wgt_wgted_L_denom");
  #ifndef NO_AD_INITIALIZE
  wgt_wgted_L_denom.initialize();
  #endif
  wgt_wgted_D_denom.allocate("wgt_wgted_D_denom");
  #ifndef NO_AD_INITIALIZE
  wgt_wgted_D_denom.initialize();
  #endif
  iter_inc_msy.allocate("iter_inc_msy");
  #ifndef NO_AD_INITIALIZE
  iter_inc_msy.initialize();
  #endif
  M.allocate(1,nages,"M");
  #ifndef NO_AD_INITIALIZE
    M.initialize();
  #endif
  M_constant.allocate(M_constant_LO,M_constant_HI,M_constant_PH,"M_constant");
  M_constant_out.allocate(1,8,"M_constant_out");
  #ifndef NO_AD_INITIALIZE
    M_constant_out.initialize();
  #endif
  smsy2msstM.allocate("smsy2msstM");
  #ifndef NO_AD_INITIALIZE
  smsy2msstM.initialize();
  #endif
  F.allocate(styr,endyr,1,nages,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  Fsum.allocate(styr,endyr,"Fsum");
  #ifndef NO_AD_INITIALIZE
    Fsum.initialize();
  #endif
  Fapex.allocate(styr,endyr,"Fapex");
  #ifndef NO_AD_INITIALIZE
    Fapex.initialize();
  #endif
  Z.allocate(styr,endyr,1,nages,"Z");
  #ifndef NO_AD_INITIALIZE
    Z.initialize();
  #endif
  log_avg_F_L_cHL.allocate(log_avg_F_cHL_LO,log_avg_F_cHL_HI,log_avg_F_cHL_PH,"log_avg_F_L_cHL");
  out_log_avg_F_L_cHL.allocate(1,8,"out_log_avg_F_L_cHL");
  #ifndef NO_AD_INITIALIZE
    out_log_avg_F_L_cHL.initialize();
  #endif
  log_dev_F_L_cHL.allocate(styr_L_cHL,endyr_L_cHL,log_dev_F_cHL_LO,log_dev_F_cHL_HI,log_dev_F_cHL_PH,"log_dev_F_L_cHL");
  log_dev_F_L_out_cHL.allocate(styr_L_cHL,endyr_L_cHL,"log_dev_F_L_out_cHL");
  #ifndef NO_AD_INITIALIZE
    log_dev_F_L_out_cHL.initialize();
  #endif
  F_cHL.allocate(styr,endyr,1,nages,"F_cHL");
  #ifndef NO_AD_INITIALIZE
    F_cHL.initialize();
  #endif
  F_cHL_out.allocate(styr,endyr,"F_cHL_out");
  #ifndef NO_AD_INITIALIZE
    F_cHL_out.initialize();
  #endif
  log_dev_F_init_cHL.allocate("log_dev_F_init_cHL");
  #ifndef NO_AD_INITIALIZE
  log_dev_F_init_cHL.initialize();
  #endif
  log_dev_F_end_cHL.allocate("log_dev_F_end_cHL");
  #ifndef NO_AD_INITIALIZE
  log_dev_F_end_cHL.initialize();
  #endif
  log_avg_F_L_cPT.allocate(log_avg_F_cPT_LO,log_avg_F_cPT_HI,log_avg_F_cPT_PH,"log_avg_F_L_cPT");
  out_log_avg_F_L_cPT.allocate(1,8,"out_log_avg_F_L_cPT");
  #ifndef NO_AD_INITIALIZE
    out_log_avg_F_L_cPT.initialize();
  #endif
  log_dev_F_L_cPT.allocate(styr_L_cPT,endyr_L_cPT,log_dev_F_cPT_LO,log_dev_F_cPT_HI,log_dev_F_cPT_PH,"log_dev_F_L_cPT");
  log_dev_F_L_out_cPT.allocate(styr_L_cPT,endyr_L_cPT,"log_dev_F_L_out_cPT");
  #ifndef NO_AD_INITIALIZE
    log_dev_F_L_out_cPT.initialize();
  #endif
  F_cPT.allocate(styr,endyr,1,nages,"F_cPT");
  #ifndef NO_AD_INITIALIZE
    F_cPT.initialize();
  #endif
  F_cPT_out.allocate(styr,endyr,"F_cPT_out");
  #ifndef NO_AD_INITIALIZE
    F_cPT_out.initialize();
  #endif
  log_dev_F_init_cPT.allocate("log_dev_F_init_cPT");
  #ifndef NO_AD_INITIALIZE
  log_dev_F_init_cPT.initialize();
  #endif
  log_dev_F_end_cPT.allocate("log_dev_F_end_cPT");
  #ifndef NO_AD_INITIALIZE
  log_dev_F_end_cPT.initialize();
  #endif
  log_avg_F_L_cTW.allocate(log_avg_F_cTW_LO,log_avg_F_cTW_HI,log_avg_F_cTW_PH,"log_avg_F_L_cTW");
  out_log_avg_F_L_cTW.allocate(1,8,"out_log_avg_F_L_cTW");
  #ifndef NO_AD_INITIALIZE
    out_log_avg_F_L_cTW.initialize();
  #endif
  log_dev_F_L_cTW.allocate(styr_L_cTW,endyr_L_cTW,log_dev_F_cTW_LO,log_dev_F_cTW_HI,log_dev_F_cTW_PH,"log_dev_F_L_cTW");
  log_dev_F_L_out_cTW.allocate(styr_L_cTW,endyr_L_cTW,"log_dev_F_L_out_cTW");
  #ifndef NO_AD_INITIALIZE
    log_dev_F_L_out_cTW.initialize();
  #endif
  F_cTW.allocate(styr,endyr,1,nages,"F_cTW");
  #ifndef NO_AD_INITIALIZE
    F_cTW.initialize();
  #endif
  F_cTW_out.allocate(styr,endyr,"F_cTW_out");
  #ifndef NO_AD_INITIALIZE
    F_cTW_out.initialize();
  #endif
  log_dev_F_init_cTW.allocate("log_dev_F_init_cTW");
  #ifndef NO_AD_INITIALIZE
  log_dev_F_init_cTW.initialize();
  #endif
  log_dev_F_end_cTW.allocate("log_dev_F_end_cTW");
  #ifndef NO_AD_INITIALIZE
  log_dev_F_end_cTW.initialize();
  #endif
  log_avg_F_L_rHB.allocate(log_avg_F_rHB_LO,log_avg_F_rHB_HI,log_avg_F_rHB_PH,"log_avg_F_L_rHB");
  out_log_avg_F_L_rHB.allocate(1,8,"out_log_avg_F_L_rHB");
  #ifndef NO_AD_INITIALIZE
    out_log_avg_F_L_rHB.initialize();
  #endif
  log_dev_F_L_rHB.allocate(styr_L_rHB,endyr_L_rHB,log_dev_F_rHB_LO,log_dev_F_rHB_HI,log_dev_F_rHB_PH,"log_dev_F_L_rHB");
  log_dev_F_L_out_rHB.allocate(styr_L_rHB,endyr_L_rHB,"log_dev_F_L_out_rHB");
  #ifndef NO_AD_INITIALIZE
    log_dev_F_L_out_rHB.initialize();
  #endif
  F_rHB.allocate(styr,endyr,1,nages,"F_rHB");
  #ifndef NO_AD_INITIALIZE
    F_rHB.initialize();
  #endif
  F_rHB_out.allocate(styr,endyr,"F_rHB_out");
  #ifndef NO_AD_INITIALIZE
    F_rHB_out.initialize();
  #endif
  log_dev_F_init_rHB.allocate("log_dev_F_init_rHB");
  #ifndef NO_AD_INITIALIZE
  log_dev_F_init_rHB.initialize();
  #endif
  log_dev_F_end_rHB.allocate("log_dev_F_end_rHB");
  #ifndef NO_AD_INITIALIZE
  log_dev_F_end_rHB.initialize();
  #endif
  log_avg_F_L_rGN.allocate(log_avg_F_rGN_LO,log_avg_F_rGN_HI,log_avg_F_rGN_PH,"log_avg_F_L_rGN");
  out_log_avg_F_L_rGN.allocate(1,8,"out_log_avg_F_L_rGN");
  #ifndef NO_AD_INITIALIZE
    out_log_avg_F_L_rGN.initialize();
  #endif
  log_dev_F_L_rGN.allocate(styr_L_rGN,endyr_L_rGN,log_dev_F_rGN_LO,log_dev_F_rGN_HI,log_dev_F_rGN_PH,"log_dev_F_L_rGN");
  log_dev_F_L_out_rGN.allocate(styr_L_rGN,endyr_L_rGN,"log_dev_F_L_out_rGN");
  #ifndef NO_AD_INITIALIZE
    log_dev_F_L_out_rGN.initialize();
  #endif
  F_rGN.allocate(styr,endyr,1,nages,"F_rGN");
  #ifndef NO_AD_INITIALIZE
    F_rGN.initialize();
  #endif
  F_rGN_out.allocate(styr,endyr,"F_rGN_out");
  #ifndef NO_AD_INITIALIZE
    F_rGN_out.initialize();
  #endif
  log_dev_F_init_rGN.allocate("log_dev_F_init_rGN");
  #ifndef NO_AD_INITIALIZE
  log_dev_F_init_rGN.initialize();
  #endif
  log_dev_F_end_rGN.allocate("log_dev_F_end_rGN");
  #ifndef NO_AD_INITIALIZE
  log_dev_F_end_rGN.initialize();
  #endif
  log_avg_F_L_rIN.allocate(log_avg_F_rIN_LO,log_avg_F_rIN_HI,log_avg_F_rIN_PH,"log_avg_F_L_rIN");
  out_log_avg_F_L_rIN.allocate(1,8,"out_log_avg_F_L_rIN");
  #ifndef NO_AD_INITIALIZE
    out_log_avg_F_L_rIN.initialize();
  #endif
  log_dev_F_L_rIN.allocate(styr_L_rIN,endyr_L_rIN,log_dev_F_rIN_LO,log_dev_F_rIN_HI,log_dev_F_rIN_PH,"log_dev_F_L_rIN");
  log_dev_F_L_out_rIN.allocate(styr_L_rIN,endyr_L_rIN,"log_dev_F_L_out_rIN");
  #ifndef NO_AD_INITIALIZE
    log_dev_F_L_out_rIN.initialize();
  #endif
  F_rIN.allocate(styr,endyr,1,nages,"F_rIN");
  #ifndef NO_AD_INITIALIZE
    F_rIN.initialize();
  #endif
  F_rIN_out.allocate(styr,endyr,"F_rIN_out");
  #ifndef NO_AD_INITIALIZE
    F_rIN_out.initialize();
  #endif
  log_dev_F_init_rIN.allocate("log_dev_F_init_rIN");
  #ifndef NO_AD_INITIALIZE
  log_dev_F_init_rIN.initialize();
  #endif
  log_dev_F_end_rIN.allocate("log_dev_F_end_rIN");
  #ifndef NO_AD_INITIALIZE
  log_dev_F_end_rIN.initialize();
  #endif
  F_init_ratio.allocate(F_init_ratio_LO,F_init_ratio_HI,F_init_ratio_PH,"F_init_ratio");
  F_init_ratio_out.allocate(1,8,"F_init_ratio_out");
  #ifndef NO_AD_INITIALIZE
    F_init_ratio_out.initialize();
  #endif
  log_avg_F_D_cGD.allocate(LO_log_avg_F_D_cGD,HI_log_avg_F_D_cGD,PH_log_avg_F_D_cGD,"log_avg_F_D_cGD");
  out_log_avg_F_D_cGD.allocate(1,8,"out_log_avg_F_D_cGD");
  #ifndef NO_AD_INITIALIZE
    out_log_avg_F_D_cGD.initialize();
  #endif
  log_dev_F_D_cGD.allocate(styr_D_cGD,endyr_D_cGD,log_dev_F_cGD_LO,log_dev_F_cGD_HI,log_dev_F_cGD_PH,"log_dev_F_D_cGD");
  log_dev_F_D_out_cGD.allocate(styr_D_cGD,endyr_D_cGD,"log_dev_F_D_out_cGD");
  #ifndef NO_AD_INITIALIZE
    log_dev_F_D_out_cGD.initialize();
  #endif
  F_cGD.allocate(styr,endyr,1,nages,"F_cGD");
  #ifndef NO_AD_INITIALIZE
    F_cGD.initialize();
  #endif
  F_cGD_out.allocate(styr,endyr,"F_cGD_out");
  #ifndef NO_AD_INITIALIZE
    F_cGD_out.initialize();
  #endif
  log_dev_F_cGD2.allocate("log_dev_F_cGD2");
  #ifndef NO_AD_INITIALIZE
  log_dev_F_cGD2.initialize();
  #endif
  log_dev_F_end_cGD.allocate("log_dev_F_end_cGD");
  #ifndef NO_AD_INITIALIZE
  log_dev_F_end_cGD.initialize();
  #endif
  log_avg_F_D_rHD.allocate(LO_log_avg_F_D_rHD,HI_log_avg_F_D_rHD,PH_log_avg_F_D_rHD,"log_avg_F_D_rHD");
  out_log_avg_F_D_rHD.allocate(1,8,"out_log_avg_F_D_rHD");
  #ifndef NO_AD_INITIALIZE
    out_log_avg_F_D_rHD.initialize();
  #endif
  log_dev_F_D_rHD.allocate(styr_D_rHD,endyr_D_rHD,log_dev_F_rHD_LO,log_dev_F_rHD_HI,log_dev_F_rHD_PH,"log_dev_F_D_rHD");
  log_dev_F_D_out_rHD.allocate(styr_D_rHD,endyr_D_rHD,"log_dev_F_D_out_rHD");
  #ifndef NO_AD_INITIALIZE
    log_dev_F_D_out_rHD.initialize();
  #endif
  F_rHD.allocate(styr,endyr,1,nages,"F_rHD");
  #ifndef NO_AD_INITIALIZE
    F_rHD.initialize();
  #endif
  F_rHD_out.allocate(styr,endyr,"F_rHD_out");
  #ifndef NO_AD_INITIALIZE
    F_rHD_out.initialize();
  #endif
  log_dev_F_end_rHD.allocate("log_dev_F_end_rHD");
  #ifndef NO_AD_INITIALIZE
  log_dev_F_end_rHD.initialize();
  #endif
  log_avg_F_D_rGD.allocate(LO_log_avg_F_D_rGD,HI_log_avg_F_D_rGD,PH_log_avg_F_D_rGD,"log_avg_F_D_rGD");
  out_log_avg_F_D_rGD.allocate(1,8,"out_log_avg_F_D_rGD");
  #ifndef NO_AD_INITIALIZE
    out_log_avg_F_D_rGD.initialize();
  #endif
  log_dev_F_D_rGD.allocate(styr_D_rGD,endyr_D_rGD,log_dev_F_rGD_LO,log_dev_F_rGD_HI,log_dev_F_rGD_PH,"log_dev_F_D_rGD");
  log_dev_F_D_out_rGD.allocate(styr_D_rGD,endyr_D_rGD,"log_dev_F_D_out_rGD");
  #ifndef NO_AD_INITIALIZE
    log_dev_F_D_out_rGD.initialize();
  #endif
  F_rGD.allocate(styr,endyr,1,nages,"F_rGD");
  #ifndef NO_AD_INITIALIZE
    F_rGD.initialize();
  #endif
  F_rGD_out.allocate(styr,endyr,"F_rGD_out");
  #ifndef NO_AD_INITIALIZE
    F_rGD_out.initialize();
  #endif
  log_dev_F_init_rGD.allocate("log_dev_F_init_rGD");
  #ifndef NO_AD_INITIALIZE
  log_dev_F_init_rGD.initialize();
  #endif
  log_dev_F_end_rGD.allocate("log_dev_F_end_rGD");
  #ifndef NO_AD_INITIALIZE
  log_dev_F_end_rGD.initialize();
  #endif
  log_avg_F_D_rID.allocate(LO_log_avg_F_D_rID,HI_log_avg_F_D_rID,PH_log_avg_F_D_rID,"log_avg_F_D_rID");
  out_log_avg_F_D_rID.allocate(1,8,"out_log_avg_F_D_rID");
  #ifndef NO_AD_INITIALIZE
    out_log_avg_F_D_rID.initialize();
  #endif
  log_dev_F_D_rID.allocate(styr_D_rID,endyr_D_rID,log_dev_F_rID_LO,log_dev_F_rID_HI,log_dev_F_rID_PH,"log_dev_F_D_rID");
  log_dev_F_D_out_rID.allocate(styr_D_rID,endyr_D_rID,"log_dev_F_D_out_rID");
  #ifndef NO_AD_INITIALIZE
    log_dev_F_D_out_rID.initialize();
  #endif
  F_rID.allocate(styr,endyr,1,nages,"F_rID");
  #ifndef NO_AD_INITIALIZE
    F_rID.initialize();
  #endif
  F_rID_out.allocate(styr,endyr,"F_rID_out");
  #ifndef NO_AD_INITIALIZE
    F_rID_out.initialize();
  #endif
  log_dev_F_init_rID.allocate("log_dev_F_init_rID");
  #ifndef NO_AD_INITIALIZE
  log_dev_F_init_rID.initialize();
  #endif
  log_dev_F_end_rID.allocate("log_dev_F_end_rID");
  #ifndef NO_AD_INITIALIZE
  log_dev_F_end_rID.initialize();
  #endif
  Dmort_cGD.allocate("Dmort_cGD");
  #ifndef NO_AD_INITIALIZE
  Dmort_cGD.initialize();
  #endif
  Dmort_rHD.allocate("Dmort_rHD");
  #ifndef NO_AD_INITIALIZE
  Dmort_rHD.initialize();
  #endif
  Dmort_rGD.allocate("Dmort_rGD");
  #ifndef NO_AD_INITIALIZE
  Dmort_rGD.initialize();
  #endif
  Dmort_rID.allocate("Dmort_rID");
  #ifndef NO_AD_INITIALIZE
  Dmort_rID.initialize();
  #endif
  N_age_spr.allocate(1,nages,"N_age_spr");
  #ifndef NO_AD_INITIALIZE
    N_age_spr.initialize();
  #endif
  N_age_spr_spawn.allocate(1,nages,"N_age_spr_spawn");
  #ifndef NO_AD_INITIALIZE
    N_age_spr_spawn.initialize();
  #endif
  L_age_spr.allocate(1,nages,"L_age_spr");
  #ifndef NO_AD_INITIALIZE
    L_age_spr.initialize();
  #endif
  T_age_spr.allocate(1,nages,"T_age_spr");
  #ifndef NO_AD_INITIALIZE
    T_age_spr.initialize();
  #endif
  Z_age_spr.allocate(1,nages,"Z_age_spr");
  #ifndef NO_AD_INITIALIZE
    Z_age_spr.initialize();
  #endif
  spr_static.allocate(styr,endyr,"spr_static");
  #ifndef NO_AD_INITIALIZE
    spr_static.initialize();
  #endif
  F_L_age_spr.allocate(1,nages,"F_L_age_spr");
  #ifndef NO_AD_INITIALIZE
    F_L_age_spr.initialize();
  #endif
  F_T_age_spr.allocate(1,nages,"F_T_age_spr");
  #ifndef NO_AD_INITIALIZE
    F_T_age_spr.initialize();
  #endif
  F_D_age_spr.allocate(1,nages,"F_D_age_spr");
  #ifndef NO_AD_INITIALIZE
    F_D_age_spr.initialize();
  #endif
  F_spr.allocate(1,n_iter_spr,"F_spr");
  #ifndef NO_AD_INITIALIZE
    F_spr.initialize();
  #endif
  spr_spr.allocate(1,n_iter_spr,"spr_spr");
  #ifndef NO_AD_INITIALIZE
    spr_spr.initialize();
  #endif
  spr_ratio.allocate(1,n_iter_spr,"spr_ratio");
  #ifndef NO_AD_INITIALIZE
    spr_ratio.initialize();
  #endif
  L_spr.allocate(1,n_iter_spr,"L_spr");
  #ifndef NO_AD_INITIALIZE
    L_spr.initialize();
  #endif
  T_spr.allocate(1,n_iter_spr,"T_spr");
  #ifndef NO_AD_INITIALIZE
    T_spr.initialize();
  #endif
  N_spr_F0.allocate(1,nages,"N_spr_F0");
  #ifndef NO_AD_INITIALIZE
    N_spr_F0.initialize();
  #endif
  N_bpr_F0.allocate(1,nages,"N_bpr_F0");
  #ifndef NO_AD_INITIALIZE
    N_bpr_F0.initialize();
  #endif
  N_spr_initial.allocate(1,nages,"N_spr_initial");
  #ifndef NO_AD_INITIALIZE
    N_spr_initial.initialize();
  #endif
  N_initial_eq.allocate(1,nages,"N_initial_eq");
  #ifndef NO_AD_INITIALIZE
    N_initial_eq.initialize();
  #endif
  F_initial.allocate(1,nages,"F_initial");
  #ifndef NO_AD_INITIALIZE
    F_initial.initialize();
  #endif
  Z_initial.allocate(1,nages,"Z_initial");
  #ifndef NO_AD_INITIALIZE
    Z_initial.initialize();
  #endif
  spr_initial.allocate("spr_initial");
  #ifndef NO_AD_INITIALIZE
  spr_initial.initialize();
  #endif
  spr_F0.allocate("spr_F0");
  #ifndef NO_AD_INITIALIZE
  spr_F0.initialize();
  #endif
  bpr_F0.allocate("bpr_F0");
  #ifndef NO_AD_INITIALIZE
  bpr_F0.initialize();
  #endif
  iter_inc_spr.allocate("iter_inc_spr");
  #ifndef NO_AD_INITIALIZE
  iter_inc_spr.initialize();
  #endif
  sdnr_lc_sBT.allocate("sdnr_lc_sBT");
  #ifndef NO_AD_INITIALIZE
  sdnr_lc_sBT.initialize();
  #endif
  sdnr_lc_sTV.allocate("sdnr_lc_sTV");
  #ifndef NO_AD_INITIALIZE
  sdnr_lc_sTV.initialize();
  #endif
  sdnr_lc_cHL.allocate("sdnr_lc_cHL");
  #ifndef NO_AD_INITIALIZE
  sdnr_lc_cHL.initialize();
  #endif
  sdnr_lc_cPT.allocate("sdnr_lc_cPT");
  #ifndef NO_AD_INITIALIZE
  sdnr_lc_cPT.initialize();
  #endif
  sdnr_lc_rHB.allocate("sdnr_lc_rHB");
  #ifndef NO_AD_INITIALIZE
  sdnr_lc_rHB.initialize();
  #endif
  sdnr_lc_rHD.allocate("sdnr_lc_rHD");
  #ifndef NO_AD_INITIALIZE
  sdnr_lc_rHD.initialize();
  #endif
  sdnr_lc_rGN.allocate("sdnr_lc_rGN");
  #ifndef NO_AD_INITIALIZE
  sdnr_lc_rGN.initialize();
  #endif
  sdnr_ac_sBT.allocate("sdnr_ac_sBT");
  #ifndef NO_AD_INITIALIZE
  sdnr_ac_sBT.initialize();
  #endif
  sdnr_ac_sTV.allocate("sdnr_ac_sTV");
  #ifndef NO_AD_INITIALIZE
  sdnr_ac_sTV.initialize();
  #endif
  sdnr_ac_cHL.allocate("sdnr_ac_cHL");
  #ifndef NO_AD_INITIALIZE
  sdnr_ac_cHL.initialize();
  #endif
  sdnr_ac_cPT.allocate("sdnr_ac_cPT");
  #ifndef NO_AD_INITIALIZE
  sdnr_ac_cPT.initialize();
  #endif
  sdnr_ac_rHB.allocate("sdnr_ac_rHB");
  #ifndef NO_AD_INITIALIZE
  sdnr_ac_rHB.initialize();
  #endif
  sdnr_I_sBT.allocate("sdnr_I_sBT");
  #ifndef NO_AD_INITIALIZE
  sdnr_I_sBT.initialize();
  #endif
  sdnr_I_sTV.allocate("sdnr_I_sTV");
  #ifndef NO_AD_INITIALIZE
  sdnr_I_sTV.initialize();
  #endif
  sdnr_I_cHL.allocate("sdnr_I_cHL");
  #ifndef NO_AD_INITIALIZE
  sdnr_I_cHL.initialize();
  #endif
  sdnr_I_rHB.allocate("sdnr_I_rHB");
  #ifndef NO_AD_INITIALIZE
  sdnr_I_rHB.initialize();
  #endif
  w_L.allocate("w_L");
  #ifndef NO_AD_INITIALIZE
  w_L.initialize();
  #endif
  w_D.allocate("w_D");
  #ifndef NO_AD_INITIALIZE
  w_D.initialize();
  #endif
  w_lenc_sBT.allocate("w_lenc_sBT");
  #ifndef NO_AD_INITIALIZE
  w_lenc_sBT.initialize();
  #endif
  w_lenc_cHL.allocate("w_lenc_cHL");
  #ifndef NO_AD_INITIALIZE
  w_lenc_cHL.initialize();
  #endif
  w_lenc_cPT.allocate("w_lenc_cPT");
  #ifndef NO_AD_INITIALIZE
  w_lenc_cPT.initialize();
  #endif
  w_lenc_rHB.allocate("w_lenc_rHB");
  #ifndef NO_AD_INITIALIZE
  w_lenc_rHB.initialize();
  #endif
  w_lenc_rHD.allocate("w_lenc_rHD");
  #ifndef NO_AD_INITIALIZE
  w_lenc_rHD.initialize();
  #endif
  w_lenc_rGN.allocate("w_lenc_rGN");
  #ifndef NO_AD_INITIALIZE
  w_lenc_rGN.initialize();
  #endif
  w_agec_sBT.allocate("w_agec_sBT");
  #ifndef NO_AD_INITIALIZE
  w_agec_sBT.initialize();
  #endif
  w_agec_sTV.allocate("w_agec_sTV");
  #ifndef NO_AD_INITIALIZE
  w_agec_sTV.initialize();
  #endif
  w_agec_cHL.allocate("w_agec_cHL");
  #ifndef NO_AD_INITIALIZE
  w_agec_cHL.initialize();
  #endif
  w_agec_cPT.allocate("w_agec_cPT");
  #ifndef NO_AD_INITIALIZE
  w_agec_cPT.initialize();
  #endif
  w_agec_rHB.allocate("w_agec_rHB");
  #ifndef NO_AD_INITIALIZE
  w_agec_rHB.initialize();
  #endif
  w_cpue_sBT.allocate("w_cpue_sBT");
  #ifndef NO_AD_INITIALIZE
  w_cpue_sBT.initialize();
  #endif
  w_cpue_sTV.allocate("w_cpue_sTV");
  #ifndef NO_AD_INITIALIZE
  w_cpue_sTV.initialize();
  #endif
  w_cpue_cHL.allocate("w_cpue_cHL");
  #ifndef NO_AD_INITIALIZE
  w_cpue_cHL.initialize();
  #endif
  w_cpue_rHB.allocate("w_cpue_rHB");
  #ifndef NO_AD_INITIALIZE
  w_cpue_rHB.initialize();
  #endif
  w_rec.allocate("w_rec");
  #ifndef NO_AD_INITIALIZE
  w_rec.initialize();
  #endif
  w_rec_early.allocate("w_rec_early");
  #ifndef NO_AD_INITIALIZE
  w_rec_early.initialize();
  #endif
  w_rec_end.allocate("w_rec_end");
  #ifndef NO_AD_INITIALIZE
  w_rec_end.initialize();
  #endif
  f_cpue_sBT.allocate("f_cpue_sBT");
  #ifndef NO_AD_INITIALIZE
  f_cpue_sBT.initialize();
  #endif
  f_cpue_sTV.allocate("f_cpue_sTV");
  #ifndef NO_AD_INITIALIZE
  f_cpue_sTV.initialize();
  #endif
  f_cpue_cHL.allocate("f_cpue_cHL");
  #ifndef NO_AD_INITIALIZE
  f_cpue_cHL.initialize();
  #endif
  f_cpue_rHB.allocate("f_cpue_rHB");
  #ifndef NO_AD_INITIALIZE
  f_cpue_rHB.initialize();
  #endif
  f_L_cHL.allocate("f_L_cHL");
  #ifndef NO_AD_INITIALIZE
  f_L_cHL.initialize();
  #endif
  f_L_cPT.allocate("f_L_cPT");
  #ifndef NO_AD_INITIALIZE
  f_L_cPT.initialize();
  #endif
  f_L_cTW.allocate("f_L_cTW");
  #ifndef NO_AD_INITIALIZE
  f_L_cTW.initialize();
  #endif
  f_L_rHB.allocate("f_L_rHB");
  #ifndef NO_AD_INITIALIZE
  f_L_rHB.initialize();
  #endif
  f_L_rGN.allocate("f_L_rGN");
  #ifndef NO_AD_INITIALIZE
  f_L_rGN.initialize();
  #endif
  f_L_rIN.allocate("f_L_rIN");
  #ifndef NO_AD_INITIALIZE
  f_L_rIN.initialize();
  #endif
  f_D_cGD.allocate("f_D_cGD");
  #ifndef NO_AD_INITIALIZE
  f_D_cGD.initialize();
  #endif
  f_D_rHD.allocate("f_D_rHD");
  #ifndef NO_AD_INITIALIZE
  f_D_rHD.initialize();
  #endif
  f_D_rGD.allocate("f_D_rGD");
  #ifndef NO_AD_INITIALIZE
  f_D_rGD.initialize();
  #endif
  f_D_rID.allocate("f_D_rID");
  #ifndef NO_AD_INITIALIZE
  f_D_rID.initialize();
  #endif
  f_lenc_sBT.allocate("f_lenc_sBT");
  #ifndef NO_AD_INITIALIZE
  f_lenc_sBT.initialize();
  #endif
  f_lenc_cHL.allocate("f_lenc_cHL");
  #ifndef NO_AD_INITIALIZE
  f_lenc_cHL.initialize();
  #endif
  f_lenc_cPT.allocate("f_lenc_cPT");
  #ifndef NO_AD_INITIALIZE
  f_lenc_cPT.initialize();
  #endif
  f_lenc_rHB.allocate("f_lenc_rHB");
  #ifndef NO_AD_INITIALIZE
  f_lenc_rHB.initialize();
  #endif
  f_lenc_rHD.allocate("f_lenc_rHD");
  #ifndef NO_AD_INITIALIZE
  f_lenc_rHD.initialize();
  #endif
  f_lenc_rGN.allocate("f_lenc_rGN");
  #ifndef NO_AD_INITIALIZE
  f_lenc_rGN.initialize();
  #endif
  f_agec_sBT.allocate("f_agec_sBT");
  #ifndef NO_AD_INITIALIZE
  f_agec_sBT.initialize();
  #endif
  f_agec_sTV.allocate("f_agec_sTV");
  #ifndef NO_AD_INITIALIZE
  f_agec_sTV.initialize();
  #endif
  f_agec_cHL.allocate("f_agec_cHL");
  #ifndef NO_AD_INITIALIZE
  f_agec_cHL.initialize();
  #endif
  f_agec_cPT.allocate("f_agec_cPT");
  #ifndef NO_AD_INITIALIZE
  f_agec_cPT.initialize();
  #endif
  f_agec_rHB.allocate("f_agec_rHB");
  #ifndef NO_AD_INITIALIZE
  f_agec_rHB.initialize();
  #endif
  f_agec_rGN.allocate("f_agec_rGN");
  #ifndef NO_AD_INITIALIZE
  f_agec_rGN.initialize();
  #endif
  f_cpue_RW_cHL.allocate("f_cpue_RW_cHL");
  #ifndef NO_AD_INITIALIZE
  f_cpue_RW_cHL.initialize();
  #endif
  f_cpue_RW_rHB.allocate("f_cpue_RW_rHB");
  #ifndef NO_AD_INITIALIZE
  f_cpue_RW_rHB.initialize();
  #endif
  f_cpue_RW_rHD.allocate("f_cpue_RW_rHD");
  #ifndef NO_AD_INITIALIZE
  f_cpue_RW_rHD.initialize();
  #endif
  f_rec_dev.allocate("f_rec_dev");
  #ifndef NO_AD_INITIALIZE
  f_rec_dev.initialize();
  #endif
  f_rec_dev_early.allocate("f_rec_dev_early");
  #ifndef NO_AD_INITIALIZE
  f_rec_dev_early.initialize();
  #endif
  f_rec_dev_end.allocate("f_rec_dev_end");
  #ifndef NO_AD_INITIALIZE
  f_rec_dev_end.initialize();
  #endif
  f_priors.allocate("f_priors");
  #ifndef NO_AD_INITIALIZE
  f_priors.initialize();
  #endif
  fval.allocate("fval");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  fval_data.allocate("fval_data");
  #ifndef NO_AD_INITIALIZE
  fval_data.initialize();
  #endif
  grad_max.allocate("grad_max");
  #ifndef NO_AD_INITIALIZE
  grad_max.initialize();
  #endif
  denom.allocate("denom");
  #ifndef NO_AD_INITIALIZE
  denom.initialize();
  #endif
  numer.allocate("numer");
  #ifndef NO_AD_INITIALIZE
  numer.initialize();
  #endif
  F_reg_proj.allocate("F_reg_proj");
  #ifndef NO_AD_INITIALIZE
  F_reg_proj.initialize();
  #endif
  F_proj.allocate(styr_proj,endyr_proj,"F_proj");
  #ifndef NO_AD_INITIALIZE
    F_proj.initialize();
  #endif
  L_knum_proj.allocate(styr_proj,endyr_proj,"L_knum_proj");
  #ifndef NO_AD_INITIALIZE
    L_knum_proj.initialize();
  #endif
  L_klb_proj.allocate(styr_proj,endyr_proj,"L_klb_proj");
  #ifndef NO_AD_INITIALIZE
    L_klb_proj.initialize();
  #endif
  B_proj.allocate(styr_proj,endyr_proj,"B_proj");
  #ifndef NO_AD_INITIALIZE
    B_proj.initialize();
  #endif
  SSB_proj.allocate(styr_proj,endyr_proj,"SSB_proj");
  #ifndef NO_AD_INITIALIZE
    SSB_proj.initialize();
  #endif
  R_proj.allocate(styr_proj,endyr_proj,"R_proj");
  #ifndef NO_AD_INITIALIZE
    R_proj.initialize();
  #endif
  FL_age_proj.allocate(1,nages,"FL_age_proj");
  #ifndef NO_AD_INITIALIZE
    FL_age_proj.initialize();
  #endif
  FD_age_proj.allocate(1,nages,"FD_age_proj");
  #ifndef NO_AD_INITIALIZE
    FD_age_proj.initialize();
  #endif
  N_proj.allocate(styr_proj,endyr_proj,1,nages,"N_proj");
  #ifndef NO_AD_INITIALIZE
    N_proj.initialize();
  #endif
  N_spawn_proj.allocate(styr_proj,endyr_proj,1,nages,"N_spawn_proj");
  #ifndef NO_AD_INITIALIZE
    N_spawn_proj.initialize();
  #endif
  Z_proj.allocate(styr_proj,endyr_proj,1,nages,"Z_proj");
  #ifndef NO_AD_INITIALIZE
    Z_proj.initialize();
  #endif
  L_age_proj.allocate(styr_proj,endyr_proj,1,nages,"L_age_proj");
  #ifndef NO_AD_INITIALIZE
    L_age_proj.initialize();
  #endif
  D_age_proj.allocate(styr_proj,endyr_proj,1,nages,"D_age_proj");
  #ifndef NO_AD_INITIALIZE
    D_age_proj.initialize();
  #endif
  D_knum_proj.allocate(styr_proj,endyr_proj,"D_knum_proj");
  #ifndef NO_AD_INITIALIZE
    D_knum_proj.initialize();
  #endif
  D_klb_proj.allocate(styr_proj,endyr_proj,"D_klb_proj");
  #ifndef NO_AD_INITIALIZE
    D_klb_proj.initialize();
  #endif
}

void model_parameters::set_runtime(void)
{
  dvector temp1("{10000, 20000, 30000, 50000, 100000, 100000, 100000;}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
  dvector temp("{1e-2, 1e-2,1e-2, 1e-4;}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
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
  F_msy(1)=0.0;  
    for (ff=2;ff<=n_iter_msy;ff++) {F_msy(ff)=F_msy(ff-1)+iter_inc_msy;}
  F_spr(1)=0.0;  
    for (ff=2;ff<=n_iter_spr;ff++){F_spr(ff)=F_spr(ff-1)+iter_inc_spr;}
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
  
 
}

void model_parameters::userfunction(void)
{
  fval =0.0;
 //cout<<"start"<<endl;
 // get_M_at_age(); 
 get_length_weight_at_age(); 
 //cout << "got length, weight, fecundity transitions" <<endl;
 get_reprod();
 get_length_at_age_dist(); 
 get_weight_at_age_landings();
 get_spr_F0();
 get_selectivity(); 
 get_mortality(); 
 get_bias_corr(); 
 get_numbers_at_age(); 
 get_landings_numbers();
 get_landings_wgt();
 get_dead_discards();
 //cout << "got discards" << endl;
 get_catchability_fcns();
 //cout << "got catchability_fcns" << endl;
 get_indices();
 //cout << "got indices" << endl;
 get_length_comps();
 //cout<< "got length comps"<< endl;
 get_age_comps();
 evaluate_objective_function();
 //cout << "objective function calculations complete" << endl;
   // cum_surv_1plus=nages_M_sum*M_constant;  
   // M=M_lorenzen*(cum_surv_1plus/sum(M_lorenzen(min_M_age,nages)));
#ifdef DEBUG
  std::cout << "DEBUG: Total gradient stack used is " << gradient_structure::get()->GRAD_STACK1->total() << " out of " << gradient_structure::get_GRADSTACK_BUFFER_SIZE() << std::endl;;
  std::cout << "DEBUG: Total dvariable address used is " << gradient_structure::get()->GRAD_LIST->total_addresses() << " out of " << gradient_structure::get_MAX_DLINKS() << std::endl;;
  std::cout << "DEBUG: Total dvariable address used is " << gradient_structure::get()->ARR_LIST1->get_max_last_offset() << " out of " << gradient_structure::get_ARRAY_MEMBLOCK_SIZE() << std::endl;;
#endif
}

void model_parameters::get_length_weight_at_age(void)
{
  //compute mean length (mm) and weight (whole) at age
    meanlen_TL=Linf*(1.0-mfexp(-K*(agebins-t0+0.5)));    //total length in mm
    wgt_g=wgtpar_a*pow(meanlen_TL,wgtpar_b);             //wgt in grams 
    wgt_kg=g2kg*wgt_g;                                   //wgt in kilograms 
    wgt_mt=g2mt*wgt_g;                                   //mt of whole wgt: g2mt converts g to mt
    wgt_klb=mt2klb*wgt_mt;                               //1000 lb of whole wgt
    wgt_lb=mt2lb*wgt_mt;                                 //1000 lb of whole wgt
    fecundity=fecpar_batches*mfexp(fecpar_a+wgt_g*fecpar_b)/fecpar_scale;    //fecundity at age, scaled
}

void model_parameters::get_reprod(void)
{
   //reprod is product of stuff going into reproductive capacity calcs
   reprod=elem_prod(wgt_klb,maturity_f);  
   reprod2=elem_prod(elem_prod(prop_f,maturity_f),wgt_klb);  
}

void model_parameters::get_length_at_age_dist(void)
{
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
}

void model_parameters::get_weight_at_age_landings(void)
{
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
}

void model_parameters::get_spr_F0(void)
{
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
}

void model_parameters::get_selectivity(void)
{
    selpar_Age0_rHD=1.0/(1.0+mfexp(-selpar_logitAge0_rHD));
    selpar_Age1_rHD=1.0/(1.0+mfexp(-selpar_logitAge1_rHD));
    selpar_Age2_rHD=1.0/(1.0+mfexp(-selpar_logitAge2_rHD));
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
}

void model_parameters::get_mortality(void)
{
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
}

void model_parameters::get_bias_corr(void)
{
  var_rec_dev=norm2(log_dev_rec(styr_rec_dev,endyr_rec_dev)-
              sum(log_dev_rec(styr_rec_dev,endyr_rec_dev))
              /(nyr_rec_dev))/(nyr_rec_dev); 
			  //cout<<"var rec dev"<<var_rec_dev<<endl;
  //if (set_BiasCor <= 0.0) {BiasCor=mfexp(var_rec_dev/2.0);}   //bias correction
  rec_sigma_sq=square(rec_sigma);
  rec_sigma_sqd2=rec_sigma_sq/2.0;
  if (set_BiasCor <= 0.0) {BiasCor=mfexp(rec_sigma_sqd2);}   //bias correction               
  else {BiasCor=set_BiasCor;}
}

void model_parameters::get_numbers_at_age(void)
{
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
  N_spr_initial(1)=1.0*mfexp(-1.0*Z_initial(1)*spawn_time_frac); //at peak spawning time;
  for (iage=2; iage<=nages; iage++)
    {
      N_spr_initial(iage)=N_spr_initial(iage-1)*
                   mfexp(-1.0*(Z_initial(iage-1)*(1.0-spawn_time_frac) + Z_initial(iage)*spawn_time_frac)); 
    }
  N_spr_initial(nages)=N_spr_initial(nages)/(1.0-mfexp(-1.0*Z_initial(nages))); //plus group
  spr_initial=sum(elem_prod(N_spr_initial,reprod));
  if (styr==styr_rec_dev) {R1=SR_eq_func(R0, steep, spr_F0, spr_initial, 1.0, SR_switch);} //without bias correction (deviation added later)
  else {R1=SR_eq_func(R0, steep, spr_F0, spr_initial, BiasCor, SR_switch);} //with bias correction                 
  if(R1<10.0) {R1=10.0;} //Avoid negative (or unreasonably low) popn sizes during search algorithm
  N_initial_eq(1)=R1;
  for (iage=2; iage<=nages; iage++)
  {
    N_initial_eq(iage)=N_initial_eq(iage-1)*
        mfexp(-1.0*(Z_initial(iage-1)));    
  }
  //plus group calculation
  N_initial_eq(nages)=N_initial_eq(nages)/(1.0-mfexp(-1.0*Z_initial(nages))); //plus group
  N(styr)(2,nages)=elem_prod(N_initial_eq(2,nages),mfexp(log_dev_Nage));
  if (styr==styr_rec_dev) {N(styr,1)=N_initial_eq(1)*mfexp(log_dev_rec(styr_rec_dev));}
  else {N(styr,1)=N_initial_eq(1);}
  N_mdyr(styr)(1,nages)=elem_prod(N(styr)(1,nages),(mfexp(-1.*(Z_initial(1,nages))*0.5))); //mid year 
  N_spawn(styr)(2,nages)=elem_prod(N(styr)(2,nages),(mfexp(-1.*(Z_initial(2,nages))*spawn_time_frac))); //peak spawning time 
  SSB(styr)=sum(elem_prod(N_spawn(styr),reprod));
  MatFemB(styr)=sum(elem_prod(N_spawn(styr),reprod2));
  B_q_DD(styr)=sum(elem_prod(N(styr)(set_q_DD_stage,nages),wgt_klb(set_q_DD_stage,nages)));
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
  rec=column(N,1);
  SdS0=SSB/S0;
}

void model_parameters::get_landings_numbers(void)
{
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
}

void model_parameters::get_landings_wgt(void)
{
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
}

void model_parameters::get_dead_discards(void)
{
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
}

void model_parameters::get_catchability_fcns(void)
{
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
}

void model_parameters::get_indices(void)
{
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
}

void model_parameters::get_length_comps(void)
{
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
}

void model_parameters::get_age_comps(void)
{
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
}

void model_parameters::get_weighted_current(void)
{
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
}

void model_parameters::get_msy(void)
{
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
}

void model_parameters::get_per_recruit_stuff(void)
{
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
}

void model_parameters::get_miscellaneous_stuff(void)
{
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
}

void model_parameters::get_projection(void)
{
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
}

void model_parameters::evaluate_objective_function(void)
{
  fval=0.0;
  fval_data=0.0;
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
  f_priors=0.0; 
  f_priors+=neg_log_prior(steep, set_steep(5), set_steep(6), set_steep(7)); 
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
}

dvar_vector model_parameters::logistic(const dvar_vector& ages, const dvariable& A50, const dvariable& slope)
{
  //ages=vector of ages, A50=age at 50% selectivity, slope=rate of increase
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=1./(1.+mfexp(-1.*slope*(ages-A50))); //logistic;  
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;
}

dvar_vector model_parameters::logistic_exponential(const dvar_vector& ages, const dvariable& A50, const dvariable& slope, const dvariable& sigma, const dvariable& joint)
{
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
}

dvar_vector model_parameters::logistic_double(const dvar_vector& ages, const dvariable& A501, const dvariable& slope1, const dvariable& A502, const dvariable& slope2)
{
  //ages=vector of ages, A50=age at 50% selectivity, slope=rate of increase, A502=age at 50% decrease additive to A501, slope2=slope of decrease
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=elem_prod( (1./(1.+mfexp(-1.*slope1*(ages-A501)))),(1.-(1./(1.+mfexp(-1.*slope2*(ages-(A501+A502)))))) );     
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;
}

dvar_vector model_parameters::logistic_joint(const dvar_vector& ages, const dvariable& A501, const dvariable& slope1, const dvariable& A502, const dvariable& slope2, const dvariable& satval, const dvariable& joint)
{
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
}

dvar_vector model_parameters::gaussian_double(const dvar_vector& ages, const dvariable& peak, const dvariable& top, const dvariable& ascwid, const dvariable& deswid, const dvariable& init, const dvariable& final)
{
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
}

dvariable model_parameters::SR_func(const dvariable& R0, const dvariable& h, const dvariable& spr_F0, const dvariable& SSB, int func)
{
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
}

dvariable model_parameters::SR_eq_func(const dvariable& R0, const dvariable& h, const dvariable& spr_F0, const dvariable& spr_F, const dvariable& BC, int func)
{
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
}

dvariable model_parameters::lk_lognormal(const dvar_vector& pred, const dvar_vector& obs, const dvar_vector& cv, const dvariable& wgt_dat)
{
  //pred=vector of predicted vals, obs=vector of observed vals, cv=vector of CVs in arithmetic space, wgt_dat=constant scaling of CVs
  //dzero is small value to avoid log(0) during search
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvar_vector var(cv.indexmin(),cv.indexmax()); //variance in log space
  var=log(1.0+square(cv/wgt_dat));   // convert cv in arithmetic space to variance in log space
  LkvalTmp=sum(0.5*elem_div(square(log(elem_div((pred+dzero),(obs+dzero)))),var) );
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;
}

dvariable model_parameters::lk_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const double& minSS, const dvariable& wgt_dat)
{
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
}

dvariable model_parameters::lk_dirichlet_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const dvariable& mbin, const double& minSS, const dvariable& log_dir_par)
{
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
}

dvariable model_parameters::neg_log_prior(dvariable pred, const double& prior, dvariable var, int pdf)
{
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
}

dvariable model_parameters::sdnr_multinomial(const double& ncomp, const dvar_vector& ages, const dvar_vector& nsamp, 
                                    const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const dvariable& wgt_dat)
{
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
}

dvariable model_parameters::sdnr_lognormal(const dvar_vector& pred, const dvar_vector& obs, const dvar_vector& cv, const dvariable& wgt_dat)
{
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
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
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
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
  time(&start);
  arrmblsize=20000000;
  gradient_structure::set_MAX_NVAR_OFFSET(1600);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(2000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(2000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(500);
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
  auto start = std::chrono::high_resolution_clock::now();
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint = defaults::iprint;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  std::cout << endl << argv[0] << " elapsed time is " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << " microseconds." << endl;
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
