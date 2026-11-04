//##  Author: NMFS, Beaufort Lab, Atlantic Fisheries Branch
//##  Analyst: Matt Vincent
//##  Species: Black Sea Bass
//##  Region: US South Atlantic
//##  SEDAR: 76 update 2025
//##  Date: 2026-04-28 08:28:18.233283


// Create a file with an R object from AD Model Builder

// open the file using the default AD Model Builder file name, and
// 6 digits of precision
open_r_file(adprogram_name + ".rdat", 6,missing);

// Example of an INFO object
open_r_info_list("info", true);
	wrt_r_item("title", "SEDAR 76 Operational Assessment");
	wrt_r_item("species", "Black sea bass");
	wrt_r_item("model", "Statistical Catch at Age");
	if(SR_switch==1)
	    {wrt_r_item("rec.model", "BH-steep");}
    if(SR_switch==2)
	    {wrt_r_item("rec.model", "Ricker-steep");}
    if(SR_switch==3)
	    {wrt_r_item("rec.model", "Mean");}
	wrt_r_item("base.run", "BSB76.tpl");
	wrt_r_item("units.length", "mm");
	wrt_r_item("units.weight", "lb");
	wrt_r_item("units.biomass", "1000 lb");
	wrt_r_item("units.ssb", "1E10 eggs");
	wrt_r_item("units.ypr", "lb");
	wrt_r_item("units.landings", "1000 lb");
	wrt_r_item("units.discards", "1000 dead fish");
    	wrt_r_item("units.numbers", "number fish");
    	wrt_r_item("units.naa", "number fish");
	wrt_r_item("units.rec", "number fish");
close_r_info_list();

// Abbreviations
open_r_info_list("abbreviations", false);
wrt_r_item("sBT", "survey blackfish trap (associated with index)");
wrt_r_item("sCT", "survey chevron trap (associated with index)");
wrt_r_item("cHL", "commercial handline (associated with landings)");
wrt_r_item("cHD", "commercial handline (associated with discards)");
wrt_r_item("cPT", "commercial pots (associated with landings)");
wrt_r_item("cPD", "commercial pots (associated with discards)");
wrt_r_item("cGN", "commercial general (associated with landings)");
wrt_r_item("cGD", "commercial general (associated with discards)");
wrt_r_item("cTW", "commercial trawl (associated with landings)");
wrt_r_item("rHB", "recreational headboat (associated with landings)");
wrt_r_item("rHD", "recreational headboat (associated with discards)");
wrt_r_item("rGN", "recreational general (associated with landings)");
wrt_r_item("rIN", "recreational general inshore (associated with landings)");
wrt_r_item("rGD", "recreational general (associated with discards)");
wrt_r_item("rID", "recreational general inshore (associated with discards)");

close_r_info_list();


// VECTOR object of parameters and estimated quantities
open_r_info_list("parms", false);
	wrt_r_item("styr", styr);
	wrt_r_item("endyr", endyr);
	wrt_r_item("styrR", styr_rec_dev);
	wrt_r_item("M.constant", M_constant);
	wrt_r_item("Linf", Linf);
	wrt_r_item("K", K);
	wrt_r_item("t0", t0);
	wrt_r_item("wgt.a", wgtpar_a);
	wrt_r_item("wgt.b", wgtpar_b);
	wrt_r_item("fec.a", fecpar_a);
	wrt_r_item("fec.b", fecpar_b);
	wrt_r_item("fec.batches", fecpar_batches);
	wrt_r_item("fec.scale", fecpar_scale);
	wrt_r_item("spawn.time", spawn_time_frac);
	wrt_r_item("Dprop.cGN.sel.D4", Dprop_cGN_sel_D4);
	wrt_r_item("Dprop.cGN.sel.cHL4", Dprop_cGN_sel_cHL4);
	wrt_r_item("Dprop.cGN.sel.cPT4", Dprop_cGN_sel_cPT4);
	wrt_r_item("q.sBT", mfexp(log_q_sBT));
	wrt_r_item("q.sTV", mfexp(log_q_sTV));
	wrt_r_item("q.cHL", mfexp(log_q_cHL));
	wrt_r_item("q.rHB", mfexp(log_q_rHB));
    //wrt_r_item("q.rHD", mfexp(log_q_rHD));
	wrt_r_item("q.rate",q_rate);
	wrt_r_item("q.beta",q_DD_beta);
	wrt_r_item("q.DD.B0.4plus", B0_q_DD);
	wrt_r_item("prop.F.L.cHL", prop_F_cHL);
	wrt_r_item("prop.F.L.cPT", prop_F_cPT);
	wrt_r_item("prop.F.L.rHB", prop_F_rHB);
	wrt_r_item("prop.F.L.rGN", prop_F_rGN);
	wrt_r_item("prop.F.L.rIN", prop_F_rIN);	
	wrt_r_item("prop.F.D.cGD", prop_F_cGD);
	wrt_r_item("prop.F.D.rHD", prop_F_rHD);
	wrt_r_item("prop.F.D.rGD", prop_F_rGD);
	wrt_r_item("prop.F.D.rID", prop_F_rID);
	wrt_r_item("F.init.ratio", F_init_ratio);
	wrt_r_item("Fend.mean", Fend_mean);

	wrt_r_item("B0", B0);
	wrt_r_item("Bstyr.B0", totB(styr)/B0);
	wrt_r_item("SSB0", S0);
	wrt_r_item("SSBstyr.SSB0", SSB(styr)/S0);
	wrt_r_item("Rstyr.R0", rec(styr)/R0);
if(SR_switch==1)
	{
    wrt_r_item("BH.biascorr",BiasCor);
	  wrt_r_item("BH.Phi0", spr_F0);
	  wrt_r_item("BH.R0", R0);
	  wrt_r_item("BH.steep", steep);
    }
    if(SR_switch==2)
	{
    wrt_r_item("Ricker.biascorr",BiasCor);
	  wrt_r_item("Ricker.Phi0", spr_F0);
	  wrt_r_item("Ricker.R0", R0);
	  wrt_r_item("Ricker.steep", steep);
    }
   if(SR_switch==3)
	{
      wrt_r_item("Mean.biascorr",BiasCor);
	  wrt_r_item("Mean.R0", R0);
	  wrt_r_item("Phi0", spr_F0);
    }		

	wrt_r_item("rec.mu.end", rec_mu);
	// wrt_r_item("earlyrec.mult", earlyrec_mult);
	wrt_r_item("R.sigma.logdevs", sigma_rec_dev);
	wrt_r_item("R.sigma.par",rec_sigma);
	wrt_r_item("R.autocorr",R_autocorr);
	wrt_r_item("R0", R0); //same as BH.R0, but used in BSR.time.plots
	wrt_r_item("R.virgin.bc", R_virgin); //bias-corrected virgin recruitment
	wrt_r_item("rec.lag", 0.0);
 	wrt_r_item("msy.klb.n", msy_n_klb_out);
	wrt_r_item("msy.knum.n", msy_n_knum_out);
	wrt_r_item("msy.klb", msy_lb_klb_out);
	wrt_r_item("msy.knum", msy_lb_knum_out);
	wrt_r_item("Fmsy.n", F_msy_n_out);
	wrt_r_item("Fmsy", F_msy_lb_out);
	wrt_r_item("SSBmsy.n", SSB_msy_n_out);
	wrt_r_item("SSBmsy", SSB_msy_lb_out);
	wrt_r_item("msst.n", (1.0-M_constant)*SSB_msy_n_out);
	wrt_r_item("msst", (1.0-M_constant)*SSB_msy_lb_out);
	wrt_r_item("Bmsy.n", B_msy_n_out);
	wrt_r_item("Bmsy", B_msy_lb_out);
	wrt_r_item("Rmsy.n", R_msy_n_out);
	wrt_r_item("Rmsy", R_msy_lb_out);
	wrt_r_item("sprmsy.n",spr_msy_n_out);
	wrt_r_item("sprmsy",spr_msy_lb_out);
	wrt_r_item("Dmsy.n.knum", D_msy_n_knum_out);
	wrt_r_item("Dmsy.knum", D_msy_lb_knum_out);
	wrt_r_item("Dmsy.n.klb", D_msy_n_klb_out);
	wrt_r_item("Dmsy.klb", D_msy_lb_klb_out);
	wrt_r_item("Fend.Fmsy.n", FdF_msy_n_end);
	wrt_r_item("Fend.Fmsy", FdF_msy_lb_end);
	wrt_r_item("Fend.Fmsy.n.mean", FdF_msy_n_end_mean);
	wrt_r_item("Fend.Fmsy.mean", FdF_msy_lb_end_mean);
	wrt_r_item("SSBend.SSBmsy.n", SdSSB_msy_n_end);
	wrt_r_item("SSBend.SSBmsy", SdSSB_msy_lb_end);
	wrt_r_item("SSBend.MSST.n", SSBend_MSST_n);
	wrt_r_item("SSBend.MSST", SSBend_MSST_lb);
	wrt_r_item("Fproxy", Fproxy_out); //For FishGraph
	wrt_r_item("SSB.Fproxy", SSB_Fproxy_out); //For FishGraph
	wrt_r_item("msst.Fproxy", (1.0-M_constant)*SSB_Fproxy_out);
	wrt_r_item("B.Fproxy", B_Fproxy_out);
	wrt_r_item("R.Fproxy", R_Fproxy_out);
	wrt_r_item("L.Fproxy.knum", L_Fproxy_knum_out);
	wrt_r_item("L.Fproxy.klb", L_Fproxy_klb_out);
	wrt_r_item("Dmort.cGD",Dmort_cGD);
	// wrt_r_item("Dmort.cPT1",Dmort_cPD1);
	// wrt_r_item("Dmort.cPT2",Dmort_cPD2);	
	// wrt_r_item("Dmort.cHD",Dmort_cHD);
	wrt_r_item("Dmort.rHD",Dmort_rHD);
	wrt_r_item("Dmort.rGD",Dmort_rGD);
close_r_info_list();

 // VECTOR object of parameters and estimated quantities
open_r_info_list("spr.brps", false);
	wrt_r_item("SPR.proxy.level", spr_proxy);
	wrt_r_item("Fproxy", Fproxy_out);
	wrt_r_item("SSB.Fproxy", SSB_Fproxy_out);
	wrt_r_item("msst.Fproxy", (1.0-M_constant)*SSB_Fproxy_out);
	wrt_r_item("B.Fproxy", B_Fproxy_out);
	wrt_r_item("R.Fproxy", R_Fproxy_out);
	wrt_r_item("L.Fproxy.knum", L_Fproxy_knum_out);
	wrt_r_item("L.Fproxy.klb", L_Fproxy_klb_out);
	wrt_r_item("Fend.Fproxy.mean", FdFproxy_end_mean);
	wrt_r_item("SSBend.SSBFproxy", SdSSB_Fproxy_end);
	wrt_r_item("SSBend.MSSTFproxy", Sdmsst_Fproxy_end);
// F30
	wrt_r_item("F30", F30_out);
	wrt_r_item("SSB.F30", SSB_F30_out);
	wrt_r_item("msst.F30", (1.0-M_constant)*SSB_F30_out);
	wrt_r_item("B.F30", B_F30_out);
	wrt_r_item("R.F30", R_F30_out);
	wrt_r_item("L.F30.knum", L_F30_knum_out);
	wrt_r_item("L.F30.klb", L_F30_klb_out);
	wrt_r_item("Fend.F30.mean", FdF30_end_mean);
	wrt_r_item("SSBend.SSBF30", SdSSB_F30_end);
	wrt_r_item("SSBend.MSSTF30", Sdmsst_F30_end);
// F40
	wrt_r_item("F40", F40_out);
	wrt_r_item("SSB.F40", SSB_F40_out);
	wrt_r_item("msst.F40", (1.0-M_constant)*SSB_F40_out);
	wrt_r_item("B.F40", B_F40_out);
	wrt_r_item("R.F40", R_F40_out);
	wrt_r_item("L.F40.knum", L_F40_knum_out);
	wrt_r_item("L.F40.klb", L_F40_klb_out);
	wrt_r_item("Fend.F40.mean", FdF40_end_mean);
	wrt_r_item("SSBend.SSBF40", SdSSB_F40_end);
	wrt_r_item("SSBend.MSSTF40", Sdmsst_F40_end);
// F50
	wrt_r_item("F50", F50_out);
	wrt_r_item("SSB.F50", SSB_F50_out);
	wrt_r_item("msst.F50", (1.0-M_constant)*SSB_F50_out);
	wrt_r_item("B.F50", B_F50_out);
	wrt_r_item("R.F50", R_F50_out);
	wrt_r_item("L.F50.knum", L_F50_knum_out);
	wrt_r_item("L.F50.klb", L_F50_klb_out);
	wrt_r_item("Fend.F50.mean", FdF50_end_mean);
	wrt_r_item("SSBend.SSBF50", SdSSB_F50_end);
	wrt_r_item("SSBend.MSSTF50", Sdmsst_F50_end);
 close_r_info_list();

// MATRIX object of parameter constraints
open_r_df("parm.cons",1,8,2);
    wrt_r_namevector(1,8);
    wrt_r_df_col("Linf",Linf_out);
    wrt_r_df_col("K",K_out);
    wrt_r_df_col("t0",t0_out);
    wrt_r_df_col("len.cv.val",len_cv_val_out);
    wrt_r_df_col("log.R0",log_R0_out);
    wrt_r_df_col("steep",steep_out);
    wrt_r_df_col("rec.sigma",rec_sigma_out);
    wrt_r_df_col("R.autocorr",R_autocorr_out);
    // wrt_r_df_col("log.dm.sTV.lc",log_dm_sTV_lc_out);
	wrt_r_df_col("log.dm.lenc.sBT",log_dm_sBT_lc_out);
	wrt_r_df_col("log.dm.lenc.cHL",log_dm_cHL_lc_out);
	wrt_r_df_col("log.dm.lenc.cPT",log_dm_cPT_lc_out);
	//wrt_r_df_col("log.dm.cTW.lc",log_dm_cTW_lc_out);
	wrt_r_df_col("log.dm.lenc.rHB",log_dm_rHB_lc_out);
	wrt_r_df_col("log.dm.lenc.rGN",log_dm_rGN_lc_out);
	wrt_r_df_col("log.dm.lenc.rHD",log_dm_rHD_lc_out);
	wrt_r_df_col("log.dm.agec.sTV",log_dm_sTV_ac_out);
	wrt_r_df_col("log.dm.agec.sBT",log_dm_sBT_ac_out);
	wrt_r_df_col("log.dm.agec.cHL",log_dm_cHL_ac_out);
	wrt_r_df_col("log.dm.agec.cPT",log_dm_cPT_ac_out);
	wrt_r_df_col("log.dm.agec.rHB",log_dm_rHB_ac_out);
	//wrt_r_df_col("log.dm.rGN.ac",log_dm_rGN_ac_out);
			
			
	wrt_r_df_col("selpar.A50.sBT", selpar_A50_sBT_out);
	wrt_r_df_col("selpar.slope.sBT", selpar_slope_sBT_out);
 	wrt_r_df_col("selpar.A502.sBT", selpar_A502_sBT_out);
	wrt_r_df_col("selpar.slope2.sBT", selpar_slope2_sBT_out);

	wrt_r_df_col("selpar.A50.sTV", selpar_A50_sTV_out);
	wrt_r_df_col("selpar.slope.sTV", selpar_slope_sTV_out);
	wrt_r_df_col("selpar.A502.sTV", selpar_A502_sTV_out);
	wrt_r_df_col("selpar.slope2.sTV", selpar_slope2_sTV_out);
	
	wrt_r_df_col("selpar.A50b02.cHL", selpar_A50b02_cHL_out);
	wrt_r_df_col("selpar.slopeb02.cHL", selpar_slopeb02_cHL_out);
	wrt_r_df_col("selpar.A50b03.cHL", selpar_A50b03_cHL_out);
	wrt_r_df_col("selpar.slopeb03.cHL", selpar_slopeb03_cHL_out);
	wrt_r_df_col("selpar.A50b04.cHL", selpar_A50b04_cHL_out);
	wrt_r_df_col("selpar.slopeb04.cHL", selpar_slopeb04_cHL_out);
    
	wrt_r_df_col("selpar.A50b02.cPT", selpar_A50b02_cPT_out);
	wrt_r_df_col("selpar.slopeb02.cPT", selpar_slopeb02_cPT_out);
	wrt_r_df_col("selpar.A50b03.cPT", selpar_A50b03_cPT_out);
	wrt_r_df_col("selpar.slopeb03.cPT", selpar_slopeb03_cPT_out);
	wrt_r_df_col("selpar.A50b04.cPT", selpar_A50b04_cPT_out);
	wrt_r_df_col("selpar.slopeb04.cPT", selpar_slopeb04_cPT_out);
    
	wrt_r_df_col("selpar.A50b01.rHB", selpar_A50b01_rHB_out);
	wrt_r_df_col("selpar.slopeb01.rHB", selpar_slopeb01_rHB_out);
	wrt_r_df_col("selpar.A50b02.rHB", selpar_A50b02_rHB_out);
	wrt_r_df_col("selpar.slopeb02.rHB", selpar_slopeb02_rHB_out);
	wrt_r_df_col("selpar.A50b03.rHB", selpar_A50b03_rHB_out);
	wrt_r_df_col("selpar.slopeb03.rHB", selpar_slopeb03_rHB_out);
	wrt_r_df_col("selpar.A50b04.rHB", selpar_A50b04_rHB_out);
	wrt_r_df_col("selpar.slopeb04.rHB", selpar_slopeb04_rHB_out);
	wrt_r_df_col("selpar.A50b05.rHB", selpar_A50b05_rHB_out);
	wrt_r_df_col("selpar.slopeb05.rHB", selpar_slopeb05_rHB_out);

  wrt_r_df_col("selpar.logitAge0.rHD",selpar_logitAge0_rHD_out);
  wrt_r_df_col("selpar.logitAge1.rHD",selpar_logitAge1_rHD_out);
  wrt_r_df_col("selpar.logitAge2.rHD",selpar_logitAge2_rHD_out);
    
	wrt_r_df_col("selpar.A50b04.rHD", selpar_A50b04_rHD_out);
	wrt_r_df_col("selpar.slopeb04.rHD", selpar_slopeb04_rHD_out);
	wrt_r_df_col("selpar.A502b04.rHD", selpar_A502b04_rHD_out);
	wrt_r_df_col("selpar.slope2b04.rHD", selpar_slope2b04_rHD_out);
	wrt_r_df_col("selpar.A50b05.rHD", selpar_A50b05_rHD_out);
	wrt_r_df_col("selpar.slopeb05.rHD", selpar_slopeb05_rHD_out);
	wrt_r_df_col("selpar.A502b05.rHD", selpar_A502b05_rHD_out);
	wrt_r_df_col("selpar.slope2b05.rHD", selpar_slope2b05_rHD_out);
	
	  wrt_r_df_col("log.q.cHL",log_q_cHL_out);
    wrt_r_df_col("log.q.rHB",log_q_rHB_out);
    wrt_r_df_col("log.q.sTV",log_q_sTV_out);
	  wrt_r_df_col("log.q.sBT",log_q_sBT_out);
    wrt_r_df_col("M.constant",M_constant_out);
    wrt_r_df_col("log.avg.F.L.cHL",out_log_avg_F_L_cHL);
	  wrt_r_df_col("log.avg.F.L.cPT",out_log_avg_F_L_cPT);
	  wrt_r_df_col("log.avg.F.L.cTW",out_log_avg_F_L_cTW);
    wrt_r_df_col("log.avg.F.L.rHB",out_log_avg_F_L_rHB);
    wrt_r_df_col("log.avg.F.L.rGN",out_log_avg_F_L_rGN);
    wrt_r_df_col("log.avg.F.L.rIN",out_log_avg_F_L_rIN);
    wrt_r_df_col("log.avg.F.D.cGD",out_log_avg_F_D_cGD);
    wrt_r_df_col("log.avg.F.D.rHD",out_log_avg_F_D_rHD);
    wrt_r_df_col("log.avg.F.D.rGD",out_log_avg_F_D_rGD);
    wrt_r_df_col("log.avg.F.D.rID",out_log_avg_F_D_rID);
    wrt_r_df_col("F.init",F_init_ratio_out);
    wrt_r_df_col("q.rate",q_rate_out);
close_r_df();


// VECTOR object of likelihood contributions
open_r_vector("like");
	wrt_r_item("lk.total", fval);            //weighted likelihood
    wrt_r_item("lk.unwgt.data", fval_data);  //likelihood of just data components
    wrt_r_item("lk.U.sBT", f_cpue_sBT);
    wrt_r_item("lk.U.sTV", f_cpue_sTV);
    wrt_r_item("lk.U.cHL", f_cpue_cHL);
    wrt_r_item("lk.U.rHB", f_cpue_rHB);
    wrt_r_item("lk.L.cHL", f_L_cHL);
    wrt_r_item("lk.L.cPT", f_L_cPT);
    wrt_r_item("lk.L.cTW", f_L_cTW);
    wrt_r_item("lk.L.rHB", f_L_rHB);
    wrt_r_item("lk.L.rGN", f_L_rGN);
    wrt_r_item("lk.L.rIN", f_L_rIN);
    wrt_r_item("lk.D.cGD", f_D_cGD);
    wrt_r_item("lk.D.rHD", f_D_rHD);
    wrt_r_item("lk.D.rGD", f_D_rGD);
    wrt_r_item("lk.D.rID", f_D_rID);
    wrt_r_item("lk.lenc.sBT", f_lenc_sBT);
    wrt_r_item("lk.lenc.cHL", f_lenc_cHL);
    wrt_r_item("lk.lenc.cPT", f_lenc_cPT);
    wrt_r_item("lk.lenc.rHB", f_lenc_rHB);
    wrt_r_item("lk.lenc.rGN", f_lenc_rGN);
    wrt_r_item("lk.lenc.rHD", f_lenc_rHD);
    wrt_r_item("lk.agec.sBT", f_agec_sBT);
    wrt_r_item("lk.agec.sTV", f_agec_sTV);
    wrt_r_item("lk.agec.cHL", f_agec_cHL);
    wrt_r_item("lk.agec.cPT", f_agec_cPT);
    wrt_r_item("lk.agec.rHB", f_agec_rHB);
    wrt_r_item("lk.priors",f_priors);
    wrt_r_item("lk.U.RW.cHL",f_cpue_RW_cHL);
    wrt_r_item("lk.U.RW.rHB",f_cpue_RW_rHB);

    wrt_r_item("lk.SRfit", f_rec_dev);
    wrt_r_item("lk.SRearly", f_rec_dev_early);
    wrt_r_item("lk.SRend", f_rec_dev_end);

	wrt_r_item("w.L", w_L);
	wrt_r_item("w.D", w_D);
	wrt_r_item("w.lc.sBT", w_lenc_sBT);
	wrt_r_item("w.lc.cHL", w_lenc_cHL);
	wrt_r_item("w.lc.cPT", w_lenc_cPT);
	wrt_r_item("w.lc.rHB", w_lenc_rHB);
	wrt_r_item("w.lc.rHD", w_lenc_rHD);
	wrt_r_item("w.lc.rGN", w_lenc_rGN);
	wrt_r_item("w.ac.sBT", w_agec_sBT);
	wrt_r_item("w.ac.sTV", w_agec_sTV);
	wrt_r_item("w.ac.cHL", w_agec_cHL);
	wrt_r_item("w.ac.cPT", w_agec_cPT);
	wrt_r_item("w.ac.rHB", w_agec_rHB);
	wrt_r_item("w.U.sBT", w_cpue_sBT);
	wrt_r_item("w.U.sTV", w_cpue_sTV);
	wrt_r_item("w.U.cHL", w_cpue_cHL);
	wrt_r_item("w.U.rHB", w_cpue_rHB);
	wrt_r_item("w.R", w_rec);
	wrt_r_item("w.R.init", w_rec_early);
	wrt_r_item("w.R.end", w_rec_end);
        wrt_r_item("gradient.max",grad_max);
close_r_vector();

// VECTOR object of parameters and estimated quantities
open_r_info_list("sel.parms",false);

	wrt_r_item("selpar.A50.sBT", selpar_A50_sBT);
	wrt_r_item("selpar.slope.sBT", selpar_slope_sBT);
	wrt_r_item("selpar.A502.sBT", selpar_A502_sBT);
	wrt_r_item("selpar.slope2.sBT", selpar_slope2_sBT);

	wrt_r_item("selpar.A50.sTV", selpar_A50_sTV);
	wrt_r_item("selpar.slope.sTV", selpar_slope_sTV);
	wrt_r_item("selpar.A502.sTV", selpar_A502_sTV);
	wrt_r_item("selpar.slope2.sTV", selpar_slope2_sTV);
	
	//wrt_r_item("selpar.A50.Vid", selpar_A50_Vid);
	//wrt_r_item("selpar.slope.Vid", selpar_slope_Vid);

	wrt_r_item("selpar.A50b02.cHL", selpar_A50b02_cHL);
	wrt_r_item("selpar.slopeb02.cHL", selpar_slopeb02_cHL);
	wrt_r_item("selpar.A50b03.cHL", selpar_A50b03_cHL);
	wrt_r_item("selpar.slopeb03.cHL", selpar_slopeb03_cHL);
	wrt_r_item("selpar.A50b04.cHL", selpar_A50b04_cHL);
	wrt_r_item("selpar.slopeb04.cHL", selpar_slopeb04_cHL);

	wrt_r_item("selpar.A50b02.cPT", selpar_A50b02_cPT);
	wrt_r_item("selpar.slopeb02.cPT", selpar_slopeb02_cPT);
	wrt_r_item("selpar.A50b03.cPT", selpar_A50b03_cPT);
	wrt_r_item("selpar.slopeb03.cPT", selpar_slopeb03_cPT);
	wrt_r_item("selpar.A50b04.cPT", selpar_A50b04_cPT);
	wrt_r_item("selpar.slopeb04.cPT", selpar_slopeb04_cPT);

	wrt_r_item("selpar.A50b01.rHB", selpar_A50b01_rHB);
	wrt_r_item("selpar.slopeb01.rHB", selpar_slopeb01_rHB);
	wrt_r_item("selpar.A50b02.rHB", selpar_A50b02_rHB);
	wrt_r_item("selpar.slopeb02.rHB", selpar_slopeb02_rHB);
	wrt_r_item("selpar.A50b03.rHB", selpar_A50b03_rHB);
	wrt_r_item("selpar.slopeb03.rHB", selpar_slopeb03_rHB);
	wrt_r_item("selpar.A50b04.rHB", selpar_A50b04_rHB);
	wrt_r_item("selpar.slopeb04.rHB", selpar_slopeb04_rHB);
	wrt_r_item("selpar.A50b05.rHB", selpar_A50b05_rHB);
	wrt_r_item("selpar.slopeb05.rHB", selpar_slopeb05_rHB);

  wrt_r_item("selpar.Age0.rHD", selpar_Age0_rHD);
  wrt_r_item("selpar.Age1.rHD", selpar_Age1_rHD);
  wrt_r_item("selpar.Age2.rHD", selpar_Age2_rHD);

	wrt_r_df_col("selpar.A50b04.rHD", selpar_A50b04_rHD_out);
	wrt_r_df_col("selpar.slopeb04.rHD", selpar_slopeb04_rHD_out);
	wrt_r_df_col("selpar.A502b04.rHD", selpar_A502b04_rHD_out);
	wrt_r_df_col("selpar.slope2b04.rHD", selpar_slope2b04_rHD_out);
	wrt_r_df_col("selpar.A50b05.rHD", selpar_A50b05_rHD_out);
	wrt_r_df_col("selpar.slopeb05.rHD", selpar_slopeb05_rHD_out);
	wrt_r_df_col("selpar.A502b05.rHD", selpar_A502b05_rHD_out);
	wrt_r_df_col("selpar.slope2b05.rHD", selpar_slope2b05_rHD_out);
	
close_r_info_list();

// DATA FRAME of time series deviation vector estimates
// names used in this object must match the names used in the "parm.tvec.cons" object
open_r_df("parm.tvec", styr, endyr, 2);
	wrt_r_namevector(styr,endyr);
	wrt_r_df_col("year", styr,endyr);
  wrt_r_df_col("log.rec.dev", log_rec_dev_out);
  wrt_r_df_col("log.dev.F.L.cHL", log_dev_F_L_out_cHL);
	wrt_r_df_col("log.dev.F.L.cPT", log_dev_F_L_out_cPT);
	wrt_r_df_col("log.dev.F.L.cTW", log_dev_F_L_out_cTW);
  wrt_r_df_col("log.dev.F.L.rHB", log_dev_F_L_out_rHB);
  wrt_r_df_col("log.dev.F.L.rGN", log_dev_F_L_out_rGN);
  wrt_r_df_col("log.dev.F.L.rIN", log_dev_F_L_out_rIN);
  wrt_r_df_col("log.dev.F.D.cGD", log_dev_F_D_out_cGD);
  wrt_r_df_col("log.dev.F.D.rHD", log_dev_F_D_out_rHD);
  wrt_r_df_col("log.dev.F.D.rGD", log_dev_F_D_out_rGD);
  wrt_r_df_col("log.dev.F.D.rID", log_dev_F_D_out_rID);
close_r_df();

// MATRIX object of deviation vector constraints
// names used in this object must match the names used in the "parm.tvec" object
open_r_df("parm.tvec.cons",1,3,2);
    wrt_r_namevector(1,3);
    wrt_r_df_col("log.rec.dev",set_log_dev_rec);
    wrt_r_df_col("log.dev.F.L.cHL",set_log_dev_F_L_cHL);
	  wrt_r_df_col("log.dev.F.L.cPT",set_log_dev_F_L_cPT);
	  wrt_r_df_col("log.dev.F.L.cTW",set_log_dev_F_L_cTW);
    wrt_r_df_col("log.dev.F.L.rHB",set_log_dev_F_L_rHB);
    wrt_r_df_col("log.dev.F.L.rGN",set_log_dev_F_L_rGN);
    wrt_r_df_col("log.dev.F.L.rIN",set_log_dev_F_L_rIN);
    wrt_r_df_col("log.dev.F.D.cGD",set_log_dev_F_D_cGD);
    wrt_r_df_col("log.dev.F.D.rHD",set_log_dev_F_D_rHD);
    wrt_r_df_col("log.dev.F.D.rGD",set_log_dev_F_D_rGD);
    wrt_r_df_col("log.dev.F.D.rID",set_log_dev_F_D_rID);
	//wrt_r_df_col("log.q.dev.rHB.RWq", set_log_dev_RWq);
	//wrt_r_df_col("log.q.dev.cHL.RWq", set_log_dev_RWq);
close_r_df();

// DATA FRAME of age vector deviation estimates
// names used in this object must match the names used in the "parm.avec.cons" object
open_r_df("parm.avec", 1, nages, 2);
	wrt_r_namevector(1,nages);
	wrt_r_df_col("age", agebins); //deviations for first age not estimated
    wrt_r_df_col("log.Nage.dev", log_Nage_dev_output);
close_r_df();

// MATRIX object of age vector deviation constraints
// names used in this object must match the names used in the "parm.avec" object
open_r_df("parm.avec.cons",1,3,2);
    wrt_r_namevector(1,3);
    wrt_r_df_col("log.Nage.dev",set_log_dev_Nage);
close_r_df();


    open_r_matrix("N.age");
    wrt_r_matrix(N, 2, 2);
    wrt_r_namevector(styr, endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("N.age.mdyr");
    wrt_r_matrix(N_mdyr, 2, 2);
    wrt_r_namevector(styr, endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("N.age.spawn");
    wrt_r_matrix(N_spawn, 2, 2);
    wrt_r_namevector(styr, endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("B.age");
    wrt_r_matrix(B, 2, 2);
    wrt_r_namevector(styr, endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("F.age");
    wrt_r_matrix(F, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("Z.age");
    wrt_r_matrix(Z, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("L.age.pred.knum");
    wrt_r_matrix(L_total_num/1000.0, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("L.age.pred.klb");
    wrt_r_matrix(L_total_klb, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();


// LIST object with annual selectivity at age by fishery

open_r_list("size.age.fishery");

    open_r_matrix("len.cHL.mm");
    wrt_r_matrix(len_cHL_mm, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("len.cPT.mm");
    wrt_r_matrix(len_cPT_mm, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("len.cTW.mm");
    wrt_r_matrix(len_cTW_mm, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("len.rHB.mm");
    wrt_r_matrix(len_rHB_mm, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("len.rGN.mm");
    wrt_r_matrix(len_rGN_mm, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();
    
    open_r_matrix("len.rIN.mm");
    wrt_r_matrix(len_rIN_mm, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();    

    open_r_matrix("len.cGD.mm");
    wrt_r_matrix(len_cGD_mm, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("len.rHD.mm");
    wrt_r_matrix(len_rHD_mm, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("len.rGD.mm");
    wrt_r_matrix(len_rGD_mm, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();
    
    open_r_matrix("len.rID.mm");
    wrt_r_matrix(len_rID_mm, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();    

    open_r_matrix("wgt.cHL.klb");
    wrt_r_matrix(wgt_cHL_klb, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("wgt.cPT.klb");
    wrt_r_matrix(wgt_cPT_klb, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("wgt.cTW.klb");
    wrt_r_matrix(wgt_cTW_klb, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("wgt.rHB.klb");
    wrt_r_matrix(wgt_rHB_klb, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("wgt.rGN.klb");
    wrt_r_matrix(wgt_rGN_klb, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();
    
    open_r_matrix("wgt.rIN.klb");
    wrt_r_matrix(wgt_rIN_klb, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();    

    open_r_matrix("wgt.cGD.klb");
    wrt_r_matrix(wgt_cGD_klb, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("wgt.rHD.klb");
    wrt_r_matrix(wgt_rHD_klb, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("wgt.rGD.klb");
    wrt_r_matrix(wgt_rGD_klb, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();
    
    open_r_matrix("wgt.rID.klb");
    wrt_r_matrix(wgt_rID_klb, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();    

close_r_list();

open_r_list("sel.age");

    wrt_r_complete_vector("sel.v.wgted.L",sel_wgted_L, agebins);
    wrt_r_complete_vector("sel.v.wgted.D",sel_wgted_D, agebins);
    wrt_r_complete_vector("sel.v.wgted.tot",sel_wgted_tot, agebins);

    // wrt_r_complete_vector("sel.v.sBT",sel_sBT_vec, agebins);
    // wrt_r_complete_vector("sel.v.sTV",sel_sTV_vec, agebins);

    open_r_matrix("sel.m.sBT");
    wrt_r_matrix(sel_sBT, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();    
    
    open_r_matrix("sel.m.sTV");
    wrt_r_matrix(sel_sTV, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();
    
    open_r_matrix("sel.m.cHL");
    wrt_r_matrix(sel_cHL, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("sel.m.cPT");
    wrt_r_matrix(sel_cPT, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("sel.m.cTW");
    wrt_r_matrix(sel_cTW, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("sel.m.rHB");
    wrt_r_matrix(sel_rHB, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("sel.m.rGN");
    wrt_r_matrix(sel_rGN, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();
    
    open_r_matrix("sel.m.rIN");
    wrt_r_matrix(sel_rIN, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();    

    open_r_matrix("sel.m.cGD");
    wrt_r_matrix(sel_cGD, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("sel.m.rHD");
    wrt_r_matrix(sel_rHD, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("sel.m.rGD");
    wrt_r_matrix(sel_rGD, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();
    
    open_r_matrix("sel.m.rID");
    wrt_r_matrix(sel_rID, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();    

close_r_list();



//LIST object with predicted and observed composition data
open_r_list("comp.mats");
    open_r_matrix("lcomp.sBT.ob");  //Excluded because the blackfish trap comp is a single year
    wrt_r_matrix(obs_lenc_sBT, 2, 2);
    wrt_r_namevector(yrs_lenc_sBT);
    wrt_r_namevector(lenbins);
    close_r_matrix();

    open_r_matrix("lcomp.sBT.pr");
    wrt_r_matrix(pred_lenc_sBT, 2, 2);
    wrt_r_namevector(yrs_lenc_sBT);
    wrt_r_namevector(lenbins);
    close_r_matrix();

    open_r_matrix("lcomp.cHL.ob");
    wrt_r_matrix(obs_lenc_cHL, 2, 2);
    wrt_r_namevector(yrs_lenc_cHL);
    wrt_r_namevector(lenbins);
    close_r_matrix();

    open_r_matrix("lcomp.cHL.pr");
    wrt_r_matrix(pred_lenc_cHL, 2, 2);
    wrt_r_namevector(yrs_lenc_cHL);
    wrt_r_namevector(lenbins);
    close_r_matrix();

    open_r_matrix("lcomp.cPT.ob");
    wrt_r_matrix(obs_lenc_cPT, 2, 2);
    wrt_r_namevector(yrs_lenc_cPT);
    wrt_r_namevector(lenbins);
    close_r_matrix();

    open_r_matrix("lcomp.cPT.pr");
    wrt_r_matrix(pred_lenc_cPT, 2, 2);
    wrt_r_namevector(yrs_lenc_cPT);
    wrt_r_namevector(lenbins);
    close_r_matrix();

    open_r_matrix("lcomp.rHB.ob");
    wrt_r_matrix(obs_lenc_rHB, 2, 2);
    wrt_r_namevector(yrs_lenc_rHB);
    wrt_r_namevector(lenbins);
    close_r_matrix();

    open_r_matrix("lcomp.rHB.pr");
    wrt_r_matrix(pred_lenc_rHB, 2, 2);
    wrt_r_namevector(yrs_lenc_rHB);
    wrt_r_namevector(lenbins);
    close_r_matrix();

    open_r_matrix("lcomp.rGN.ob");
    wrt_r_matrix(obs_lenc_rGN, 2, 2);
    wrt_r_namevector(yrs_lenc_rGN);
    wrt_r_namevector(lenbins);
    close_r_matrix();

    open_r_matrix("lcomp.rGN.pr");
    wrt_r_matrix(pred_lenc_rGN, 2, 2);
    wrt_r_namevector(yrs_lenc_rGN);
    wrt_r_namevector(lenbins);
    close_r_matrix();

    open_r_matrix("lcomp.rHD.ob");
    wrt_r_matrix(obs_lenc_rHD, 2, 2);
    wrt_r_namevector(yrs_lenc_rHD);
    wrt_r_namevector(lenbins);
    close_r_matrix();

    open_r_matrix("lcomp.rHD.pr");
    wrt_r_matrix(pred_lenc_rHD, 2, 2);
    wrt_r_namevector(yrs_lenc_rHD);
    wrt_r_namevector(lenbins);
    close_r_matrix();
    //-----------------------------------

    open_r_matrix("acomp.sBT.ob");
    wrt_r_matrix(obs_agec_sBT, 2, 2);
    wrt_r_namevector(yrs_agec_sBT);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("acomp.sBT.pr");
    wrt_r_matrix(pred_agec_sBT, 2, 2);
    wrt_r_namevector(yrs_agec_sBT);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("acomp.sTV.ob");
    wrt_r_matrix(obs_agec_sTV, 2, 2);
    wrt_r_namevector(yrs_agec_sTV);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("acomp.sTV.pr");
    wrt_r_matrix(pred_agec_sTV, 2, 2);
    wrt_r_namevector(yrs_agec_sTV);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("acomp.cHL.ob");
    wrt_r_matrix(obs_agec_cHL, 2, 2);
    wrt_r_namevector(yrs_agec_cHL);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("acomp.cHL.pr");
    wrt_r_matrix(pred_agec_cHL, 2, 2);
    wrt_r_namevector(yrs_agec_cHL);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("acomp.cPT.ob");
    wrt_r_matrix(obs_agec_cPT, 2, 2);
    wrt_r_namevector(yrs_agec_cPT);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("acomp.cPT.pr");
    wrt_r_matrix(pred_agec_cPT, 2, 2);
    wrt_r_namevector(yrs_agec_cPT);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("acomp.rHB.ob");
    wrt_r_matrix(obs_agec_rHB, 2, 2);
    wrt_r_namevector(yrs_agec_rHB);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("acomp.rHB.pr");
    wrt_r_matrix(pred_agec_rHB, 2, 2);
    wrt_r_namevector(yrs_agec_rHB);
    wrt_r_namevector(agebins);
    close_r_matrix();

    //open_r_matrix("acomp.rGN.ob");
    //wrt_r_matrix(obs_rGN_agec, 2, 2);
    //wrt_r_namevector(yrs_rGN_agec);
    //wrt_r_namevector(agebins);
    //close_r_matrix();
    //
    //open_r_matrix("acomp.rGN.pr");
    //wrt_r_matrix(pred_agec_rGN, 2, 2);
    //wrt_r_namevector(yrs_rGN_agec);
    //wrt_r_namevector(agebins);
    //close_r_matrix();
close_r_list();

// VECTOR object of SDNR calculations
open_r_vector("sdnr");
    wrt_r_item("sdnr.U.cHL", sdnr_I_cHL);
	//wrt_r_item("sdnr.U.Vid", sdnr_I_Vid);
    wrt_r_item("sdnr.U.rHB", sdnr_I_rHB);
    wrt_r_item("sdnr.U.sTV", sdnr_I_sTV);
	wrt_r_item("sdnr.U.sBT", sdnr_I_sBT);
	 
    //wrt_r_item("sdnr.lc.sBT", sdnr_lc_sBT);
    //wrt_r_item("sdnr.lc.sTV", sdnr_lc_sTV);	
	//wrt_r_item("sdnr.lc.cHL", sdnr_lc_cHL);
	//wrt_r_item("sdnr.lc.cPT", sdnr_lc_cPT);	
    //wrt_r_item("sdnr.lc.rHB", sdnr_lc_rHB);
	//wrt_r_item("sdnr.lc.rHD", sdnr_lc_rHD);
	//wrt_r_item("sdnr.lc.rGN", sdnr_lc_rGN);
	//
    //wrt_r_item("sdnr.ac.sBT", sdnr_ac_sBT);
    //wrt_r_item("sdnr.ac.sTV", sdnr_ac_sTV);
	//wrt_r_item("sdnr.ac.cHL", sdnr_ac_cHL);
	//wrt_r_item("sdnr.ac.rHB", sdnr_ac_rHB);
	//wrt_r_item("sdnr.ac.Vid", sdnr_ac_Vid);
	
 close_r_vector();

// DATA FRAME of time series
open_r_df("t.series", styr, endyr, 2);
	wrt_r_namevector(styr,endyr);
	wrt_r_df_col("year", styr,endyr);
	//wrt_r_df_col("F.Fmsy.n", FdF_msy_n);
	wrt_r_df_col("F.Fmsy", FdF_msy_lb);
	wrt_r_df_col("F.full", Fapex);
	wrt_r_df_col("F.L.cHL", F_cHL_out);
	wrt_r_df_col("F.L.cPT", F_cPT_out);
	wrt_r_df_col("F.L.cTW", F_cTW_out);
	wrt_r_df_col("F.L.rHB", F_rHB_out);
	wrt_r_df_col("F.L.rGN", F_rGN_out);
	wrt_r_df_col("F.L.rIN", F_rIN_out);
	wrt_r_df_col("F.D.cGD", F_cGD_out);
	wrt_r_df_col("F.D.rHD", F_rHD_out);
	wrt_r_df_col("F.D.rGD", F_rGD_out);
	wrt_r_df_col("F.D.rID", F_rID_out);
	wrt_r_df_col("Fsum", Fsum);
    wrt_r_df_col("N", totN);   //abundance at start of year
    wrt_r_df_col("recruits", rec);
    wrt_r_df_col("logR.dev", log_rec_dev_output); //places zeros in yrs deviations not estimated
    wrt_r_df_col("SSB", SSB);
    wrt_r_df_col("SSB.SSBmsy.n", SdSSB_msy_n);
    wrt_r_df_col("SSB.SSBmsy", SdSSB_msy_lb);
    wrt_r_df_col("SSB.msst.n", SdSSB_msy_n/(1.0-M_constant));
    wrt_r_df_col("SSB.msst", SdSSB_msy_lb/(1.0-M_constant));
    wrt_r_df_col("FdFproxy", FdFproxy);
    wrt_r_df_col("SSB.SSBFproxy", SdSSB_Fproxy);
    wrt_r_df_col("SSB.msstFproxy", Sdmsst_Fproxy);
    wrt_r_df_col("B", totB);
    wrt_r_df_col("B.B0", totB/B0);
    wrt_r_df_col("MatFemB", MatFemB);
    wrt_r_df_col("SPR.static", spr_static);

    wrt_r_df_col("total.L.klb", L_total_klb_yr);
    wrt_r_df_col("total.L.knum", L_total_knum_yr);
    wrt_r_df_col("total.D.klb",D_total_klb_yr);
    wrt_r_df_col("total.D.knum",D_total_knum_yr);

    wrt_r_df_col("U.sBT.ob", obs_cpue_sBT);
    wrt_r_df_col("U.sBT.pr", pred_cpue_sBT);
    wrt_r_df_col("cv.U.sBT", obs_cv_cpue_sBT);
    wrt_r_df_col("U.sTV.ob", obs_cpue_sTV); // use observed values, for which negative values were replaced with NA inside the tpl
    wrt_r_df_col("U.sTV.pr", pred_cpue_sTV);
    wrt_r_df_col("cv.U.sTV", obs_cv_cpue_corrected_sTV);
	//wrt_r_df_col("U.Vid.ob", obs_Vid_cpue);
    //wrt_r_df_col("U.Vid.pr", pred_cpue_Vid);
    //wrt_r_df_col("cv.U.Vid", Vid_cpue_cv);
    wrt_r_df_col("U.cHL.ob", obs_cpue_cHL);
    wrt_r_df_col("U.cHL.pr", pred_cpue_cHL);
    wrt_r_df_col("cv.U.cHL", obs_cv_cpue_cHL);
    wrt_r_df_col("U.rHB.ob", obs_cpue_rHB);
    wrt_r_df_col("U.rHB.pr", pred_cpue_rHB);
    wrt_r_df_col("cv.U.rHB", obs_cv_cpue_rHB);
    //wrt_r_df_col("U.rHD.ob", obs_rHD_cpue);
    //wrt_r_df_col("U.rHD.pr", pred_cpue_rHD);
    //wrt_r_df_col("cv.U.rHD", rHD_cpue_cv);

    wrt_r_df_col("q.cHL", q_cHL);
    wrt_r_df_col("q.rHB", q_rHB);
    
    wrt_r_df_col("q.cHL.rate.mult",q_rate_fcn_cHL);
    wrt_r_df_col("q.rHB.rate.mult",q_rate_fcn_rHB);
    //wrt_r_df_col("q.rHD.rate.mult",q_rate_fcn_rHD);
    wrt_r_df_col("q.DD.mult", q_DD_fcn);
    wrt_r_df_col("q.DD.B.4plus", B_q_DD);
    wrt_r_df_col("q.cHL.RW.log.dev",q_RW_log_dev_cHL);
    wrt_r_df_col("q.rHB.RW.log.dev",q_RW_log_dev_rHB);
    //wrt_r_df_col("q.rHD.RW.log.dev",q_RW_log_dev_rHD);

    wrt_r_df_col("L.cHL.ob", obs_L_cHL);
    wrt_r_df_col("L.cHL.pr", pred_L_cHL_klb);
    wrt_r_df_col("cv.L.cHL", obs_cv_L_cHL);
    
    wrt_r_df_col("L.cPT.ob", obs_L_cPT);
    wrt_r_df_col("L.cPT.pr", pred_L_cPT_klb);
    wrt_r_df_col("cv.L.cPT", obs_cv_L_cPT);
    
    wrt_r_df_col("L.cTW.ob", obs_L_cTW);
    wrt_r_df_col("L.cTW.pr", pred_L_cTW_klb);
    wrt_r_df_col("cv.L.cTW", obs_cv_L_cTW);
    
    wrt_r_df_col("L.rHB.ob", obs_L_rHB);
    wrt_r_df_col("L.rHB.pr", pred_L_rHB_knum);
    wrt_r_df_col("cv.L.rHB", obs_cv_L_rHB);
    
    wrt_r_df_col("L.rGN.ob", obs_L_rGN);
    wrt_r_df_col("L.rGN.pr", pred_L_rGN_knum);
    wrt_r_df_col("cv.L.rGN", obs_cv_L_rGN);
    
    wrt_r_df_col("L.rIN.ob", obs_L_rIN);
    wrt_r_df_col("L.rIN.pr", pred_L_rIN_knum);
    wrt_r_df_col("cv.L.rIN", obs_cv_L_rIN);    

    // wrt_r_df_col("D.cHD.ob", obs_D_cHD);
    // wrt_r_df_col("D.cHD.pr", pred_D_cHD_knum);
    // 
    // wrt_r_df_col("D.cPD.ob", obs_D_cPD);
    // wrt_r_df_col("D.cPD.pr", pred_D_cPD_knum);
    
    wrt_r_df_col("D.cGD.ob", obs_D_cGD);
    wrt_r_df_col("D.cGD.pr", pred_D_cGD_knum);
    wrt_r_df_col("cv.D.cGD", obs_cv_D_cGD);
    
    wrt_r_df_col("D.rHD.ob", obs_D_rHD);
    wrt_r_df_col("D.rHD.pr", pred_D_rHD_knum);
    wrt_r_df_col("cv.D.rHD", obs_cv_D_rHD);
    
    wrt_r_df_col("D.rGD.ob", obs_D_rGD);
    wrt_r_df_col("D.rGD.pr", pred_D_rGD_knum);
    wrt_r_df_col("cv.D.rGD", obs_cv_D_rGD);
    
    wrt_r_df_col("D.rID.ob", obs_D_rID);
    wrt_r_df_col("D.rID.pr", pred_D_rID_knum);
    wrt_r_df_col("cv.D.rID", obs_cv_D_rID);

    //comp sample sizes
    wrt_r_df_col("lcomp.sBT.n", nsamp_allyr_lenc_sBT);
    wrt_r_df_col("lcomp.cHL.n", nsamp_allyr_lenc_cHL);
    wrt_r_df_col("lcomp.cPT.n", nsamp_allyr_lenc_cPT);
    wrt_r_df_col("lcomp.rHB.n", nsamp_allyr_lenc_rHB);
    wrt_r_df_col("lcomp.rGN.n", nsamp_allyr_lenc_rGN);
    wrt_r_df_col("lcomp.rHD.n", nsamp_allyr_lenc_rHD);

    wrt_r_df_col("acomp.sBT.n", nsamp_allyr_agec_sBT);
    wrt_r_df_col("acomp.sTV.n", nsamp_allyr_agec_sTV);
    wrt_r_df_col("acomp.cHL.n", nsamp_allyr_agec_cHL);
    wrt_r_df_col("acomp.cPT.n", nsamp_allyr_agec_cPT);
    wrt_r_df_col("acomp.rHB.n", nsamp_allyr_agec_rHB);
    //wrt_r_df_col("acomp.rGN.n", nsamp_allyr_agec_rGN);

    wrt_r_df_col("lcomp.sBT.nfish", nfish_allyr_lenc_sBT);
    wrt_r_df_col("lcomp.cHL.nfish", nfish_allyr_lenc_cHL);
    wrt_r_df_col("lcomp.cPT.nfish", nfish_allyr_lenc_cPT);
    wrt_r_df_col("lcomp.rHB.nfish", nfish_allyr_lenc_rHB);
    wrt_r_df_col("lcomp.rGN.nfish", nfish_allyr_lenc_rGN);
    wrt_r_df_col("lcomp.rHD.nfish", nfish_allyr_lenc_rHD);

    wrt_r_df_col("acomp.sBT.nfish", nfish_allyr_agec_sBT);
    wrt_r_df_col("acomp.sTV.nfish", nfish_allyr_agec_sTV);
    wrt_r_df_col("acomp.cHL.nfish", nfish_allyr_agec_cHL);
    wrt_r_df_col("acomp.cPT.nfish", nfish_allyr_agec_cPT);
    wrt_r_df_col("acomp.rHB.nfish", nfish_allyr_agec_rHB);

    wrt_r_df_col("lcomp.sBT.neff", out_neff_allyr_lenc_sBT);
    wrt_r_df_col("lcomp.cHL.neff", out_neff_allyr_lenc_cHL);
    wrt_r_df_col("lcomp.cPT.neff", out_neff_allyr_lenc_cPT);
    wrt_r_df_col("lcomp.rHB.neff", out_neff_allyr_lenc_rHB);
    wrt_r_df_col("lcomp.rGN.neff", out_neff_allyr_lenc_rGN);
    wrt_r_df_col("lcomp.rHD.neff", out_neff_allyr_lenc_rHD);

    wrt_r_df_col("acomp.sBT.neff", out_neff_allyr_agec_sBT);
    wrt_r_df_col("acomp.sTV.neff", out_neff_allyr_agec_sTV);
    wrt_r_df_col("acomp.cHL.neff", out_neff_allyr_agec_cHL);
    wrt_r_df_col("acomp.cPT.neff", out_neff_allyr_agec_cPT);
    wrt_r_df_col("acomp.rHB.neff", out_neff_allyr_agec_rHB);
close_r_df();

// DATA FRAME of L and D time series by fishery
open_r_df("LD.pr.tseries", styr, (endyr+1), 2);
	wrt_r_namevector(styr,(endyr+1));
	wrt_r_df_col("year", styr,(endyr+1));

    wrt_r_df_col("L.cHL.klb", pred_L_cHL_klb);
    wrt_r_df_col("L.cHL.knum", pred_L_cHL_knum);
    wrt_r_df_col("L.cPT.klb", pred_L_cPT_klb);
    wrt_r_df_col("L.cPT.knum", pred_L_cPT_knum);
    wrt_r_df_col("L.cTW.klb", pred_L_cTW_klb);
    wrt_r_df_col("L.cTW.knum", pred_L_cTW_knum);
    wrt_r_df_col("L.rHB.klb", pred_L_rHB_klb);
    wrt_r_df_col("L.rHB.knum", pred_L_rHB_knum);
    wrt_r_df_col("L.rGN.klb", pred_L_rGN_klb);
    wrt_r_df_col("L.rGN.knum", pred_L_rGN_knum);
    wrt_r_df_col("L.rIN.klb", pred_L_rIN_klb);
    wrt_r_df_col("L.rIN.knum", pred_L_rIN_knum);    

    wrt_r_df_col("D.cGD.klb", pred_D_cGD_klb);
    wrt_r_df_col("D.cGD.knum", pred_D_cGD_knum);
    wrt_r_df_col("D.rHD.klb", pred_D_rHD_klb);
    wrt_r_df_col("D.rHD.knum", pred_D_rHD_knum);
    wrt_r_df_col("D.rGD.klb", pred_D_rGD_klb);
    wrt_r_df_col("D.rGD.knum", pred_D_rGD_knum);
    wrt_r_df_col("D.rID.klb", pred_D_rID_klb);
    wrt_r_df_col("D.rID.knum", pred_D_rID_knum);    
    
close_r_df();


open_r_df("a.series", 1, nages, 2);
	wrt_r_namevector(agebins);
	wrt_r_df_col("age", agebins);
	wrt_r_df_col("length", meanlen_TL);
	wrt_r_df_col("length.cv", len_cv);
	wrt_r_df_col("length.sd", len_sd);
	wrt_r_df_col("weight", wgt_lb);     //for FishGraph
	wrt_r_df_col("wgt.klb", wgt_klb);
	wrt_r_df_col("wgt.mt", wgt_mt);
	wrt_r_df_col("wgt.wgted.L.klb", wgt_wgted_L_klb);
	wrt_r_df_col("wgt.wgted.D.klb", wgt_wgted_D_klb);
	wrt_r_df_col("prop.female", prop_f);
	wrt_r_df_col("mat.female", maturity_f);
	wrt_r_df_col("fecundity.scaled", fecundity);
	wrt_r_df_col("reprod", reprod);
	wrt_r_df_col("reprod.MatFemB", reprod2);
	wrt_r_df_col("M", M);
	wrt_r_df_col("F.initial", F_initial);
	wrt_r_df_col("Z.initial", Z_initial);
	wrt_r_df_col("Nage.eq.init",N_initial_eq);
	wrt_r_df_col("log.Nage2plus.init.dev",log_Nage_dev_output);
close_r_df();

open_r_df("eq.series", 1, n_iter_msy, 2);
	wrt_r_namevector(1,n_iter_msy);
	wrt_r_df_col("F.eq", F_msy);
	wrt_r_df_col("spr.eq", spr_msy);
	wrt_r_df_col("R.eq", R_eq);
	wrt_r_df_col("SSB.eq", SSB_eq);
	wrt_r_df_col("B.eq", B_eq);
	wrt_r_df_col("L.eq.klb", L_eq_klb);
	wrt_r_df_col("L.eq.knum", L_eq_knum);
	wrt_r_df_col("D.eq.knum", D_eq_knum);
	wrt_r_df_col("D.eq.klb", D_eq_klb);
	wrt_r_df_col("T.eq.knum", T_eq_knum);
	wrt_r_df_col("T.eq.klb", T_eq_klb);
close_r_df();

open_r_df("pr.series", 1, n_iter_spr, 2);
	wrt_r_namevector(1,n_iter_spr);
	wrt_r_df_col("F.spr", F_spr);
	wrt_r_df_col("spr", spr_spr);
	wrt_r_df_col("SPR", spr_ratio);
	wrt_r_df_col("ypr.lb", L_spr);
	wrt_r_df_col("total.lb", T_spr);
close_r_df();


open_r_list("CLD.est.mats");

    open_r_matrix("Lw.cHL");
        wrt_r_matrix(L_cHL_klb, 1,1);
    close_r_matrix();

    open_r_matrix("Lw.cPT");
        wrt_r_matrix(L_cPT_klb, 1,1);
    close_r_matrix();

    open_r_matrix("Lw.cTW");
        wrt_r_matrix(L_cTW_klb, 1,1);
    close_r_matrix();

    open_r_matrix("Lw.rHB");
        wrt_r_matrix(L_rHB_klb, 1,1);
    close_r_matrix();

    open_r_matrix("Lw.rGN");
        wrt_r_matrix(L_rGN_klb, 1,1);
    close_r_matrix();
    
    open_r_matrix("Lw.rIN");
    wrt_r_matrix(L_rIN_klb, 1,1);
    close_r_matrix();    

    open_r_matrix("Lw.total");
        wrt_r_matrix(L_total_klb, 1,1);
    close_r_matrix();

    open_r_matrix("Dw.cGD");
        wrt_r_matrix(D_cGD_klb, 1,1);
    close_r_matrix();

    open_r_matrix("Dw.rHD");
        wrt_r_matrix(D_rHD_klb, 1,1);
    close_r_matrix();

    open_r_matrix("Dw.rGD");
        wrt_r_matrix(D_rGD_klb, 1,1);
    close_r_matrix();
    
    open_r_matrix("Dw.rID");
    wrt_r_matrix(D_rID_klb, 1,1);
    close_r_matrix();    

    open_r_matrix("Dw.total");
        wrt_r_matrix(D_total_klb, 1,1);
    close_r_matrix();


    open_r_matrix("Ln.cHL");
        wrt_r_matrix(L_cHL_num, 1,1);
    close_r_matrix();

    open_r_matrix("Ln.cPT");
        wrt_r_matrix(L_cPT_num, 1,1);
    close_r_matrix();

    open_r_matrix("Ln.cTW");
        wrt_r_matrix(L_cTW_num, 1,1);
    close_r_matrix();

    open_r_matrix("Ln.rHB");
        wrt_r_matrix(L_rHB_num, 1,1);
    close_r_matrix();

    open_r_matrix("Ln.rGN");
        wrt_r_matrix(L_rGN_num, 1,1);
    close_r_matrix();
    
    open_r_matrix("Ln.rIN");
    wrt_r_matrix(L_rIN_num, 1,1);
    close_r_matrix();    

    open_r_matrix("Ln.total");
        wrt_r_matrix(L_total_num, 1,1);
    close_r_matrix();

    open_r_matrix("Dn.cGD");
        wrt_r_matrix(D_cGD_num,1,1);
    close_r_matrix();

    open_r_matrix("Dn.rHD");
        wrt_r_matrix(D_rHD_num, 1,1);
    close_r_matrix();

    open_r_matrix("Dn.rGD");
        wrt_r_matrix(D_rGD_num, 1,1);
    close_r_matrix();
    
    open_r_matrix("Dn.rID");
    wrt_r_matrix(D_rID_num, 1,1);
    close_r_matrix();    

    open_r_matrix("Dn.total");
        wrt_r_matrix(D_total_num, 1,1);
    close_r_matrix();

close_r_list();

  //LIST object of age error matrix
open_r_list("age.error");

    open_r_matrix("error.mat");
    wrt_r_matrix(age_error, 2, 2);
    wrt_r_namevector(agebins);
    wrt_r_namevector(agebins);
    close_r_matrix();
close_r_list();

// DATA FRAME of projection time series
open_r_df("proj.t.series", styr_proj, endyr_proj, 2);
  wrt_r_namevector(styr_proj,endyr_proj);
  wrt_r_df_col("year", styr_proj,endyr_proj);
  wrt_r_df_col("F.proj", F_proj);
  wrt_r_df_col("SSB.proj", SSB_proj);
  wrt_r_df_col("B.proj", B_proj);
  wrt_r_df_col("R.proj", R_proj);
  wrt_r_df_col("L.knum.proj", L_knum_proj);
  wrt_r_df_col("L.klb.proj", L_klb_proj);
 close_r_df();

 open_r_matrix("proj.N.age");
    wrt_r_matrix(N_proj, 2, 2);
    wrt_r_namevector(styr_proj, endyr_proj);
    wrt_r_namevector(agebins);
 close_r_matrix();


close_r_file();
