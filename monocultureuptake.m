function [model]=monocultureuptake(model)
%carbon substrates
%in excess-nutrient state uncomment and use the following bounds
% model=changeRxnBounds(model,{'EX_xyl_D(e)'},-10,'l');
% model=changeRxnBounds(model,{'EX_glc_D(e)'},-30,'l');

%in minimal-nutrient state uncomment and use the following bounds
% model=changeRxnBounds(model,{'EX_xyl_D(e)'},-1,'l');
% model=changeRxnBounds(model,{'EX_glc_D(e)'},-1,'l');

%aminoacids
model=changeRxnBounds(model,{'EX_ala_L(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_arg_L(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_asn_L(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_asp_L(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_cys_L(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_gln_L(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_glu_L(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_gly(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_his_L(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_ile_L(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_leu_L(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_lys_L(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_met_L(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_phe_L(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_pro_L(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_ser_L(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_thr_L(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_trp_L(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_tyr_L(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_val_L(e)'},-1,'l');

%other trace elements, vitamins, MRS media components 
model=changeRxnBounds(model,{'EX_k(e)'},-0.2069,'l');
model=changeRxnBounds(model,{'EX_mg2(e)'},-0.0073,'l');
model=changeRxnBounds(model,{'EX_mn2(e)'},-0.0027,'l');
model=changeRxnBounds(model,{'EX_na1(e)'},-0.5488,'l');
model=changeRxnBounds(model,{'EX_nac(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_nh3(e)'},-0.16,'l');
model=changeRxnBounds(model,{'EX_ni2(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_nmn(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_no2(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_no3(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_zn2(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_co2(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_btn(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_ca2(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_cbl1(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_cbl2(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_cit(e)'},-0.08,'l');
model=changeRxnBounds(model,{'EX_cl(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_cobalt2(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_cu2(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_fe2(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_fe3(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_fol(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_gua(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_ins(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_h(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_h2(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_h2o(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_hxan(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_pi(e)'},-0.1034,'l');
model=changeRxnBounds(model,{'EX_pnto_R(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_pydam(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_pydx(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_pydxn(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_ribflv(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_sheme(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_pheme(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_so4(e)'},-0.0099,'l');
model=changeRxnBounds(model,{'EX_spmd(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_thm(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_thymd(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_ptrc(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_ura(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_urea(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_xan(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_malt(u)'},-0.0421,'l');
model=changeRxnBounds(model,{'EX_pb(u)'},-1,'l');
model=changeRxnBounds(model,{'EX_cd2(u)'},-1,'l');

%model-specific metabolites
model=changeRxnBounds(model,{'EX_12dgr180(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_adocbl(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_Lcyst(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_26dap_M(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_2dmmq8(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_cgly(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_mqn8(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_mqn7(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_q8(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_ocdca(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_orn(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_3mop(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_glyphe(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_gthox(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_gthrd(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_4hbz(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_4abz(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_crn(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_csn(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_chol(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_ade(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_hdcea(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_ocdcea(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_ocdcya(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_ocdctr(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_lac_L(e)'},0,'b');

% sink reactions for some models
model=changeRxnBounds(model,{'sink_PGPm1[c]'},-1000,'l');
model=changeRxnBounds(model,{'sink_gthrd(c)'},-1,'l');
model=changeRxnBounds(model,{'sink_mqn7(c)'},-1000,'l');
model=changeRxnBounds(model,{'sink_mqn8(c)'},-1000,'l');
model=changeRxnBounds(model,{'sink_q8(c)'},-1000,'l');
model=changeRxnBounds(model,{'sink_sheme(c)'},-1000,'l');
model=changeRxnBounds(model,{'sink_s(c)'},-1000,'l');
model=changeRxnBounds(model,{'sink_dmbzid(c)'},-1000,'l');
model=changeRxnBounds(model,{'EX_sink1'},-1000,'l');
model=changeRxnBounds(model,{'EX_sink2'},-1000,'l');
model=changeRxnBounds(model,{'EX_sink3'},-1000,'l');
model = model;
