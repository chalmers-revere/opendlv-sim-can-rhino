//
// File: L_f_L_f_h_ang_Ccode.cpp
//
// MATLAB Coder version            : 3.2
// C/C++ source code generated on  : 25-Jan-2019 20:24:04
//

// Include Files
#include "rt_nonfinite.h"
#include "L_f_L_f_h_ang_Ccode.h"

// Function Declarations
static double rt_powd_snf(double u0, double u1);

// Function Definitions

//
// Arguments    : double u0
//                double u1
// Return Type  : double
//
static double rt_powd_snf(double u0, double u1)
{
  double y;
  /*double d0;
  double d1;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    d0 = std::abs(u0);
    d1 = std::abs(u1);
    if (rtIsInf(u1)) {
      if (d0 == 1.0) {
        y = rtNaN;
      } else if (d0 > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d1 == 0.0) {
      y = 1.0;
    } else if (d1 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = std::sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > std::floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }*/

  y = pow(u0, u1);
  return y;
}

//
// yp_dot = 1;
//  psi_dot =1;
//  epsi =1;
//  ey= 1;
//  s = 1;
//
//
//  %%states of ob:
//   pos_ob_x = 10;  pos_ob_y= 1;
//   vel_ob_x = 1; vel_ob_y= 1;
//   acc_ob_x = 1; acc_ob_y= 1;
//
//
//  %safety distance:
//  Ds = 1;
//
//  %the maximum acc:
//  a_m = 4;
//
//
//  %the longitudinal velocity, assume it is constant:
//  xp_dot= 1;
//
//
//  %constants:
//
//
//      a = 1.68;
//      b = 1.715;
//      cf =  3.4812e+05;
//      cr = 3.5537e+05;
//      m = 9840;
//      Iz = 41340;
//      psi_dot_com = 0;
// Arguments    : double yp_dot
//                double psi_dot
//                double epsi
//                double ey
//                double s
//                double pos_ob_x
//                double pos_ob_y
//                double vel_ob_x
//                double vel_ob_y
//                double acc_ob_x
//                double acc_ob_y
//                double Ds
//                double a_m
//                double xp_dot
//                double a
//                double b
//                double cf
//                double cr
//                double m
//                double Iz
//                double psi_dot_com
//                double steer
//                double acc
//                double out[8]
// Return Type  : void
//
void L_f_L_f_h_ang_Ccode(double yp_dot, double psi_dot, double epsi, double ey,
  double s, double pos_ob_x, double pos_ob_y, double vel_ob_x, double vel_ob_y,
  double acc_ob_x, double acc_ob_y, double Ds, double, double xp_dot, double a,
  double b, double cf, double cr, double m, double Iz, double psi_dot_com,
  double steer, double, double out[8])
{
  double b_a;
  double c_a;
  double d_a;
  double e_a;
  double f_a;
  double g_a;
  double h_a;
  double i_a;
  double j_a;
  double k_a;
  double l_a;
  double m_a;
  double n_a;
  double o_a;
  double p_a;
  double q_a;
  double r_a;
  double s_a;
  double t_a;
  double u_a;
  double v_a;
  double w_a;
  double x_a;
  double y_a;
  double ab_a;
  double bb_a;
  double cb_a;
  double db_a;
  double eb_a;
  double fb_a;
  double gb_a;
  double hb_a;
  double ib_a;
  double jb_a;
  double kb_a;
  double lb_a;
  double mb_a;
  double nb_a;
  double ob_a;
  double pb_a;
  double qb_a;
  double rb_a;
  double sb_a;
  double tb_a;
  double ub_a;
  double vb_a;
  double wb_a;
  double xb_a;
  double yb_a;
  double ac_a;
  double bc_a;
  double cc_a;
  double dc_a;
  double ec_a;
  double fc_a;
  double gc_a;
  double hc_a;
  double ic_a;
  double jc_a;
  double kc_a;
  double lc_a;
  double mc_a;
  double nc_a;
  double oc_a;
  double pc_a;
  double qc_a;
  double rc_a;
  double sc_a;
  double tc_a;
  double uc_a;
  double vc_a;
  double wc_a;
  double xc_a;
  double yc_a;
  double ad_a;
  double bd_a;
  double cd_a;
  double dd_a;
  double ed_a;
  double fd_a;
  double gd_a;
  double hd_a;
  double id_a;
  double jd_a;
  double kd_a;
  double ld_a;
  double md_a;
  double nd_a;
  double od_a;
  double pd_a;
  double qd_a;
  double rd_a;
  double sd_a;
  double td_a;
  double ud_a;
  double vd_a;
  double wd_a;
  double xd_a;
  double yd_a;
  double ae_a;
  double be_a;
  double ce_a;
  double de_a;
  double ee_a;
  double fe_a;
  double ge_a;
  double he_a;
  double ie_a;
  double je_a;
  double ke_a;
  double le_a;
  double me_a;
  double ne_a;
  double oe_a;
  double pe_a;
  double qe_a;
  double re_a;
  double se_a;
  double te_a;
  double ue_a;
  double ve_a;
  double we_a;
  double xe_a;
  double ye_a;
  double af_a;
  double bf_a;
  double cf_a;
  double df_a;
  double ef_a;
  double ff_a;
  double gf_a;
  double hf_a;
  double if_a;
  double jf_a;
  double kf_a;
  double lf_a;
  double mf_a;
  double nf_a;
  double of_a;
  double pf_a;
  double qf_a;
  double rf_a;
  double sf_a;
  double tf_a;
  double uf_a;
  double vf_a;
  double wf_a;
  double xf_a;
  double yf_a;
  double ag_a;
  double bg_a;
  double cg_a;
  double dg_a;
  double eg_a;
  double fg_a;
  double gg_a;
  double hg_a;
  double ig_a;
  double jg_a;
  double kg_a;
  double lg_a;
  double mg_a;
  double ng_a;
  double og_a;
  double pg_a;
  double qg_a;
  double rg_a;
  double sg_a;
  double tg_a;
  double ug_a;
  double vg_a;
  double wg_a;
  double xg_a;
  double yg_a;
  double ah_a;
  double bh_a;
  double ch_a;
  double dh_a;
  double eh_a;
  double fh_a;
  double gh_a;
  double hh_a;
  double ih_a;
  double jh_a;
  double kh_a;
  double lh_a;
  double mh_a;
  double nh_a;
  double oh_a;
  double ph_a;
  double qh_a;
  double rh_a;
  double sh_a;
  double th_a;
  double uh_a;
  double vh_a;
  double wh_a;
  double xh_a;
  double yh_a;
  double ai_a;
  double bi_a;
  double ci_a;
  double di_a;
  double ei_a;
  double fi_a;
  double gi_a;
  double hi_a;
  double ii_a;
  double ji_a;
  double ki_a;
  double li_a;
  double mi_a;
  double ni_a;
  double oi_a;
  double pi_a;
  double qi_a;
  double ri_a;
  double si_a;
  double ti_a;
  double ui_a;
  double vi_a;
  double wi_a;
  double xi_a;
  double yi_a;
  double aj_a;
  double bj_a;
  double cj_a;
  double dj_a;
  double ej_a;
  double fj_a;
  double gj_a;
  double hj_a;
  double ij_a;
  double jj_a;
  double kj_a;
  double lj_a;
  double mj_a;
  double nj_a;
  double oj_a;
  double pj_a;
  double qj_a;
  double rj_a;
  double sj_a;
  double tj_a;
  double uj_a;
  double vj_a;
  double wj_a;
  double xj_a;
  double yj_a;
  double ak_a;
  double bk_a;
  double ck_a;
  double dk_a;
  double ek_a;
  double fk_a;
  double gk_a;
  double hk_a;
  double ik_a;
  double jk_a;
  double kk_a;
  double lk_a;
  double mk_a;
  double nk_a;
  double ok_a;
  double pk_a;
  double qk_a;
  double rk_a;
  double sk_a;
  double tk_a;
  double uk_a;
  double vk_a;
  double wk_a;
  double xk_a;
  double yk_a;
  double al_a;
  double bl_a;
  double cl_a;
  double dl_a;
  double el_a;
  double fl_a;
  double gl_a;
  double hl_a;
  double il_a;
  double jl_a;
  double kl_a;
  double ll_a;
  double ml_a;
  double nl_a;
  double ol_a;
  double pl_a;
  double ql_a;
  double rl_a;
  double sl_a;
  double tl_a;
  double ul_a;
  double vl_a;
  double wl_a;
  double xl_a;
  double yl_a;
  double am_a;
  double bm_a;
  double cm_a;
  double dm_a;
  double em_a;
  double fm_a;
  double gm_a;
  double hm_a;
  double im_a;
  double jm_a;
  double km_a;
  double lm_a;
  double mm_a;
  double nm_a;
  double om_a;
  double pm_a;
  double qm_a;
  double rm_a;
  double sm_a;
  double tm_a;
  double um_a;
  double vm_a;
  double wm_a;
  double xm_a;
  double ym_a;
  double an_a;
  double bn_a;
  double cn_a;
  double dn_a;
  double en_a;
  double fn_a;
  double gn_a;
  double hn_a;
  double in_a;
  double jn_a;
  double kn_a;
  double ln_a;
  double mn_a;
  double nn_a;
  double on_a;
  double pn_a;
  double qn_a;
  double rn_a;
  double sn_a;
  double tn_a;
  double un_a;
  double vn_a;
  double wn_a;
  double xn_a;
  double yn_a;
  double ao_a;
  double bo_a;
  double co_a;
  double do_a;
  double eo_a;
  double fo_a;
  double go_a;
  double ho_a;
  double io_a;
  double jo_a;
  double ko_a;
  double lo_a;
  double mo_a;
  double no_a;
  double oo_a;
  double po_a;
  double qo_a;
  double ro_a;
  double so_a;
  double to_a;
  double uo_a;
  double vo_a;
  double wo_a;
  double xo_a;
  double yo_a;
  double ap_a;
  double bp_a;
  double cp_a;
  double dp_a;
  double ep_a;
  double fp_a;
  double gp_a;
  double hp_a;
  double ip_a;
  double jp_a;
  double kp_a;
  double lp_a;
  double mp_a;
  double np_a;
  double op_a;
  double pp_a;
  double qp_a;
  double rp_a;
  double sp_a;
  double tp_a;
  double up_a;
  double vp_a;
  double wp_a;
  double xp_a;
  double yp_a;
  double aq_a;
  double bq_a;
  double cq_a;
  double dq_a;
  double eq_a;
  double fq_a;
  double gq_a;
  double hq_a;
  double iq_a;
  double jq_a;
  double kq_a;
  double lq_a;
  double mq_a;
  double nq_a;
  double oq_a;
  double pq_a;
  double qq_a;
  double rq_a;
  double sq_a;
  double tq_a;
  double uq_a;
  double vq_a;
  double wq_a;
  double xq_a;
  double yq_a;
  double ar_a;
  double br_a;
  double cr_a;
  double dr_a;
  double er_a;
  double fr_a;
  double gr_a;
  double hr_a;
  double ir_a;
  double jr_a;
  double kr_a;
  double lr_a;
  double mr_a;
  double nr_a;
  double or_a;
  double pr_a;
  double qr_a;
  double rr_a;
  double sr_a;
  double tr_a;
  double ur_a;
  double vr_a;
  double wr_a;
  double xr_a;
  double yr_a;
  double as_a;
  double bs_a;
  double cs_a;
  double ds_a;
  double es_a;
  double fs_a;
  double gs_a;
  double hs_a;
  double is_a;
  double js_a;
  double ks_a;
  double ls_a;
  double ms_a;
  double ns_a;
  double os_a;
  double ps_a;
  double qs_a;
  double rs_a;
  double ss_a;
  double ts_a;
  double us_a;
  double vs_a;
  double ws_a;
  double xs_a;
  double ys_a;
  double at_a;
  double bt_a;
  double ct_a;
  double dt_a;
  double et_a;
  double ft_a;
  double gt_a;
  double ht_a;
  double it_a;
  double jt_a;
  double kt_a;
  double lt_a;
  double mt_a;
  double nt_a;
  double ot_a;
  double pt_a;
  double qt_a;
  double rt_a;
  double st_a;
  double tt_a;
  double ut_a;
  double vt_a;
  double wt_a;
  double xt_a;
  double yt_a;
  double au_a;
  double bu_a;
  double cu_a;
  double du_a;
  double eu_a;
  double fu_a;
  double gu_a;
  double hu_a;
  double iu_a;
  double ju_a;
  double ku_a;
  double lu_a;
  double mu_a;
  double nu_a;
  double ou_a;
  double pu_a;
  double qu_a;
  double ru_a;
  double su_a;
  double tu_a;
  double uu_a;
  double vu_a;
  double wu_a;
  double xu_a;
  double yu_a;
  double av_a;
  double bv_a;
  double cv_a;
  double dv_a;
  double ev_a;
  double fv_a;
  double gv_a;
  double hv_a;
  double iv_a;
  double jv_a;
  double kv_a;
  double lv_a;
  double mv_a;
  double nv_a;
  double ov_a;
  double pv_a;
  double qv_a;
  double rv_a;
  double sv_a;
  double tv_a;
  double uv_a;
  double vv_a;
  double wv_a;
  double xv_a;
  double yv_a;
  double aw_a;
  double bw_a;
  double cw_a;
  double dw_a;
  double ew_a;
  double fw_a;
  double gw_a;
  double hw_a;
  double iw_a;
  double jw_a;
  double kw_a;
  double lw_a;
  double mw_a;
  double nw_a;
  double ow_a;
  double pw_a;
  double qw_a;
  double rw_a;
  double sw_a;
  double tw_a;
  double uw_a;
  double vw_a;
  double ww_a;
  double xw_a;
  double yw_a;
  double ax_a;
  double bx_a;
  double cx_a;
  double dx_a;
  double ex_a;
  double fx_a;
  double gx_a;
  double hx_a;
  double ix_a;
  double jx_a;
  double kx_a;
  double lx_a;
  double mx_a;
  double nx_a;
  double ox_a;
  double px_a;
  double qx_a;
  double rx_a;
  double sx_a;
  double tx_a;
  double ux_a;
  double vx_a;
  double wx_a;
  double xx_a;
  double yx_a;
  double ay_a;
  double by_a;
  double cy_a;
  double dy_a;
  double ey_a;
  double fy_a;
  double gy_a;
  double hy_a;
  double iy_a;
  double jy_a;
  double ky_a;
  double ly_a;
  double my_a;
  double ny_a;
  double oy_a;
  double py_a;
  double qy_a;
  double ry_a;
  double sy_a;
  double ty_a;
  double uy_a;
  double vy_a;
  double wy_a;
  double xy_a;
  double yy_a;
  double aab_a;
  double bab_a;
  double cab_a;
  double dab_a;
  double eab_a;
  double fab_a;
  double gab_a;
  double hab_a;
  double iab_a;
  double jab_a;
  double kab_a;
  double lab_a;
  double mab_a;
  double nab_a;
  double oab_a;
  double pab_a;
  double qab_a;
  double rab_a;
  double sab_a;
  double tab_a;
  double uab_a;
  double vab_a;
  double wab_a;
  double xab_a;
  double yab_a;
  double abb_a;
  double bbb_a;
  double cbb_a;
  double dbb_a;
  double ebb_a;
  double fbb_a;
  double gbb_a;
  double hbb_a;
  double ibb_a;
  double jbb_a;
  double kbb_a;
  double lbb_a;
  double mbb_a;
  double nbb_a;
  double obb_a;
  double pbb_a;
  double qbb_a;
  double rbb_a;
  double sbb_a;
  double tbb_a;
  double ubb_a;
  double vbb_a;
  double wbb_a;
  double xbb_a;
  double ybb_a;
  double acb_a;
  double bcb_a;
  double ccb_a;
  double dcb_a;
  double ecb_a;
  double fcb_a;
  double gcb_a;
  double hcb_a;
  double icb_a;
  double jcb_a;
  double kcb_a;
  double lcb_a;
  double mcb_a;
  double ncb_a;
  double ocb_a;
  double pcb_a;
  double qcb_a;
  double rcb_a;
  double scb_a;
  double tcb_a;
  double ucb_a;
  double vcb_a;
  double wcb_a;
  double xcb_a;
  double ycb_a;
  double adb_a;
  double bdb_a;
  double cdb_a;
  double ddb_a;
  double edb_a;
  double fdb_a;
  double gdb_a;
  double hdb_a;
  double idb_a;
  double jdb_a;
  double kdb_a;
  double ldb_a;
  double mdb_a;
  double ndb_a;
  double odb_a;
  double pdb_a;
  double qdb_a;
  double rdb_a;
  double sdb_a;
  double tdb_a;
  double udb_a;
  double vdb_a;
  double wdb_a;
  double xdb_a;
  double ydb_a;
  double aeb_a;
  double beb_a;
  double ceb_a;
  double deb_a;
  double eeb_a;
  double feb_a;
  double geb_a;
  double heb_a;
  double ieb_a;
  double jeb_a;
  double keb_a;
  double leb_a;
  double meb_a;
  double neb_a;
  double oeb_a;
  double peb_a;
  double qeb_a;
  double reb_a;
  double seb_a;
  double teb_a;
  double ueb_a;
  double veb_a;
  double web_a;
  double xeb_a;
  double yeb_a;
  double afb_a;
  double bfb_a;
  double cfb_a;
  double dfb_a;
  double efb_a;
  double ffb_a;
  double gfb_a;
  double hfb_a;
  double ifb_a;
  double jfb_a;
  double kfb_a;
  double lfb_a;
  double mfb_a;
  double nfb_a;
  double ofb_a;
  double pfb_a;
  double qfb_a;
  double rfb_a;
  double sfb_a;
  double tfb_a;
  double ufb_a;
  double vfb_a;
  double wfb_a;
  double xfb_a;
  double yfb_a;
  double agb_a;
  double bgb_a;
  double cgb_a;
  double dgb_a;
  double egb_a;
  double fgb_a;
  double ggb_a;
  double hgb_a;
  double igb_a;
  double jgb_a;
  double kgb_a;
  double lgb_a;
  double mgb_a;
  double ngb_a;
  double ogb_a;
  double pgb_a;
  double qgb_a;
  double rgb_a;
  double sgb_a;
  double tgb_a;
  double ugb_a;
  double vgb_a;
  double wgb_a;
  double xgb_a;
  double ygb_a;
  double ahb_a;
  double bhb_a;
  double chb_a;
  double dhb_a;
  double ehb_a;
  double fhb_a;
  double ghb_a;
  double hhb_a;
  double ihb_a;
  double jhb_a;
  double khb_a;
  double lhb_a;
  double mhb_a;
  double nhb_a;
  double ohb_a;
  double phb_a;
  double qhb_a;
  double rhb_a;
  double shb_a;
  double thb_a;
  double uhb_a;
  double vhb_a;
  double whb_a;
  double xhb_a;
  double yhb_a;
  double aib_a;
  double bib_a;
  double cib_a;
  double dib_a;
  double eib_a;
  double fib_a;
  double gib_a;
  double hib_a;
  double iib_a;
  double jib_a;
  double kib_a;
  double lib_a;
  double mib_a;
  double nib_a;
  double oib_a;
  double pib_a;
  double qib_a;
  double rib_a;
  double sib_a;
  double tib_a;
  double uib_a;
  double vib_a;
  double wib_a;
  double xib_a;
  double yib_a;
  double ajb_a;
  double bjb_a;
  double cjb_a;
  double djb_a;
  double ejb_a;
  double fjb_a;
  double gjb_a;
  double hjb_a;
  double ijb_a;
  double jjb_a;
  double kjb_a;
  double ljb_a;
  double mjb_a;
  double njb_a;
  double ojb_a;
  double pjb_a;
  double qjb_a;
  double rjb_a;
  double sjb_a;
  double tjb_a;
  double ujb_a;
  double vjb_a;
  double wjb_a;
  double xjb_a;
  double yjb_a;
  double akb_a;
  double bkb_a;
  double ckb_a;
  double dkb_a;
  double ekb_a;
  double fkb_a;
  double gkb_a;
  double hkb_a;
  double ikb_a;
  double jkb_a;
  double kkb_a;
  double lkb_a;
  double mkb_a;
  double nkb_a;
  double okb_a;
  double pkb_a;
  double qkb_a;
  double rkb_a;
  double skb_a;
  double tkb_a;
  double ukb_a;
  double vkb_a;
  double wkb_a;
  double xkb_a;
  double ykb_a;
  double alb_a;
  double blb_a;
  double clb_a;
  double dlb_a;
  double elb_a;
  double flb_a;
  double glb_a;
  double hlb_a;
  double ilb_a;
  double jlb_a;
  double klb_a;
  double llb_a;
  double mlb_a;
  double nlb_a;
  double olb_a;
  double plb_a;
  double qlb_a;
  double rlb_a;
  double slb_a;
  double tlb_a;
  double ulb_a;
  double vlb_a;
  double wlb_a;
  double xlb_a;
  double ylb_a;
  double amb_a;
  double bmb_a;
  double cmb_a;
  double dmb_a;
  double emb_a;
  double fmb_a;
  double gmb_a;
  double hmb_a;
  double imb_a;
  double jmb_a;
  double kmb_a;
  double lmb_a;
  double mmb_a;
  double nmb_a;
  double omb_a;
  double pmb_a;
  double qmb_a;
  double rmb_a;
  double smb_a;
  double tmb_a;
  double umb_a;
  double vmb_a;
  double wmb_a;
  double xmb_a;
  double ymb_a;
  double anb_a;
  double bnb_a;
  double cnb_a;
  double dnb_a;
  double enb_a;
  double fnb_a;
  double gnb_a;
  double hnb_a;
  double inb_a;
  double jnb_a;
  double knb_a;
  double lnb_a;
  double mnb_a;
  double nnb_a;
  double onb_a;
  double pnb_a;
  double qnb_a;
  double rnb_a;
  double snb_a;
  double tnb_a;
  double unb_a;
  double vnb_a;
  double wnb_a;
  double xnb_a;
  double ynb_a;
  double aob_a;
  double bob_a;
  double cob_a;
  double dob_a;
  double eob_a;
  double fob_a;
  double gob_a;
  double hob_a;
  double iob_a;
  double job_a;
  double kob_a;
  double lob_a;
  double mob_a;
  double nob_a;
  double oob_a;
  double pob_a;
  double qob_a;
  double rob_a;
  double sob_a;
  double tob_a;
  double uob_a;
  double vob_a;
  double wob_a;
  double xob_a;
  double yob_a;
  double apb_a;
  double bpb_a;
  double cpb_a;
  double dpb_a;
  double epb_a;
  double fpb_a;
  double gpb_a;
  double hpb_a;
  double ipb_a;
  double jpb_a;
  double kpb_a;
  double lpb_a;
  double mpb_a;
  double npb_a;
  double opb_a;
  double ppb_a;
  double qpb_a;
  double rpb_a;
  double spb_a;
  double tpb_a;
  double upb_a;
  double vpb_a;
  double wpb_a;
  double xpb_a;
  double ypb_a;
  double aqb_a;
  double bqb_a;
  double cqb_a;
  double dqb_a;
  double eqb_a;
  double fqb_a;
  double gqb_a;
  double hqb_a;
  double iqb_a;
  double jqb_a;
  double kqb_a;
  double lqb_a;
  double mqb_a;
  double nqb_a;
  double oqb_a;
  double pqb_a;
  double qqb_a;
  double rqb_a;
  double sqb_a;
  double tqb_a;
  double uqb_a;
  double vqb_a;
  double wqb_a;
  double xqb_a;
  double yqb_a;
  double arb_a;
  double brb_a;
  double crb_a;
  double drb_a;
  double erb_a;
  double frb_a;
  double grb_a;
  double hrb_a;
  double irb_a;
  double jrb_a;
  double krb_a;
  double lrb_a;
  double mrb_a;
  double nrb_a;
  double orb_a;
  double prb_a;
  double qrb_a;
  double rrb_a;
  double srb_a;
  double trb_a;
  double urb_a;
  double vrb_a;
  double wrb_a;
  double xrb_a;
  double yrb_a;
  double asb_a;
  double bsb_a;
  double csb_a;
  double dsb_a;
  double esb_a;
  double fsb_a;
  double gsb_a;
  double hsb_a;
  double isb_a;
  double jsb_a;
  double ksb_a;
  double lsb_a;
  double msb_a;
  double nsb_a;
  double osb_a;
  double psb_a;
  double qsb_a;
  double rsb_a;
  double ssb_a;
  double tsb_a;
  double usb_a;
  double vsb_a;
  double wsb_a;
  double xsb_a;
  double ysb_a;
  double atb_a;
  double btb_a;
  double ctb_a;
  double dtb_a;
  double etb_a;
  double ftb_a;
  double gtb_a;
  double htb_a;
  double itb_a;
  double jtb_a;
  double ktb_a;
  double ltb_a;
  double mtb_a;
  double ntb_a;
  double otb_a;
  double ptb_a;
  double qtb_a;
  double rtb_a;
  double stb_a;
  double ttb_a;
  double utb_a;
  double vtb_a;
  double wtb_a;
  double xtb_a;
  double ytb_a;
  double aub_a;
  double bub_a;
  double cub_a;
  double dub_a;
  double eub_a;
  double fub_a;
  double gub_a;
  double hub_a;
  double iub_a;
  double jub_a;
  double kub_a;
  double lub_a;
  double mub_a;
  double nub_a;
  double oub_a;
  double pub_a;
  double qub_a;
  double rub_a;
  double sub_a;
  double tub_a;
  double uub_a;
  double vub_a;
  double wub_a;
  double xub_a;
  double yub_a;
  double avb_a;
  double bvb_a;
  double cvb_a;
  double dvb_a;
  double evb_a;
  double fvb_a;
  double gvb_a;
  double hvb_a;
  double ivb_a;
  double jvb_a;
  double kvb_a;
  double lvb_a;
  double mvb_a;
  double nvb_a;
  double ovb_a;
  double pvb_a;
  double qvb_a;
  double rvb_a;
  double svb_a;
  double tvb_a;
  double uvb_a;
  double vvb_a;
  double wvb_a;
  double xvb_a;
  double yvb_a;
  double awb_a;
  double bwb_a;
  double cwb_a;
  double dwb_a;
  double ewb_a;
  double fwb_a;
  double gwb_a;
  double hwb_a;
  double iwb_a;
  double jwb_a;
  double kwb_a;
  double lwb_a;
  double mwb_a;
  double nwb_a;
  double owb_a;
  double pwb_a;
  double qwb_a;
  double rwb_a;
  double swb_a;
  double twb_a;
  double uwb_a;
  double vwb_a;
  double wwb_a;
  double xwb_a;
  double ywb_a;
  double axb_a;
  double bxb_a;
  double cxb_a;
  double dxb_a;
  double exb_a;
  double fxb_a;
  double gxb_a;
  double hxb_a;
  double ixb_a;
  double jxb_a;
  double kxb_a;
  double lxb_a;
  double mxb_a;
  double nxb_a;
  double oxb_a;
  double pxb_a;
  double qxb_a;
  double rxb_a;
  double sxb_a;
  double txb_a;
  double uxb_a;
  double vxb_a;
  double wxb_a;
  double xxb_a;
  double yxb_a;
  double ayb_a;
  double byb_a;
  double cyb_a;
  double dyb_a;
  double eyb_a;

  // p_x_p_Lgh_baru is 1-by-5
  b_a = pos_ob_x - s;
  c_a = ey - pos_ob_y;
  d_a = pos_ob_x - s;
  e_a = ey - pos_ob_y;
  f_a = pos_ob_x - s;
  g_a = ey - pos_ob_y;
  h_a = pos_ob_x - s;
  i_a = ey - pos_ob_y;
  j_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  k_a = pos_ob_x - s;
  l_a = ey - pos_ob_y;
  m_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  n_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  o_a = pos_ob_x - s;
  p_a = ey - pos_ob_y;
  q_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  r_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  s_a = pos_ob_x - s;
  t_a = ey - pos_ob_y;
  u_a = pos_ob_x - s;
  v_a = ey - pos_ob_y;
  w_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  x_a = pos_ob_x - s;
  y_a = ey - pos_ob_y;
  ab_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  bb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  cb_a = pos_ob_x - s;
  db_a = ey - pos_ob_y;
  eb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  fb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  gb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  hb_a = pos_ob_x - s;
  ib_a = ey - pos_ob_y;
  jb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  kb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  lb_a = pos_ob_x - s;
  mb_a = ey - pos_ob_y;
  nb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ob_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  pb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  qb_a = pos_ob_x - s;
  rb_a = ey - pos_ob_y;
  sb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  tb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  ub_a = pos_ob_x - s;
  vb_a = ey - pos_ob_y;
  wb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  xb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  yb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  ac_a = pos_ob_x - s;
  bc_a = ey - pos_ob_y;
  cc_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  dc_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  ec_a = pos_ob_x - s;
  fc_a = ey - pos_ob_y;
  gc_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  hc_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  ic_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  jc_a = pos_ob_x - s;
  kc_a = ey - pos_ob_y;
  lc_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  mc_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  nc_a = pos_ob_x - s;
  oc_a = ey - pos_ob_y;
  pc_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  qc_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  rc_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  sc_a = pos_ob_x - s;
  tc_a = ey - pos_ob_y;
  uc_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  vc_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  wc_a = pos_ob_x - s;
  xc_a = ey - pos_ob_y;
  yc_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ad_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  bd_a = pos_ob_x - s;
  cd_a = ey - pos_ob_y;
  dd_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ed_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  fd_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  gd_a = pos_ob_x - s;
  hd_a = ey - pos_ob_y;
  gd_a = gd_a * gd_a + hd_a * hd_a;
  hd_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  id_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  jd_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  kd_a = pos_ob_x - s;
  ld_a = ey - pos_ob_y;
  md_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  nd_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  od_a = pos_ob_x - s;
  pd_a = ey - pos_ob_y;
  qd_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  rd_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  sd_a = pos_ob_x - s;
  td_a = ey - pos_ob_y;
  ud_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  vd_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  wd_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  xd_a = pos_ob_x - s;
  yd_a = ey - pos_ob_y;
  xd_a = xd_a * xd_a + yd_a * yd_a;
  yd_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ae_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  be_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  ce_a = pos_ob_x - s;
  de_a = ey - pos_ob_y;
  ee_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  fe_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  ge_a = pos_ob_x - s;
  he_a = ey - pos_ob_y;
  ie_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  je_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  ke_a = pos_ob_x - s;
  le_a = ey - pos_ob_y;
  me_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ne_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  oe_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  pe_a = pos_ob_x - s;
  qe_a = ey - pos_ob_y;
  pe_a = pe_a * pe_a + qe_a * qe_a;
  qe_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  re_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  se_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  te_a = pos_ob_x - s;
  ue_a = ey - pos_ob_y;
  ve_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  we_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  xe_a = pos_ob_x - s;
  ye_a = ey - pos_ob_y;
  af_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  bf_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  cf_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  df_a = pos_ob_x - s;
  ef_a = ey - pos_ob_y;
  ff_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  gf_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  hf_a = pos_ob_x - s;
  if_a = ey - pos_ob_y;
  jf_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  kf_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  lf_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  mf_a = pos_ob_x - s;
  nf_a = ey - pos_ob_y;
  of_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  pf_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  qf_a = pos_ob_x - s;
  rf_a = ey - pos_ob_y;
  sf_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  tf_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  uf_a = pos_ob_x - s;
  vf_a = ey - pos_ob_y;
  wf_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  xf_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  yf_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  ag_a = pos_ob_x - s;
  bg_a = ey - pos_ob_y;
  ag_a = ag_a * ag_a + bg_a * bg_a;
  bg_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  cg_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  dg_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  eg_a = pos_ob_x - s;
  fg_a = ey - pos_ob_y;
  gg_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  hg_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  ig_a = pos_ob_x - s;
  jg_a = ey - pos_ob_y;
  kg_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  lg_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  mg_a = pos_ob_x - s;
  ng_a = ey - pos_ob_y;
  og_a = pos_ob_x - s;
  pg_a = ey - pos_ob_y;
  qg_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  rg_a = pos_ob_x - s;
  sg_a = ey - pos_ob_y;
  tg_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ug_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  vg_a = pos_ob_x - s;
  wg_a = ey - pos_ob_y;
  xg_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  yg_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  ah_a = pos_ob_x - s;
  bh_a = ey - pos_ob_y;
  ch_a = pos_ob_x - s;
  dh_a = ey - pos_ob_y;
  eh_a = pos_ob_x - s;
  fh_a = ey - pos_ob_y;
  gh_a = pos_ob_x - s;
  hh_a = ey - pos_ob_y;
  ih_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  jh_a = pos_ob_x - s;
  kh_a = ey - pos_ob_y;
  lh_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  mh_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  nh_a = pos_ob_x - s;
  oh_a = ey - pos_ob_y;
  ph_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  qh_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  rh_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  sh_a = pos_ob_x - s;
  th_a = ey - pos_ob_y;
  uh_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  vh_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  wh_a = pos_ob_x - s;
  xh_a = ey - pos_ob_y;
  yh_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ai_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  bi_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  ci_a = pos_ob_x - s;
  di_a = ey - pos_ob_y;
  ei_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  fi_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  gi_a = pos_ob_x - s;
  hi_a = ey - pos_ob_y;
  ii_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ji_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  ki_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  li_a = pos_ob_x - s;
  mi_a = ey - pos_ob_y;
  ni_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  oi_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  pi_a = pos_ob_x - s;
  qi_a = ey - pos_ob_y;
  ri_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  si_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  ti_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  ui_a = pos_ob_x - s;
  vi_a = ey - pos_ob_y;
  wi_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  xi_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  yi_a = pos_ob_x - s;
  aj_a = ey - pos_ob_y;
  bj_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  cj_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  dj_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  ej_a = pos_ob_x - s;
  fj_a = ey - pos_ob_y;
  ej_a = ej_a * ej_a + fj_a * fj_a;
  fj_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  gj_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  hj_a = pos_ob_x - s;
  ij_a = ey - pos_ob_y;
  jj_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  kj_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  lj_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  mj_a = pos_ob_x - s;
  nj_a = ey - pos_ob_y;
  oj_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  pj_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  qj_a = pos_ob_x - s;
  rj_a = ey - pos_ob_y;
  sj_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  tj_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  uj_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  vj_a = pos_ob_x - s;
  wj_a = ey - pos_ob_y;
  vj_a = vj_a * vj_a + wj_a * wj_a;
  wj_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  xj_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  yj_a = pos_ob_x - s;
  ak_a = ey - pos_ob_y;
  bk_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ck_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  dk_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  ek_a = pos_ob_x - s;
  fk_a = ey - pos_ob_y;
  gk_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  hk_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  ik_a = pos_ob_x - s;
  jk_a = ey - pos_ob_y;
  kk_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  lk_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  mk_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  nk_a = pos_ob_x - s;
  ok_a = ey - pos_ob_y;
  nk_a = nk_a * nk_a + ok_a * ok_a;
  ok_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  pk_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  qk_a = pos_ob_x - s;
  rk_a = ey - pos_ob_y;
  sk_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  tk_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  uk_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  vk_a = pos_ob_x - s;
  wk_a = ey - pos_ob_y;
  xk_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  yk_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  al_a = pos_ob_x - s;
  bl_a = ey - pos_ob_y;
  cl_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  dl_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  el_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  fl_a = pos_ob_x - s;
  gl_a = ey - pos_ob_y;
  hl_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  il_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  jl_a = pos_ob_x - s;
  kl_a = ey - pos_ob_y;
  ll_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ml_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  nl_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  ol_a = pos_ob_x - s;
  pl_a = ey - pos_ob_y;
  ql_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  rl_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  sl_a = pos_ob_x - s;
  tl_a = ey - pos_ob_y;
  ul_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  vl_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  wl_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  xl_a = pos_ob_x - s;
  yl_a = ey - pos_ob_y;
  am_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  bm_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  cm_a = pos_ob_x - s;
  dm_a = ey - pos_ob_y;
  em_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  fm_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  gm_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  hm_a = pos_ob_x - s;
  im_a = ey - pos_ob_y;
  hm_a = hm_a * hm_a + im_a * im_a;
  im_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  jm_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  km_a = pos_ob_x - s;
  lm_a = ey - pos_ob_y;
  mm_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  nm_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  om_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  pm_a = pos_ob_x - s;
  qm_a = ey - pos_ob_y;
  rm_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  sm_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  tm_a = pos_ob_x - s;
  um_a = ey - pos_ob_y;
  vm_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  wm_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  xm_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  ym_a = pos_ob_x - s;
  an_a = ey - pos_ob_y;
  bn_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  cn_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  dn_a = pos_ob_x - s;
  en_a = ey - pos_ob_y;
  fn_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  gn_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  hn_a = pos_ob_x - s;
  in_a = ey - pos_ob_y;
  jn_a = pos_ob_x - s;
  kn_a = ey - pos_ob_y;
  ln_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  mn_a = pos_ob_x - s;
  nn_a = ey - pos_ob_y;
  on_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  pn_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  qn_a = pos_ob_x - s;
  rn_a = ey - pos_ob_y;
  sn_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  tn_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  un_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  vn_a = pos_ob_x - s;
  wn_a = ey - pos_ob_y;
  xn_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  yn_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  ao_a = pos_ob_x - s;
  bo_a = ey - pos_ob_y;
  co_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  do_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  eo_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  fo_a = pos_ob_x - s;
  go_a = ey - pos_ob_y;
  ho_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  io_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  jo_a = pos_ob_x - s;
  ko_a = ey - pos_ob_y;
  lo_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  mo_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  no_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  oo_a = pos_ob_x - s;
  po_a = ey - pos_ob_y;
  qo_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ro_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  so_a = pos_ob_x - s;
  to_a = ey - pos_ob_y;
  uo_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  vo_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  wo_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  xo_a = pos_ob_x - s;
  yo_a = ey - pos_ob_y;
  ap_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  bp_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  cp_a = pos_ob_x - s;
  dp_a = ey - pos_ob_y;
  ep_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  fp_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  gp_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  hp_a = pos_ob_x - s;
  ip_a = ey - pos_ob_y;
  jp_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  kp_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  lp_a = pos_ob_x - s;
  mp_a = ey - pos_ob_y;
  np_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  op_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  pp_a = pos_ob_x - s;
  qp_a = ey - pos_ob_y;
  rp_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  sp_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  tp_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  up_a = pos_ob_x - s;
  vp_a = ey - pos_ob_y;
  wp_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  xp_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  wp_a = wp_a * wp_a + xp_a * xp_a;
  xp_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  yp_a = pos_ob_x - s;
  aq_a = ey - pos_ob_y;
  bq_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  cq_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  dq_a = pos_ob_x - s;
  eq_a = ey - pos_ob_y;
  fq_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  gq_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  hq_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  iq_a = pos_ob_x - s;
  jq_a = ey - pos_ob_y;
  kq_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  lq_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  mq_a = pos_ob_x - s;
  nq_a = ey - pos_ob_y;
  oq_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  pq_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  qq_a = pos_ob_x - s;
  rq_a = ey - pos_ob_y;
  sq_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  tq_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  uq_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  vq_a = pos_ob_x - s;
  wq_a = ey - pos_ob_y;
  xq_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  yq_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  xq_a = xq_a * xq_a + yq_a * yq_a;
  yq_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  ar_a = pos_ob_x - s;
  br_a = ey - pos_ob_y;
  cr_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  dr_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  er_a = pos_ob_x - s;
  fr_a = ey - pos_ob_y;
  gr_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  hr_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  ir_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  jr_a = pos_ob_x - s;
  kr_a = ey - pos_ob_y;
  lr_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  mr_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  nr_a = pos_ob_x - s;
  or_a = ey - pos_ob_y;
  pr_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  qr_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  rr_a = pos_ob_x - s;
  sr_a = ey - pos_ob_y;
  tr_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ur_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  vr_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  wr_a = pos_ob_x - s;
  xr_a = ey - pos_ob_y;
  yr_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  as_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  yr_a = yr_a * yr_a + as_a * as_a;
  as_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  bs_a = pos_ob_x - s;
  cs_a = ey - pos_ob_y;
  ds_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  es_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  fs_a = pos_ob_x - s;
  gs_a = ey - pos_ob_y;
  hs_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  is_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  js_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  ks_a = pos_ob_x - s;
  ls_a = ey - pos_ob_y;
  ms_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ns_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  os_a = pos_ob_x - s;
  ps_a = ey - pos_ob_y;
  qs_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  rs_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  ss_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  ts_a = pos_ob_x - s;
  us_a = ey - pos_ob_y;
  vs_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ws_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  xs_a = pos_ob_x - s;
  ys_a = ey - pos_ob_y;
  at_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  bt_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  ct_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  dt_a = pos_ob_x - s;
  et_a = ey - pos_ob_y;
  ft_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  gt_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  ht_a = pos_ob_x - s;
  it_a = ey - pos_ob_y;
  jt_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  kt_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  lt_a = pos_ob_x - s;
  mt_a = ey - pos_ob_y;
  nt_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ot_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  pt_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  qt_a = pos_ob_x - s;
  rt_a = ey - pos_ob_y;
  st_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  tt_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  st_a = st_a * st_a + tt_a * tt_a;
  tt_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  ut_a = pos_ob_x - s;
  vt_a = ey - pos_ob_y;
  wt_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  xt_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  yt_a = pos_ob_x - s;
  au_a = ey - pos_ob_y;
  bu_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  cu_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  du_a = pos_ob_x - s;
  eu_a = ey - pos_ob_y;
  fu_a = pos_ob_x - s;
  gu_a = ey - pos_ob_y;
  hu_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  iu_a = pos_ob_x - s;
  ju_a = ey - pos_ob_y;
  ku_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  lu_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  mu_a = pos_ob_x - s;
  nu_a = ey - pos_ob_y;
  ou_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  pu_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  qu_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  ru_a = pos_ob_x - s;
  su_a = ey - pos_ob_y;
  tu_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  uu_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  vu_a = pos_ob_x - s;
  wu_a = ey - pos_ob_y;
  xu_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  yu_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  av_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  bv_a = pos_ob_x - s;
  cv_a = ey - pos_ob_y;
  dv_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ev_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  fv_a = pos_ob_x - s;
  gv_a = ey - pos_ob_y;
  hv_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  iv_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  jv_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  kv_a = pos_ob_x - s;
  lv_a = ey - pos_ob_y;
  mv_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  nv_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  ov_a = pos_ob_x - s;
  pv_a = ey - pos_ob_y;
  qv_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  rv_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  sv_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  tv_a = pos_ob_x - s;
  uv_a = ey - pos_ob_y;
  vv_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  wv_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  xv_a = pos_ob_x - s;
  yv_a = ey - pos_ob_y;
  aw_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  bw_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  cw_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  dw_a = pos_ob_x - s;
  ew_a = ey - pos_ob_y;
  fw_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  gw_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  hw_a = pos_ob_x - s;
  iw_a = ey - pos_ob_y;
  jw_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  kw_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  lw_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  mw_a = pos_ob_x - s;
  nw_a = ey - pos_ob_y;
  ow_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  pw_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  qw_a = pos_ob_x - s;
  rw_a = ey - pos_ob_y;
  sw_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  tw_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  uw_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  vw_a = pos_ob_x - s;
  ww_a = ey - pos_ob_y;
  xw_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  yw_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  ax_a = pos_ob_x - s;
  bx_a = ey - pos_ob_y;
  cx_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  dx_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  ex_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  fx_a = pos_ob_x - s;
  gx_a = ey - pos_ob_y;
  hx_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ix_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  jx_a = pos_ob_x - s;
  kx_a = ey - pos_ob_y;
  lx_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  mx_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  nx_a = pos_ob_x - s;
  ox_a = ey - pos_ob_y;
  px_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  qx_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  rx_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  sx_a = pos_ob_x - s;
  tx_a = ey - pos_ob_y;
  ux_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  vx_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  ux_a = ux_a * ux_a + vx_a * vx_a;
  vx_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  wx_a = pos_ob_x - s;
  xx_a = ey - pos_ob_y;
  yx_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ay_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  by_a = pos_ob_x - s;
  cy_a = ey - pos_ob_y;
  dy_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ey_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  fy_a = pos_ob_x - s;
  gy_a = ey - pos_ob_y;
  hy_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  iy_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  jy_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  ky_a = pos_ob_x - s;
  ly_a = ey - pos_ob_y;
  my_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ny_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  my_a = my_a * my_a + ny_a * ny_a;
  ny_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  oy_a = pos_ob_x - s;
  py_a = ey - pos_ob_y;
  qy_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ry_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  sy_a = pos_ob_x - s;
  ty_a = ey - pos_ob_y;
  uy_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  vy_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  wy_a = pos_ob_x - s;
  xy_a = ey - pos_ob_y;
  yy_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  aab_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  bab_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  cab_a = pos_ob_x - s;
  dab_a = ey - pos_ob_y;
  eab_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  fab_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  eab_a = eab_a * eab_a + fab_a * fab_a;
  fab_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  gab_a = pos_ob_x - s;
  hab_a = ey - pos_ob_y;
  iab_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  jab_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  kab_a = pos_ob_x - s;
  lab_a = ey - pos_ob_y;
  mab_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  nab_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  oab_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  pab_a = pos_ob_x - s;
  qab_a = ey - pos_ob_y;
  rab_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  sab_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  tab_a = pos_ob_x - s;
  uab_a = ey - pos_ob_y;
  vab_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  wab_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  xab_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  yab_a = pos_ob_x - s;
  abb_a = ey - pos_ob_y;
  bbb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  cbb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  dbb_a = pos_ob_x - s;
  ebb_a = ey - pos_ob_y;
  fbb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  gbb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  hbb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  ibb_a = pos_ob_x - s;
  jbb_a = ey - pos_ob_y;
  kbb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  lbb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  mbb_a = pos_ob_x - s;
  nbb_a = ey - pos_ob_y;
  obb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  pbb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  qbb_a = pos_ob_x - s;
  rbb_a = ey - pos_ob_y;
  sbb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  tbb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  ubb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  vbb_a = pos_ob_x - s;
  wbb_a = ey - pos_ob_y;
  xbb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ybb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  xbb_a = xbb_a * xbb_a + ybb_a * ybb_a;
  ybb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  acb_a = pos_ob_x - s;
  bcb_a = ey - pos_ob_y;
  ccb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  dcb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  ecb_a = pos_ob_x - s;
  fcb_a = ey - pos_ob_y;
  gcb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  hcb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  icb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  jcb_a = pos_ob_x - s;
  kcb_a = ey - pos_ob_y;
  lcb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  mcb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  ncb_a = pos_ob_x - s;
  ocb_a = ey - pos_ob_y;
  pcb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  qcb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  rcb_a = pos_ob_x - s;
  scb_a = ey - pos_ob_y;
  tcb_a = pos_ob_x - s;
  ucb_a = ey - pos_ob_y;
  vcb_a = pos_ob_x - s;
  wcb_a = ey - pos_ob_y;
  xcb_a = pos_ob_x - s;
  ycb_a = ey - pos_ob_y;
  adb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  bdb_a = pos_ob_x - s;
  cdb_a = ey - pos_ob_y;
  ddb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  edb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  fdb_a = pos_ob_x - s;
  gdb_a = ey - pos_ob_y;
  hdb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  idb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  jdb_a = pos_ob_x - s;
  kdb_a = ey - pos_ob_y;
  ldb_a = pos_ob_x - s;
  mdb_a = ey - pos_ob_y;
  ndb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  odb_a = pos_ob_x - s;
  pdb_a = ey - pos_ob_y;
  qdb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  rdb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  sdb_a = pos_ob_x - s;
  tdb_a = ey - pos_ob_y;
  udb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  vdb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  wdb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  xdb_a = pos_ob_x - s;
  ydb_a = ey - pos_ob_y;
  xdb_a = xdb_a * xdb_a + ydb_a * ydb_a;
  ydb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  aeb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  beb_a = pos_ob_x - s;
  ceb_a = ey - pos_ob_y;
  deb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  eeb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  feb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  geb_a = pos_ob_x - s;
  heb_a = ey - pos_ob_y;
  ieb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  jeb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  keb_a = pos_ob_x - s;
  leb_a = ey - pos_ob_y;
  meb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  neb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  oeb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  peb_a = pos_ob_x - s;
  qeb_a = ey - pos_ob_y;
  reb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  seb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  teb_a = pos_ob_x - s;
  ueb_a = ey - pos_ob_y;
  veb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  web_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  xeb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  yeb_a = pos_ob_x - s;
  afb_a = ey - pos_ob_y;
  bfb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  cfb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  dfb_a = pos_ob_x - s;
  efb_a = ey - pos_ob_y;
  ffb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  gfb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  hfb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  ifb_a = pos_ob_x - s;
  jfb_a = ey - pos_ob_y;
  kfb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  lfb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  mfb_a = pos_ob_x - s;
  nfb_a = ey - pos_ob_y;
  ofb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  pfb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  qfb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  rfb_a = pos_ob_x - s;
  sfb_a = ey - pos_ob_y;
  tfb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ufb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  vfb_a = pos_ob_x - s;
  wfb_a = ey - pos_ob_y;
  xfb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  yfb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  agb_a = pos_ob_x - s;
  bgb_a = ey - pos_ob_y;
  cgb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  dgb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  egb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  fgb_a = pos_ob_x - s;
  ggb_a = ey - pos_ob_y;
  hgb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  igb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  hgb_a = hgb_a * hgb_a + igb_a * igb_a;
  igb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  jgb_a = pos_ob_x - s;
  kgb_a = ey - pos_ob_y;
  lgb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  mgb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  ngb_a = pos_ob_x - s;
  ogb_a = ey - pos_ob_y;
  pgb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  qgb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  rgb_a = pos_ob_x - s;
  sgb_a = ey - pos_ob_y;
  tgb_a = pos_ob_x - s;
  ugb_a = ey - pos_ob_y;
  vgb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  wgb_a = pos_ob_x - s;
  xgb_a = ey - pos_ob_y;
  ygb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ahb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  bhb_a = pos_ob_x - s;
  chb_a = ey - pos_ob_y;
  dhb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ehb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  fhb_a = pos_ob_x - s;
  ghb_a = ey - pos_ob_y;
  hhb_a = pos_ob_x - s;
  ihb_a = ey - pos_ob_y;
  jhb_a = pos_ob_x - s;
  khb_a = ey - pos_ob_y;
  lhb_a = pos_ob_x - s;
  mhb_a = ey - pos_ob_y;
  nhb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  ohb_a = pos_ob_x - s;
  phb_a = ey - pos_ob_y;
  qhb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  rhb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  shb_a = pos_ob_x - s;
  thb_a = ey - pos_ob_y;
  uhb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  vhb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  whb_a = pos_ob_x - s;
  xhb_a = ey - pos_ob_y;
  yhb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  aib_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  bib_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  cib_a = pos_ob_x - s;
  dib_a = ey - pos_ob_y;
  cib_a = cib_a * cib_a + dib_a * dib_a;
  dib_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  eib_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  fib_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  gib_a = pos_ob_x - s;
  hib_a = ey - pos_ob_y;
  iib_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  jib_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  kib_a = pos_ob_x - s;
  lib_a = ey - pos_ob_y;
  mib_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  nib_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  oib_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  pib_a = pos_ob_x - s;
  qib_a = ey - pos_ob_y;
  rib_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  sib_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  tib_a = pos_ob_x - s;
  uib_a = ey - pos_ob_y;
  vib_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  wib_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  xib_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  yib_a = pos_ob_x - s;
  ajb_a = ey - pos_ob_y;
  bjb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  cjb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  djb_a = pos_ob_x - s;
  ejb_a = ey - pos_ob_y;
  fjb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  gjb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  hjb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  ijb_a = pos_ob_x - s;
  jjb_a = ey - pos_ob_y;
  kjb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ljb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  mjb_a = pos_ob_x - s;
  njb_a = ey - pos_ob_y;
  ojb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  pjb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  qjb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  rjb_a = pos_ob_x - s;
  sjb_a = ey - pos_ob_y;
  tjb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ujb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  vjb_a = pos_ob_x - s;
  wjb_a = ey - pos_ob_y;
  xjb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  yjb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  akb_a = pos_ob_x - s;
  bkb_a = ey - pos_ob_y;
  ckb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  dkb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  ekb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  fkb_a = pos_ob_x - s;
  gkb_a = ey - pos_ob_y;
  hkb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ikb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  hkb_a = hkb_a * hkb_a + ikb_a * ikb_a;
  ikb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  jkb_a = pos_ob_x - s;
  kkb_a = ey - pos_ob_y;
  lkb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  mkb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  nkb_a = pos_ob_x - s;
  okb_a = ey - pos_ob_y;
  pkb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  qkb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  rkb_a = pos_ob_x - s;
  skb_a = ey - pos_ob_y;
  tkb_a = pos_ob_x - s;
  ukb_a = ey - pos_ob_y;
  vkb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  wkb_a = pos_ob_x - s;
  xkb_a = ey - pos_ob_y;
  ykb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  alb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  blb_a = pos_ob_x - s;
  clb_a = ey - pos_ob_y;
  dlb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  elb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  flb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  glb_a = pos_ob_x - s;
  hlb_a = ey - pos_ob_y;
  ilb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  jlb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  klb_a = pos_ob_x - s;
  llb_a = ey - pos_ob_y;
  mlb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  nlb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  olb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  plb_a = pos_ob_x - s;
  qlb_a = ey - pos_ob_y;
  rlb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  slb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  tlb_a = pos_ob_x - s;
  ulb_a = ey - pos_ob_y;
  vlb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  wlb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  xlb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  ylb_a = pos_ob_x - s;
  amb_a = ey - pos_ob_y;
  bmb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  cmb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  dmb_a = pos_ob_x - s;
  emb_a = ey - pos_ob_y;
  fmb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  gmb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  hmb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  imb_a = pos_ob_x - s;
  jmb_a = ey - pos_ob_y;
  kmb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  lmb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  mmb_a = pos_ob_x - s;
  nmb_a = ey - pos_ob_y;
  omb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  pmb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  qmb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  rmb_a = pos_ob_x - s;
  smb_a = ey - pos_ob_y;
  tmb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  umb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  vmb_a = pos_ob_x - s;
  wmb_a = ey - pos_ob_y;
  xmb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ymb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  anb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  bnb_a = pos_ob_x - s;
  cnb_a = ey - pos_ob_y;
  dnb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  enb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  fnb_a = pos_ob_x - s;
  gnb_a = ey - pos_ob_y;
  hnb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  inb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  jnb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  knb_a = pos_ob_x - s;
  lnb_a = ey - pos_ob_y;
  mnb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  nnb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  onb_a = pos_ob_x - s;
  pnb_a = ey - pos_ob_y;
  qnb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  rnb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  snb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  tnb_a = pos_ob_x - s;
  unb_a = ey - pos_ob_y;
  vnb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  wnb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  xnb_a = pos_ob_x - s;
  ynb_a = ey - pos_ob_y;
  aob_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  bob_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  cob_a = pos_ob_x - s;
  dob_a = ey - pos_ob_y;
  eob_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  fob_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  gob_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  hob_a = pos_ob_x - s;
  iob_a = ey - pos_ob_y;
  job_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  kob_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  job_a = job_a * job_a + kob_a * kob_a;
  kob_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  lob_a = pos_ob_x - s;
  mob_a = ey - pos_ob_y;
  nob_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  oob_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  pob_a = pos_ob_x - s;
  qob_a = ey - pos_ob_y;
  rob_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  sob_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  tob_a = pos_ob_x - s;
  uob_a = ey - pos_ob_y;
  vob_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  wob_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  xob_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  yob_a = pos_ob_x - s;
  apb_a = ey - pos_ob_y;
  bpb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  cpb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  bpb_a = bpb_a * bpb_a + cpb_a * cpb_a;
  cpb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  dpb_a = pos_ob_x - s;
  epb_a = ey - pos_ob_y;
  fpb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  gpb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  hpb_a = pos_ob_x - s;
  ipb_a = ey - pos_ob_y;
  jpb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  kpb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  lpb_a = pos_ob_x - s;
  mpb_a = ey - pos_ob_y;
  npb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  opb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  ppb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  qpb_a = pos_ob_x - s;
  rpb_a = ey - pos_ob_y;
  spb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  tpb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  spb_a = spb_a * spb_a + tpb_a * tpb_a;
  tpb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  upb_a = pos_ob_x - s;
  vpb_a = ey - pos_ob_y;
  wpb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  xpb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  ypb_a = pos_ob_x - s;
  aqb_a = ey - pos_ob_y;
  bqb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  cqb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  dqb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  eqb_a = pos_ob_x - s;
  fqb_a = ey - pos_ob_y;
  gqb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  hqb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  iqb_a = pos_ob_x - s;
  jqb_a = ey - pos_ob_y;
  kqb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  lqb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  mqb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  nqb_a = pos_ob_x - s;
  oqb_a = ey - pos_ob_y;
  pqb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  qqb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  rqb_a = pos_ob_x - s;
  sqb_a = ey - pos_ob_y;
  tqb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  uqb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  vqb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  wqb_a = pos_ob_x - s;
  xqb_a = ey - pos_ob_y;
  yqb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  arb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  brb_a = pos_ob_x - s;
  crb_a = ey - pos_ob_y;
  drb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  erb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  frb_a = pos_ob_x - s;
  grb_a = ey - pos_ob_y;
  hrb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  irb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  jrb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  krb_a = pos_ob_x - s;
  lrb_a = ey - pos_ob_y;
  mrb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  nrb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  mrb_a = mrb_a * mrb_a + nrb_a * nrb_a;
  nrb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  orb_a = pos_ob_x - s;
  prb_a = ey - pos_ob_y;
  qrb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  rrb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  srb_a = pos_ob_x - s;
  trb_a = ey - pos_ob_y;
  urb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  vrb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  wrb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  xrb_a = pos_ob_x - s;
  yrb_a = ey - pos_ob_y;
  asb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  bsb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  csb_a = pos_ob_x - s;
  dsb_a = ey - pos_ob_y;
  esb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  fsb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  gsb_a = pos_ob_x - s;
  hsb_a = ey - pos_ob_y;
  isb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  jsb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  ksb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  lsb_a = pos_ob_x - s;
  msb_a = ey - pos_ob_y;
  nsb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  osb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  nsb_a = nsb_a * nsb_a + osb_a * osb_a;
  osb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  psb_a = pos_ob_x - s;
  qsb_a = ey - pos_ob_y;
  rsb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ssb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  tsb_a = pos_ob_x - s;
  usb_a = ey - pos_ob_y;
  vsb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  wsb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  xsb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  ysb_a = pos_ob_x - s;
  atb_a = ey - pos_ob_y;
  btb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ctb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  dtb_a = pos_ob_x - s;
  etb_a = ey - pos_ob_y;
  ftb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  gtb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  htb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  itb_a = pos_ob_x - s;
  jtb_a = ey - pos_ob_y;
  ktb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ltb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  mtb_a = pos_ob_x - s;
  ntb_a = ey - pos_ob_y;
  otb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ptb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  qtb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  rtb_a = pos_ob_x - s;
  stb_a = ey - pos_ob_y;
  ttb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  utb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  vtb_a = pos_ob_x - s;
  wtb_a = ey - pos_ob_y;
  xtb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ytb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  aub_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  bub_a = pos_ob_x - s;
  cub_a = ey - pos_ob_y;
  dub_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  eub_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  fub_a = pos_ob_x - s;
  gub_a = ey - pos_ob_y;
  hub_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  iub_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  jub_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  kub_a = pos_ob_x - s;
  lub_a = ey - pos_ob_y;
  mub_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  nub_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  oub_a = pos_ob_x - s;
  pub_a = ey - pos_ob_y;
  qub_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  rub_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  sub_a = pos_ob_x - s;
  tub_a = ey - pos_ob_y;
  uub_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  vub_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  wub_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  xub_a = pos_ob_x - s;
  yub_a = ey - pos_ob_y;
  avb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  bvb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  avb_a = avb_a * avb_a + bvb_a * bvb_a;
  bvb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  cvb_a = pos_ob_x - s;
  dvb_a = ey - pos_ob_y;
  evb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  fvb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  gvb_a = pos_ob_x - s;
  hvb_a = ey - pos_ob_y;
  ivb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  jvb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  kvb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  lvb_a = pos_ob_x - s;
  mvb_a = ey - pos_ob_y;
  nvb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ovb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  pvb_a = pos_ob_x - s;
  qvb_a = ey - pos_ob_y;
  rvb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  svb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  tvb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  uvb_a = pos_ob_x - s;
  vvb_a = ey - pos_ob_y;
  uvb_a = uvb_a * uvb_a + vvb_a * vvb_a;
  vvb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  wvb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  xvb_a = pos_ob_x - s;
  yvb_a = ey - pos_ob_y;
  awb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  bwb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  cwb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  dwb_a = pos_ob_x - s;
  ewb_a = ey - pos_ob_y;
  fwb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  gwb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  hwb_a = pos_ob_x - s;
  iwb_a = ey - pos_ob_y;
  jwb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  kwb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  lwb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  mwb_a = pos_ob_x - s;
  nwb_a = ey - pos_ob_y;
  owb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  pwb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  qwb_a = pos_ob_x - s;
  rwb_a = ey - pos_ob_y;
  swb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  twb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  uwb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  vwb_a = pos_ob_x - s;
  wwb_a = ey - pos_ob_y;
  xwb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ywb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  axb_a = pos_ob_x - s;
  bxb_a = ey - pos_ob_y;
  cxb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  dxb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  exb_a = pos_ob_x - s;
  fxb_a = ey - pos_ob_y;
  gxb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  hxb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  ixb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  jxb_a = pos_ob_x - s;
  kxb_a = ey - pos_ob_y;
  jxb_a = jxb_a * jxb_a + kxb_a * kxb_a;
  kxb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  lxb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  mxb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  nxb_a = pos_ob_x - s;
  oxb_a = ey - pos_ob_y;
  pxb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  qxb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  rxb_a = pos_ob_x - s;
  sxb_a = ey - pos_ob_y;
  txb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  uxb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  vxb_a = (pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi));
  wxb_a = pos_ob_x - s;
  xxb_a = ey - pos_ob_y;
  yxb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  ayb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  byb_a = pos_ob_x - s;
  cyb_a = ey - pos_ob_y;
  dyb_a = (vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi);
  eyb_a = (yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi);
  out[0] = ((((xp_dot * std::cos(epsi) - yp_dot * std::sin(epsi)) *
              (((((((((((((((pow(Ds, 3.0) * (2.0 * pos_ob_x - 2.0 * s) *
    (((((((ey * yp_dot * std::cos(epsi) - pos_ob_x * xp_dot * std::cos(epsi)) -
    pos_ob_y * yp_dot * std::cos(epsi)) + s * xp_dot * std::cos(epsi)) + ey *
    xp_dot * std::sin(epsi)) - pos_ob_y * xp_dot * std::sin(epsi)) + pos_ob_x *
    yp_dot * std::sin(epsi)) - s * yp_dot * std::sin(epsi)) / (2.0 * pow
    (1.0 - Ds * Ds / (b_a * b_a + c_a * c_a), 1.5) * pow(d_a * d_a + e_a
    * e_a, 3.5)) + Ds * (ey - pos_ob_y) * (((((((ey * xp_dot * std::cos(epsi) -
    pos_ob_y * xp_dot * std::cos(epsi)) + pos_ob_x * yp_dot * std::cos(epsi)) -
    s * yp_dot * std::cos(epsi)) - ey * yp_dot * std::sin(epsi)) + pos_ob_x *
    xp_dot * std::sin(epsi)) + pos_ob_y * yp_dot * std::sin(epsi)) - s * xp_dot *
    std::sin(epsi)) / (std::sqrt(1.0 - Ds * Ds / (f_a * f_a + g_a * g_a)) *
                       pow(h_a * h_a + i_a * i_a, 2.5))) + (yp_dot * std::
    cos(epsi) + xp_dot * std::sin(epsi)) * (((((((((((ey * vel_ob_x + pos_ob_x *
    vel_ob_y) - pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos
    (epsi)) + pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos
    (epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) -
    pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) +
    s * xp_dot * std::sin(epsi)) / (std::sqrt(1.0 - j_a * j_a / ((k_a * k_a +
    l_a * l_a) * (m_a * m_a + n_a * n_a))) * pow(o_a * o_a + p_a * p_a,
    1.5) * std::sqrt(q_a * q_a + r_a * r_a))) + Ds * (2.0 * pos_ob_x - 2.0 * s) *
    (((((((ey * yp_dot * std::cos(epsi) - pos_ob_x * xp_dot * std::cos(epsi)) -
    pos_ob_y * yp_dot * std::cos(epsi)) + s * xp_dot * std::cos(epsi)) + ey *
    xp_dot * std::sin(epsi)) - pos_ob_y * xp_dot * std::sin(epsi)) + pos_ob_x *
      yp_dot * std::sin(epsi)) - s * yp_dot * std::sin(epsi)) / (std::sqrt(1.0 -
    Ds * Ds / (s_a * s_a + t_a * t_a)) * pow(u_a * u_a + v_a * v_a, 2.5)))
    - (psi_dot - psi_dot_com) * ((yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot *
    std::sin(epsi)) * (((((xp_dot * xp_dot + yp_dot * yp_dot) - vel_ob_x *
    xp_dot * std::cos(epsi)) - vel_ob_y * yp_dot * std::cos(epsi)) - vel_ob_y *
                        xp_dot * std::sin(epsi)) + vel_ob_x * yp_dot * std::sin
                       (epsi)) / (std::sqrt(1.0 - w_a * w_a / ((x_a * x_a + y_a *
    y_a) * (ab_a * ab_a + bb_a * bb_a))) * std::sqrt(cb_a * cb_a + db_a * db_a) *
    pow(eb_a * eb_a + fb_a * fb_a, 1.5))) - (ey - pos_ob_y) * (xp_dot *
    std::cos(epsi) - yp_dot * std::sin(epsi)) * ((yp_dot * std::cos(epsi) -
    vel_ob_y) + xp_dot * std::sin(epsi)) / (std::sqrt(1.0 - gb_a * gb_a / ((hb_a
    * hb_a + ib_a * ib_a) * (jb_a * jb_a + kb_a * kb_a))) * pow(lb_a *
    lb_a + mb_a * mb_a, 1.5) * std::sqrt(nb_a * nb_a + ob_a * ob_a))) -
                        (pos_ob_x - s) * (yp_dot * std::cos(epsi) + xp_dot * std::
    sin(epsi)) * ((yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi))
                        / (std::sqrt(1.0 - pb_a * pb_a / ((qb_a * qb_a + rb_a *
    rb_a) * (sb_a * sb_a + tb_a * tb_a))) * pow(ub_a * ub_a + vb_a *
    vb_a, 1.5) * std::sqrt(wb_a * wb_a + xb_a * xb_a))) - (2.0 * pos_ob_x - 2.0 *
    s) * (psi_dot - psi_dot_com) * (((((xp_dot * xp_dot + yp_dot * yp_dot) -
    vel_ob_x * xp_dot * std::cos(epsi)) - vel_ob_y * yp_dot * std::cos(epsi)) -
    vel_ob_y * xp_dot * std::sin(epsi)) + vel_ob_x * yp_dot * std::sin(epsi)) *
                       (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) -
    pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) +
    pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) +
    s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x *
    xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot *
                        std::sin(epsi)) / (2.0 * std::sqrt(1.0 - yb_a * yb_a /
    ((ac_a * ac_a + bc_a * bc_a) * (cc_a * cc_a + dc_a * dc_a))) * pow
    (ec_a * ec_a + fc_a * fc_a, 1.5) * pow(gc_a * gc_a + hc_a * hc_a,
    1.5))) - 3.0 * (2.0 * pos_ob_x - 2.0 * s) * (ey - pos_ob_y) * (xp_dot * std::
    cos(epsi) - yp_dot * std::sin(epsi)) * (((((((((((ey * vel_ob_x + pos_ob_x *
    vel_ob_y) - pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos
    (epsi)) + pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos
    (epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) -
    pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) +
    s * xp_dot * std::sin(epsi)) / (2.0 * std::sqrt(1.0 - ic_a * ic_a / ((jc_a *
    jc_a + kc_a * kc_a) * (lc_a * lc_a + mc_a * mc_a))) * pow(nc_a *
    nc_a + oc_a * oc_a, 2.5) * std::sqrt(pc_a * pc_a + qc_a * qc_a))) - 3.0 *
                     (2.0 * pos_ob_x - 2.0 * s) * (pos_ob_x - s) * (yp_dot * std::
    cos(epsi) + xp_dot * std::sin(epsi)) * (((((((((((ey * vel_ob_x + pos_ob_x *
    vel_ob_y) - pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos
    (epsi)) + pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos
    (epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) -
    pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) +
    s * xp_dot * std::sin(epsi)) / (2.0 * std::sqrt(1.0 - rc_a * rc_a / ((sc_a *
    sc_a + tc_a * tc_a) * (uc_a * uc_a + vc_a * vc_a))) * pow(wc_a *
    wc_a + xc_a * xc_a, 2.5) * std::sqrt(yc_a * yc_a + ad_a * ad_a))) + (2.0 *
    ((pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin
                       (epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) -
    vel_ob_y) + xp_dot * std::sin(epsi))) * ((vel_ob_x - xp_dot * std::cos(epsi))
    + yp_dot * std::sin(epsi)) / ((bd_a * bd_a + cd_a * cd_a) * (dd_a * dd_a +
    ed_a * ed_a)) - (2.0 * pos_ob_x - 2.0 * s) * (fd_a * fd_a) / (gd_a * gd_a *
    (hd_a * hd_a + id_a * id_a))) * (psi_dot - psi_dot_com) * (((((xp_dot *
    xp_dot + yp_dot * yp_dot) - vel_ob_x * xp_dot * std::cos(epsi)) - vel_ob_y *
    yp_dot * std::cos(epsi)) - vel_ob_y * xp_dot * std::sin(epsi)) + vel_ob_x *
    yp_dot * std::sin(epsi)) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) -
    pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) +
    pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) +
    s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x *
    xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot *
    std::sin(epsi)) / (2.0 * pow(1.0 - jd_a * jd_a / ((kd_a * kd_a +
    ld_a * ld_a) * (md_a * md_a + nd_a * nd_a)), 1.5) * std::sqrt(od_a * od_a +
    pd_a * pd_a) * pow(qd_a * qd_a + rd_a * rd_a, 1.5))) + (2.0 *
    ((pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin
                       (epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) -
    vel_ob_y) + xp_dot * std::sin(epsi))) * ((vel_ob_x - xp_dot * std::cos(epsi))
    + yp_dot * std::sin(epsi)) / ((sd_a * sd_a + td_a * td_a) * (ud_a * ud_a +
    vd_a * vd_a)) - (2.0 * pos_ob_x - 2.0 * s) * (wd_a * wd_a) / (xd_a * xd_a *
    (yd_a * yd_a + ae_a * ae_a))) * (ey - pos_ob_y) * (xp_dot * std::cos(epsi) -
    yp_dot * std::sin(epsi)) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) -
    pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) +
    pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) +
    s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x *
    xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot *
    std::sin(epsi)) / (2.0 * pow(1.0 - be_a * be_a / ((ce_a * ce_a +
    de_a * de_a) * (ee_a * ee_a + fe_a * fe_a)), 1.5) * pow(ge_a * ge_a
    + he_a * he_a, 1.5) * std::sqrt(ie_a * ie_a + je_a * je_a))) + (2.0 *
    ((pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin
                       (epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) -
    vel_ob_y) + xp_dot * std::sin(epsi))) * ((vel_ob_x - xp_dot * std::cos(epsi))
    + yp_dot * std::sin(epsi)) / ((ke_a * ke_a + le_a * le_a) * (me_a * me_a +
    ne_a * ne_a)) - (2.0 * pos_ob_x - 2.0 * s) * (oe_a * oe_a) / (pe_a * pe_a *
    (qe_a * qe_a + re_a * re_a))) * (pos_ob_x - s) * (yp_dot * std::cos(epsi) +
    xp_dot * std::sin(epsi)) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) -
    pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) +
    pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) +
    s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x *
    xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot *
    std::sin(epsi)) / (2.0 * pow(1.0 - se_a * se_a / ((te_a * te_a +
    ue_a * ue_a) * (ve_a * ve_a + we_a * we_a)), 1.5) * pow(xe_a * xe_a
    + ye_a * ye_a, 1.5) * std::sqrt(af_a * af_a + bf_a * bf_a))) - ((vel_ob_x *
    std::cos(epsi) - xp_dot) + vel_ob_y * std::sin(epsi)) * ((yp_dot * std::cos
    (epsi) - vel_ob_y) + xp_dot * std::sin(epsi)) * ((((2.0 * cf * yp_dot + 2.0 *
    cr * yp_dot) + m * psi_dot * (xp_dot * xp_dot)) + 2.0 * a * cf * psi_dot) -
    2.0 * b * cr * psi_dot) / (m * xp_dot * std::sqrt(1.0 - cf_a * cf_a / ((df_a
    * df_a + ef_a * ef_a) * (ff_a * ff_a + gf_a * gf_a))) * std::sqrt(hf_a *
    hf_a + if_a * if_a) * pow(jf_a * jf_a + kf_a * kf_a, 1.5))) - (2.0 *
    pos_ob_x - 2.0 * s) * ((vel_ob_x * std::cos(epsi) - xp_dot) + vel_ob_y * std::
    sin(epsi)) * ((((2.0 * cf * yp_dot + 2.0 * cr * yp_dot) + m * psi_dot *
                    (xp_dot * xp_dot)) + 2.0 * a * cf * psi_dot) - 2.0 * b * cr *
                  psi_dot) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) -
    pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) +
    pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) +
    s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x *
    xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot *
    std::sin(epsi)) / (2.0 * m * xp_dot * std::sqrt(1.0 - lf_a * lf_a / ((mf_a *
    mf_a + nf_a * nf_a) * (of_a * of_a + pf_a * pf_a))) * pow(qf_a *
    qf_a + rf_a * rf_a, 1.5) * pow(sf_a * sf_a + tf_a * tf_a, 1.5))) +
               (2.0 * ((pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) +
    yp_dot * std::sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) -
    vel_ob_y) + xp_dot * std::sin(epsi))) * ((vel_ob_x - xp_dot * std::cos(epsi))
    + yp_dot * std::sin(epsi)) / ((uf_a * uf_a + vf_a * vf_a) * (wf_a * wf_a +
    xf_a * xf_a)) - (2.0 * pos_ob_x - 2.0 * s) * (yf_a * yf_a) / (ag_a * ag_a *
    (bg_a * bg_a + cg_a * cg_a))) * ((vel_ob_x * std::cos(epsi) - xp_dot) +
    vel_ob_y * std::sin(epsi)) * ((((2.0 * cf * yp_dot + 2.0 * cr * yp_dot) + m *
    psi_dot * (xp_dot * xp_dot)) + 2.0 * a * cf * psi_dot) - 2.0 * b * cr *
    psi_dot) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y *
    vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y *
                      xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos
                     (epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::
                   sin(epsi)) - pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y *
                 yp_dot * std::sin(epsi)) + s * xp_dot * std::sin(epsi)) / (2.0 *
    m * xp_dot * pow(1.0 - dg_a * dg_a / ((eg_a * eg_a + fg_a * fg_a) *
    (gg_a * gg_a + hg_a * hg_a)), 1.5) * std::sqrt(ig_a * ig_a + jg_a * jg_a) *
    pow(kg_a * kg_a + lg_a * lg_a, 1.5))) + (yp_dot * std::cos(epsi) +
    xp_dot * std::sin(epsi)) * (((((((((((((((Ds * (pos_ob_x - s) * (((((((ey *
    xp_dot * std::cos(epsi) - pos_ob_y * xp_dot * std::cos(epsi)) + pos_ob_x *
    yp_dot * std::cos(epsi)) - s * yp_dot * std::cos(epsi)) - ey * yp_dot * std::
    sin(epsi)) + pos_ob_x * xp_dot * std::sin(epsi)) + pos_ob_y * yp_dot * std::
    sin(epsi)) - s * xp_dot * std::sin(epsi)) / (std::sqrt(1.0 - Ds * Ds / (mg_a
    * mg_a + ng_a * ng_a)) * pow(og_a * og_a + pg_a * pg_a, 2.5)) -
    (xp_dot * std::cos(epsi) - yp_dot * std::sin(epsi)) * (((((((((((ey *
    vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey *
    xp_dot * std::cos(epsi)) + pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x *
    yp_dot * std::cos(epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::
    sin(epsi)) - pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::
    sin(epsi)) + s * xp_dot * std::sin(epsi)) / (std::sqrt(1.0 - qg_a * qg_a /
    ((rg_a * rg_a + sg_a * sg_a) * (tg_a * tg_a + ug_a * ug_a))) * pow
    (vg_a * vg_a + wg_a * wg_a, 1.5) * std::sqrt(xg_a * xg_a + yg_a * yg_a))) -
    pow(Ds, 3.0) * (2.0 * ey - 2.0 * pos_ob_y) * (((((((ey * yp_dot *
    std::cos(epsi) - pos_ob_x * xp_dot * std::cos(epsi)) - pos_ob_y * yp_dot *
    std::cos(epsi)) + s * xp_dot * std::cos(epsi)) + ey * xp_dot * std::sin(epsi))
    - pos_ob_y * xp_dot * std::sin(epsi)) + pos_ob_x * yp_dot * std::sin(epsi))
    - s * yp_dot * std::sin(epsi)) / (2.0 * pow(1.0 - Ds * Ds / (ah_a *
    ah_a + bh_a * bh_a), 1.5) * pow(ch_a * ch_a + dh_a * dh_a, 3.5))) -
    Ds * (2.0 * ey - 2.0 * pos_ob_y) * (((((((ey * yp_dot * std::cos(epsi) -
    pos_ob_x * xp_dot * std::cos(epsi)) - pos_ob_y * yp_dot * std::cos(epsi)) +
    s * xp_dot * std::cos(epsi)) + ey * xp_dot * std::sin(epsi)) - pos_ob_y *
    xp_dot * std::sin(epsi)) + pos_ob_x * yp_dot * std::sin(epsi)) - s * yp_dot *
    std::sin(epsi)) / (std::sqrt(1.0 - Ds * Ds / (eh_a * eh_a + fh_a * fh_a)) *
                       pow(gh_a * gh_a + hh_a * hh_a, 2.5))) - (ey -
    pos_ob_y) * (xp_dot * std::cos(epsi) - yp_dot * std::sin(epsi)) * ((vel_ob_x
    - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi)) / (std::sqrt(1.0 -
    ih_a * ih_a / ((jh_a * jh_a + kh_a * kh_a) * (lh_a * lh_a + mh_a * mh_a))) *
    pow(nh_a * nh_a + oh_a * oh_a, 1.5) * std::sqrt(ph_a * ph_a + qh_a *
    qh_a))) - (pos_ob_x - s) * (yp_dot * std::cos(epsi) + xp_dot * std::sin(epsi))
    * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi)) / (std::
    sqrt(1.0 - rh_a * rh_a / ((sh_a * sh_a + th_a * th_a) * (uh_a * uh_a + vh_a *
    vh_a))) * pow(wh_a * wh_a + xh_a * xh_a, 1.5) * std::sqrt(yh_a *
    yh_a + ai_a * ai_a))) - (psi_dot - psi_dot_com) * ((vel_ob_x - xp_dot * std::
    cos(epsi)) + yp_dot * std::sin(epsi)) * (((((xp_dot * xp_dot + yp_dot *
    yp_dot) - vel_ob_x * xp_dot * std::cos(epsi)) - vel_ob_y * yp_dot * std::cos
    (epsi)) - vel_ob_y * xp_dot * std::sin(epsi)) + vel_ob_x * yp_dot * std::sin
    (epsi)) / (std::sqrt(1.0 - bi_a * bi_a / ((ci_a * ci_a + di_a * di_a) *
    (ei_a * ei_a + fi_a * fi_a))) * std::sqrt(gi_a * gi_a + hi_a * hi_a) *
               pow(ii_a * ii_a + ji_a * ji_a, 1.5))) + 3.0 * (ey -
    pos_ob_y) * (xp_dot * std::cos(epsi) - yp_dot * std::sin(epsi)) * (2.0 * ey
    - 2.0 * pos_ob_y) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) -
    pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) +
    pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) +
    s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x *
    xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot *
    std::sin(epsi)) / (2.0 * std::sqrt(1.0 - ki_a * ki_a / ((li_a * li_a + mi_a *
    mi_a) * (ni_a * ni_a + oi_a * oi_a))) * pow(pi_a * pi_a + qi_a *
    qi_a, 2.5) * std::sqrt(ri_a * ri_a + si_a * si_a))) + 3.0 * (pos_ob_x - s) *
    (yp_dot * std::cos(epsi) + xp_dot * std::sin(epsi)) * (2.0 * ey - 2.0 *
    pos_ob_y) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y *
    vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y *
                       xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos
                      (epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::
                    sin(epsi)) - pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y *
                  yp_dot * std::sin(epsi)) + s * xp_dot * std::sin(epsi)) / (2.0
    * std::sqrt(1.0 - ti_a * ti_a / ((ui_a * ui_a + vi_a * vi_a) * (wi_a * wi_a
    + xi_a * xi_a))) * pow(yi_a * yi_a + aj_a * aj_a, 2.5) * std::sqrt
    (bj_a * bj_a + cj_a * cj_a))) + (psi_dot - psi_dot_com) * (dj_a * dj_a *
    (2.0 * ey - 2.0 * pos_ob_y) / (ej_a * ej_a * (fj_a * fj_a + gj_a * gj_a)) -
    2.0 * ((pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::
    sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi))) * ((yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot *
    std::sin(epsi)) / ((hj_a * hj_a + ij_a * ij_a) * (jj_a * jj_a + kj_a * kj_a)))
    * (((((xp_dot * xp_dot + yp_dot * yp_dot) - vel_ob_x * xp_dot * std::cos
          (epsi)) - vel_ob_y * yp_dot * std::cos(epsi)) - vel_ob_y * xp_dot *
        std::sin(epsi)) + vel_ob_x * yp_dot * std::sin(epsi)) * (((((((((((ey *
    vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey *
    xp_dot * std::cos(epsi)) + pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x *
    yp_dot * std::cos(epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::
    sin(epsi)) - pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::
    sin(epsi)) + s * xp_dot * std::sin(epsi)) / (2.0 * pow(1.0 - lj_a *
    lj_a / ((mj_a * mj_a + nj_a * nj_a) * (oj_a * oj_a + pj_a * pj_a)), 1.5) *
    std::sqrt(qj_a * qj_a + rj_a * rj_a) * pow(sj_a * sj_a + tj_a * tj_a,
    1.5))) + (ey - pos_ob_y) * (xp_dot * std::cos(epsi) - yp_dot * std::sin(epsi))
    * (uj_a * uj_a * (2.0 * ey - 2.0 * pos_ob_y) / (vj_a * vj_a * (wj_a * wj_a +
    xj_a * xj_a)) - 2.0 * ((pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi))
    + yp_dot * std::sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) -
    vel_ob_y) + xp_dot * std::sin(epsi))) * ((yp_dot * std::cos(epsi) - vel_ob_y)
    + xp_dot * std::sin(epsi)) / ((yj_a * yj_a + ak_a * ak_a) * (bk_a * bk_a +
    ck_a * ck_a))) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y *
    vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y *
    xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) + s * yp_dot *
    std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x * xp_dot * std::
                        sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s *
                      xp_dot * std::sin(epsi)) / (2.0 * pow(1.0 - dk_a *
    dk_a / ((ek_a * ek_a + fk_a * fk_a) * (gk_a * gk_a + hk_a * hk_a)), 1.5) *
    pow(ik_a * ik_a + jk_a * jk_a, 1.5) * std::sqrt(kk_a * kk_a + lk_a *
    lk_a))) + (pos_ob_x - s) * (yp_dot * std::cos(epsi) + xp_dot * std::sin(epsi))
    * (mk_a * mk_a * (2.0 * ey - 2.0 * pos_ob_y) / (nk_a * nk_a * (ok_a * ok_a +
    pk_a * pk_a)) - 2.0 * ((pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi))
    + yp_dot * std::sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) -
    vel_ob_y) + xp_dot * std::sin(epsi))) * ((yp_dot * std::cos(epsi) - vel_ob_y)
    + xp_dot * std::sin(epsi)) / ((qk_a * qk_a + rk_a * rk_a) * (sk_a * sk_a +
    tk_a * tk_a))) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y *
    vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y *
    xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) + s * yp_dot *
    std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x * xp_dot * std::
                        sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s *
                      xp_dot * std::sin(epsi)) / (2.0 * pow(1.0 - uk_a *
    uk_a / ((vk_a * vk_a + wk_a * wk_a) * (xk_a * xk_a + yk_a * yk_a)), 1.5) *
    pow(al_a * al_a + bl_a * bl_a, 1.5) * std::sqrt(cl_a * cl_a + dl_a *
    dl_a))) + (psi_dot - psi_dot_com) * (2.0 * ey - 2.0 * pos_ob_y) *
    (((((xp_dot * xp_dot + yp_dot * yp_dot) - vel_ob_x * xp_dot * std::cos(epsi))
       - vel_ob_y * yp_dot * std::cos(epsi)) - vel_ob_y * xp_dot * std::sin(epsi))
     + vel_ob_x * yp_dot * std::sin(epsi)) * (((((((((((ey * vel_ob_x + pos_ob_x
    * vel_ob_y) - pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos
    (epsi)) + pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos
    (epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) -
    pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) +
    s * xp_dot * std::sin(epsi)) / (2.0 * std::sqrt(1.0 - el_a * el_a / ((fl_a *
    fl_a + gl_a * gl_a) * (hl_a * hl_a + il_a * il_a))) * pow(jl_a *
    jl_a + kl_a * kl_a, 1.5) * pow(ll_a * ll_a + ml_a * ml_a, 1.5))) -
    ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi)) *
    ((vel_ob_x * std::cos(epsi) - xp_dot) + vel_ob_y * std::sin(epsi)) * ((((2.0
    * cf * yp_dot + 2.0 * cr * yp_dot) + m * psi_dot * (xp_dot * xp_dot)) + 2.0 *
    a * cf * psi_dot) - 2.0 * b * cr * psi_dot) / (m * xp_dot * std::sqrt(1.0 -
    nl_a * nl_a / ((ol_a * ol_a + pl_a * pl_a) * (ql_a * ql_a + rl_a * rl_a))) *
    std::sqrt(sl_a * sl_a + tl_a * tl_a) * pow(ul_a * ul_a + vl_a * vl_a,
    1.5))) + (2.0 * ey - 2.0 * pos_ob_y) * ((vel_ob_x * std::cos(epsi) - xp_dot)
    + vel_ob_y * std::sin(epsi)) * ((((2.0 * cf * yp_dot + 2.0 * cr * yp_dot) +
    m * psi_dot * (xp_dot * xp_dot)) + 2.0 * a * cf * psi_dot) - 2.0 * b * cr *
    psi_dot) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y *
    vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y *
                      xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos
                     (epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::
                   sin(epsi)) - pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y *
                 yp_dot * std::sin(epsi)) + s * xp_dot * std::sin(epsi)) / (2.0 *
    m * xp_dot * std::sqrt(1.0 - wl_a * wl_a / ((xl_a * xl_a + yl_a * yl_a) *
    (am_a * am_a + bm_a * bm_a))) * pow(cm_a * cm_a + dm_a * dm_a, 1.5) *
    pow(em_a * em_a + fm_a * fm_a, 1.5))) + (gm_a * gm_a * (2.0 * ey -
    2.0 * pos_ob_y) / (hm_a * hm_a * (im_a * im_a + jm_a * jm_a)) - 2.0 *
    ((pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin
                       (epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) -
    vel_ob_y) + xp_dot * std::sin(epsi))) * ((yp_dot * std::cos(epsi) - vel_ob_y)
    + xp_dot * std::sin(epsi)) / ((km_a * km_a + lm_a * lm_a) * (mm_a * mm_a +
    nm_a * nm_a))) * ((vel_ob_x * std::cos(epsi) - xp_dot) + vel_ob_y * std::sin
                      (epsi)) * ((((2.0 * cf * yp_dot + 2.0 * cr * yp_dot) + m *
    psi_dot * (xp_dot * xp_dot)) + 2.0 * a * cf * psi_dot) - 2.0 * b * cr *
    psi_dot) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y *
    vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y *
                      xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos
                     (epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::
                   sin(epsi)) - pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y *
                 yp_dot * std::sin(epsi)) + s * xp_dot * std::sin(epsi)) / (2.0 *
    m * xp_dot * pow(1.0 - om_a * om_a / ((pm_a * pm_a + qm_a * qm_a) *
    (rm_a * rm_a + sm_a * sm_a)), 1.5) * std::sqrt(tm_a * tm_a + um_a * um_a) *
    pow(vm_a * vm_a + wm_a * wm_a, 1.5)))) - (psi_dot - psi_dot_com) *
             (((((((((((((((((ey - pos_ob_y) * (xp_dot * std::cos(epsi) - yp_dot
    * std::sin(epsi)) * (((((((ey * yp_dot * std::cos(epsi) - pos_ob_x * xp_dot *
    std::cos(epsi)) - pos_ob_y * yp_dot * std::cos(epsi)) + s * xp_dot * std::
    cos(epsi)) + ey * xp_dot * std::sin(epsi)) - pos_ob_y * xp_dot * std::sin
    (epsi)) + pos_ob_x * yp_dot * std::sin(epsi)) - s * yp_dot * std::sin(epsi))
    / (std::sqrt(1.0 - xm_a * xm_a / ((ym_a * ym_a + an_a * an_a) * (bn_a * bn_a
    + cn_a * cn_a))) * pow(dn_a * dn_a + en_a * en_a, 1.5) * std::sqrt
       (fn_a * fn_a + gn_a * gn_a)) - Ds * (((((((ey * xp_dot * std::cos(epsi) -
    pos_ob_y * xp_dot * std::cos(epsi)) + pos_ob_x * yp_dot * std::cos(epsi)) -
    s * yp_dot * std::cos(epsi)) - ey * yp_dot * std::sin(epsi)) + pos_ob_x *
    xp_dot * std::sin(epsi)) + pos_ob_y * yp_dot * std::sin(epsi)) - s * xp_dot *
    std::sin(epsi)) / (std::sqrt(1.0 - Ds * Ds / (hn_a * hn_a + in_a * in_a)) *
                       pow(jn_a * jn_a + kn_a * kn_a, 1.5))) + (pos_ob_x
    - s) * (yp_dot * std::cos(epsi) + xp_dot * std::sin(epsi)) * (((((((ey *
    yp_dot * std::cos(epsi) - pos_ob_x * xp_dot * std::cos(epsi)) - pos_ob_y *
    yp_dot * std::cos(epsi)) + s * xp_dot * std::cos(epsi)) + ey * xp_dot * std::
    sin(epsi)) - pos_ob_y * xp_dot * std::sin(epsi)) + pos_ob_x * yp_dot * std::
    sin(epsi)) - s * yp_dot * std::sin(epsi)) / (std::sqrt(1.0 - ln_a * ln_a /
    ((mn_a * mn_a + nn_a * nn_a) * (on_a * on_a + pn_a * pn_a))) * pow
    (qn_a * qn_a + rn_a * rn_a, 1.5) * std::sqrt(sn_a * sn_a + tn_a * tn_a))) -
    (ey - pos_ob_y) * (yp_dot * std::cos(epsi) + xp_dot * std::sin(epsi)) *
    (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y * vel_ob_x) - s *
    vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y * xp_dot * std::cos
    (epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) + s * yp_dot * std::cos(epsi))
    + ey * yp_dot * std::sin(epsi)) - pos_ob_x * xp_dot * std::sin(epsi)) -
      pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot * std::sin(epsi)) / (std::
    sqrt(1.0 - un_a * un_a / ((vn_a * vn_a + wn_a * wn_a) * (xn_a * xn_a + yn_a *
    yn_a))) * pow(ao_a * ao_a + bo_a * bo_a, 1.5) * std::sqrt(co_a *
    co_a + do_a * do_a))) + (pos_ob_x - s) * (xp_dot * std::cos(epsi) - yp_dot *
    std::sin(epsi)) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y
    * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y *
    xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) + s * yp_dot *
    std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x * xp_dot * std::
    sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot * std::sin
                       (epsi)) / (std::sqrt(1.0 - eo_a * eo_a / ((fo_a * fo_a +
    go_a * go_a) * (ho_a * ho_a + io_a * io_a))) * pow(jo_a * jo_a +
    ko_a * ko_a, 1.5) * std::sqrt(lo_a * lo_a + mo_a * mo_a))) + (psi_dot -
    psi_dot_com) * (((vel_ob_x * yp_dot * std::cos(epsi) - vel_ob_y * xp_dot *
                      std::cos(epsi)) + vel_ob_x * xp_dot * std::sin(epsi)) +
                    vel_ob_y * yp_dot * std::sin(epsi)) * (((((((((((ey *
    vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey *
    xp_dot * std::cos(epsi)) + pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x *
    yp_dot * std::cos(epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::
    sin(epsi)) - pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::
    sin(epsi)) + s * xp_dot * std::sin(epsi)) / (std::sqrt(1.0 - no_a * no_a /
    ((oo_a * oo_a + po_a * po_a) * (qo_a * qo_a + ro_a * ro_a))) * std::sqrt
    (so_a * so_a + to_a * to_a) * pow(uo_a * uo_a + vo_a * vo_a, 1.5)))
                        + (psi_dot - psi_dot_com) * (((((xp_dot * xp_dot +
    yp_dot * yp_dot) - vel_ob_x * xp_dot * std::cos(epsi)) - vel_ob_y * yp_dot *
    std::cos(epsi)) - vel_ob_y * xp_dot * std::sin(epsi)) + vel_ob_x * yp_dot *
    std::sin(epsi)) * (((((((ey * yp_dot * std::cos(epsi) - pos_ob_x * xp_dot *
    std::cos(epsi)) - pos_ob_y * yp_dot * std::cos(epsi)) + s * xp_dot * std::
    cos(epsi)) + ey * xp_dot * std::sin(epsi)) - pos_ob_y * xp_dot * std::sin
    (epsi)) + pos_ob_x * yp_dot * std::sin(epsi)) - s * yp_dot * std::sin(epsi))
                        / (std::sqrt(1.0 - wo_a * wo_a / ((xo_a * xo_a + yo_a *
    yo_a) * (ap_a * ap_a + bp_a * bp_a))) * std::sqrt(cp_a * cp_a + dp_a * dp_a)
    * pow(ep_a * ep_a + fp_a * fp_a, 1.5))) - (pos_ob_x - s) * (yp_dot *
    std::cos(epsi) + xp_dot * std::sin(epsi)) * (((2.0 * vel_ob_x * yp_dot * std::
    cos(epsi) - 2.0 * vel_ob_y * xp_dot * std::cos(epsi)) + 2.0 * vel_ob_x *
    xp_dot * std::sin(epsi)) + 2.0 * vel_ob_y * yp_dot * std::sin(epsi)) *
                       (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) -
    pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) +
    pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) +
    s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x *
    xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot *
                        std::sin(epsi)) / (2.0 * std::sqrt(1.0 - gp_a * gp_a /
    ((hp_a * hp_a + ip_a * ip_a) * (jp_a * jp_a + kp_a * kp_a))) * pow
    (lp_a * lp_a + mp_a * mp_a, 1.5) * pow(np_a * np_a + op_a * op_a,
    1.5))) + (pos_ob_x - s) * (yp_dot * std::cos(epsi) + xp_dot * std::sin(epsi))
                      * (2.0 * ((pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos
    (epsi)) + yp_dot * std::sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos
    (epsi) - vel_ob_y) + xp_dot * std::sin(epsi))) * ((ey - pos_ob_y) * (xp_dot *
    std::cos(epsi) - yp_dot * std::sin(epsi)) + (pos_ob_x - s) * (yp_dot * std::
    cos(epsi) + xp_dot * std::sin(epsi))) / ((pp_a * pp_a + qp_a * qp_a) * (rp_a
    * rp_a + sp_a * sp_a)) - tp_a * tp_a * (((2.0 * vel_ob_x * yp_dot * std::cos
    (epsi) - 2.0 * vel_ob_y * xp_dot * std::cos(epsi)) + 2.0 * vel_ob_x * xp_dot
    * std::sin(epsi)) + 2.0 * vel_ob_y * yp_dot * std::sin(epsi)) / ((up_a *
    up_a + vp_a * vp_a) * (wp_a * wp_a))) * (((((((((((ey * vel_ob_x + pos_ob_x *
    vel_ob_y) - pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos
    (epsi)) + pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos
    (epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) -
    pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) +
    s * xp_dot * std::sin(epsi)) / (2.0 * pow(1.0 - xp_a * xp_a / ((yp_a
    * yp_a + aq_a * aq_a) * (bq_a * bq_a + cq_a * cq_a)), 1.5) * pow
    (dq_a * dq_a + eq_a * eq_a, 1.5) * std::sqrt(fq_a * fq_a + gq_a * gq_a))) -
                     3.0 * (psi_dot - psi_dot_com) * (((2.0 * vel_ob_x * yp_dot *
    std::cos(epsi) - 2.0 * vel_ob_y * xp_dot * std::cos(epsi)) + 2.0 * vel_ob_x *
    xp_dot * std::sin(epsi)) + 2.0 * vel_ob_y * yp_dot * std::sin(epsi)) *
                     (((((xp_dot * xp_dot + yp_dot * yp_dot) - vel_ob_x * xp_dot
    * std::cos(epsi)) - vel_ob_y * yp_dot * std::cos(epsi)) - vel_ob_y * xp_dot *
                       std::sin(epsi)) + vel_ob_x * yp_dot * std::sin(epsi)) *
                     (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y *
    vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y *
    xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) + s * yp_dot *
    std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x * xp_dot * std::
                        sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s *
                      xp_dot * std::sin(epsi)) / (2.0 * std::sqrt(1.0 - hq_a *
    hq_a / ((iq_a * iq_a + jq_a * jq_a) * (kq_a * kq_a + lq_a * lq_a))) * std::
    sqrt(mq_a * mq_a + nq_a * nq_a) * pow(oq_a * oq_a + pq_a * pq_a, 2.5)))
                    + (psi_dot - psi_dot_com) * (2.0 * ((pos_ob_x - s) *
    ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi)) + (ey -
    pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi)))
    * ((ey - pos_ob_y) * (xp_dot * std::cos(epsi) - yp_dot * std::sin(epsi)) +
       (pos_ob_x - s) * (yp_dot * std::cos(epsi) + xp_dot * std::sin(epsi))) /
    ((qq_a * qq_a + rq_a * rq_a) * (sq_a * sq_a + tq_a * tq_a)) - uq_a * uq_a *
    (((2.0 * vel_ob_x * yp_dot * std::cos(epsi) - 2.0 * vel_ob_y * xp_dot * std::
       cos(epsi)) + 2.0 * vel_ob_x * xp_dot * std::sin(epsi)) + 2.0 * vel_ob_y *
     yp_dot * std::sin(epsi)) / ((vq_a * vq_a + wq_a * wq_a) * (xq_a * xq_a))) *
                    (((((xp_dot * xp_dot + yp_dot * yp_dot) - vel_ob_x * xp_dot *
                        std::cos(epsi)) - vel_ob_y * yp_dot * std::cos(epsi)) -
                      vel_ob_y * xp_dot * std::sin(epsi)) + vel_ob_x * yp_dot *
                     std::sin(epsi)) * (((((((((((ey * vel_ob_x + pos_ob_x *
    vel_ob_y) - pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos
    (epsi)) + pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos
    (epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) -
    pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) +
    s * xp_dot * std::sin(epsi)) / (2.0 * pow(1.0 - yq_a * yq_a / ((ar_a
    * ar_a + br_a * br_a) * (cr_a * cr_a + dr_a * dr_a)), 1.5) * std::sqrt(er_a *
    er_a + fr_a * fr_a) * pow(gr_a * gr_a + hr_a * hr_a, 1.5))) - (ey -
    pos_ob_y) * (xp_dot * std::cos(epsi) - yp_dot * std::sin(epsi)) * (((2.0 *
    vel_ob_x * yp_dot * std::cos(epsi) - 2.0 * vel_ob_y * xp_dot * std::cos(epsi))
    + 2.0 * vel_ob_x * xp_dot * std::sin(epsi)) + 2.0 * vel_ob_y * yp_dot * std::
    sin(epsi)) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y *
    vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y *
                        xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos
                       (epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot *
                     std::sin(epsi)) - pos_ob_x * xp_dot * std::sin(epsi)) -
                   pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot * std::sin
                  (epsi)) / (2.0 * std::sqrt(1.0 - ir_a * ir_a / ((jr_a * jr_a +
    kr_a * kr_a) * (lr_a * lr_a + mr_a * mr_a))) * pow(nr_a * nr_a +
    or_a * or_a, 1.5) * pow(pr_a * pr_a + qr_a * qr_a, 1.5))) + (ey -
    pos_ob_y) * (xp_dot * std::cos(epsi) - yp_dot * std::sin(epsi)) * (2.0 *
    ((pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin
                       (epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) -
    vel_ob_y) + xp_dot * std::sin(epsi))) * ((ey - pos_ob_y) * (xp_dot * std::
    cos(epsi) - yp_dot * std::sin(epsi)) + (pos_ob_x - s) * (yp_dot * std::cos
    (epsi) + xp_dot * std::sin(epsi))) / ((rr_a * rr_a + sr_a * sr_a) * (tr_a *
    tr_a + ur_a * ur_a)) - vr_a * vr_a * (((2.0 * vel_ob_x * yp_dot * std::cos
    (epsi) - 2.0 * vel_ob_y * xp_dot * std::cos(epsi)) + 2.0 * vel_ob_x * xp_dot
    * std::sin(epsi)) + 2.0 * vel_ob_y * yp_dot * std::sin(epsi)) / ((wr_a *
    wr_a + xr_a * xr_a) * (yr_a * yr_a))) * (((((((((((ey * vel_ob_x + pos_ob_x *
    vel_ob_y) - pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos
    (epsi)) + pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos
    (epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) -
    pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) +
    s * xp_dot * std::sin(epsi)) / (2.0 * pow(1.0 - as_a * as_a / ((bs_a
    * bs_a + cs_a * cs_a) * (ds_a * ds_a + es_a * es_a)), 1.5) * pow
    (fs_a * fs_a + gs_a * gs_a, 1.5) * std::sqrt(hs_a * hs_a + is_a * is_a))) +
                 ((vel_ob_x * std::cos(epsi) - xp_dot) + vel_ob_y * std::sin
                  (epsi)) * ((((2.0 * cf * yp_dot + 2.0 * cr * yp_dot) + m *
    psi_dot * (xp_dot * xp_dot)) + 2.0 * a * cf * psi_dot) - 2.0 * b * cr *
    psi_dot) * (((((((ey * yp_dot * std::cos(epsi) - pos_ob_x * xp_dot * std::
                      cos(epsi)) - pos_ob_y * yp_dot * std::cos(epsi)) + s *
                    xp_dot * std::cos(epsi)) + ey * xp_dot * std::sin(epsi)) -
                  pos_ob_y * xp_dot * std::sin(epsi)) + pos_ob_x * yp_dot * std::
                 sin(epsi)) - s * yp_dot * std::sin(epsi)) / (m * xp_dot * std::
    sqrt(1.0 - js_a * js_a / ((ks_a * ks_a + ls_a * ls_a) * (ms_a * ms_a + ns_a *
    ns_a))) * std::sqrt(os_a * os_a + ps_a * ps_a) * pow(qs_a * qs_a +
    rs_a * rs_a, 1.5))) + (vel_ob_y * std::cos(epsi) - vel_ob_x * std::sin(epsi))
                * ((((2.0 * cf * yp_dot + 2.0 * cr * yp_dot) + m * psi_dot *
                     (xp_dot * xp_dot)) + 2.0 * a * cf * psi_dot) - 2.0 * b * cr
                   * psi_dot) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y)
    - pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) +
    pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) +
    s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x *
    xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot *
    std::sin(epsi)) / (m * xp_dot * std::sqrt(1.0 - ss_a * ss_a / ((ts_a * ts_a
    + us_a * us_a) * (vs_a * vs_a + ws_a * ws_a))) * std::sqrt(xs_a * xs_a +
    ys_a * ys_a) * pow(at_a * at_a + bt_a * bt_a, 1.5))) - 3.0 *
               ((vel_ob_x * std::cos(epsi) - xp_dot) + vel_ob_y * std::sin(epsi))
               * (((2.0 * vel_ob_x * yp_dot * std::cos(epsi) - 2.0 * vel_ob_y *
                    xp_dot * std::cos(epsi)) + 2.0 * vel_ob_x * xp_dot * std::
                   sin(epsi)) + 2.0 * vel_ob_y * yp_dot * std::sin(epsi)) *
               ((((2.0 * cf * yp_dot + 2.0 * cr * yp_dot) + m * psi_dot *
                  (xp_dot * xp_dot)) + 2.0 * a * cf * psi_dot) - 2.0 * b * cr *
                psi_dot) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) -
    pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) +
    pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) +
    s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x *
    xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot *
    std::sin(epsi)) / (2.0 * m * xp_dot * std::sqrt(1.0 - ct_a * ct_a / ((dt_a *
    dt_a + et_a * et_a) * (ft_a * ft_a + gt_a * gt_a))) * std::sqrt(ht_a * ht_a
    + it_a * it_a) * pow(jt_a * jt_a + kt_a * kt_a, 2.5))) + (2.0 *
    ((pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin
                       (epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) -
    vel_ob_y) + xp_dot * std::sin(epsi))) * ((ey - pos_ob_y) * (xp_dot * std::
    cos(epsi) - yp_dot * std::sin(epsi)) + (pos_ob_x - s) * (yp_dot * std::cos
    (epsi) + xp_dot * std::sin(epsi))) / ((lt_a * lt_a + mt_a * mt_a) * (nt_a *
    nt_a + ot_a * ot_a)) - pt_a * pt_a * (((2.0 * vel_ob_x * yp_dot * std::cos
    (epsi) - 2.0 * vel_ob_y * xp_dot * std::cos(epsi)) + 2.0 * vel_ob_x * xp_dot
    * std::sin(epsi)) + 2.0 * vel_ob_y * yp_dot * std::sin(epsi)) / ((qt_a *
    qt_a + rt_a * rt_a) * (st_a * st_a))) * ((vel_ob_x * std::cos(epsi) - xp_dot)
    + vel_ob_y * std::sin(epsi)) * ((((2.0 * cf * yp_dot + 2.0 * cr * yp_dot) +
    m * psi_dot * (xp_dot * xp_dot)) + 2.0 * a * cf * psi_dot) - 2.0 * b * cr *
    psi_dot) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y *
    vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y *
                      xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos
                     (epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::
                   sin(epsi)) - pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y *
                 yp_dot * std::sin(epsi)) + s * xp_dot * std::sin(epsi)) / (2.0 *
    m * xp_dot * pow(1.0 - tt_a * tt_a / ((ut_a * ut_a + vt_a * vt_a) *
    (wt_a * wt_a + xt_a * xt_a)), 1.5) * std::sqrt(yt_a * yt_a + au_a * au_a) *
    pow(bu_a * bu_a + cu_a * cu_a, 1.5)))) - ((((2.0 * cf * yp_dot + 2.0
    * cr * yp_dot) + m * psi_dot * (xp_dot * xp_dot)) + 2.0 * a * cf * psi_dot)
             - 2.0 * b * cr * psi_dot) * ((((((((((((((((Ds * (((ey * std::cos
    (epsi) - pos_ob_y * std::cos(epsi)) + pos_ob_x * std::sin(epsi)) - s * std::
    sin(epsi)) / (std::sqrt(1.0 - Ds * Ds / (du_a * du_a + eu_a * eu_a)) *
                  pow(fu_a * fu_a + gu_a * gu_a, 1.5)) - std::cos(epsi) *
    (pos_ob_x - s) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y *
    vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y *
    xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) + s * yp_dot *
    std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x * xp_dot * std::
    sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot * std::sin
                      (epsi)) / (std::sqrt(1.0 - hu_a * hu_a / ((iu_a * iu_a +
    ju_a * ju_a) * (ku_a * ku_a + lu_a * lu_a))) * pow(mu_a * mu_a +
    nu_a * nu_a, 1.5) * std::sqrt(ou_a * ou_a + pu_a * pu_a))) + std::sin(epsi) *
    (ey - pos_ob_y) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y
    * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y *
    xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) + s * yp_dot *
    std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x * xp_dot * std::
    sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot * std::sin
                       (epsi)) / (std::sqrt(1.0 - qu_a * qu_a / ((ru_a * ru_a +
    su_a * su_a) * (tu_a * tu_a + uu_a * uu_a))) * pow(vu_a * vu_a +
    wu_a * wu_a, 1.5) * std::sqrt(xu_a * xu_a + yu_a * yu_a))) + (psi_dot -
    psi_dot_com) * (((pos_ob_x * std::cos(epsi) - s * std::cos(epsi)) - ey * std::
                     sin(epsi)) + pos_ob_y * std::sin(epsi)) * (((((xp_dot *
    xp_dot + yp_dot * yp_dot) - vel_ob_x * xp_dot * std::cos(epsi)) - vel_ob_y *
    yp_dot * std::cos(epsi)) - vel_ob_y * xp_dot * std::sin(epsi)) + vel_ob_x *
    yp_dot * std::sin(epsi)) / (std::sqrt(1.0 - av_a * av_a / ((bv_a * bv_a +
    cv_a * cv_a) * (dv_a * dv_a + ev_a * ev_a))) * std::sqrt(fv_a * fv_a + gv_a *
    gv_a) * pow(hv_a * hv_a + iv_a * iv_a, 1.5))) + (ey - pos_ob_y) *
    (xp_dot * std::cos(epsi) - yp_dot * std::sin(epsi)) * (((pos_ob_x * std::cos
    (epsi) - s * std::cos(epsi)) - ey * std::sin(epsi)) + pos_ob_y * std::sin
    (epsi)) / (std::sqrt(1.0 - jv_a * jv_a / ((kv_a * kv_a + lv_a * lv_a) *
    (mv_a * mv_a + nv_a * nv_a))) * pow(ov_a * ov_a + pv_a * pv_a, 1.5) *
               std::sqrt(qv_a * qv_a + rv_a * rv_a))) + (pos_ob_x - s) * (yp_dot
    * std::cos(epsi) + xp_dot * std::sin(epsi)) * (((pos_ob_x * std::cos(epsi) -
    s * std::cos(epsi)) - ey * std::sin(epsi)) + pos_ob_y * std::sin(epsi)) /
    (std::sqrt(1.0 - sv_a * sv_a / ((tv_a * tv_a + uv_a * uv_a) * (vv_a * vv_a +
    wv_a * wv_a))) * pow(xv_a * xv_a + yv_a * yv_a, 1.5) * std::sqrt
     (aw_a * aw_a + bw_a * bw_a))) - (psi_dot - psi_dot_com) * ((2.0 * yp_dot -
    vel_ob_y * std::cos(epsi)) + vel_ob_x * std::sin(epsi)) * (((((((((((ey *
    vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey *
    xp_dot * std::cos(epsi)) + pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x *
    yp_dot * std::cos(epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::
    sin(epsi)) - pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::
    sin(epsi)) + s * xp_dot * std::sin(epsi)) / (std::sqrt(1.0 - cw_a * cw_a /
    ((dw_a * dw_a + ew_a * ew_a) * (fw_a * fw_a + gw_a * gw_a))) * std::sqrt
    (hw_a * hw_a + iw_a * iw_a) * pow(jw_a * jw_a + kw_a * kw_a, 1.5)))
    + 3.0 * (psi_dot - psi_dot_com) * ((2.0 * yp_dot - 2.0 * vel_ob_y * std::cos
    (epsi)) + 2.0 * vel_ob_x * std::sin(epsi)) * (((((xp_dot * xp_dot + yp_dot *
    yp_dot) - vel_ob_x * xp_dot * std::cos(epsi)) - vel_ob_y * yp_dot * std::cos
    (epsi)) - vel_ob_y * xp_dot * std::sin(epsi)) + vel_ob_x * yp_dot * std::sin
    (epsi)) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y *
    vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y *
                     xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos
                    (epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::
                  sin(epsi)) - pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y *
                yp_dot * std::sin(epsi)) + s * xp_dot * std::sin(epsi)) / (2.0 *
    std::sqrt(1.0 - lw_a * lw_a / ((mw_a * mw_a + nw_a * nw_a) * (ow_a * ow_a +
    pw_a * pw_a))) * std::sqrt(qw_a * qw_a + rw_a * rw_a) * pow(sw_a *
    sw_a + tw_a * tw_a, 2.5))) + (ey - pos_ob_y) * (xp_dot * std::cos(epsi) -
    yp_dot * std::sin(epsi)) * ((2.0 * yp_dot - 2.0 * vel_ob_y * std::cos(epsi))
    + 2.0 * vel_ob_x * std::sin(epsi)) * (((((((((((ey * vel_ob_x + pos_ob_x *
    vel_ob_y) - pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos
    (epsi)) + pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos
    (epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) -
    pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) +
    s * xp_dot * std::sin(epsi)) / (2.0 * std::sqrt(1.0 - uw_a * uw_a / ((vw_a *
    vw_a + ww_a * ww_a) * (xw_a * xw_a + yw_a * yw_a))) * pow(ax_a *
    ax_a + bx_a * bx_a, 1.5) * pow(cx_a * cx_a + dx_a * dx_a, 1.5))) +
    (pos_ob_x - s) * (yp_dot * std::cos(epsi) + xp_dot * std::sin(epsi)) * ((2.0
    * yp_dot - 2.0 * vel_ob_y * std::cos(epsi)) + 2.0 * vel_ob_x * std::sin(epsi))
    * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y * vel_ob_x) - s
               * vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y * xp_dot *
             std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) + s * yp_dot *
           std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x * xp_dot *
         std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot *
       std::sin(epsi)) / (2.0 * std::sqrt(1.0 - ex_a * ex_a / ((fx_a * fx_a +
    gx_a * gx_a) * (hx_a * hx_a + ix_a * ix_a))) * pow(jx_a * jx_a +
    kx_a * kx_a, 1.5) * pow(lx_a * lx_a + mx_a * mx_a, 1.5))) - (psi_dot
    - psi_dot_com) * (2.0 * ((pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos
    (epsi)) + yp_dot * std::sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos
    (epsi) - vel_ob_y) + xp_dot * std::sin(epsi))) * (((ey * std::cos(epsi) -
    pos_ob_y * std::cos(epsi)) + pos_ob_x * std::sin(epsi)) - s * std::sin(epsi))
                      / ((nx_a * nx_a + ox_a * ox_a) * (px_a * px_a + qx_a *
    qx_a)) - rx_a * rx_a * ((2.0 * yp_dot - 2.0 * vel_ob_y * std::cos(epsi)) +
    2.0 * vel_ob_x * std::sin(epsi)) / ((sx_a * sx_a + tx_a * tx_a) * (ux_a *
    ux_a))) * (((((xp_dot * xp_dot + yp_dot * yp_dot) - vel_ob_x * xp_dot * std::
                  cos(epsi)) - vel_ob_y * yp_dot * std::cos(epsi)) - vel_ob_y *
                xp_dot * std::sin(epsi)) + vel_ob_x * yp_dot * std::sin(epsi)) *
    (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y * vel_ob_x) - s *
             vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y * xp_dot * std::
           cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) + s * yp_dot * std::
         cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x * xp_dot * std::
       sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot * std::sin
     (epsi)) / (2.0 * pow(1.0 - vx_a * vx_a / ((wx_a * wx_a + xx_a *
    xx_a) * (yx_a * yx_a + ay_a * ay_a)), 1.5) * std::sqrt(by_a * by_a + cy_a *
    cy_a) * pow(dy_a * dy_a + ey_a * ey_a, 1.5))) - (ey - pos_ob_y) *
    (2.0 * ((pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot *
    std::sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi))) * (((ey * std::cos(epsi) - pos_ob_y * std::cos
    (epsi)) + pos_ob_x * std::sin(epsi)) - s * std::sin(epsi)) / ((fy_a * fy_a +
    gy_a * gy_a) * (hy_a * hy_a + iy_a * iy_a)) - jy_a * jy_a * ((2.0 * yp_dot -
    2.0 * vel_ob_y * std::cos(epsi)) + 2.0 * vel_ob_x * std::sin(epsi)) / ((ky_a
    * ky_a + ly_a * ly_a) * (my_a * my_a))) * (xp_dot * std::cos(epsi) - yp_dot *
    std::sin(epsi)) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y
    * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y *
    xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) + s * yp_dot *
    std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x * xp_dot * std::
    sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot * std::sin
                       (epsi)) / (2.0 * pow(1.0 - ny_a * ny_a / ((oy_a *
    oy_a + py_a * py_a) * (qy_a * qy_a + ry_a * ry_a)), 1.5) * pow(sy_a *
    sy_a + ty_a * ty_a, 1.5) * std::sqrt(uy_a * uy_a + vy_a * vy_a))) -
    (pos_ob_x - s) * (2.0 * ((pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos
    (epsi)) + yp_dot * std::sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos
    (epsi) - vel_ob_y) + xp_dot * std::sin(epsi))) * (((ey * std::cos(epsi) -
    pos_ob_y * std::cos(epsi)) + pos_ob_x * std::sin(epsi)) - s * std::sin(epsi))
                      / ((wy_a * wy_a + xy_a * xy_a) * (yy_a * yy_a + aab_a *
    aab_a)) - bab_a * bab_a * ((2.0 * yp_dot - 2.0 * vel_ob_y * std::cos(epsi))
    + 2.0 * vel_ob_x * std::sin(epsi)) / ((cab_a * cab_a + dab_a * dab_a) *
    (eab_a * eab_a))) * (yp_dot * std::cos(epsi) + xp_dot * std::sin(epsi)) *
    (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y * vel_ob_x) - s *
             vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y * xp_dot * std::
           cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) + s * yp_dot * std::
         cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x * xp_dot * std::
       sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot * std::sin
     (epsi)) / (2.0 * pow(1.0 - fab_a * fab_a / ((gab_a * gab_a + hab_a *
    hab_a) * (iab_a * iab_a + jab_a * jab_a)), 1.5) * pow(kab_a * kab_a
    + lab_a * lab_a, 1.5) * std::sqrt(mab_a * mab_a + nab_a * nab_a))) +
    ((vel_ob_x * std::cos(epsi) - xp_dot) + vel_ob_y * std::sin(epsi)) *
    (((pos_ob_x * std::cos(epsi) - s * std::cos(epsi)) - ey * std::sin(epsi)) +
     pos_ob_y * std::sin(epsi)) * ((((2.0 * cf * yp_dot + 2.0 * cr * yp_dot) + m
    * psi_dot * (xp_dot * xp_dot)) + 2.0 * a * cf * psi_dot) - 2.0 * b * cr *
    psi_dot) / (m * xp_dot * std::sqrt(1.0 - oab_a * oab_a / ((pab_a * pab_a +
    qab_a * qab_a) * (rab_a * rab_a + sab_a * sab_a))) * std::sqrt(tab_a * tab_a
    + uab_a * uab_a) * pow(vab_a * vab_a + wab_a * wab_a, 1.5))) - (2.0 *
    cf + 2.0 * cr) * ((vel_ob_x * std::cos(epsi) - xp_dot) + vel_ob_y * std::sin
                      (epsi)) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y)
    - pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) +
    pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) +
    s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x *
    xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot *
    std::sin(epsi)) / (m * xp_dot * std::sqrt(1.0 - xab_a * xab_a / ((yab_a *
    yab_a + abb_a * abb_a) * (bbb_a * bbb_a + cbb_a * cbb_a))) * std::sqrt(dbb_a
    * dbb_a + ebb_a * ebb_a) * pow(fbb_a * fbb_a + gbb_a * gbb_a, 1.5)))
              + 3.0 * ((vel_ob_x * std::cos(epsi) - xp_dot) + vel_ob_y * std::
                       sin(epsi)) * ((2.0 * yp_dot - 2.0 * vel_ob_y * std::cos
    (epsi)) + 2.0 * vel_ob_x * std::sin(epsi)) * ((((2.0 * cf * yp_dot + 2.0 *
    cr * yp_dot) + m * psi_dot * (xp_dot * xp_dot)) + 2.0 * a * cf * psi_dot) -
    2.0 * b * cr * psi_dot) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) -
    pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) +
    pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) +
    s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x *
    xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot *
    std::sin(epsi)) / (2.0 * m * xp_dot * std::sqrt(1.0 - hbb_a * hbb_a /
    ((ibb_a * ibb_a + jbb_a * jbb_a) * (kbb_a * kbb_a + lbb_a * lbb_a))) * std::
                       sqrt(mbb_a * mbb_a + nbb_a * nbb_a) * pow(obb_a *
    obb_a + pbb_a * pbb_a, 2.5))) - (2.0 * ((pos_ob_x - s) * ((vel_ob_x - xp_dot
    * std::cos(epsi)) + yp_dot * std::sin(epsi)) + (ey - pos_ob_y) * ((yp_dot *
    std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi))) * (((ey * std::cos
    (epsi) - pos_ob_y * std::cos(epsi)) + pos_ob_x * std::sin(epsi)) - s * std::
    sin(epsi)) / ((qbb_a * qbb_a + rbb_a * rbb_a) * (sbb_a * sbb_a + tbb_a *
    tbb_a)) - ubb_a * ubb_a * ((2.0 * yp_dot - 2.0 * vel_ob_y * std::cos(epsi))
    + 2.0 * vel_ob_x * std::sin(epsi)) / ((vbb_a * vbb_a + wbb_a * wbb_a) *
    (xbb_a * xbb_a))) * ((vel_ob_x * std::cos(epsi) - xp_dot) + vel_ob_y * std::
              sin(epsi)) * ((((2.0 * cf * yp_dot + 2.0 * cr * yp_dot) + m *
    psi_dot * (xp_dot * xp_dot)) + 2.0 * a * cf * psi_dot) - 2.0 * b * cr *
              psi_dot) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) -
    pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) +
    pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) +
    s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x *
    xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot *
              std::sin(epsi)) / (2.0 * m * xp_dot * pow(1.0 - ybb_a *
    ybb_a / ((acb_a * acb_a + bcb_a * bcb_a) * (ccb_a * ccb_a + dcb_a * dcb_a)),
    1.5) * std::sqrt(ecb_a * ecb_a + fcb_a * fcb_a) * pow(gcb_a * gcb_a
    + hcb_a * hcb_a, 1.5))) / (m * xp_dot)) + 2.0 * (((a * cf * yp_dot - b * cr *
    yp_dot) + a * a * cf * psi_dot) + b * b * cr * psi_dot) * ((((((((m * xp_dot
    * (yp_dot * yp_dot) - 2.0 * a * cf * xp_dot) + 2.0 * b * cr * xp_dot) + 2.0 *
    a * cf * vel_ob_x * std::cos(epsi)) - 2.0 * b * cr * vel_ob_x * std::cos
    (epsi)) + 2.0 * a * cf * vel_ob_y * std::sin(epsi)) - 2.0 * b * cr *
    vel_ob_y * std::sin(epsi)) + m * vel_ob_x * xp_dot * yp_dot * std::sin(epsi))
    - m * vel_ob_y * xp_dot * yp_dot * std::cos(epsi)) * (((((((((((ey *
    vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey *
    xp_dot * std::cos(epsi)) + pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x *
    yp_dot * std::cos(epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::
    sin(epsi)) - pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::
    sin(epsi)) + s * xp_dot * std::sin(epsi)) / (Iz * m * (xp_dot * xp_dot) *
    std::sqrt(1.0 - icb_a * icb_a / ((jcb_a * jcb_a + kcb_a * kcb_a) * (lcb_a *
    lcb_a + mcb_a * mcb_a))) * std::sqrt(ncb_a * ncb_a + ocb_a * ocb_a) *
    pow(pcb_a * pcb_a + qcb_a * qcb_a, 1.5));
  out[1] = (((yp_dot * std::cos(epsi) + xp_dot * std::sin(epsi)) * ((((((3.0 *
    Ds * (2.0 * ey - 2.0 * pos_ob_y) * (((ey * vel_ob_y - pos_ob_x * vel_ob_x) -
    pos_ob_y * vel_ob_y) + s * vel_ob_x) / (2.0 * std::sqrt(1.0 - Ds * Ds /
    (rcb_a * rcb_a + scb_a * scb_a)) * pow(tcb_a * tcb_a + ucb_a * ucb_a,
    2.5)) - Ds * vel_ob_y / (std::sqrt(1.0 - Ds * Ds / (vcb_a * vcb_a + wcb_a *
    wcb_a)) * pow(xcb_a * xcb_a + ycb_a * ycb_a, 1.5))) + (((((((((((ey *
    vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey *
    xp_dot * std::cos(epsi)) + pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x *
    yp_dot * std::cos(epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::
    sin(epsi)) - pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::
    sin(epsi)) + s * xp_dot * std::sin(epsi)) * (((((((((((((((((((vel_ob_x *
    (vel_ob_y * vel_ob_y) + vel_ob_x * (xp_dot * xp_dot)) + vel_ob_x * (yp_dot *
    yp_dot)) + pow(vel_ob_x, 3.0)) - 2.0 * (vel_ob_x * vel_ob_x) *
    xp_dot * std::cos(epsi)) + 2.0 * (vel_ob_x * vel_ob_x) * yp_dot * std::sin
    (epsi)) + 2.0 * acc_ob_x * ey * vel_ob_y) - 2.0 * acc_ob_y * ey * vel_ob_x)
    - 2.0 * acc_ob_x * pos_ob_y * vel_ob_y) + 2.0 * acc_ob_y * pos_ob_y *
    vel_ob_x) + 2.0 * acc_ob_y * ey * xp_dot * std::cos(epsi)) - 2.0 * acc_ob_x *
    ey * yp_dot * std::cos(epsi)) - 2.0 * acc_ob_y * pos_ob_y * xp_dot * std::
    cos(epsi)) + 2.0 * acc_ob_x * pos_ob_y * yp_dot * std::cos(epsi)) - 2.0 *
    acc_ob_x * ey * xp_dot * std::sin(epsi)) - 2.0 * acc_ob_y * ey * yp_dot *
    std::sin(epsi)) + 2.0 * acc_ob_x * pos_ob_y * xp_dot * std::sin(epsi)) + 2.0
    * acc_ob_y * pos_ob_y * yp_dot * std::sin(epsi)) - 2.0 * vel_ob_x * vel_ob_y
    * yp_dot * std::cos(epsi)) - 2.0 * vel_ob_x * vel_ob_y * xp_dot * std::sin
    (epsi)) / (std::sqrt(1.0 - adb_a * adb_a / ((bdb_a * bdb_a + cdb_a * cdb_a) *
    (ddb_a * ddb_a + edb_a * edb_a))) * pow(fdb_a * fdb_a + gdb_a *
    gdb_a, 1.5) * pow(hdb_a * hdb_a + idb_a * idb_a, 1.5))) +
    pow(Ds, 3.0) * (2.0 * ey - 2.0 * pos_ob_y) * (((ey * vel_ob_y -
    pos_ob_x * vel_ob_x) - pos_ob_y * vel_ob_y) + s * vel_ob_x) / (2.0 *
    pow(1.0 - Ds * Ds / (jdb_a * jdb_a + kdb_a * kdb_a), 1.5) *
    pow(ldb_a * ldb_a + mdb_a * mdb_a, 3.5))) - ((vel_ob_x - xp_dot *
    std::cos(epsi)) + yp_dot * std::sin(epsi)) *
    (((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((pos_ob_y *
    pow(vel_ob_x, 3.0) - pos_ob_x * pow(vel_ob_y, 3.0)) - ey *
    pow(vel_ob_x, 3.0)) + s * pow(vel_ob_y, 3.0)) - acc_ob_x *
    (pos_ob_x * pos_ob_x) * vel_ob_y) + acc_ob_y * (pos_ob_x * pos_ob_x) *
    vel_ob_x) - acc_ob_x * (pos_ob_y * pos_ob_y) * vel_ob_y) + acc_ob_y *
    (pos_ob_y * pos_ob_y) * vel_ob_x) - acc_ob_x * (s * s) * vel_ob_y) +
    acc_ob_y * (s * s) * vel_ob_x) - ey * vel_ob_x * (vel_ob_y * vel_ob_y)) - ey
    * vel_ob_x * (xp_dot * xp_dot)) - ey * vel_ob_x * (yp_dot * yp_dot)) -
    pos_ob_x * (vel_ob_x * vel_ob_x) * vel_ob_y) + pos_ob_y * vel_ob_x *
    (vel_ob_y * vel_ob_y)) - pos_ob_x * vel_ob_y * (xp_dot * xp_dot)) + pos_ob_y
    * vel_ob_x * (xp_dot * xp_dot)) + s * (vel_ob_x * vel_ob_x) * vel_ob_y) -
    pos_ob_x * vel_ob_y * (yp_dot * yp_dot)) + pos_ob_y * vel_ob_x * (yp_dot *
    yp_dot)) + s * vel_ob_y * (xp_dot * xp_dot)) + s * vel_ob_y * (yp_dot *
    yp_dot)) - acc_ob_x * (ey * ey) * vel_ob_y) + acc_ob_y * (ey * ey) *
    vel_ob_x) + 2.0 * acc_ob_x * ey * pos_ob_y * vel_ob_y) - 2.0 * acc_ob_y * ey
    * pos_ob_y * vel_ob_x) + 2.0 * acc_ob_x * pos_ob_x * s * vel_ob_y) - 2.0 *
    acc_ob_y * pos_ob_x * s * vel_ob_x) - acc_ob_y * (ey * ey) * xp_dot * std::
    cos(epsi)) + acc_ob_x * (ey * ey) * yp_dot * std::cos(epsi)) - acc_ob_y *
    (pos_ob_x * pos_ob_x) * xp_dot * std::cos(epsi)) - acc_ob_y * (pos_ob_y *
    pos_ob_y) * xp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_x * pos_ob_x) *
    yp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_y * pos_ob_y) * yp_dot * std::
    cos(epsi)) - acc_ob_y * (s * s) * xp_dot * std::cos(epsi)) + acc_ob_x * (s *
    s) * yp_dot * std::cos(epsi)) + acc_ob_x * (ey * ey) * xp_dot * std::sin
    (epsi)) + acc_ob_y * (ey * ey) * yp_dot * std::sin(epsi)) + 2.0 * ey *
    (vel_ob_x * vel_ob_x) * xp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_x *
    pos_ob_x) * xp_dot * std::sin(epsi)) + acc_ob_x * (pos_ob_y * pos_ob_y) *
    xp_dot * std::sin(epsi)) + acc_ob_y * (pos_ob_x * pos_ob_x) * yp_dot * std::
    sin(epsi)) + acc_ob_y * (pos_ob_y * pos_ob_y) * yp_dot * std::sin(epsi)) +
    acc_ob_x * (s * s) * xp_dot * std::sin(epsi)) + acc_ob_y * (s * s) * yp_dot *
    std::sin(epsi)) - 2.0 * pos_ob_y * (vel_ob_x * vel_ob_x) * xp_dot * std::cos
    (epsi)) + 2.0 * pos_ob_x * (vel_ob_y * vel_ob_y) * yp_dot * std::cos(epsi))
    - 2.0 * s * (vel_ob_y * vel_ob_y) * yp_dot * std::cos(epsi)) - 2.0 * ey *
    (vel_ob_x * vel_ob_x) * yp_dot * std::sin(epsi)) + 2.0 * pos_ob_x *
    (vel_ob_y * vel_ob_y) * xp_dot * std::sin(epsi)) + 2.0 * pos_ob_y *
    (vel_ob_x * vel_ob_x) * yp_dot * std::sin(epsi)) - 2.0 * s * (vel_ob_y *
    vel_ob_y) * xp_dot * std::sin(epsi)) - 2.0 * s * vel_ob_x * vel_ob_y *
    xp_dot * std::cos(epsi)) + 2.0 * ey * vel_ob_x * vel_ob_y * xp_dot * std::
    sin(epsi)) - 2.0 * pos_ob_y * vel_ob_x * vel_ob_y * xp_dot * std::sin(epsi))
                 - 2.0 * pos_ob_x * vel_ob_x * vel_ob_y * yp_dot * std::sin(epsi))
                + 2.0 * s * vel_ob_x * vel_ob_y * yp_dot * std::sin(epsi)) + 2.0
               * acc_ob_y * ey * pos_ob_y * xp_dot * std::cos(epsi)) - 2.0 *
              acc_ob_x * ey * pos_ob_y * yp_dot * std::cos(epsi)) + 2.0 *
             acc_ob_y * pos_ob_x * s * xp_dot * std::cos(epsi)) - 2.0 * acc_ob_x
            * pos_ob_x * s * yp_dot * std::cos(epsi)) - 2.0 * acc_ob_x * ey *
           pos_ob_y * xp_dot * std::sin(epsi)) - 2.0 * acc_ob_y * ey * pos_ob_y *
          yp_dot * std::sin(epsi)) + 2.0 * ey * vel_ob_x * vel_ob_y * yp_dot *
         std::cos(epsi)) - 2.0 * acc_ob_x * pos_ob_x * s * xp_dot * std::sin
        (epsi)) - 2.0 * acc_ob_y * pos_ob_x * s * yp_dot * std::sin(epsi)) + 2.0
      * pos_ob_x * vel_ob_x * vel_ob_y * xp_dot * std::cos(epsi)) - 2.0 *
     pos_ob_y * vel_ob_x * vel_ob_y * yp_dot * std::cos(epsi)) / (std::sqrt(1.0
    - ndb_a * ndb_a / ((odb_a * odb_a + pdb_a * pdb_a) * (qdb_a * qdb_a + rdb_a *
    rdb_a))) * pow(sdb_a * sdb_a + tdb_a * tdb_a, 1.5) * pow
    (udb_a * udb_a + vdb_a * vdb_a, 1.5))) + (wdb_a * wdb_a * (2.0 * ey - 2.0 *
    pos_ob_y) / (xdb_a * xdb_a * (ydb_a * ydb_a + aeb_a * aeb_a)) - 2.0 *
    ((pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin
                       (epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) -
    vel_ob_y) + xp_dot * std::sin(epsi))) * ((yp_dot * std::cos(epsi) - vel_ob_y)
    + xp_dot * std::sin(epsi)) / ((beb_a * beb_a + ceb_a * ceb_a) * (deb_a *
    deb_a + eeb_a * eeb_a))) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) -
    pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) +
    pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) +
    s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x *
    xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot *
    std::sin(epsi)) *
    (((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((pos_ob_y *
    pow(vel_ob_x, 3.0) - pos_ob_x * pow(vel_ob_y, 3.0)) - ey *
    pow(vel_ob_x, 3.0)) + s * pow(vel_ob_y, 3.0)) - acc_ob_x *
    (pos_ob_x * pos_ob_x) * vel_ob_y) + acc_ob_y * (pos_ob_x * pos_ob_x) *
    vel_ob_x) - acc_ob_x * (pos_ob_y * pos_ob_y) * vel_ob_y) + acc_ob_y *
    (pos_ob_y * pos_ob_y) * vel_ob_x) - acc_ob_x * (s * s) * vel_ob_y) +
    acc_ob_y * (s * s) * vel_ob_x) - ey * vel_ob_x * (vel_ob_y * vel_ob_y)) - ey
    * vel_ob_x * (xp_dot * xp_dot)) - ey * vel_ob_x * (yp_dot * yp_dot)) -
    pos_ob_x * (vel_ob_x * vel_ob_x) * vel_ob_y) + pos_ob_y * vel_ob_x *
    (vel_ob_y * vel_ob_y)) - pos_ob_x * vel_ob_y * (xp_dot * xp_dot)) + pos_ob_y
    * vel_ob_x * (xp_dot * xp_dot)) + s * (vel_ob_x * vel_ob_x) * vel_ob_y) -
    pos_ob_x * vel_ob_y * (yp_dot * yp_dot)) + pos_ob_y * vel_ob_x * (yp_dot *
    yp_dot)) + s * vel_ob_y * (xp_dot * xp_dot)) + s * vel_ob_y * (yp_dot *
    yp_dot)) - acc_ob_x * (ey * ey) * vel_ob_y) + acc_ob_y * (ey * ey) *
    vel_ob_x) + 2.0 * acc_ob_x * ey * pos_ob_y * vel_ob_y) - 2.0 * acc_ob_y * ey
    * pos_ob_y * vel_ob_x) + 2.0 * acc_ob_x * pos_ob_x * s * vel_ob_y) - 2.0 *
    acc_ob_y * pos_ob_x * s * vel_ob_x) - acc_ob_y * (ey * ey) * xp_dot * std::
    cos(epsi)) + acc_ob_x * (ey * ey) * yp_dot * std::cos(epsi)) - acc_ob_y *
    (pos_ob_x * pos_ob_x) * xp_dot * std::cos(epsi)) - acc_ob_y * (pos_ob_y *
    pos_ob_y) * xp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_x * pos_ob_x) *
    yp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_y * pos_ob_y) * yp_dot * std::
    cos(epsi)) - acc_ob_y * (s * s) * xp_dot * std::cos(epsi)) + acc_ob_x * (s *
    s) * yp_dot * std::cos(epsi)) + acc_ob_x * (ey * ey) * xp_dot * std::sin
    (epsi)) + acc_ob_y * (ey * ey) * yp_dot * std::sin(epsi)) + 2.0 * ey *
    (vel_ob_x * vel_ob_x) * xp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_x *
    pos_ob_x) * xp_dot * std::sin(epsi)) + acc_ob_x * (pos_ob_y * pos_ob_y) *
    xp_dot * std::sin(epsi)) + acc_ob_y * (pos_ob_x * pos_ob_x) * yp_dot * std::
    sin(epsi)) + acc_ob_y * (pos_ob_y * pos_ob_y) * yp_dot * std::sin(epsi)) +
    acc_ob_x * (s * s) * xp_dot * std::sin(epsi)) + acc_ob_y * (s * s) * yp_dot *
    std::sin(epsi)) - 2.0 * pos_ob_y * (vel_ob_x * vel_ob_x) * xp_dot * std::cos
    (epsi)) + 2.0 * pos_ob_x * (vel_ob_y * vel_ob_y) * yp_dot * std::cos(epsi))
    - 2.0 * s * (vel_ob_y * vel_ob_y) * yp_dot * std::cos(epsi)) - 2.0 * ey *
    (vel_ob_x * vel_ob_x) * yp_dot * std::sin(epsi)) + 2.0 * pos_ob_x *
    (vel_ob_y * vel_ob_y) * xp_dot * std::sin(epsi)) + 2.0 * pos_ob_y *
    (vel_ob_x * vel_ob_x) * yp_dot * std::sin(epsi)) - 2.0 * s * (vel_ob_y *
    vel_ob_y) * xp_dot * std::sin(epsi)) - 2.0 * s * vel_ob_x * vel_ob_y *
    xp_dot * std::cos(epsi)) + 2.0 * ey * vel_ob_x * vel_ob_y * xp_dot * std::
                   sin(epsi)) - 2.0 * pos_ob_y * vel_ob_x * vel_ob_y * xp_dot *
                  std::sin(epsi)) - 2.0 * pos_ob_x * vel_ob_x * vel_ob_y *
                 yp_dot * std::sin(epsi)) + 2.0 * s * vel_ob_x * vel_ob_y *
                yp_dot * std::sin(epsi)) + 2.0 * acc_ob_y * ey * pos_ob_y *
               xp_dot * std::cos(epsi)) - 2.0 * acc_ob_x * ey * pos_ob_y *
              yp_dot * std::cos(epsi)) + 2.0 * acc_ob_y * pos_ob_x * s * xp_dot *
             std::cos(epsi)) - 2.0 * acc_ob_x * pos_ob_x * s * yp_dot * std::cos
            (epsi)) - 2.0 * acc_ob_x * ey * pos_ob_y * xp_dot * std::sin(epsi))
          - 2.0 * acc_ob_y * ey * pos_ob_y * yp_dot * std::sin(epsi)) + 2.0 * ey
         * vel_ob_x * vel_ob_y * yp_dot * std::cos(epsi)) - 2.0 * acc_ob_x *
        pos_ob_x * s * xp_dot * std::sin(epsi)) - 2.0 * acc_ob_y * pos_ob_x * s *
       yp_dot * std::sin(epsi)) + 2.0 * pos_ob_x * vel_ob_x * vel_ob_y * xp_dot *
      std::cos(epsi)) - 2.0 * pos_ob_y * vel_ob_x * vel_ob_y * yp_dot * std::cos
     (epsi)) / (2.0 * pow(1.0 - feb_a * feb_a / ((geb_a * geb_a + heb_a *
    heb_a) * (ieb_a * ieb_a + jeb_a * jeb_a)), 1.5) * pow(keb_a * keb_a
    + leb_a * leb_a, 1.5) * pow(meb_a * meb_a + neb_a * neb_a, 1.5))) +
              3.0 * (2.0 * ey - 2.0 * pos_ob_y) * (((((((((((ey * vel_ob_x +
    pos_ob_x * vel_ob_y) - pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot *
    std::cos(epsi)) + pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot *
    std::cos(epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi))
    - pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi))
    + s * xp_dot * std::sin(epsi)) *
              (((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
    pos_ob_y * pow(vel_ob_x, 3.0) - pos_ob_x * pow(vel_ob_y, 3.0))
    - ey * pow(vel_ob_x, 3.0)) + s * pow(vel_ob_y, 3.0)) -
    acc_ob_x * (pos_ob_x * pos_ob_x) * vel_ob_y) + acc_ob_y * (pos_ob_x *
    pos_ob_x) * vel_ob_x) - acc_ob_x * (pos_ob_y * pos_ob_y) * vel_ob_y) +
    acc_ob_y * (pos_ob_y * pos_ob_y) * vel_ob_x) - acc_ob_x * (s * s) * vel_ob_y)
    + acc_ob_y * (s * s) * vel_ob_x) - ey * vel_ob_x * (vel_ob_y * vel_ob_y)) -
    ey * vel_ob_x * (xp_dot * xp_dot)) - ey * vel_ob_x * (yp_dot * yp_dot)) -
    pos_ob_x * (vel_ob_x * vel_ob_x) * vel_ob_y) + pos_ob_y * vel_ob_x *
    (vel_ob_y * vel_ob_y)) - pos_ob_x * vel_ob_y * (xp_dot * xp_dot)) + pos_ob_y
    * vel_ob_x * (xp_dot * xp_dot)) + s * (vel_ob_x * vel_ob_x) * vel_ob_y) -
    pos_ob_x * vel_ob_y * (yp_dot * yp_dot)) + pos_ob_y * vel_ob_x * (yp_dot *
    yp_dot)) + s * vel_ob_y * (xp_dot * xp_dot)) + s * vel_ob_y * (yp_dot *
    yp_dot)) - acc_ob_x * (ey * ey) * vel_ob_y) + acc_ob_y * (ey * ey) *
    vel_ob_x) + 2.0 * acc_ob_x * ey * pos_ob_y * vel_ob_y) - 2.0 * acc_ob_y * ey
    * pos_ob_y * vel_ob_x) + 2.0 * acc_ob_x * pos_ob_x * s * vel_ob_y) - 2.0 *
    acc_ob_y * pos_ob_x * s * vel_ob_x) - acc_ob_y * (ey * ey) * xp_dot * std::
    cos(epsi)) + acc_ob_x * (ey * ey) * yp_dot * std::cos(epsi)) - acc_ob_y *
    (pos_ob_x * pos_ob_x) * xp_dot * std::cos(epsi)) - acc_ob_y * (pos_ob_y *
    pos_ob_y) * xp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_x * pos_ob_x) *
    yp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_y * pos_ob_y) * yp_dot * std::
    cos(epsi)) - acc_ob_y * (s * s) * xp_dot * std::cos(epsi)) + acc_ob_x * (s *
    s) * yp_dot * std::cos(epsi)) + acc_ob_x * (ey * ey) * xp_dot * std::sin
    (epsi)) + acc_ob_y * (ey * ey) * yp_dot * std::sin(epsi)) + 2.0 * ey *
    (vel_ob_x * vel_ob_x) * xp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_x *
    pos_ob_x) * xp_dot * std::sin(epsi)) + acc_ob_x * (pos_ob_y * pos_ob_y) *
    xp_dot * std::sin(epsi)) + acc_ob_y * (pos_ob_x * pos_ob_x) * yp_dot * std::
    sin(epsi)) + acc_ob_y * (pos_ob_y * pos_ob_y) * yp_dot * std::sin(epsi)) +
    acc_ob_x * (s * s) * xp_dot * std::sin(epsi)) + acc_ob_y * (s * s) * yp_dot *
    std::sin(epsi)) - 2.0 * pos_ob_y * (vel_ob_x * vel_ob_x) * xp_dot * std::cos
    (epsi)) + 2.0 * pos_ob_x * (vel_ob_y * vel_ob_y) * yp_dot * std::cos(epsi))
    - 2.0 * s * (vel_ob_y * vel_ob_y) * yp_dot * std::cos(epsi)) - 2.0 * ey *
    (vel_ob_x * vel_ob_x) * yp_dot * std::sin(epsi)) + 2.0 * pos_ob_x *
    (vel_ob_y * vel_ob_y) * xp_dot * std::sin(epsi)) + 2.0 * pos_ob_y *
    (vel_ob_x * vel_ob_x) * yp_dot * std::sin(epsi)) - 2.0 * s * (vel_ob_y *
    vel_ob_y) * xp_dot * std::sin(epsi)) - 2.0 * s * vel_ob_x * vel_ob_y *
    xp_dot * std::cos(epsi)) + 2.0 * ey * vel_ob_x * vel_ob_y * xp_dot * std::
    sin(epsi)) - 2.0 * pos_ob_y * vel_ob_x * vel_ob_y * xp_dot * std::sin(epsi))
    - 2.0 * pos_ob_x * vel_ob_x * vel_ob_y * yp_dot * std::sin(epsi)) + 2.0 * s *
    vel_ob_x * vel_ob_y * yp_dot * std::sin(epsi)) + 2.0 * acc_ob_y * ey *
    pos_ob_y * xp_dot * std::cos(epsi)) - 2.0 * acc_ob_x * ey * pos_ob_y *
                        yp_dot * std::cos(epsi)) + 2.0 * acc_ob_y * pos_ob_x * s
                       * xp_dot * std::cos(epsi)) - 2.0 * acc_ob_x * pos_ob_x *
                      s * yp_dot * std::cos(epsi)) - 2.0 * acc_ob_x * ey *
                     pos_ob_y * xp_dot * std::sin(epsi)) - 2.0 * acc_ob_y * ey *
                    pos_ob_y * yp_dot * std::sin(epsi)) + 2.0 * ey * vel_ob_x *
                   vel_ob_y * yp_dot * std::cos(epsi)) - 2.0 * acc_ob_x *
                  pos_ob_x * s * xp_dot * std::sin(epsi)) - 2.0 * acc_ob_y *
                 pos_ob_x * s * yp_dot * std::sin(epsi)) + 2.0 * pos_ob_x *
                vel_ob_x * vel_ob_y * xp_dot * std::cos(epsi)) - 2.0 * pos_ob_y *
               vel_ob_x * vel_ob_y * yp_dot * std::cos(epsi)) / (2.0 * std::sqrt
    (1.0 - oeb_a * oeb_a / ((peb_a * peb_a + qeb_a * qeb_a) * (reb_a * reb_a +
    seb_a * seb_a))) * pow(teb_a * teb_a + ueb_a * ueb_a, 2.5) *
    pow(veb_a * veb_a + web_a * web_a, 1.5))) - (psi_dot - psi_dot_com) *
             ((((((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y *
    vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y *
                       xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos
                      (epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::
                    sin(epsi)) - pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y *
                  yp_dot * std::sin(epsi)) + s * xp_dot * std::sin(epsi)) *
                (((((((((((((((((((((((((((((((((((((((acc_ob_x * (ey * ey) *
    xp_dot * std::cos(epsi) + acc_ob_y * (ey * ey) * yp_dot * std::cos(epsi)) +
    acc_ob_x * (pos_ob_x * pos_ob_x) * xp_dot * std::cos(epsi)) + acc_ob_x *
    (pos_ob_y * pos_ob_y) * xp_dot * std::cos(epsi)) + acc_ob_y * (pos_ob_x *
    pos_ob_x) * yp_dot * std::cos(epsi)) + acc_ob_y * (pos_ob_y * pos_ob_y) *
    yp_dot * std::cos(epsi)) + acc_ob_x * (s * s) * xp_dot * std::cos(epsi)) +
    acc_ob_y * (s * s) * yp_dot * std::cos(epsi)) + acc_ob_y * (ey * ey) *
    xp_dot * std::sin(epsi)) - acc_ob_x * (ey * ey) * yp_dot * std::sin(epsi)) -
    2.0 * ey * (vel_ob_x * vel_ob_x) * yp_dot * std::cos(epsi)) + acc_ob_y *
    (pos_ob_x * pos_ob_x) * xp_dot * std::sin(epsi)) + acc_ob_y * (pos_ob_y *
    pos_ob_y) * xp_dot * std::sin(epsi)) - acc_ob_x * (pos_ob_x * pos_ob_x) *
    yp_dot * std::sin(epsi)) - acc_ob_x * (pos_ob_y * pos_ob_y) * yp_dot * std::
    sin(epsi)) + acc_ob_y * (s * s) * xp_dot * std::sin(epsi)) - acc_ob_x * (s *
    s) * yp_dot * std::sin(epsi)) + 2.0 * pos_ob_x * (vel_ob_y * vel_ob_y) *
    xp_dot * std::cos(epsi)) + 2.0 * pos_ob_y * (vel_ob_x * vel_ob_x) * yp_dot *
    std::cos(epsi)) - 2.0 * s * (vel_ob_y * vel_ob_y) * xp_dot * std::cos(epsi))
    - 2.0 * ey * (vel_ob_x * vel_ob_x) * xp_dot * std::sin(epsi)) + 2.0 *
    pos_ob_y * (vel_ob_x * vel_ob_x) * xp_dot * std::sin(epsi)) - 2.0 * pos_ob_x
    * (vel_ob_y * vel_ob_y) * yp_dot * std::sin(epsi)) + 2.0 * s * (vel_ob_y *
    vel_ob_y) * yp_dot * std::sin(epsi)) + 2.0 * s * vel_ob_x * vel_ob_y *
    yp_dot * std::cos(epsi)) - 2.0 * ey * vel_ob_x * vel_ob_y * yp_dot * std::
    sin(epsi)) - 2.0 * pos_ob_x * vel_ob_x * vel_ob_y * xp_dot * std::sin(epsi))
    + 2.0 * pos_ob_y * vel_ob_x * vel_ob_y * yp_dot * std::sin(epsi)) + 2.0 * s *
    vel_ob_x * vel_ob_y * xp_dot * std::sin(epsi)) - 2.0 * acc_ob_x * ey *
    pos_ob_y * xp_dot * std::cos(epsi)) - 2.0 * acc_ob_y * ey * pos_ob_y *
    yp_dot * std::cos(epsi)) - 2.0 * acc_ob_x * pos_ob_x * s * xp_dot * std::cos
    (epsi)) - 2.0 * acc_ob_y * pos_ob_x * s * yp_dot * std::cos(epsi)) - 2.0 *
                       acc_ob_y * ey * pos_ob_y * xp_dot * std::sin(epsi)) + 2.0
                      * acc_ob_x * ey * pos_ob_y * yp_dot * std::sin(epsi)) +
                     2.0 * ey * vel_ob_x * vel_ob_y * xp_dot * std::cos(epsi)) -
                    2.0 * acc_ob_y * pos_ob_x * s * xp_dot * std::sin(epsi)) +
                   2.0 * acc_ob_x * pos_ob_x * s * yp_dot * std::sin(epsi)) -
                  2.0 * pos_ob_y * vel_ob_x * vel_ob_y * xp_dot * std::cos(epsi))
                 - 2.0 * pos_ob_x * vel_ob_x * vel_ob_y * yp_dot * std::cos(epsi))
                / (std::sqrt(1.0 - xeb_a * xeb_a / ((yeb_a * yeb_a + afb_a *
    afb_a) * (bfb_a * bfb_a + cfb_a * cfb_a))) * pow(dfb_a * dfb_a +
    efb_a * efb_a, 1.5) * pow(ffb_a * ffb_a + gfb_a * gfb_a, 1.5)) +
                (((((((ey * yp_dot * std::cos(epsi) - pos_ob_x * xp_dot * std::
                       cos(epsi)) - pos_ob_y * yp_dot * std::cos(epsi)) + s *
                     xp_dot * std::cos(epsi)) + ey * xp_dot * std::sin(epsi)) -
                   pos_ob_y * xp_dot * std::sin(epsi)) + pos_ob_x * yp_dot * std::
                  sin(epsi)) - s * yp_dot * std::sin(epsi)) *
                (((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
    ((pos_ob_y * pow(vel_ob_x, 3.0) - pos_ob_x * pow(vel_ob_y,
    3.0)) - ey * pow(vel_ob_x, 3.0)) + s * pow(vel_ob_y, 3.0)) -
    acc_ob_x * (pos_ob_x * pos_ob_x) * vel_ob_y) + acc_ob_y * (pos_ob_x *
    pos_ob_x) * vel_ob_x) - acc_ob_x * (pos_ob_y * pos_ob_y) * vel_ob_y) +
    acc_ob_y * (pos_ob_y * pos_ob_y) * vel_ob_x) - acc_ob_x * (s * s) * vel_ob_y)
    + acc_ob_y * (s * s) * vel_ob_x) - ey * vel_ob_x * (vel_ob_y * vel_ob_y)) -
    ey * vel_ob_x * (xp_dot * xp_dot)) - ey * vel_ob_x * (yp_dot * yp_dot)) -
    pos_ob_x * (vel_ob_x * vel_ob_x) * vel_ob_y) + pos_ob_y * vel_ob_x *
    (vel_ob_y * vel_ob_y)) - pos_ob_x * vel_ob_y * (xp_dot * xp_dot)) + pos_ob_y
    * vel_ob_x * (xp_dot * xp_dot)) + s * (vel_ob_x * vel_ob_x) * vel_ob_y) -
    pos_ob_x * vel_ob_y * (yp_dot * yp_dot)) + pos_ob_y * vel_ob_x * (yp_dot *
    yp_dot)) + s * vel_ob_y * (xp_dot * xp_dot)) + s * vel_ob_y * (yp_dot *
    yp_dot)) - acc_ob_x * (ey * ey) * vel_ob_y) + acc_ob_y * (ey * ey) *
    vel_ob_x) + 2.0 * acc_ob_x * ey * pos_ob_y * vel_ob_y) - 2.0 * acc_ob_y * ey
    * pos_ob_y * vel_ob_x) + 2.0 * acc_ob_x * pos_ob_x * s * vel_ob_y) - 2.0 *
    acc_ob_y * pos_ob_x * s * vel_ob_x) - acc_ob_y * (ey * ey) * xp_dot * std::
    cos(epsi)) + acc_ob_x * (ey * ey) * yp_dot * std::cos(epsi)) - acc_ob_y *
    (pos_ob_x * pos_ob_x) * xp_dot * std::cos(epsi)) - acc_ob_y * (pos_ob_y *
    pos_ob_y) * xp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_x * pos_ob_x) *
    yp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_y * pos_ob_y) * yp_dot * std::
    cos(epsi)) - acc_ob_y * (s * s) * xp_dot * std::cos(epsi)) + acc_ob_x * (s *
    s) * yp_dot * std::cos(epsi)) + acc_ob_x * (ey * ey) * xp_dot * std::sin
    (epsi)) + acc_ob_y * (ey * ey) * yp_dot * std::sin(epsi)) + 2.0 * ey *
    (vel_ob_x * vel_ob_x) * xp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_x *
    pos_ob_x) * xp_dot * std::sin(epsi)) + acc_ob_x * (pos_ob_y * pos_ob_y) *
    xp_dot * std::sin(epsi)) + acc_ob_y * (pos_ob_x * pos_ob_x) * yp_dot * std::
    sin(epsi)) + acc_ob_y * (pos_ob_y * pos_ob_y) * yp_dot * std::sin(epsi)) +
    acc_ob_x * (s * s) * xp_dot * std::sin(epsi)) + acc_ob_y * (s * s) * yp_dot *
    std::sin(epsi)) - 2.0 * pos_ob_y * (vel_ob_x * vel_ob_x) * xp_dot * std::cos
    (epsi)) + 2.0 * pos_ob_x * (vel_ob_y * vel_ob_y) * yp_dot * std::cos(epsi))
    - 2.0 * s * (vel_ob_y * vel_ob_y) * yp_dot * std::cos(epsi)) - 2.0 * ey *
    (vel_ob_x * vel_ob_x) * yp_dot * std::sin(epsi)) + 2.0 * pos_ob_x *
    (vel_ob_y * vel_ob_y) * xp_dot * std::sin(epsi)) + 2.0 * pos_ob_y *
    (vel_ob_x * vel_ob_x) * yp_dot * std::sin(epsi)) - 2.0 * s * (vel_ob_y *
    vel_ob_y) * xp_dot * std::sin(epsi)) - 2.0 * s * vel_ob_x * vel_ob_y *
    xp_dot * std::cos(epsi)) + 2.0 * ey * vel_ob_x * vel_ob_y * xp_dot * std::
    sin(epsi)) - 2.0 * pos_ob_y * vel_ob_x * vel_ob_y * xp_dot * std::sin(epsi))
    - 2.0 * pos_ob_x * vel_ob_x * vel_ob_y * yp_dot * std::sin(epsi)) + 2.0 * s *
    vel_ob_x * vel_ob_y * yp_dot * std::sin(epsi)) + 2.0 * acc_ob_y * ey *
    pos_ob_y * xp_dot * std::cos(epsi)) - 2.0 * acc_ob_x * ey * pos_ob_y *
    yp_dot * std::cos(epsi)) + 2.0 * acc_ob_y * pos_ob_x * s * xp_dot * std::cos
    (epsi)) - 2.0 * acc_ob_x * pos_ob_x * s * yp_dot * std::cos(epsi)) - 2.0 *
                       acc_ob_x * ey * pos_ob_y * xp_dot * std::sin(epsi)) - 2.0
                      * acc_ob_y * ey * pos_ob_y * yp_dot * std::sin(epsi)) +
                     2.0 * ey * vel_ob_x * vel_ob_y * yp_dot * std::cos(epsi)) -
                    2.0 * acc_ob_x * pos_ob_x * s * xp_dot * std::sin(epsi)) -
                   2.0 * acc_ob_y * pos_ob_x * s * yp_dot * std::sin(epsi)) +
                  2.0 * pos_ob_x * vel_ob_x * vel_ob_y * xp_dot * std::cos(epsi))
                 - 2.0 * pos_ob_y * vel_ob_x * vel_ob_y * yp_dot * std::cos(epsi))
                / (std::sqrt(1.0 - hfb_a * hfb_a / ((ifb_a * ifb_a + jfb_a *
    jfb_a) * (kfb_a * kfb_a + lfb_a * lfb_a))) * pow(mfb_a * mfb_a +
    nfb_a * nfb_a, 1.5) * pow(ofb_a * ofb_a + pfb_a * pfb_a, 1.5))) -
               3.0 * (((2.0 * vel_ob_x * yp_dot * std::cos(epsi) - 2.0 *
                        vel_ob_y * xp_dot * std::cos(epsi)) + 2.0 * vel_ob_x *
                       xp_dot * std::sin(epsi)) + 2.0 * vel_ob_y * yp_dot * std::
                      sin(epsi)) * (((((((((((ey * vel_ob_x + pos_ob_x *
    vel_ob_y) - pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos
    (epsi)) + pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos
    (epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) -
    pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) +
    s * xp_dot * std::sin(epsi)) *
               ((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
    (pos_ob_y * pow(vel_ob_x, 3.0) - pos_ob_x * pow(vel_ob_y,
    3.0)) - ey * pow(vel_ob_x, 3.0)) + s * pow(vel_ob_y, 3.0)) -
    acc_ob_x * (pos_ob_x * pos_ob_x) * vel_ob_y) + acc_ob_y * (pos_ob_x *
    pos_ob_x) * vel_ob_x) - acc_ob_x * (pos_ob_y * pos_ob_y) * vel_ob_y) +
    acc_ob_y * (pos_ob_y * pos_ob_y) * vel_ob_x) - acc_ob_x * (s * s) * vel_ob_y)
    + acc_ob_y * (s * s) * vel_ob_x) - ey * vel_ob_x * (vel_ob_y * vel_ob_y)) -
    ey * vel_ob_x * (xp_dot * xp_dot)) - ey * vel_ob_x * (yp_dot * yp_dot)) -
    pos_ob_x * (vel_ob_x * vel_ob_x) * vel_ob_y) + pos_ob_y * vel_ob_x *
    (vel_ob_y * vel_ob_y)) - pos_ob_x * vel_ob_y * (xp_dot * xp_dot)) + pos_ob_y
    * vel_ob_x * (xp_dot * xp_dot)) + s * (vel_ob_x * vel_ob_x) * vel_ob_y) -
    pos_ob_x * vel_ob_y * (yp_dot * yp_dot)) + pos_ob_y * vel_ob_x * (yp_dot *
    yp_dot)) + s * vel_ob_y * (xp_dot * xp_dot)) + s * vel_ob_y * (yp_dot *
    yp_dot)) - acc_ob_x * (ey * ey) * vel_ob_y) + acc_ob_y * (ey * ey) *
    vel_ob_x) + 2.0 * acc_ob_x * ey * pos_ob_y * vel_ob_y) - 2.0 * acc_ob_y * ey
    * pos_ob_y * vel_ob_x) + 2.0 * acc_ob_x * pos_ob_x * s * vel_ob_y) - 2.0 *
    acc_ob_y * pos_ob_x * s * vel_ob_x) - acc_ob_y * (ey * ey) * xp_dot * std::
    cos(epsi)) + acc_ob_x * (ey * ey) * yp_dot * std::cos(epsi)) - acc_ob_y *
    (pos_ob_x * pos_ob_x) * xp_dot * std::cos(epsi)) - acc_ob_y * (pos_ob_y *
    pos_ob_y) * xp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_x * pos_ob_x) *
    yp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_y * pos_ob_y) * yp_dot * std::
    cos(epsi)) - acc_ob_y * (s * s) * xp_dot * std::cos(epsi)) + acc_ob_x * (s *
    s) * yp_dot * std::cos(epsi)) + acc_ob_x * (ey * ey) * xp_dot * std::sin
    (epsi)) + acc_ob_y * (ey * ey) * yp_dot * std::sin(epsi)) + 2.0 * ey *
    (vel_ob_x * vel_ob_x) * xp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_x *
    pos_ob_x) * xp_dot * std::sin(epsi)) + acc_ob_x * (pos_ob_y * pos_ob_y) *
    xp_dot * std::sin(epsi)) + acc_ob_y * (pos_ob_x * pos_ob_x) * yp_dot * std::
    sin(epsi)) + acc_ob_y * (pos_ob_y * pos_ob_y) * yp_dot * std::sin(epsi)) +
    acc_ob_x * (s * s) * xp_dot * std::sin(epsi)) + acc_ob_y * (s * s) * yp_dot *
    std::sin(epsi)) - 2.0 * pos_ob_y * (vel_ob_x * vel_ob_x) * xp_dot * std::cos
    (epsi)) + 2.0 * pos_ob_x * (vel_ob_y * vel_ob_y) * yp_dot * std::cos(epsi))
    - 2.0 * s * (vel_ob_y * vel_ob_y) * yp_dot * std::cos(epsi)) - 2.0 * ey *
    (vel_ob_x * vel_ob_x) * yp_dot * std::sin(epsi)) + 2.0 * pos_ob_x *
    (vel_ob_y * vel_ob_y) * xp_dot * std::sin(epsi)) + 2.0 * pos_ob_y *
    (vel_ob_x * vel_ob_x) * yp_dot * std::sin(epsi)) - 2.0 * s * (vel_ob_y *
    vel_ob_y) * xp_dot * std::sin(epsi)) - 2.0 * s * vel_ob_x * vel_ob_y *
    xp_dot * std::cos(epsi)) + 2.0 * ey * vel_ob_x * vel_ob_y * xp_dot * std::
    sin(epsi)) - 2.0 * pos_ob_y * vel_ob_x * vel_ob_y * xp_dot * std::sin(epsi))
    - 2.0 * pos_ob_x * vel_ob_x * vel_ob_y * yp_dot * std::sin(epsi)) + 2.0 * s *
    vel_ob_x * vel_ob_y * yp_dot * std::sin(epsi)) + 2.0 * acc_ob_y * ey *
    pos_ob_y * xp_dot * std::cos(epsi)) - 2.0 * acc_ob_x * ey * pos_ob_y *
    yp_dot * std::cos(epsi)) + 2.0 * acc_ob_y * pos_ob_x * s * xp_dot * std::cos
                        (epsi)) - 2.0 * acc_ob_x * pos_ob_x * s * yp_dot * std::
                       cos(epsi)) - 2.0 * acc_ob_x * ey * pos_ob_y * xp_dot *
                      std::sin(epsi)) - 2.0 * acc_ob_y * ey * pos_ob_y * yp_dot *
                     std::sin(epsi)) + 2.0 * ey * vel_ob_x * vel_ob_y * yp_dot *
                    std::cos(epsi)) - 2.0 * acc_ob_x * pos_ob_x * s * xp_dot *
                   std::sin(epsi)) - 2.0 * acc_ob_y * pos_ob_x * s * yp_dot *
                  std::sin(epsi)) + 2.0 * pos_ob_x * vel_ob_x * vel_ob_y *
                 xp_dot * std::cos(epsi)) - 2.0 * pos_ob_y * vel_ob_x * vel_ob_y
                * yp_dot * std::cos(epsi)) / (2.0 * std::sqrt(1.0 - qfb_a *
    qfb_a / ((rfb_a * rfb_a + sfb_a * sfb_a) * (tfb_a * tfb_a + ufb_a * ufb_a)))
    * pow(vfb_a * vfb_a + wfb_a * wfb_a, 1.5) * pow(xfb_a *
    xfb_a + yfb_a * yfb_a, 2.5))) + (2.0 * ((pos_ob_x - s) * ((vel_ob_x - xp_dot
    * std::cos(epsi)) + yp_dot * std::sin(epsi)) + (ey - pos_ob_y) * ((yp_dot *
    std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi))) * ((ey - pos_ob_y) *
    (xp_dot * std::cos(epsi) - yp_dot * std::sin(epsi)) + (pos_ob_x - s) *
    (yp_dot * std::cos(epsi) + xp_dot * std::sin(epsi))) / ((agb_a * agb_a +
    bgb_a * bgb_a) * (cgb_a * cgb_a + dgb_a * dgb_a)) - egb_a * egb_a * (((2.0 *
    vel_ob_x * yp_dot * std::cos(epsi) - 2.0 * vel_ob_y * xp_dot * std::cos(epsi))
    + 2.0 * vel_ob_x * xp_dot * std::sin(epsi)) + 2.0 * vel_ob_y * yp_dot * std::
    sin(epsi)) / ((fgb_a * fgb_a + ggb_a * ggb_a) * (hgb_a * hgb_a))) *
              (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y *
                        vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi))
                     + pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot *
                    std::cos(epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot
                  * std::sin(epsi)) - pos_ob_x * xp_dot * std::sin(epsi)) -
                pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot * std::sin(epsi))
              * (((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
    ((pos_ob_y * pow(vel_ob_x, 3.0) - pos_ob_x * pow(vel_ob_y,
    3.0)) - ey * pow(vel_ob_x, 3.0)) + s * pow(vel_ob_y, 3.0)) -
    acc_ob_x * (pos_ob_x * pos_ob_x) * vel_ob_y) + acc_ob_y * (pos_ob_x *
    pos_ob_x) * vel_ob_x) - acc_ob_x * (pos_ob_y * pos_ob_y) * vel_ob_y) +
    acc_ob_y * (pos_ob_y * pos_ob_y) * vel_ob_x) - acc_ob_x * (s * s) * vel_ob_y)
    + acc_ob_y * (s * s) * vel_ob_x) - ey * vel_ob_x * (vel_ob_y * vel_ob_y)) -
    ey * vel_ob_x * (xp_dot * xp_dot)) - ey * vel_ob_x * (yp_dot * yp_dot)) -
    pos_ob_x * (vel_ob_x * vel_ob_x) * vel_ob_y) + pos_ob_y * vel_ob_x *
    (vel_ob_y * vel_ob_y)) - pos_ob_x * vel_ob_y * (xp_dot * xp_dot)) + pos_ob_y
    * vel_ob_x * (xp_dot * xp_dot)) + s * (vel_ob_x * vel_ob_x) * vel_ob_y) -
    pos_ob_x * vel_ob_y * (yp_dot * yp_dot)) + pos_ob_y * vel_ob_x * (yp_dot *
    yp_dot)) + s * vel_ob_y * (xp_dot * xp_dot)) + s * vel_ob_y * (yp_dot *
    yp_dot)) - acc_ob_x * (ey * ey) * vel_ob_y) + acc_ob_y * (ey * ey) *
    vel_ob_x) + 2.0 * acc_ob_x * ey * pos_ob_y * vel_ob_y) - 2.0 * acc_ob_y * ey
    * pos_ob_y * vel_ob_x) + 2.0 * acc_ob_x * pos_ob_x * s * vel_ob_y) - 2.0 *
    acc_ob_y * pos_ob_x * s * vel_ob_x) - acc_ob_y * (ey * ey) * xp_dot * std::
    cos(epsi)) + acc_ob_x * (ey * ey) * yp_dot * std::cos(epsi)) - acc_ob_y *
    (pos_ob_x * pos_ob_x) * xp_dot * std::cos(epsi)) - acc_ob_y * (pos_ob_y *
    pos_ob_y) * xp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_x * pos_ob_x) *
    yp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_y * pos_ob_y) * yp_dot * std::
    cos(epsi)) - acc_ob_y * (s * s) * xp_dot * std::cos(epsi)) + acc_ob_x * (s *
    s) * yp_dot * std::cos(epsi)) + acc_ob_x * (ey * ey) * xp_dot * std::sin
    (epsi)) + acc_ob_y * (ey * ey) * yp_dot * std::sin(epsi)) + 2.0 * ey *
    (vel_ob_x * vel_ob_x) * xp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_x *
    pos_ob_x) * xp_dot * std::sin(epsi)) + acc_ob_x * (pos_ob_y * pos_ob_y) *
    xp_dot * std::sin(epsi)) + acc_ob_y * (pos_ob_x * pos_ob_x) * yp_dot * std::
    sin(epsi)) + acc_ob_y * (pos_ob_y * pos_ob_y) * yp_dot * std::sin(epsi)) +
    acc_ob_x * (s * s) * xp_dot * std::sin(epsi)) + acc_ob_y * (s * s) * yp_dot *
    std::sin(epsi)) - 2.0 * pos_ob_y * (vel_ob_x * vel_ob_x) * xp_dot * std::cos
    (epsi)) + 2.0 * pos_ob_x * (vel_ob_y * vel_ob_y) * yp_dot * std::cos(epsi))
    - 2.0 * s * (vel_ob_y * vel_ob_y) * yp_dot * std::cos(epsi)) - 2.0 * ey *
    (vel_ob_x * vel_ob_x) * yp_dot * std::sin(epsi)) + 2.0 * pos_ob_x *
    (vel_ob_y * vel_ob_y) * xp_dot * std::sin(epsi)) + 2.0 * pos_ob_y *
    (vel_ob_x * vel_ob_x) * yp_dot * std::sin(epsi)) - 2.0 * s * (vel_ob_y *
    vel_ob_y) * xp_dot * std::sin(epsi)) - 2.0 * s * vel_ob_x * vel_ob_y *
    xp_dot * std::cos(epsi)) + 2.0 * ey * vel_ob_x * vel_ob_y * xp_dot * std::
    sin(epsi)) - 2.0 * pos_ob_y * vel_ob_x * vel_ob_y * xp_dot * std::sin(epsi))
    - 2.0 * pos_ob_x * vel_ob_x * vel_ob_y * yp_dot * std::sin(epsi)) + 2.0 * s *
    vel_ob_x * vel_ob_y * yp_dot * std::sin(epsi)) + 2.0 * acc_ob_y * ey *
    pos_ob_y * xp_dot * std::cos(epsi)) - 2.0 * acc_ob_x * ey * pos_ob_y *
    yp_dot * std::cos(epsi)) + 2.0 * acc_ob_y * pos_ob_x * s * xp_dot * std::cos
    (epsi)) - 2.0 * acc_ob_x * pos_ob_x * s * yp_dot * std::cos(epsi)) - 2.0 *
                       acc_ob_x * ey * pos_ob_y * xp_dot * std::sin(epsi)) - 2.0
                      * acc_ob_y * ey * pos_ob_y * yp_dot * std::sin(epsi)) +
                     2.0 * ey * vel_ob_x * vel_ob_y * yp_dot * std::cos(epsi)) -
                    2.0 * acc_ob_x * pos_ob_x * s * xp_dot * std::sin(epsi)) -
                   2.0 * acc_ob_y * pos_ob_x * s * yp_dot * std::sin(epsi)) +
                  2.0 * pos_ob_x * vel_ob_x * vel_ob_y * xp_dot * std::cos(epsi))
                 - 2.0 * pos_ob_y * vel_ob_x * vel_ob_y * yp_dot * std::cos(epsi))
              / (2.0 * pow(1.0 - igb_a * igb_a / ((jgb_a * jgb_a + kgb_a
    * kgb_a) * (lgb_a * lgb_a + mgb_a * mgb_a)), 1.5) * pow(ngb_a *
    ngb_a + ogb_a * ogb_a, 1.5) * pow(pgb_a * pgb_a + qgb_a * qgb_a, 1.5))))
            - (xp_dot * std::cos(epsi) - yp_dot * std::sin(epsi)) * ((((((Ds *
    vel_ob_x / (std::sqrt(1.0 - Ds * Ds / (rgb_a * rgb_a + sgb_a * sgb_a)) *
                pow(tgb_a * tgb_a + ugb_a * ugb_a, 1.5)) + ((yp_dot *
    std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi)) *
    (((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((pos_ob_y *
    pow(vel_ob_x, 3.0) - pos_ob_x * pow(vel_ob_y, 3.0)) - ey *
    pow(vel_ob_x, 3.0)) + s * pow(vel_ob_y, 3.0)) - acc_ob_x *
    (pos_ob_x * pos_ob_x) * vel_ob_y) + acc_ob_y * (pos_ob_x * pos_ob_x) *
    vel_ob_x) - acc_ob_x * (pos_ob_y * pos_ob_y) * vel_ob_y) + acc_ob_y *
    (pos_ob_y * pos_ob_y) * vel_ob_x) - acc_ob_x * (s * s) * vel_ob_y) +
    acc_ob_y * (s * s) * vel_ob_x) - ey * vel_ob_x * (vel_ob_y * vel_ob_y)) - ey
    * vel_ob_x * (xp_dot * xp_dot)) - ey * vel_ob_x * (yp_dot * yp_dot)) -
    pos_ob_x * (vel_ob_x * vel_ob_x) * vel_ob_y) + pos_ob_y * vel_ob_x *
    (vel_ob_y * vel_ob_y)) - pos_ob_x * vel_ob_y * (xp_dot * xp_dot)) + pos_ob_y
    * vel_ob_x * (xp_dot * xp_dot)) + s * (vel_ob_x * vel_ob_x) * vel_ob_y) -
    pos_ob_x * vel_ob_y * (yp_dot * yp_dot)) + pos_ob_y * vel_ob_x * (yp_dot *
    yp_dot)) + s * vel_ob_y * (xp_dot * xp_dot)) + s * vel_ob_y * (yp_dot *
    yp_dot)) - acc_ob_x * (ey * ey) * vel_ob_y) + acc_ob_y * (ey * ey) *
    vel_ob_x) + 2.0 * acc_ob_x * ey * pos_ob_y * vel_ob_y) - 2.0 * acc_ob_y * ey
    * pos_ob_y * vel_ob_x) + 2.0 * acc_ob_x * pos_ob_x * s * vel_ob_y) - 2.0 *
    acc_ob_y * pos_ob_x * s * vel_ob_x) - acc_ob_y * (ey * ey) * xp_dot * std::
    cos(epsi)) + acc_ob_x * (ey * ey) * yp_dot * std::cos(epsi)) - acc_ob_y *
    (pos_ob_x * pos_ob_x) * xp_dot * std::cos(epsi)) - acc_ob_y * (pos_ob_y *
    pos_ob_y) * xp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_x * pos_ob_x) *
    yp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_y * pos_ob_y) * yp_dot * std::
    cos(epsi)) - acc_ob_y * (s * s) * xp_dot * std::cos(epsi)) + acc_ob_x * (s *
    s) * yp_dot * std::cos(epsi)) + acc_ob_x * (ey * ey) * xp_dot * std::sin
    (epsi)) + acc_ob_y * (ey * ey) * yp_dot * std::sin(epsi)) + 2.0 * ey *
    (vel_ob_x * vel_ob_x) * xp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_x *
    pos_ob_x) * xp_dot * std::sin(epsi)) + acc_ob_x * (pos_ob_y * pos_ob_y) *
    xp_dot * std::sin(epsi)) + acc_ob_y * (pos_ob_x * pos_ob_x) * yp_dot * std::
    sin(epsi)) + acc_ob_y * (pos_ob_y * pos_ob_y) * yp_dot * std::sin(epsi)) +
    acc_ob_x * (s * s) * xp_dot * std::sin(epsi)) + acc_ob_y * (s * s) * yp_dot *
    std::sin(epsi)) - 2.0 * pos_ob_y * (vel_ob_x * vel_ob_x) * xp_dot * std::cos
    (epsi)) + 2.0 * pos_ob_x * (vel_ob_y * vel_ob_y) * yp_dot * std::cos(epsi))
    - 2.0 * s * (vel_ob_y * vel_ob_y) * yp_dot * std::cos(epsi)) - 2.0 * ey *
    (vel_ob_x * vel_ob_x) * yp_dot * std::sin(epsi)) + 2.0 * pos_ob_x *
    (vel_ob_y * vel_ob_y) * xp_dot * std::sin(epsi)) + 2.0 * pos_ob_y *
    (vel_ob_x * vel_ob_x) * yp_dot * std::sin(epsi)) - 2.0 * s * (vel_ob_y *
    vel_ob_y) * xp_dot * std::sin(epsi)) - 2.0 * s * vel_ob_x * vel_ob_y *
    xp_dot * std::cos(epsi)) + 2.0 * ey * vel_ob_x * vel_ob_y * xp_dot * std::
    sin(epsi)) - 2.0 * pos_ob_y * vel_ob_x * vel_ob_y * xp_dot * std::sin(epsi))
    - 2.0 * pos_ob_x * vel_ob_x * vel_ob_y * yp_dot * std::sin(epsi)) + 2.0 * s *
                vel_ob_x * vel_ob_y * yp_dot * std::sin(epsi)) + 2.0 * acc_ob_y *
               ey * pos_ob_y * xp_dot * std::cos(epsi)) - 2.0 * acc_ob_x * ey *
              pos_ob_y * yp_dot * std::cos(epsi)) + 2.0 * acc_ob_y * pos_ob_x *
             s * xp_dot * std::cos(epsi)) - 2.0 * acc_ob_x * pos_ob_x * s *
            yp_dot * std::cos(epsi)) - 2.0 * acc_ob_x * ey * pos_ob_y * xp_dot *
           std::sin(epsi)) - 2.0 * acc_ob_y * ey * pos_ob_y * yp_dot * std::sin
          (epsi)) + 2.0 * ey * vel_ob_x * vel_ob_y * yp_dot * std::cos(epsi)) -
        2.0 * acc_ob_x * pos_ob_x * s * xp_dot * std::sin(epsi)) - 2.0 *
       acc_ob_y * pos_ob_x * s * yp_dot * std::sin(epsi)) + 2.0 * pos_ob_x *
      vel_ob_x * vel_ob_y * xp_dot * std::cos(epsi)) - 2.0 * pos_ob_y * vel_ob_x
     * vel_ob_y * yp_dot * std::cos(epsi)) / (std::sqrt(1.0 - vgb_a * vgb_a /
    ((wgb_a * wgb_a + xgb_a * xgb_a) * (ygb_a * ygb_a + ahb_a * ahb_a))) *
    pow(bhb_a * bhb_a + chb_a * chb_a, 1.5) * pow(dhb_a * dhb_a
    + ehb_a * ehb_a, 1.5))) + 3.0 * Ds * (2.0 * pos_ob_x - 2.0 * s) * (((ey *
    vel_ob_y - pos_ob_x * vel_ob_x) - pos_ob_y * vel_ob_y) + s * vel_ob_x) /
    (2.0 * std::sqrt(1.0 - Ds * Ds / (fhb_a * fhb_a + ghb_a * ghb_a)) *
     pow(hhb_a * hhb_a + ihb_a * ihb_a, 2.5))) + pow(Ds, 3.0) *
    (2.0 * pos_ob_x - 2.0 * s) * (((ey * vel_ob_y - pos_ob_x * vel_ob_x) -
    pos_ob_y * vel_ob_y) + s * vel_ob_x) / (2.0 * pow(1.0 - Ds * Ds /
    (jhb_a * jhb_a + khb_a * khb_a), 1.5) * pow(lhb_a * lhb_a + mhb_a *
    mhb_a, 3.5))) + (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y *
    vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y *
    xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) + s * yp_dot *
    std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x * xp_dot * std::
                       sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s *
                     xp_dot * std::sin(epsi)) * (((((((((((((((((((vel_ob_x *
    vel_ob_x * vel_ob_y + vel_ob_y * (xp_dot * xp_dot)) + vel_ob_y * (yp_dot *
    yp_dot)) + pow(vel_ob_y, 3.0)) - 2.0 * (vel_ob_y * vel_ob_y) *
    yp_dot * std::cos(epsi)) - 2.0 * (vel_ob_y * vel_ob_y) * xp_dot * std::sin
    (epsi)) + 2.0 * acc_ob_x * pos_ob_x * vel_ob_y) - 2.0 * acc_ob_y * pos_ob_x *
    vel_ob_x) - 2.0 * acc_ob_x * s * vel_ob_y) + 2.0 * acc_ob_y * s * vel_ob_x)
    + 2.0 * acc_ob_y * pos_ob_x * xp_dot * std::cos(epsi)) - 2.0 * acc_ob_x *
    pos_ob_x * yp_dot * std::cos(epsi)) - 2.0 * acc_ob_y * s * xp_dot * std::cos
    (epsi)) + 2.0 * acc_ob_x * s * yp_dot * std::cos(epsi)) - 2.0 * acc_ob_x *
    pos_ob_x * xp_dot * std::sin(epsi)) - 2.0 * acc_ob_y * pos_ob_x * yp_dot *
    std::sin(epsi)) + 2.0 * acc_ob_x * s * xp_dot * std::sin(epsi)) + 2.0 *
    acc_ob_y * s * yp_dot * std::sin(epsi)) - 2.0 * vel_ob_x * vel_ob_y * xp_dot
    * std::cos(epsi)) + 2.0 * vel_ob_x * vel_ob_y * yp_dot * std::sin(epsi)) /
    (std::sqrt(1.0 - nhb_a * nhb_a / ((ohb_a * ohb_a + phb_a * phb_a) * (qhb_a *
    qhb_a + rhb_a * rhb_a))) * pow(shb_a * shb_a + thb_a * thb_a, 1.5) *
     pow(uhb_a * uhb_a + vhb_a * vhb_a, 1.5))) - (2.0 * ((pos_ob_x - s) *
    ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi)) + (ey -
    pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi)))
    * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi)) / ((whb_a
    * whb_a + xhb_a * xhb_a) * (yhb_a * yhb_a + aib_a * aib_a)) - (2.0 *
    pos_ob_x - 2.0 * s) * (bib_a * bib_a) / (cib_a * cib_a * (dib_a * dib_a +
    eib_a * eib_a))) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) -
    pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) +
    pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) +
    s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x *
    xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot *
                        std::sin(epsi)) *
              (((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
    pos_ob_y * pow(vel_ob_x, 3.0) - pos_ob_x * pow(vel_ob_y, 3.0))
    - ey * pow(vel_ob_x, 3.0)) + s * pow(vel_ob_y, 3.0)) -
    acc_ob_x * (pos_ob_x * pos_ob_x) * vel_ob_y) + acc_ob_y * (pos_ob_x *
    pos_ob_x) * vel_ob_x) - acc_ob_x * (pos_ob_y * pos_ob_y) * vel_ob_y) +
    acc_ob_y * (pos_ob_y * pos_ob_y) * vel_ob_x) - acc_ob_x * (s * s) * vel_ob_y)
    + acc_ob_y * (s * s) * vel_ob_x) - ey * vel_ob_x * (vel_ob_y * vel_ob_y)) -
    ey * vel_ob_x * (xp_dot * xp_dot)) - ey * vel_ob_x * (yp_dot * yp_dot)) -
    pos_ob_x * (vel_ob_x * vel_ob_x) * vel_ob_y) + pos_ob_y * vel_ob_x *
    (vel_ob_y * vel_ob_y)) - pos_ob_x * vel_ob_y * (xp_dot * xp_dot)) + pos_ob_y
    * vel_ob_x * (xp_dot * xp_dot)) + s * (vel_ob_x * vel_ob_x) * vel_ob_y) -
    pos_ob_x * vel_ob_y * (yp_dot * yp_dot)) + pos_ob_y * vel_ob_x * (yp_dot *
    yp_dot)) + s * vel_ob_y * (xp_dot * xp_dot)) + s * vel_ob_y * (yp_dot *
    yp_dot)) - acc_ob_x * (ey * ey) * vel_ob_y) + acc_ob_y * (ey * ey) *
    vel_ob_x) + 2.0 * acc_ob_x * ey * pos_ob_y * vel_ob_y) - 2.0 * acc_ob_y * ey
    * pos_ob_y * vel_ob_x) + 2.0 * acc_ob_x * pos_ob_x * s * vel_ob_y) - 2.0 *
    acc_ob_y * pos_ob_x * s * vel_ob_x) - acc_ob_y * (ey * ey) * xp_dot * std::
    cos(epsi)) + acc_ob_x * (ey * ey) * yp_dot * std::cos(epsi)) - acc_ob_y *
    (pos_ob_x * pos_ob_x) * xp_dot * std::cos(epsi)) - acc_ob_y * (pos_ob_y *
    pos_ob_y) * xp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_x * pos_ob_x) *
    yp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_y * pos_ob_y) * yp_dot * std::
    cos(epsi)) - acc_ob_y * (s * s) * xp_dot * std::cos(epsi)) + acc_ob_x * (s *
    s) * yp_dot * std::cos(epsi)) + acc_ob_x * (ey * ey) * xp_dot * std::sin
    (epsi)) + acc_ob_y * (ey * ey) * yp_dot * std::sin(epsi)) + 2.0 * ey *
    (vel_ob_x * vel_ob_x) * xp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_x *
    pos_ob_x) * xp_dot * std::sin(epsi)) + acc_ob_x * (pos_ob_y * pos_ob_y) *
    xp_dot * std::sin(epsi)) + acc_ob_y * (pos_ob_x * pos_ob_x) * yp_dot * std::
    sin(epsi)) + acc_ob_y * (pos_ob_y * pos_ob_y) * yp_dot * std::sin(epsi)) +
    acc_ob_x * (s * s) * xp_dot * std::sin(epsi)) + acc_ob_y * (s * s) * yp_dot *
    std::sin(epsi)) - 2.0 * pos_ob_y * (vel_ob_x * vel_ob_x) * xp_dot * std::cos
    (epsi)) + 2.0 * pos_ob_x * (vel_ob_y * vel_ob_y) * yp_dot * std::cos(epsi))
    - 2.0 * s * (vel_ob_y * vel_ob_y) * yp_dot * std::cos(epsi)) - 2.0 * ey *
    (vel_ob_x * vel_ob_x) * yp_dot * std::sin(epsi)) + 2.0 * pos_ob_x *
    (vel_ob_y * vel_ob_y) * xp_dot * std::sin(epsi)) + 2.0 * pos_ob_y *
    (vel_ob_x * vel_ob_x) * yp_dot * std::sin(epsi)) - 2.0 * s * (vel_ob_y *
    vel_ob_y) * xp_dot * std::sin(epsi)) - 2.0 * s * vel_ob_x * vel_ob_y *
    xp_dot * std::cos(epsi)) + 2.0 * ey * vel_ob_x * vel_ob_y * xp_dot * std::
    sin(epsi)) - 2.0 * pos_ob_y * vel_ob_x * vel_ob_y * xp_dot * std::sin(epsi))
    - 2.0 * pos_ob_x * vel_ob_x * vel_ob_y * yp_dot * std::sin(epsi)) + 2.0 * s *
    vel_ob_x * vel_ob_y * yp_dot * std::sin(epsi)) + 2.0 * acc_ob_y * ey *
    pos_ob_y * xp_dot * std::cos(epsi)) - 2.0 * acc_ob_x * ey * pos_ob_y *
                        yp_dot * std::cos(epsi)) + 2.0 * acc_ob_y * pos_ob_x * s
                       * xp_dot * std::cos(epsi)) - 2.0 * acc_ob_x * pos_ob_x *
                      s * yp_dot * std::cos(epsi)) - 2.0 * acc_ob_x * ey *
                     pos_ob_y * xp_dot * std::sin(epsi)) - 2.0 * acc_ob_y * ey *
                    pos_ob_y * yp_dot * std::sin(epsi)) + 2.0 * ey * vel_ob_x *
                   vel_ob_y * yp_dot * std::cos(epsi)) - 2.0 * acc_ob_x *
                  pos_ob_x * s * xp_dot * std::sin(epsi)) - 2.0 * acc_ob_y *
                 pos_ob_x * s * yp_dot * std::sin(epsi)) + 2.0 * pos_ob_x *
                vel_ob_x * vel_ob_y * xp_dot * std::cos(epsi)) - 2.0 * pos_ob_y *
               vel_ob_x * vel_ob_y * yp_dot * std::cos(epsi)) / (2.0 *
    pow(1.0 - fib_a * fib_a / ((gib_a * gib_a + hib_a * hib_a) * (iib_a *
    iib_a + jib_a * jib_a)), 1.5) * pow(kib_a * kib_a + lib_a * lib_a,
    1.5) * pow(mib_a * mib_a + nib_a * nib_a, 1.5))) + 3.0 * (2.0 *
              pos_ob_x - 2.0 * s) * (((((((((((ey * vel_ob_x + pos_ob_x *
    vel_ob_y) - pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos
    (epsi)) + pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos
    (epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) -
    pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) +
              s * xp_dot * std::sin(epsi)) *
             (((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((pos_ob_y
    * pow(vel_ob_x, 3.0) - pos_ob_x * pow(vel_ob_y, 3.0)) - ey *
    pow(vel_ob_x, 3.0)) + s * pow(vel_ob_y, 3.0)) - acc_ob_x *
    (pos_ob_x * pos_ob_x) * vel_ob_y) + acc_ob_y * (pos_ob_x * pos_ob_x) *
    vel_ob_x) - acc_ob_x * (pos_ob_y * pos_ob_y) * vel_ob_y) + acc_ob_y *
    (pos_ob_y * pos_ob_y) * vel_ob_x) - acc_ob_x * (s * s) * vel_ob_y) +
    acc_ob_y * (s * s) * vel_ob_x) - ey * vel_ob_x * (vel_ob_y * vel_ob_y)) - ey
    * vel_ob_x * (xp_dot * xp_dot)) - ey * vel_ob_x * (yp_dot * yp_dot)) -
    pos_ob_x * (vel_ob_x * vel_ob_x) * vel_ob_y) + pos_ob_y * vel_ob_x *
    (vel_ob_y * vel_ob_y)) - pos_ob_x * vel_ob_y * (xp_dot * xp_dot)) + pos_ob_y
    * vel_ob_x * (xp_dot * xp_dot)) + s * (vel_ob_x * vel_ob_x) * vel_ob_y) -
    pos_ob_x * vel_ob_y * (yp_dot * yp_dot)) + pos_ob_y * vel_ob_x * (yp_dot *
    yp_dot)) + s * vel_ob_y * (xp_dot * xp_dot)) + s * vel_ob_y * (yp_dot *
    yp_dot)) - acc_ob_x * (ey * ey) * vel_ob_y) + acc_ob_y * (ey * ey) *
    vel_ob_x) + 2.0 * acc_ob_x * ey * pos_ob_y * vel_ob_y) - 2.0 * acc_ob_y * ey
    * pos_ob_y * vel_ob_x) + 2.0 * acc_ob_x * pos_ob_x * s * vel_ob_y) - 2.0 *
    acc_ob_y * pos_ob_x * s * vel_ob_x) - acc_ob_y * (ey * ey) * xp_dot * std::
    cos(epsi)) + acc_ob_x * (ey * ey) * yp_dot * std::cos(epsi)) - acc_ob_y *
    (pos_ob_x * pos_ob_x) * xp_dot * std::cos(epsi)) - acc_ob_y * (pos_ob_y *
    pos_ob_y) * xp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_x * pos_ob_x) *
    yp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_y * pos_ob_y) * yp_dot * std::
    cos(epsi)) - acc_ob_y * (s * s) * xp_dot * std::cos(epsi)) + acc_ob_x * (s *
    s) * yp_dot * std::cos(epsi)) + acc_ob_x * (ey * ey) * xp_dot * std::sin
    (epsi)) + acc_ob_y * (ey * ey) * yp_dot * std::sin(epsi)) + 2.0 * ey *
    (vel_ob_x * vel_ob_x) * xp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_x *
    pos_ob_x) * xp_dot * std::sin(epsi)) + acc_ob_x * (pos_ob_y * pos_ob_y) *
    xp_dot * std::sin(epsi)) + acc_ob_y * (pos_ob_x * pos_ob_x) * yp_dot * std::
    sin(epsi)) + acc_ob_y * (pos_ob_y * pos_ob_y) * yp_dot * std::sin(epsi)) +
    acc_ob_x * (s * s) * xp_dot * std::sin(epsi)) + acc_ob_y * (s * s) * yp_dot *
    std::sin(epsi)) - 2.0 * pos_ob_y * (vel_ob_x * vel_ob_x) * xp_dot * std::cos
    (epsi)) + 2.0 * pos_ob_x * (vel_ob_y * vel_ob_y) * yp_dot * std::cos(epsi))
    - 2.0 * s * (vel_ob_y * vel_ob_y) * yp_dot * std::cos(epsi)) - 2.0 * ey *
    (vel_ob_x * vel_ob_x) * yp_dot * std::sin(epsi)) + 2.0 * pos_ob_x *
    (vel_ob_y * vel_ob_y) * xp_dot * std::sin(epsi)) + 2.0 * pos_ob_y *
    (vel_ob_x * vel_ob_x) * yp_dot * std::sin(epsi)) - 2.0 * s * (vel_ob_y *
    vel_ob_y) * xp_dot * std::sin(epsi)) - 2.0 * s * vel_ob_x * vel_ob_y *
    xp_dot * std::cos(epsi)) + 2.0 * ey * vel_ob_x * vel_ob_y * xp_dot * std::
    sin(epsi)) - 2.0 * pos_ob_y * vel_ob_x * vel_ob_y * xp_dot * std::sin(epsi))
    - 2.0 * pos_ob_x * vel_ob_x * vel_ob_y * yp_dot * std::sin(epsi)) + 2.0 * s *
    vel_ob_x * vel_ob_y * yp_dot * std::sin(epsi)) + 2.0 * acc_ob_y * ey *
                        pos_ob_y * xp_dot * std::cos(epsi)) - 2.0 * acc_ob_x *
                       ey * pos_ob_y * yp_dot * std::cos(epsi)) + 2.0 * acc_ob_y
                      * pos_ob_x * s * xp_dot * std::cos(epsi)) - 2.0 * acc_ob_x
                     * pos_ob_x * s * yp_dot * std::cos(epsi)) - 2.0 * acc_ob_x *
                    ey * pos_ob_y * xp_dot * std::sin(epsi)) - 2.0 * acc_ob_y *
                   ey * pos_ob_y * yp_dot * std::sin(epsi)) + 2.0 * ey *
                  vel_ob_x * vel_ob_y * yp_dot * std::cos(epsi)) - 2.0 *
                 acc_ob_x * pos_ob_x * s * xp_dot * std::sin(epsi)) - 2.0 *
                acc_ob_y * pos_ob_x * s * yp_dot * std::sin(epsi)) + 2.0 *
               pos_ob_x * vel_ob_x * vel_ob_y * xp_dot * std::cos(epsi)) - 2.0 *
              pos_ob_y * vel_ob_x * vel_ob_y * yp_dot * std::cos(epsi)) / (2.0 *
              std::sqrt(1.0 - oib_a * oib_a / ((pib_a * pib_a + qib_a * qib_a) *
    (rib_a * rib_a + sib_a * sib_a))) * pow(tib_a * tib_a + uib_a *
    uib_a, 2.5) * pow(vib_a * vib_a + wib_a * wib_a, 1.5)))) -
    ((((((pos_ob_x * std::cos(epsi) - s * std::cos(epsi)) - ey * std::sin(epsi))
        + pos_ob_y * std::sin(epsi)) *
       (((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((pos_ob_y
    * pow(vel_ob_x, 3.0) - pos_ob_x * pow(vel_ob_y, 3.0)) - ey *
    pow(vel_ob_x, 3.0)) + s * pow(vel_ob_y, 3.0)) - acc_ob_x *
    (pos_ob_x * pos_ob_x) * vel_ob_y) + acc_ob_y * (pos_ob_x * pos_ob_x) *
    vel_ob_x) - acc_ob_x * (pos_ob_y * pos_ob_y) * vel_ob_y) + acc_ob_y *
    (pos_ob_y * pos_ob_y) * vel_ob_x) - acc_ob_x * (s * s) * vel_ob_y) +
    acc_ob_y * (s * s) * vel_ob_x) - ey * vel_ob_x * (vel_ob_y * vel_ob_y)) - ey
    * vel_ob_x * (xp_dot * xp_dot)) - ey * vel_ob_x * (yp_dot * yp_dot)) -
    pos_ob_x * (vel_ob_x * vel_ob_x) * vel_ob_y) + pos_ob_y * vel_ob_x *
    (vel_ob_y * vel_ob_y)) - pos_ob_x * vel_ob_y * (xp_dot * xp_dot)) + pos_ob_y
    * vel_ob_x * (xp_dot * xp_dot)) + s * (vel_ob_x * vel_ob_x) * vel_ob_y) -
    pos_ob_x * vel_ob_y * (yp_dot * yp_dot)) + pos_ob_y * vel_ob_x * (yp_dot *
    yp_dot)) + s * vel_ob_y * (xp_dot * xp_dot)) + s * vel_ob_y * (yp_dot *
    yp_dot)) - acc_ob_x * (ey * ey) * vel_ob_y) + acc_ob_y * (ey * ey) *
    vel_ob_x) + 2.0 * acc_ob_x * ey * pos_ob_y * vel_ob_y) - 2.0 * acc_ob_y * ey
    * pos_ob_y * vel_ob_x) + 2.0 * acc_ob_x * pos_ob_x * s * vel_ob_y) - 2.0 *
    acc_ob_y * pos_ob_x * s * vel_ob_x) - acc_ob_y * (ey * ey) * xp_dot * std::
    cos(epsi)) + acc_ob_x * (ey * ey) * yp_dot * std::cos(epsi)) - acc_ob_y *
    (pos_ob_x * pos_ob_x) * xp_dot * std::cos(epsi)) - acc_ob_y * (pos_ob_y *
    pos_ob_y) * xp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_x * pos_ob_x) *
    yp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_y * pos_ob_y) * yp_dot * std::
    cos(epsi)) - acc_ob_y * (s * s) * xp_dot * std::cos(epsi)) + acc_ob_x * (s *
    s) * yp_dot * std::cos(epsi)) + acc_ob_x * (ey * ey) * xp_dot * std::sin
    (epsi)) + acc_ob_y * (ey * ey) * yp_dot * std::sin(epsi)) + 2.0 * ey *
    (vel_ob_x * vel_ob_x) * xp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_x *
    pos_ob_x) * xp_dot * std::sin(epsi)) + acc_ob_x * (pos_ob_y * pos_ob_y) *
    xp_dot * std::sin(epsi)) + acc_ob_y * (pos_ob_x * pos_ob_x) * yp_dot * std::
    sin(epsi)) + acc_ob_y * (pos_ob_y * pos_ob_y) * yp_dot * std::sin(epsi)) +
    acc_ob_x * (s * s) * xp_dot * std::sin(epsi)) + acc_ob_y * (s * s) * yp_dot *
    std::sin(epsi)) - 2.0 * pos_ob_y * (vel_ob_x * vel_ob_x) * xp_dot * std::cos
    (epsi)) + 2.0 * pos_ob_x * (vel_ob_y * vel_ob_y) * yp_dot * std::cos(epsi))
    - 2.0 * s * (vel_ob_y * vel_ob_y) * yp_dot * std::cos(epsi)) - 2.0 * ey *
    (vel_ob_x * vel_ob_x) * yp_dot * std::sin(epsi)) + 2.0 * pos_ob_x *
    (vel_ob_y * vel_ob_y) * xp_dot * std::sin(epsi)) + 2.0 * pos_ob_y *
    (vel_ob_x * vel_ob_x) * yp_dot * std::sin(epsi)) - 2.0 * s * (vel_ob_y *
    vel_ob_y) * xp_dot * std::sin(epsi)) - 2.0 * s * vel_ob_x * vel_ob_y *
                       xp_dot * std::cos(epsi)) + 2.0 * ey * vel_ob_x * vel_ob_y
                      * xp_dot * std::sin(epsi)) - 2.0 * pos_ob_y * vel_ob_x *
                     vel_ob_y * xp_dot * std::sin(epsi)) - 2.0 * pos_ob_x *
                    vel_ob_x * vel_ob_y * yp_dot * std::sin(epsi)) + 2.0 * s *
                   vel_ob_x * vel_ob_y * yp_dot * std::sin(epsi)) + 2.0 *
                  acc_ob_y * ey * pos_ob_y * xp_dot * std::cos(epsi)) - 2.0 *
                 acc_ob_x * ey * pos_ob_y * yp_dot * std::cos(epsi)) + 2.0 *
                acc_ob_y * pos_ob_x * s * xp_dot * std::cos(epsi)) - 2.0 *
               acc_ob_x * pos_ob_x * s * yp_dot * std::cos(epsi)) - 2.0 *
              acc_ob_x * ey * pos_ob_y * xp_dot * std::sin(epsi)) - 2.0 *
             acc_ob_y * ey * pos_ob_y * yp_dot * std::sin(epsi)) + 2.0 * ey *
            vel_ob_x * vel_ob_y * yp_dot * std::cos(epsi)) - 2.0 * acc_ob_x *
           pos_ob_x * s * xp_dot * std::sin(epsi)) - 2.0 * acc_ob_y * pos_ob_x *
          s * yp_dot * std::sin(epsi)) + 2.0 * pos_ob_x * vel_ob_x * vel_ob_y *
         xp_dot * std::cos(epsi)) - 2.0 * pos_ob_y * vel_ob_x * vel_ob_y *
        yp_dot * std::cos(epsi)) / (std::sqrt(1.0 - xib_a * xib_a / ((yib_a *
           yib_a + ajb_a * ajb_a) * (bjb_a * bjb_a + cjb_a * cjb_a))) *
        pow(djb_a * djb_a + ejb_a * ejb_a, 1.5) * pow(fjb_a *
         fjb_a + gjb_a * gjb_a, 1.5)) - (((((((((((ey * vel_ob_x + pos_ob_x *
    vel_ob_y) - pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos
    (epsi)) + pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos
             (epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin
           (epsi)) - pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot *
         std::sin(epsi)) + s * xp_dot * std::sin(epsi)) *
       (((((((((((((((((((((((acc_ob_x * (ey * ey) * std::cos(epsi) + acc_ob_x *
    (pos_ob_x * pos_ob_x) * std::cos(epsi)) + acc_ob_x * (pos_ob_y * pos_ob_y) *
    std::cos(epsi)) + acc_ob_x * (s * s) * std::cos(epsi)) + acc_ob_y * (ey * ey)
    * std::sin(epsi)) + acc_ob_y * (pos_ob_x * pos_ob_x) * std::sin(epsi)) +
    acc_ob_y * (pos_ob_y * pos_ob_y) * std::sin(epsi)) + acc_ob_y * (s * s) *
    std::sin(epsi)) + 2.0 * pos_ob_x * (vel_ob_y * vel_ob_y) * std::cos(epsi)) -
                      2.0 * s * (vel_ob_y * vel_ob_y) * std::cos(epsi)) - 2.0 *
                     ey * (vel_ob_x * vel_ob_x) * std::sin(epsi)) + 2.0 *
                    pos_ob_y * (vel_ob_x * vel_ob_x) * std::sin(epsi)) - 2.0 *
                   ey * vel_ob_x * yp_dot) - 2.0 * pos_ob_x * vel_ob_y * yp_dot)
                 + 2.0 * pos_ob_y * vel_ob_x * yp_dot) + 2.0 * s * vel_ob_y *
                yp_dot) - 2.0 * acc_ob_x * ey * pos_ob_y * std::cos(epsi)) - 2.0
              * acc_ob_x * pos_ob_x * s * std::cos(epsi)) - 2.0 * acc_ob_y * ey *
             pos_ob_y * std::sin(epsi)) + 2.0 * ey * vel_ob_x * vel_ob_y * std::
            cos(epsi)) - 2.0 * acc_ob_y * pos_ob_x * s * std::sin(epsi)) - 2.0 *
          pos_ob_y * vel_ob_x * vel_ob_y * std::cos(epsi)) - 2.0 * pos_ob_x *
         vel_ob_x * vel_ob_y * std::sin(epsi)) + 2.0 * s * vel_ob_x * vel_ob_y *
        std::sin(epsi)) / (std::sqrt(1.0 - hjb_a * hjb_a / ((ijb_a * ijb_a +
           jjb_a * jjb_a) * (kjb_a * kjb_a + ljb_a * ljb_a))) * pow
        (mjb_a * mjb_a + njb_a * njb_a, 1.5) * pow(ojb_a * ojb_a + pjb_a
         * pjb_a, 1.5))) + 3.0 * ((2.0 * yp_dot - 2.0 * vel_ob_y * std::cos(epsi))
       + 2.0 * vel_ob_x * std::sin(epsi)) * (((((((((((ey * vel_ob_x + pos_ob_x *
    vel_ob_y) - pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos
              (epsi)) + pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot *
            std::cos(epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::
          sin(epsi)) - pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot *
        std::sin(epsi)) + s * xp_dot * std::sin(epsi)) *
      (((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((pos_ob_y
    * pow(vel_ob_x, 3.0) - pos_ob_x * pow(vel_ob_y, 3.0)) - ey *
    pow(vel_ob_x, 3.0)) + s * pow(vel_ob_y, 3.0)) - acc_ob_x *
    (pos_ob_x * pos_ob_x) * vel_ob_y) + acc_ob_y * (pos_ob_x * pos_ob_x) *
    vel_ob_x) - acc_ob_x * (pos_ob_y * pos_ob_y) * vel_ob_y) + acc_ob_y *
    (pos_ob_y * pos_ob_y) * vel_ob_x) - acc_ob_x * (s * s) * vel_ob_y) +
    acc_ob_y * (s * s) * vel_ob_x) - ey * vel_ob_x * (vel_ob_y * vel_ob_y)) - ey
    * vel_ob_x * (xp_dot * xp_dot)) - ey * vel_ob_x * (yp_dot * yp_dot)) -
    pos_ob_x * (vel_ob_x * vel_ob_x) * vel_ob_y) + pos_ob_y * vel_ob_x *
    (vel_ob_y * vel_ob_y)) - pos_ob_x * vel_ob_y * (xp_dot * xp_dot)) + pos_ob_y
    * vel_ob_x * (xp_dot * xp_dot)) + s * (vel_ob_x * vel_ob_x) * vel_ob_y) -
    pos_ob_x * vel_ob_y * (yp_dot * yp_dot)) + pos_ob_y * vel_ob_x * (yp_dot *
    yp_dot)) + s * vel_ob_y * (xp_dot * xp_dot)) + s * vel_ob_y * (yp_dot *
    yp_dot)) - acc_ob_x * (ey * ey) * vel_ob_y) + acc_ob_y * (ey * ey) *
    vel_ob_x) + 2.0 * acc_ob_x * ey * pos_ob_y * vel_ob_y) - 2.0 * acc_ob_y * ey
    * pos_ob_y * vel_ob_x) + 2.0 * acc_ob_x * pos_ob_x * s * vel_ob_y) - 2.0 *
    acc_ob_y * pos_ob_x * s * vel_ob_x) - acc_ob_y * (ey * ey) * xp_dot * std::
    cos(epsi)) + acc_ob_x * (ey * ey) * yp_dot * std::cos(epsi)) - acc_ob_y *
    (pos_ob_x * pos_ob_x) * xp_dot * std::cos(epsi)) - acc_ob_y * (pos_ob_y *
    pos_ob_y) * xp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_x * pos_ob_x) *
    yp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_y * pos_ob_y) * yp_dot * std::
    cos(epsi)) - acc_ob_y * (s * s) * xp_dot * std::cos(epsi)) + acc_ob_x * (s *
    s) * yp_dot * std::cos(epsi)) + acc_ob_x * (ey * ey) * xp_dot * std::sin
    (epsi)) + acc_ob_y * (ey * ey) * yp_dot * std::sin(epsi)) + 2.0 * ey *
    (vel_ob_x * vel_ob_x) * xp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_x *
    pos_ob_x) * xp_dot * std::sin(epsi)) + acc_ob_x * (pos_ob_y * pos_ob_y) *
    xp_dot * std::sin(epsi)) + acc_ob_y * (pos_ob_x * pos_ob_x) * yp_dot * std::
    sin(epsi)) + acc_ob_y * (pos_ob_y * pos_ob_y) * yp_dot * std::sin(epsi)) +
    acc_ob_x * (s * s) * xp_dot * std::sin(epsi)) + acc_ob_y * (s * s) * yp_dot *
    std::sin(epsi)) - 2.0 * pos_ob_y * (vel_ob_x * vel_ob_x) * xp_dot * std::cos
    (epsi)) + 2.0 * pos_ob_x * (vel_ob_y * vel_ob_y) * yp_dot * std::cos(epsi))
    - 2.0 * s * (vel_ob_y * vel_ob_y) * yp_dot * std::cos(epsi)) - 2.0 * ey *
    (vel_ob_x * vel_ob_x) * yp_dot * std::sin(epsi)) + 2.0 * pos_ob_x *
    (vel_ob_y * vel_ob_y) * xp_dot * std::sin(epsi)) + 2.0 * pos_ob_y *
    (vel_ob_x * vel_ob_x) * yp_dot * std::sin(epsi)) - 2.0 * s * (vel_ob_y *
    vel_ob_y) * xp_dot * std::sin(epsi)) - 2.0 * s * vel_ob_x * vel_ob_y *
                      xp_dot * std::cos(epsi)) + 2.0 * ey * vel_ob_x * vel_ob_y *
                     xp_dot * std::sin(epsi)) - 2.0 * pos_ob_y * vel_ob_x *
                    vel_ob_y * xp_dot * std::sin(epsi)) - 2.0 * pos_ob_x *
                   vel_ob_x * vel_ob_y * yp_dot * std::sin(epsi)) + 2.0 * s *
                  vel_ob_x * vel_ob_y * yp_dot * std::sin(epsi)) + 2.0 *
                 acc_ob_y * ey * pos_ob_y * xp_dot * std::cos(epsi)) - 2.0 *
                acc_ob_x * ey * pos_ob_y * yp_dot * std::cos(epsi)) + 2.0 *
               acc_ob_y * pos_ob_x * s * xp_dot * std::cos(epsi)) - 2.0 *
              acc_ob_x * pos_ob_x * s * yp_dot * std::cos(epsi)) - 2.0 *
             acc_ob_x * ey * pos_ob_y * xp_dot * std::sin(epsi)) - 2.0 *
            acc_ob_y * ey * pos_ob_y * yp_dot * std::sin(epsi)) + 2.0 * ey *
           vel_ob_x * vel_ob_y * yp_dot * std::cos(epsi)) - 2.0 * acc_ob_x *
          pos_ob_x * s * xp_dot * std::sin(epsi)) - 2.0 * acc_ob_y * pos_ob_x *
         s * yp_dot * std::sin(epsi)) + 2.0 * pos_ob_x * vel_ob_x * vel_ob_y *
        xp_dot * std::cos(epsi)) - 2.0 * pos_ob_y * vel_ob_x * vel_ob_y * yp_dot
       * std::cos(epsi)) / (2.0 * std::sqrt(1.0 - qjb_a * qjb_a / ((rjb_a *
          rjb_a + sjb_a * sjb_a) * (tjb_a * tjb_a + ujb_a * ujb_a))) *
       pow(vjb_a * vjb_a + wjb_a * wjb_a, 1.5) * pow(xjb_a *
        xjb_a + yjb_a * yjb_a, 2.5))) - (2.0 * ((pos_ob_x - s) * ((vel_ob_x -
         xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi)) + (ey - pos_ob_y) *
       ((yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi))) * (((ey
         * std::cos(epsi) - pos_ob_y * std::cos(epsi)) + pos_ob_x * std::sin
        (epsi)) - s * std::sin(epsi)) / ((akb_a * akb_a + bkb_a * bkb_a) *
       (ckb_a * ckb_a + dkb_a * dkb_a)) - ekb_a * ekb_a * ((2.0 * yp_dot - 2.0 *
        vel_ob_y * std::cos(epsi)) + 2.0 * vel_ob_x * std::sin(epsi)) / ((fkb_a *
        fkb_a + gkb_a * gkb_a) * (hkb_a * hkb_a))) * (((((((((((ey * vel_ob_x +
    pos_ob_x * vel_ob_y) - pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot *
             std::cos(epsi)) + pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x *
           yp_dot * std::cos(epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot
         * std::sin(epsi)) - pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y *
       yp_dot * std::sin(epsi)) + s * xp_dot * std::sin(epsi)) *
     (((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((pos_ob_y
    * pow(vel_ob_x, 3.0) - pos_ob_x * pow(vel_ob_y, 3.0)) - ey *
    pow(vel_ob_x, 3.0)) + s * pow(vel_ob_y, 3.0)) - acc_ob_x *
    (pos_ob_x * pos_ob_x) * vel_ob_y) + acc_ob_y * (pos_ob_x * pos_ob_x) *
    vel_ob_x) - acc_ob_x * (pos_ob_y * pos_ob_y) * vel_ob_y) + acc_ob_y *
    (pos_ob_y * pos_ob_y) * vel_ob_x) - acc_ob_x * (s * s) * vel_ob_y) +
    acc_ob_y * (s * s) * vel_ob_x) - ey * vel_ob_x * (vel_ob_y * vel_ob_y)) - ey
    * vel_ob_x * (xp_dot * xp_dot)) - ey * vel_ob_x * (yp_dot * yp_dot)) -
    pos_ob_x * (vel_ob_x * vel_ob_x) * vel_ob_y) + pos_ob_y * vel_ob_x *
    (vel_ob_y * vel_ob_y)) - pos_ob_x * vel_ob_y * (xp_dot * xp_dot)) + pos_ob_y
    * vel_ob_x * (xp_dot * xp_dot)) + s * (vel_ob_x * vel_ob_x) * vel_ob_y) -
    pos_ob_x * vel_ob_y * (yp_dot * yp_dot)) + pos_ob_y * vel_ob_x * (yp_dot *
    yp_dot)) + s * vel_ob_y * (xp_dot * xp_dot)) + s * vel_ob_y * (yp_dot *
    yp_dot)) - acc_ob_x * (ey * ey) * vel_ob_y) + acc_ob_y * (ey * ey) *
    vel_ob_x) + 2.0 * acc_ob_x * ey * pos_ob_y * vel_ob_y) - 2.0 * acc_ob_y * ey
    * pos_ob_y * vel_ob_x) + 2.0 * acc_ob_x * pos_ob_x * s * vel_ob_y) - 2.0 *
    acc_ob_y * pos_ob_x * s * vel_ob_x) - acc_ob_y * (ey * ey) * xp_dot * std::
    cos(epsi)) + acc_ob_x * (ey * ey) * yp_dot * std::cos(epsi)) - acc_ob_y *
    (pos_ob_x * pos_ob_x) * xp_dot * std::cos(epsi)) - acc_ob_y * (pos_ob_y *
    pos_ob_y) * xp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_x * pos_ob_x) *
    yp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_y * pos_ob_y) * yp_dot * std::
    cos(epsi)) - acc_ob_y * (s * s) * xp_dot * std::cos(epsi)) + acc_ob_x * (s *
    s) * yp_dot * std::cos(epsi)) + acc_ob_x * (ey * ey) * xp_dot * std::sin
    (epsi)) + acc_ob_y * (ey * ey) * yp_dot * std::sin(epsi)) + 2.0 * ey *
    (vel_ob_x * vel_ob_x) * xp_dot * std::cos(epsi)) + acc_ob_x * (pos_ob_x *
    pos_ob_x) * xp_dot * std::sin(epsi)) + acc_ob_x * (pos_ob_y * pos_ob_y) *
    xp_dot * std::sin(epsi)) + acc_ob_y * (pos_ob_x * pos_ob_x) * yp_dot * std::
    sin(epsi)) + acc_ob_y * (pos_ob_y * pos_ob_y) * yp_dot * std::sin(epsi)) +
    acc_ob_x * (s * s) * xp_dot * std::sin(epsi)) + acc_ob_y * (s * s) * yp_dot *
    std::sin(epsi)) - 2.0 * pos_ob_y * (vel_ob_x * vel_ob_x) * xp_dot * std::cos
    (epsi)) + 2.0 * pos_ob_x * (vel_ob_y * vel_ob_y) * yp_dot * std::cos(epsi))
    - 2.0 * s * (vel_ob_y * vel_ob_y) * yp_dot * std::cos(epsi)) - 2.0 * ey *
    (vel_ob_x * vel_ob_x) * yp_dot * std::sin(epsi)) + 2.0 * pos_ob_x *
    (vel_ob_y * vel_ob_y) * xp_dot * std::sin(epsi)) + 2.0 * pos_ob_y *
                       (vel_ob_x * vel_ob_x) * yp_dot * std::sin(epsi)) - 2.0 *
                      s * (vel_ob_y * vel_ob_y) * xp_dot * std::sin(epsi)) - 2.0
                     * s * vel_ob_x * vel_ob_y * xp_dot * std::cos(epsi)) + 2.0 *
                    ey * vel_ob_x * vel_ob_y * xp_dot * std::sin(epsi)) - 2.0 *
                   pos_ob_y * vel_ob_x * vel_ob_y * xp_dot * std::sin(epsi)) -
                  2.0 * pos_ob_x * vel_ob_x * vel_ob_y * yp_dot * std::sin(epsi))
                 + 2.0 * s * vel_ob_x * vel_ob_y * yp_dot * std::sin(epsi)) +
                2.0 * acc_ob_y * ey * pos_ob_y * xp_dot * std::cos(epsi)) - 2.0 *
               acc_ob_x * ey * pos_ob_y * yp_dot * std::cos(epsi)) + 2.0 *
              acc_ob_y * pos_ob_x * s * xp_dot * std::cos(epsi)) - 2.0 *
             acc_ob_x * pos_ob_x * s * yp_dot * std::cos(epsi)) - 2.0 * acc_ob_x
            * ey * pos_ob_y * xp_dot * std::sin(epsi)) - 2.0 * acc_ob_y * ey *
           pos_ob_y * yp_dot * std::sin(epsi)) + 2.0 * ey * vel_ob_x * vel_ob_y *
          yp_dot * std::cos(epsi)) - 2.0 * acc_ob_x * pos_ob_x * s * xp_dot *
         std::sin(epsi)) - 2.0 * acc_ob_y * pos_ob_x * s * yp_dot * std::sin
        (epsi)) + 2.0 * pos_ob_x * vel_ob_x * vel_ob_y * xp_dot * std::cos(epsi))
      - 2.0 * pos_ob_y * vel_ob_x * vel_ob_y * yp_dot * std::cos(epsi)) / (2.0 *
      pow(1.0 - ikb_a * ikb_a / ((jkb_a * jkb_a + kkb_a * kkb_a) *
        (lkb_a * lkb_a + mkb_a * mkb_a)), 1.5) * pow(nkb_a * nkb_a +
       okb_a * okb_a, 1.5) * pow(pkb_a * pkb_a + qkb_a * qkb_a, 1.5))) *
    ((((2.0 * cf * yp_dot + 2.0 * cr * yp_dot) + m * psi_dot * (xp_dot * xp_dot))
      + 2.0 * a * cf * psi_dot) - 2.0 * b * cr * psi_dot) / (m * xp_dot);
  out[2] = 2.0 * cf * ((((((((((((((((Ds * (((ey * std::cos(epsi) - pos_ob_y *
    std::cos(epsi)) + pos_ob_x * std::sin(epsi)) - s * std::sin(epsi)) / (std::
    sqrt(1.0 - Ds * Ds / (rkb_a * rkb_a + skb_a * skb_a)) * pow(tkb_a *
    tkb_a + ukb_a * ukb_a, 1.5)) - std::cos(epsi) * (pos_ob_x - s) *
    (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y * vel_ob_x) - s *
    vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y * xp_dot * std::cos
    (epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) + s * yp_dot * std::cos(epsi))
    + ey * yp_dot * std::sin(epsi)) - pos_ob_x * xp_dot * std::sin(epsi)) -
      pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot * std::sin(epsi)) / (std::
    sqrt(1.0 - vkb_a * vkb_a / ((wkb_a * wkb_a + xkb_a * xkb_a) * (ykb_a * ykb_a
    + alb_a * alb_a))) * pow(blb_a * blb_a + clb_a * clb_a, 1.5) * std::
    sqrt(dlb_a * dlb_a + elb_a * elb_a))) + std::sin(epsi) * (ey - pos_ob_y) *
    (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y * vel_ob_x) - s *
    vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y * xp_dot * std::cos
    (epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) + s * yp_dot * std::cos(epsi))
        + ey * yp_dot * std::sin(epsi)) - pos_ob_x * xp_dot * std::sin(epsi)) -
      pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot * std::sin(epsi)) / (std::
    sqrt(1.0 - flb_a * flb_a / ((glb_a * glb_a + hlb_a * hlb_a) * (ilb_a * ilb_a
    + jlb_a * jlb_a))) * pow(klb_a * klb_a + llb_a * llb_a, 1.5) * std::
    sqrt(mlb_a * mlb_a + nlb_a * nlb_a))) + (psi_dot - psi_dot_com) *
    (((pos_ob_x * std::cos(epsi) - s * std::cos(epsi)) - ey * std::sin(epsi)) +
     pos_ob_y * std::sin(epsi)) * (((((xp_dot * xp_dot + yp_dot * yp_dot) -
    vel_ob_x * xp_dot * std::cos(epsi)) - vel_ob_y * yp_dot * std::cos(epsi)) -
    vel_ob_y * xp_dot * std::sin(epsi)) + vel_ob_x * yp_dot * std::sin(epsi)) /
    (std::sqrt(1.0 - olb_a * olb_a / ((plb_a * plb_a + qlb_a * qlb_a) * (rlb_a *
    rlb_a + slb_a * slb_a))) * std::sqrt(tlb_a * tlb_a + ulb_a * ulb_a) *
     pow(vlb_a * vlb_a + wlb_a * wlb_a, 1.5))) + (ey - pos_ob_y) *
    (xp_dot * std::cos(epsi) - yp_dot * std::sin(epsi)) * (((pos_ob_x * std::cos
    (epsi) - s * std::cos(epsi)) - ey * std::sin(epsi)) + pos_ob_y * std::sin
    (epsi)) / (std::sqrt(1.0 - xlb_a * xlb_a / ((ylb_a * ylb_a + amb_a * amb_a) *
    (bmb_a * bmb_a + cmb_a * cmb_a))) * pow(dmb_a * dmb_a + emb_a *
    emb_a, 1.5) * std::sqrt(fmb_a * fmb_a + gmb_a * gmb_a))) + (pos_ob_x - s) *
    (yp_dot * std::cos(epsi) + xp_dot * std::sin(epsi)) * (((pos_ob_x * std::cos
    (epsi) - s * std::cos(epsi)) - ey * std::sin(epsi)) + pos_ob_y * std::sin
    (epsi)) / (std::sqrt(1.0 - hmb_a * hmb_a / ((imb_a * imb_a + jmb_a * jmb_a) *
    (kmb_a * kmb_a + lmb_a * lmb_a))) * pow(mmb_a * mmb_a + nmb_a *
    nmb_a, 1.5) * std::sqrt(omb_a * omb_a + pmb_a * pmb_a))) - (psi_dot -
    psi_dot_com) * ((2.0 * yp_dot - vel_ob_y * std::cos(epsi)) + vel_ob_x * std::
                    sin(epsi)) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y)
    - pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) +
    pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) +
    s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x *
    xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot *
    std::sin(epsi)) / (std::sqrt(1.0 - qmb_a * qmb_a / ((rmb_a * rmb_a + smb_a *
    smb_a) * (tmb_a * tmb_a + umb_a * umb_a))) * std::sqrt(vmb_a * vmb_a + wmb_a
    * wmb_a) * pow(xmb_a * xmb_a + ymb_a * ymb_a, 1.5))) + 3.0 *
    (psi_dot - psi_dot_com) * ((2.0 * yp_dot - 2.0 * vel_ob_y * std::cos(epsi))
    + 2.0 * vel_ob_x * std::sin(epsi)) * (((((xp_dot * xp_dot + yp_dot * yp_dot)
    - vel_ob_x * xp_dot * std::cos(epsi)) - vel_ob_y * yp_dot * std::cos(epsi))
    - vel_ob_y * xp_dot * std::sin(epsi)) + vel_ob_x * yp_dot * std::sin(epsi)) *
    (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y * vel_ob_x) - s *
             vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y * xp_dot * std::
           cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) + s * yp_dot * std::
         cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x * xp_dot * std::
       sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot * std::sin
     (epsi)) / (2.0 * std::sqrt(1.0 - anb_a * anb_a / ((bnb_a * bnb_a + cnb_a *
    cnb_a) * (dnb_a * dnb_a + enb_a * enb_a))) * std::sqrt(fnb_a * fnb_a + gnb_a
    * gnb_a) * pow(hnb_a * hnb_a + inb_a * inb_a, 2.5))) + (ey -
    pos_ob_y) * (xp_dot * std::cos(epsi) - yp_dot * std::sin(epsi)) * ((2.0 *
    yp_dot - 2.0 * vel_ob_y * std::cos(epsi)) + 2.0 * vel_ob_x * std::sin(epsi))
    * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y * vel_ob_x) - s
               * vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y * xp_dot *
             std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) + s * yp_dot *
           std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x * xp_dot *
         std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot *
       std::sin(epsi)) / (2.0 * std::sqrt(1.0 - jnb_a * jnb_a / ((knb_a * knb_a
    + lnb_a * lnb_a) * (mnb_a * mnb_a + nnb_a * nnb_a))) * pow(onb_a *
    onb_a + pnb_a * pnb_a, 1.5) * pow(qnb_a * qnb_a + rnb_a * rnb_a, 1.5)))
    + (pos_ob_x - s) * (yp_dot * std::cos(epsi) + xp_dot * std::sin(epsi)) *
    ((2.0 * yp_dot - 2.0 * vel_ob_y * std::cos(epsi)) + 2.0 * vel_ob_x * std::
     sin(epsi)) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y *
    vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y *
    xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) + s * yp_dot *
                       std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) -
                     pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot *
                    std::sin(epsi)) + s * xp_dot * std::sin(epsi)) / (2.0 * std::
    sqrt(1.0 - snb_a * snb_a / ((tnb_a * tnb_a + unb_a * unb_a) * (vnb_a * vnb_a
    + wnb_a * wnb_a))) * pow(xnb_a * xnb_a + ynb_a * ynb_a, 1.5) *
    pow(aob_a * aob_a + bob_a * bob_a, 1.5))) - (psi_dot - psi_dot_com) *
    (2.0 * ((pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot *
    std::sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y) +
    xp_dot * std::sin(epsi))) * (((ey * std::cos(epsi) - pos_ob_y * std::cos
    (epsi)) + pos_ob_x * std::sin(epsi)) - s * std::sin(epsi)) / ((cob_a * cob_a
    + dob_a * dob_a) * (eob_a * eob_a + fob_a * fob_a)) - gob_a * gob_a * ((2.0 *
    yp_dot - 2.0 * vel_ob_y * std::cos(epsi)) + 2.0 * vel_ob_x * std::sin(epsi))
     / ((hob_a * hob_a + iob_a * iob_a) * (job_a * job_a))) * (((((xp_dot *
    xp_dot + yp_dot * yp_dot) - vel_ob_x * xp_dot * std::cos(epsi)) - vel_ob_y *
    yp_dot * std::cos(epsi)) - vel_ob_y * xp_dot * std::sin(epsi)) + vel_ob_x *
    yp_dot * std::sin(epsi)) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) -
    pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) +
    pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) +
    s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x *
    xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot *
    std::sin(epsi)) / (2.0 * pow(1.0 - kob_a * kob_a / ((lob_a * lob_a +
    mob_a * mob_a) * (nob_a * nob_a + oob_a * oob_a)), 1.5) * std::sqrt(pob_a *
    pob_a + qob_a * qob_a) * pow(rob_a * rob_a + sob_a * sob_a, 1.5))) -
    (ey - pos_ob_y) * (2.0 * ((pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos
    (epsi)) + yp_dot * std::sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos
    (epsi) - vel_ob_y) + xp_dot * std::sin(epsi))) * (((ey * std::cos(epsi) -
    pos_ob_y * std::cos(epsi)) + pos_ob_x * std::sin(epsi)) - s * std::sin(epsi))
                       / ((tob_a * tob_a + uob_a * uob_a) * (vob_a * vob_a +
    wob_a * wob_a)) - xob_a * xob_a * ((2.0 * yp_dot - 2.0 * vel_ob_y * std::cos
    (epsi)) + 2.0 * vel_ob_x * std::sin(epsi)) / ((yob_a * yob_a + apb_a * apb_a)
    * (bpb_a * bpb_a))) * (xp_dot * std::cos(epsi) - yp_dot * std::sin(epsi)) *
    (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y * vel_ob_x) - s *
             vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y * xp_dot * std::
           cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) + s * yp_dot * std::
         cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x * xp_dot * std::
       sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot * std::sin
     (epsi)) / (2.0 * pow(1.0 - cpb_a * cpb_a / ((dpb_a * dpb_a + epb_a *
    epb_a) * (fpb_a * fpb_a + gpb_a * gpb_a)), 1.5) * pow(hpb_a * hpb_a
    + ipb_a * ipb_a, 1.5) * std::sqrt(jpb_a * jpb_a + kpb_a * kpb_a))) -
    (pos_ob_x - s) * (2.0 * ((pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos
    (epsi)) + yp_dot * std::sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos
    (epsi) - vel_ob_y) + xp_dot * std::sin(epsi))) * (((ey * std::cos(epsi) -
    pos_ob_y * std::cos(epsi)) + pos_ob_x * std::sin(epsi)) - s * std::sin(epsi))
                      / ((lpb_a * lpb_a + mpb_a * mpb_a) * (npb_a * npb_a +
    opb_a * opb_a)) - ppb_a * ppb_a * ((2.0 * yp_dot - 2.0 * vel_ob_y * std::cos
    (epsi)) + 2.0 * vel_ob_x * std::sin(epsi)) / ((qpb_a * qpb_a + rpb_a * rpb_a)
    * (spb_a * spb_a))) * (yp_dot * std::cos(epsi) + xp_dot * std::sin(epsi)) *
    (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y * vel_ob_x) - s *
             vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y * xp_dot * std::
           cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) + s * yp_dot * std::
         cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x * xp_dot * std::
       sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot * std::sin
     (epsi)) / (2.0 * pow(1.0 - tpb_a * tpb_a / ((upb_a * upb_a + vpb_a *
    vpb_a) * (wpb_a * wpb_a + xpb_a * xpb_a)), 1.5) * pow(ypb_a * ypb_a
    + aqb_a * aqb_a, 1.5) * std::sqrt(bqb_a * bqb_a + cqb_a * cqb_a))) +
    ((vel_ob_x * std::cos(epsi) - xp_dot) + vel_ob_y * std::sin(epsi)) *
    (((pos_ob_x * std::cos(epsi) - s * std::cos(epsi)) - ey * std::sin(epsi)) +
     pos_ob_y * std::sin(epsi)) * ((((m * psi_dot * (xp_dot * xp_dot) + 2.0 * cf
    * yp_dot) + 2.0 * cr * yp_dot) + 2.0 * a * cf * psi_dot) - 2.0 * b * cr *
    psi_dot) / (m * xp_dot * std::sqrt(1.0 - dqb_a * dqb_a / ((eqb_a * eqb_a +
    fqb_a * fqb_a) * (gqb_a * gqb_a + hqb_a * hqb_a))) * std::sqrt(iqb_a * iqb_a
    + jqb_a * jqb_a) * pow(kqb_a * kqb_a + lqb_a * lqb_a, 1.5))) - (2.0 *
    cf + 2.0 * cr) * ((vel_ob_x * std::cos(epsi) - xp_dot) + vel_ob_y * std::sin
                      (epsi)) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y)
    - pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) +
    pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) +
    s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x *
    xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot *
    std::sin(epsi)) / (m * xp_dot * std::sqrt(1.0 - mqb_a * mqb_a / ((nqb_a *
    nqb_a + oqb_a * oqb_a) * (pqb_a * pqb_a + qqb_a * qqb_a))) * std::sqrt(rqb_a
    * rqb_a + sqb_a * sqb_a) * pow(tqb_a * tqb_a + uqb_a * uqb_a, 1.5)))
                        + 3.0 * ((vel_ob_x * std::cos(epsi) - xp_dot) + vel_ob_y
    * std::sin(epsi)) * ((2.0 * yp_dot - 2.0 * vel_ob_y * std::cos(epsi)) + 2.0 *
    vel_ob_x * std::sin(epsi)) * ((((m * psi_dot * (xp_dot * xp_dot) + 2.0 * cf *
    yp_dot) + 2.0 * cr * yp_dot) + 2.0 * a * cf * psi_dot) - 2.0 * b * cr *
    psi_dot) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y *
    vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y *
                      xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos
                     (epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::
                   sin(epsi)) - pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y *
                 yp_dot * std::sin(epsi)) + s * xp_dot * std::sin(epsi)) / (2.0 *
    m * xp_dot * std::sqrt(1.0 - vqb_a * vqb_a / ((wqb_a * wqb_a + xqb_a * xqb_a)
    * (yqb_a * yqb_a + arb_a * arb_a))) * std::sqrt(brb_a * brb_a + crb_a *
    crb_a) * pow(drb_a * drb_a + erb_a * erb_a, 2.5))) - (2.0 *
    ((pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin
                       (epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) -
    vel_ob_y) + xp_dot * std::sin(epsi))) * (((ey * std::cos(epsi) - pos_ob_y *
    std::cos(epsi)) + pos_ob_x * std::sin(epsi)) - s * std::sin(epsi)) / ((frb_a
    * frb_a + grb_a * grb_a) * (hrb_a * hrb_a + irb_a * irb_a)) - jrb_a * jrb_a *
    ((2.0 * yp_dot - 2.0 * vel_ob_y * std::cos(epsi)) + 2.0 * vel_ob_x * std::
     sin(epsi)) / ((krb_a * krb_a + lrb_a * lrb_a) * (mrb_a * mrb_a))) *
                       ((vel_ob_x * std::cos(epsi) - xp_dot) + vel_ob_y * std::
                        sin(epsi)) * ((((m * psi_dot * (xp_dot * xp_dot) + 2.0 *
    cf * yp_dot) + 2.0 * cr * yp_dot) + 2.0 * a * cf * psi_dot) - 2.0 * b * cr *
    psi_dot) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y *
    vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y *
                      xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos
                     (epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::
                   sin(epsi)) - pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y *
                 yp_dot * std::sin(epsi)) + s * xp_dot * std::sin(epsi)) / (2.0 *
    m * xp_dot * pow(1.0 - nrb_a * nrb_a / ((orb_a * orb_a + prb_a *
    prb_a) * (qrb_a * qrb_a + rrb_a * rrb_a)), 1.5) * std::sqrt(srb_a * srb_a +
    trb_a * trb_a) * pow(urb_a * urb_a + vrb_a * vrb_a, 1.5))) / m - 2.0
    * a * cf * ((((((((m * xp_dot * (yp_dot * yp_dot) - 2.0 * a * cf * xp_dot) +
                      2.0 * b * cr * xp_dot) + 2.0 * a * cf * vel_ob_x * std::
                     cos(epsi)) - 2.0 * b * cr * vel_ob_x * std::cos(epsi)) +
                   2.0 * a * cf * vel_ob_y * std::sin(epsi)) - 2.0 * b * cr *
                  vel_ob_y * std::sin(epsi)) + m * vel_ob_x * xp_dot * yp_dot *
                 std::sin(epsi)) - m * vel_ob_y * xp_dot * yp_dot * std::cos
                (epsi)) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) -
    pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) +
    pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) +
    s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x *
    xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot *
    std::sin(epsi)) / (Iz * m * xp_dot * std::sqrt(1.0 - wrb_a * wrb_a / ((xrb_a
    * xrb_a + yrb_a * yrb_a) * (asb_a * asb_a + bsb_a * bsb_a))) * std::sqrt
                       (csb_a * csb_a + dsb_a * dsb_a) * pow(esb_a *
    esb_a + fsb_a * fsb_a, 1.5));
  out[3] = (cf * steer * (2.0 * ((pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos
    (epsi)) + yp_dot * std::sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos
    (epsi) - vel_ob_y) + xp_dot * std::sin(epsi))) * (((ey * std::cos(epsi) -
    pos_ob_y * std::cos(epsi)) + pos_ob_x * std::sin(epsi)) - s * std::sin(epsi))
             / ((gsb_a * gsb_a + hsb_a * hsb_a) * (isb_a * isb_a + jsb_a * jsb_a))
             - ksb_a * ksb_a * ((2.0 * yp_dot - 2.0 * vel_ob_y * std::cos(epsi))
              + 2.0 * vel_ob_x * std::sin(epsi)) / ((lsb_a * lsb_a + msb_a *
    msb_a) * (nsb_a * nsb_a))) * ((vel_ob_x * std::cos(epsi) - xp_dot) +
             vel_ob_y * std::sin(epsi)) * (((((((((((ey * vel_ob_x + pos_ob_x *
    vel_ob_y) - pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos
    (epsi)) + pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos
    (epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) -
    pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) +
             s * xp_dot * std::sin(epsi)) / (m * pow(1.0 - osb_a * osb_a
              / ((psb_a * psb_a + qsb_a * qsb_a) * (rsb_a * rsb_a + ssb_a *
    ssb_a)), 1.5) * std::sqrt(tsb_a * tsb_a + usb_a * usb_a) * pow(vsb_a
              * vsb_a + wsb_a * wsb_a, 1.5)) - 3.0 * cf * steer * ((vel_ob_x *
              std::cos(epsi) - xp_dot) + vel_ob_y * std::sin(epsi)) * ((2.0 *
              yp_dot - 2.0 * vel_ob_y * std::cos(epsi)) + 2.0 * vel_ob_x * std::
             sin(epsi)) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) -
    pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) +
    pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) +
    s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x *
    xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot *
             std::sin(epsi)) / (m * std::sqrt(1.0 - xsb_a * xsb_a / ((ysb_a *
    ysb_a + atb_a * atb_a) * (btb_a * btb_a + ctb_a * ctb_a))) * std::sqrt(dtb_a
              * dtb_a + etb_a * etb_a) * pow(ftb_a * ftb_a + gtb_a *
              gtb_a, 2.5))) - 2.0 * cf * steer * ((vel_ob_x * std::cos(epsi) -
    xp_dot) + vel_ob_y * std::sin(epsi)) * (((pos_ob_x * std::cos(epsi) - s *
    std::cos(epsi)) - ey * std::sin(epsi)) + pos_ob_y * std::sin(epsi)) / (m *
    std::sqrt(1.0 - htb_a * htb_a / ((itb_a * itb_a + jtb_a * jtb_a) * (ktb_a *
    ktb_a + ltb_a * ltb_a))) * std::sqrt(mtb_a * mtb_a + ntb_a * ntb_a) *
    pow(otb_a * otb_a + ptb_a * ptb_a, 1.5));
  out[4] = 0.0;
  out[5] = ((2.0 * cf * steer * (vel_ob_y * std::cos(epsi) - vel_ob_x * std::sin
              (epsi)) * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) -
    pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi)) +
    pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) +
    s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x *
    xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot *
              std::sin(epsi)) / (m * std::sqrt(1.0 - qtb_a * qtb_a / ((rtb_a *
    rtb_a + stb_a * stb_a) * (ttb_a * ttb_a + utb_a * utb_a))) * std::sqrt(vtb_a
    * vtb_a + wtb_a * wtb_a) * pow(xtb_a * xtb_a + ytb_a * ytb_a, 1.5))
             + 2.0 * cf * steer * ((vel_ob_x * std::cos(epsi) - xp_dot) +
              vel_ob_y * std::sin(epsi)) * (((((((ey * yp_dot * std::cos(epsi) -
    pos_ob_x * xp_dot * std::cos(epsi)) - pos_ob_y * yp_dot * std::cos(epsi)) +
    s * xp_dot * std::cos(epsi)) + ey * xp_dot * std::sin(epsi)) - pos_ob_y *
    xp_dot * std::sin(epsi)) + pos_ob_x * yp_dot * std::sin(epsi)) - s * yp_dot *
              std::sin(epsi)) / (m * std::sqrt(1.0 - aub_a * aub_a / ((bub_a *
    bub_a + cub_a * cub_a) * (dub_a * dub_a + eub_a * eub_a))) * std::sqrt(fub_a
    * fub_a + gub_a * gub_a) * pow(hub_a * hub_a + iub_a * iub_a, 1.5)))
            - 3.0 * cf * steer * ((vel_ob_x * std::cos(epsi) - xp_dot) +
             vel_ob_y * std::sin(epsi)) * (((2.0 * vel_ob_x * yp_dot * std::cos
    (epsi) - 2.0 * vel_ob_y * xp_dot * std::cos(epsi)) + 2.0 * vel_ob_x * xp_dot
              * std::sin(epsi)) + 2.0 * vel_ob_y * yp_dot * std::sin(epsi)) *
            (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y *
                      vel_ob_x) - s * vel_ob_y) - ey * xp_dot * std::cos(epsi))
                   + pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot *
                  std::cos(epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot *
                std::sin(epsi)) - pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y
              * yp_dot * std::sin(epsi)) + s * xp_dot * std::sin(epsi)) / (m *
             std::sqrt(1.0 - jub_a * jub_a / ((kub_a * kub_a + lub_a * lub_a) *
    (mub_a * mub_a + nub_a * nub_a))) * std::sqrt(oub_a * oub_a + pub_a * pub_a)
             * pow(qub_a * qub_a + rub_a * rub_a, 2.5))) + cf * steer *
    (2.0 * ((pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot *
       std::sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) - vel_ob_y)
       + xp_dot * std::sin(epsi))) * ((ey - pos_ob_y) * (xp_dot * std::cos(epsi)
       - yp_dot * std::sin(epsi)) + (pos_ob_x - s) * (yp_dot * std::cos(epsi) +
       xp_dot * std::sin(epsi))) / ((sub_a * sub_a + tub_a * tub_a) * (uub_a *
       uub_a + vub_a * vub_a)) - wub_a * wub_a * (((2.0 * vel_ob_x * yp_dot *
        std::cos(epsi) - 2.0 * vel_ob_y * xp_dot * std::cos(epsi)) + 2.0 *
       vel_ob_x * xp_dot * std::sin(epsi)) + 2.0 * vel_ob_y * yp_dot * std::sin
      (epsi)) / ((xub_a * xub_a + yub_a * yub_a) * (avb_a * avb_a))) *
    ((vel_ob_x * std::cos(epsi) - xp_dot) + vel_ob_y * std::sin(epsi)) *
    (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y * vel_ob_x) - s *
             vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y * xp_dot * std::
           cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) + s * yp_dot * std::
         cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x * xp_dot * std::
       sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot * std::sin
     (epsi)) / (m * pow(1.0 - bvb_a * bvb_a / ((cvb_a * cvb_a + dvb_a *
    dvb_a) * (evb_a * evb_a + fvb_a * fvb_a)), 1.5) * std::sqrt(gvb_a * gvb_a +
    hvb_a * hvb_a) * pow(ivb_a * ivb_a + jvb_a * jvb_a, 1.5));
  out[6] = (2.0 * cf * steer * ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot *
             std::sin(epsi)) * ((vel_ob_x * std::cos(epsi) - xp_dot) + vel_ob_y *
             std::sin(epsi)) / (m * std::sqrt(1.0 - kvb_a * kvb_a / ((lvb_a *
    lvb_a + mvb_a * mvb_a) * (nvb_a * nvb_a + ovb_a * ovb_a))) * std::sqrt(pvb_a
              * pvb_a + qvb_a * qvb_a) * pow(rvb_a * rvb_a + svb_a *
              svb_a, 1.5)) - cf * steer * (tvb_a * tvb_a * (2.0 * ey - 2.0 *
              pos_ob_y) / (uvb_a * uvb_a * (vvb_a * vvb_a + wvb_a * wvb_a)) -
             2.0 * ((pos_ob_x - s) * ((vel_ob_x - xp_dot * std::cos(epsi)) +
    yp_dot * std::sin(epsi)) + (ey - pos_ob_y) * ((yp_dot * std::cos(epsi) -
    vel_ob_y) + xp_dot * std::sin(epsi))) * ((yp_dot * std::cos(epsi) - vel_ob_y)
              + xp_dot * std::sin(epsi)) / ((xvb_a * xvb_a + yvb_a * yvb_a) *
              (awb_a * awb_a + bwb_a * bwb_a))) * ((vel_ob_x * std::cos(epsi) -
              xp_dot) + vel_ob_y * std::sin(epsi)) * (((((((((((ey * vel_ob_x +
    pos_ob_x * vel_ob_y) - pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot *
    std::cos(epsi)) + pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot *
    std::cos(epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi))
    - pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi))
             + s * xp_dot * std::sin(epsi)) / (m * pow(1.0 - cwb_a *
              cwb_a / ((dwb_a * dwb_a + ewb_a * ewb_a) * (fwb_a * fwb_a + gwb_a *
    gwb_a)), 1.5) * std::sqrt(hwb_a * hwb_a + iwb_a * iwb_a) * pow(jwb_a
              * jwb_a + kwb_a * kwb_a, 1.5))) - cf * steer * (2.0 * ey - 2.0 *
    pos_ob_y) * ((vel_ob_x * std::cos(epsi) - xp_dot) + vel_ob_y * std::sin(epsi))
    * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y * vel_ob_x) - s
               * vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y * xp_dot *
             std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) + s * yp_dot *
           std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x * xp_dot *
         std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot *
       std::sin(epsi)) / (m * std::sqrt(1.0 - lwb_a * lwb_a / ((mwb_a * mwb_a +
    nwb_a * nwb_a) * (owb_a * owb_a + pwb_a * pwb_a))) * pow(qwb_a *
    qwb_a + rwb_a * rwb_a, 1.5) * pow(swb_a * swb_a + twb_a * twb_a, 1.5));
  out[7] = (2.0 * cf * steer * ((vel_ob_x * std::cos(epsi) - xp_dot) + vel_ob_y *
             std::sin(epsi)) * ((yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot *
             std::sin(epsi)) / (m * std::sqrt(1.0 - uwb_a * uwb_a / ((vwb_a *
    vwb_a + wwb_a * wwb_a) * (xwb_a * xwb_a + ywb_a * ywb_a))) * std::sqrt(axb_a
              * axb_a + bxb_a * bxb_a) * pow(cxb_a * cxb_a + dxb_a *
              dxb_a, 1.5)) - cf * steer * (2.0 * ((pos_ob_x - s) * ((vel_ob_x -
    xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi)) + (ey - pos_ob_y) *
              ((yp_dot * std::cos(epsi) - vel_ob_y) + xp_dot * std::sin(epsi))) *
             ((vel_ob_x - xp_dot * std::cos(epsi)) + yp_dot * std::sin(epsi)) /
             ((exb_a * exb_a + fxb_a * fxb_a) * (gxb_a * gxb_a + hxb_a * hxb_a))
             - (2.0 * pos_ob_x - 2.0 * s) * (ixb_a * ixb_a) / (jxb_a * jxb_a *
              (kxb_a * kxb_a + lxb_a * lxb_a))) * ((vel_ob_x * std::cos(epsi) -
              xp_dot) + vel_ob_y * std::sin(epsi)) * (((((((((((ey * vel_ob_x +
    pos_ob_x * vel_ob_y) - pos_ob_y * vel_ob_x) - s * vel_ob_y) - ey * xp_dot *
    std::cos(epsi)) + pos_ob_y * xp_dot * std::cos(epsi)) - pos_ob_x * yp_dot *
    std::cos(epsi)) + s * yp_dot * std::cos(epsi)) + ey * yp_dot * std::sin(epsi))
    - pos_ob_x * xp_dot * std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi))
             + s * xp_dot * std::sin(epsi)) / (m * pow(1.0 - mxb_a *
              mxb_a / ((nxb_a * nxb_a + oxb_a * oxb_a) * (pxb_a * pxb_a + qxb_a *
    qxb_a)), 1.5) * std::sqrt(rxb_a * rxb_a + sxb_a * sxb_a) * pow(txb_a
              * txb_a + uxb_a * uxb_a, 1.5))) + cf * steer * (2.0 * pos_ob_x -
    2.0 * s) * ((vel_ob_x * std::cos(epsi) - xp_dot) + vel_ob_y * std::sin(epsi))
    * (((((((((((ey * vel_ob_x + pos_ob_x * vel_ob_y) - pos_ob_y * vel_ob_x) - s
               * vel_ob_y) - ey * xp_dot * std::cos(epsi)) + pos_ob_y * xp_dot *
             std::cos(epsi)) - pos_ob_x * yp_dot * std::cos(epsi)) + s * yp_dot *
           std::cos(epsi)) + ey * yp_dot * std::sin(epsi)) - pos_ob_x * xp_dot *
         std::sin(epsi)) - pos_ob_y * yp_dot * std::sin(epsi)) + s * xp_dot *
       std::sin(epsi)) / (m * std::sqrt(1.0 - vxb_a * vxb_a / ((wxb_a * wxb_a +
    xxb_a * xxb_a) * (yxb_a * yxb_a + ayb_a * ayb_a))) * pow(byb_a *
    byb_a + cyb_a * cyb_a, 1.5) * pow(dyb_a * dyb_a + eyb_a * eyb_a, 1.5));
}

//
// File trailer for L_f_L_f_h_ang_Ccode.cpp
//
// [EOF]
//
