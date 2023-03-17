function [alpha_el,CN_e,alpha_ed,CD_e,alpha_em,CM_e] = experimental_results(airfoil,k,M,alphabase,A_alpha,H,phi)

if strcmp(airfoil,'NACA0012') && k == 0.1 && alphabase == deg2rad(5) && A_alpha == deg2rad(10) && H==0 && phi==0 && M<=0.3
    load('experimental_data/CN_NACA0012.mat');
    alpha_el = CN_NACA0012(:,1);
    CN_e = CN_NACA0012(:,2);
    load('experimental_data/CD_NACA0012.mat');
    alpha_ed = CD_NACA0012(:,1);
    CD_e = CD_NACA0012(:,2);
    load('experimental_data/CM_NACA0012_510.mat');
    alpha_em = CM_NACA0012_510(:,1);
    CM_e = CM_NACA0012_510(:,2);
elseif strcmp(airfoil,'HH02') && k == 0.1 && alphabase == deg2rad(5) && A_alpha == deg2rad(10) && H==0 && phi==0 && M<=0.3
    load('experimental_data/CN_HH02.mat');
    alpha_el = CN_HH02(:,1);
    CN_e = CN_HH02(:,2);
    load('experimental_data/CD_HH02.mat');
    alpha_ed = CD_HH02(:,1);
    CD_e = CD_HH02(:,2);
    alpha_em = [];
    CM_e = [];
elseif strcmp(airfoil,'SC1095') && k == 0.1 && alphabase == deg2rad(5) && A_alpha == deg2rad(10) && H==0 && phi==0 && M<=0.3
    load('experimental_data/CN_SC1095.mat');
    alpha_el = CN_SC1095(:,1);
    CN_e = CN_SC1095(:,2);
    load('experimental_data/CD_SC1095.mat');
    alpha_ed = CD_SC1095(:,1);
    CD_e = CD_SC1095(:,2);
    alpha_em = [];
    CM_e = [];
elseif strcmp(airfoil,'NACA0012') && k == 0.1 && alphabase == deg2rad(10) && A_alpha == deg2rad(10) && H==0 && phi==0 && M<=0.3
    load('experimental_data/CN_NACA0012_1010.mat');
    alpha_el = CN_NACA0012_1010(:,1);
    CN_e = CN_NACA0012_1010(:,2);
    alpha_em = [];
    CM_e = [];
    load('experimental_data/CD_NACA0012_1010.mat');
    alpha_ed = CD_NACA0012_1010(:,1);
    CD_e = CD_NACA0012_1010(:,2);
    load('experimental_data/CM_NACA0012_1010.mat');
    alpha_em = CM_NACA0012_1010(:,1);
    CM_e = CM_NACA0012_1010(:,2);
elseif strcmp(airfoil,'HH02') && k == 0.1 && alphabase == deg2rad(10) && A_alpha == deg2rad(10) && H==0 && phi==0 && M<=0.3
    load('experimental_data/CN_HH02_1010.mat');
    alpha_el = CN_HH02_1010(:,1);
    CN_e = CN_HH02_1010(:,2);
    load('experimental_data/CD_HH02_1010.mat');
    alpha_ed = CD_HH02_1010(:,1);
    CD_e = CD_HH02_1010(:,2);
    alpha_em = [];
    CM_e = [];
elseif strcmp(airfoil,'SC1095') && k == 0.1 && alphabase == deg2rad(10) && A_alpha == deg2rad(10) && H==0 && phi==0 && M<=0.3
    load('experimental_data/CN_SC1095_1010.mat');
    alpha_el = CN_SC1095_1010(:,1);
    CN_e = CN_SC1095_1010(:,2);
    load('experimental_data/CD_SC1095_1010.mat');
    alpha_ed = CD_SC1095_1010(:,1);
    CD_e = CD_SC1095_1010(:,2);
    alpha_em = [];
    CM_e = [];
elseif strcmp(airfoil,'NACA0012') && k == 0.1 && alphabase == deg2rad(15) && A_alpha == deg2rad(10) && H==0 && phi==0 && M<=0.3
    load('experimental_data/CN_NACA0012_1510.mat');
    alpha_el = CN_NACA0012_1510(:,1);
    CN_e = CN_NACA0012_1510(:,2);
    load('experimental_data/CD_NACA0012_1510.mat');
    alpha_ed = CD_NACA0012_1510(:,1);
    CD_e = CD_NACA0012_1510(:,2);
    load('experimental_data/CM_NACA0012_1510.mat');
    alpha_em = CM_NACA0012_1510(:,1);
    CM_e = CM_NACA0012_1510(:,2);
elseif strcmp(airfoil,'HH02') && k == 0.1 && alphabase == deg2rad(15) && A_alpha == deg2rad(10) && H==0 && phi==0 && M<=0.3
    load('experimental_data/CN_HH02_1510.mat');
    alpha_el = CN_HH02_1510(:,1);
    CN_e = CN_HH02_1510(:,2);
    alpha_em = [];
    CM_e = [];
    load('experimental_data/CD_HH02_1510.mat');
    alpha_ed = CD_HH02_1510(:,1);
    CD_e = CD_HH02_1510(:,2);
    alpha_em = [];
    CM_e = [];
elseif strcmp(airfoil,'SC1095') && k == 0.1 && alphabase == deg2rad(15) && A_alpha == deg2rad(10) && H==0 && phi==0 && M<=0.3
    load('experimental_data/CN_SC1095_1510.mat');
    alpha_el = CN_SC1095_1510(:,1);
    CN_e = CN_SC1095_1510(:,2);
    load('experimental_data/CD_SC1095_1510.mat');
    alpha_ed = CD_SC1095_1510(:,1);
    CD_e = CD_SC1095_1510(:,2);
    alpha_em = [];
    CM_e = [];
elseif strcmp(airfoil,'NACA0012') && k == 0.074 && alphabase == deg2rad(2.1) && A_alpha == deg2rad(8.2) && H==0 && phi==0 && M<=0.3
    load('experimental_data/CN_NACA0012_28.mat');
    alpha_el = CN_NACA0012_28(:,1);
    CN_e = CN_NACA0012_28(:,2);
    alpha_ed = [];
    CD_e = [];
    alpha_em = [];
    CM_e = [];
elseif strcmp(airfoil,'NACA0012') && k == 0.075 && alphabase == deg2rad(10.3) && A_alpha == deg2rad(8.1) && H==0 && phi==0 && M<=0.3
    load('experimental_data/CN_NACA0012_108.mat');
    alpha_el = CN_NACA0012_108(:,1);
    CN_e = CN_NACA0012_108(:,2);
    alpha_ed = [];
    CD_e = [];
    alpha_em = [];
    CM_e = [];
elseif strcmp(airfoil,'NACA0012') && k == 0.1 && alphabase == deg2rad(10) && A_alpha == deg2rad(15) && H==0 && phi==0 && M<=0.3
    load('experimental_data/Cl_3D_DES_Re_135_10_15sinwt_shift_4_1cycle.mat');
    alpha_el = AOA_ds_Cl_3D_DES_Re_135_10_15sinwt_shift_4_1cycle;
    CN_e = Cl_3D_DES_Re_135_10_15sinwt_shift_4_1cycle;
    alpha_ed = [];
    CD_e = [];
    alpha_em = [];
    CM_e = [];
elseif strcmp(airfoil,'NACA0012') && k == 0.1 && alphabase == deg2rad(0) && A_alpha == deg2rad(10) && H==0 && phi==0 && M<=0.3
    load('experimental_data/Cl_3D_DES_Re_135_0_10sinwt_shift_2_1cycle.mat');
    alpha_el = AOA_ds_Cl_3D_DES_Re_135_0_10sinwt_shift_2_1cycle;
    CN_e = Cl_3D_DES_Re_135_0_10sinwt_shift_2_1cycle;
    alpha_ed = [];
    CD_e = [];
    alpha_em = [];
    CM_e = [];
elseif strcmp(airfoil,'NACA0012') && k == 0.1 && alphabase == deg2rad(0) && A_alpha == deg2rad(15) && H==0 && phi==0 && M<=0.3
    load('experimental_data/Cl_3D_DES_Re_135_0_15sinwt_shift_2_1cycle.mat');
    alpha_el = AOA_ds_Cl_3D_DES_Re_135_0_15sinwt_shift_2_1cycle;
    CN_e = Cl_3D_DES_Re_135_0_15sinwt_shift_2_1cycle;
    alpha_ed = [];
    CD_e = [];
    alpha_em = [];
    CM_e = [];
elseif strcmp(airfoil,'NACA0012') && k == 0.051 && alphabase == deg2rad(8.99) && A_alpha == deg2rad(4.91) && H==0 && phi==0 && M<=0.3
    load('experimental_data/NASA Data/frame_7019.mat',...
        'alpha_exp_cl','alpha_exp_cd','alpha_exp_cm','cl_exp','cd_exp','cm_exp');
    alpha_el = alpha_exp_cl;
    CN_e = cl_exp;
    alpha_ed = alpha_exp_cd;
    CD_e = cd_exp;
    alpha_em = alpha_exp_cm;
    CM_e = cm_exp;
elseif strcmp(airfoil,'NACA0012') && k == 0.100 && alphabase == deg2rad(8.99) && A_alpha == deg2rad(4.91) && H==0 && phi==0 && M<=0.3
    load('experimental_data/NASA Data/frame_7021.mat',...
        'alpha_exp_cl','alpha_exp_cd','alpha_exp_cm','cl_exp','cd_exp','cm_exp');
    alpha_el = alpha_exp_cl;
    CN_e = cl_exp;
    alpha_ed = alpha_exp_cd;
    CD_e = cd_exp;
    alpha_em = alpha_exp_cm;
    CM_e = cm_exp;
elseif strcmp(airfoil,'NACA0012') && k == 0.201 && alphabase == deg2rad(8.99) && A_alpha == deg2rad(4.91) && H==0 && phi==0 && M<=0.3
    load('experimental_data/NASA Data/frame_7023.mat',...
        'alpha_exp_cl','alpha_exp_cd','alpha_exp_cm','cl_exp','cd_exp','cm_exp');
    alpha_el = alpha_exp_cl;
    CN_e = cl_exp;
    alpha_ed = alpha_exp_cd;
    CD_e = cd_exp;
    alpha_em = alpha_exp_cm;
    CM_e = cm_exp;
elseif strcmp(airfoil,'NACA0012') && k == 0.150 && alphabase == deg2rad(8.99) && A_alpha == deg2rad(4.91) && H==0 && phi==0 && M<=0.3
    load('experimental_data/NASA Data/frame_7101.mat',...
        'alpha_exp_cl','alpha_exp_cd','alpha_exp_cm','cl_exp','cd_exp','cm_exp');
    alpha_el = alpha_exp_cl;
    CN_e = cl_exp;
    alpha_ed = alpha_exp_cd;
    CD_e = cd_exp;
    alpha_em = alpha_exp_cm;
    CM_e = cm_exp;
elseif strcmp(airfoil,'NACA0012') && k == 0.025 && alphabase == deg2rad(8.97) && A_alpha == deg2rad(4.90) && H==0 && phi==0 && M>=0.3 && M<0.45
    load('experimental_data/NASA Data/frame_7104.mat',...
        'alpha_exp_cl','alpha_exp_cd','alpha_exp_cm','cl_exp','cd_exp','cm_exp');
    alpha_el = alpha_exp_cl;
    CN_e = cl_exp;
    alpha_ed = alpha_exp_cd;
    CD_e = cd_exp;
    alpha_em = alpha_exp_cm;
    CM_e = cm_exp;
elseif strcmp(airfoil,'NACA0012') && k == 0.025 && alphabase == deg2rad(7.97) && A_alpha == deg2rad(4.90) && H==0 && phi==0 && M>=0.3 && M<0.45
    load('experimental_data/NASA Data/frame_7108.mat',...
        'alpha_exp_cl','alpha_exp_cd','alpha_exp_cm','cl_exp','cd_exp','cm_exp');
    alpha_el = alpha_exp_cl;
    CN_e = cl_exp;
    alpha_ed = alpha_exp_cd;
    CD_e = cd_exp;
    alpha_em = alpha_exp_cm;
    CM_e = cm_exp;
elseif strcmp(airfoil,'NACA0012') && k == 0.100 && alphabase == deg2rad(7.95) && A_alpha == deg2rad(4.90) && H==0 && phi==0 && M>=0.3 && M<0.45
    load('experimental_data/NASA Data/frame_7110.mat',...
        'alpha_exp_cl','alpha_exp_cd','alpha_exp_cm','cl_exp','cd_exp','cm_exp');
    alpha_el = alpha_exp_cl;
    CN_e = cl_exp;
    alpha_ed = alpha_exp_cd;
    CD_e = cd_exp;
    alpha_em = alpha_exp_cm;
    CM_e = cm_exp;
elseif strcmp(airfoil,'NACA0012') && k == 0.199 && alphabase == deg2rad(7.97) && A_alpha == deg2rad(4.91) && H==0 && phi==0 && M>=0.3 && M<0.45
    load('experimental_data/NASA Data/frame_7111.mat',...
        'alpha_exp_cl','alpha_exp_cd','alpha_exp_cm','cl_exp','cd_exp','cm_exp');
    alpha_el = alpha_exp_cl;
    CN_e = cl_exp;
    alpha_ed = alpha_exp_cd;
    CD_e = cd_exp;
    alpha_em = alpha_exp_cm;
    CM_e = cm_exp;
elseif strcmp(airfoil,'NACA0012') && k == 0.025 && alphabase == deg2rad(9.97) && A_alpha == deg2rad(4.90) && H==0 && phi==0 && M>=0.3 && M<0.45
    load('experimental_data/NASA Data/frame_7112.mat',...
        'alpha_exp_cl','alpha_exp_cd','alpha_exp_cm','cl_exp','cd_exp','cm_exp');
    alpha_el = alpha_exp_cl;
    CN_e = cl_exp;
    alpha_ed = alpha_exp_cd;
    CD_e = cd_exp;
    alpha_em = alpha_exp_cm;
    CM_e = cm_exp;
elseif strcmp(airfoil,'NACA0012') && k == 0.099 && alphabase == deg2rad(9.95) && A_alpha == deg2rad(4.91) && H==0 && phi==0 && M>=0.3 && M<0.45
    load('experimental_data/NASA Data/frame_7113.mat',...
        'alpha_exp_cl','alpha_exp_cd','alpha_exp_cm','cl_exp','cd_exp','cm_exp');
    alpha_el = alpha_exp_cl;
    CN_e = cl_exp;
    alpha_ed = alpha_exp_cd;
    CD_e = cd_exp;
    alpha_em = alpha_exp_cm;
    CM_e = cm_exp;
elseif strcmp(airfoil,'NACA0012') && k == 0.199 && alphabase == deg2rad(9.91) && A_alpha == deg2rad(4.93) && H==0 && phi==0 && M>=0.3 && M<0.45
    load('experimental_data/NASA Data/frame_7114.mat',...
        'alpha_exp_cl','alpha_exp_cd','alpha_exp_cm','cl_exp','cd_exp','cm_exp');
    alpha_el = alpha_exp_cl;
    CN_e = cl_exp;
    alpha_ed = alpha_exp_cd;
    CD_e = cd_exp;
    alpha_em = alpha_exp_cm;
    CM_e = cm_exp;
elseif strcmp(airfoil,'NACA0012') && k == 0.025 && alphabase == deg2rad(10.92) && A_alpha == deg2rad(4.90) && H==0 && phi==0 && M>=0.3 && M<0.45
    load('experimental_data/NASA Data/frame_7117.mat',...
        'alpha_exp_cl','alpha_exp_cd','alpha_exp_cm','cl_exp','cd_exp','cm_exp');
    alpha_el = alpha_exp_cl;
    CN_e = cl_exp;
    alpha_ed = alpha_exp_cd;
    CD_e = cd_exp;
    alpha_em = alpha_exp_cm;
    CM_e = cm_exp;
elseif strcmp(airfoil,'NACA0012') && k == 0.050 && alphabase == deg2rad(10.92) && A_alpha == deg2rad(4.90) && H==0 && phi==0 && M>=0.3 && M<0.45
    load('experimental_data/NASA Data/frame_7118.mat',...
        'alpha_exp_cl','alpha_exp_cd','alpha_exp_cm','cl_exp','cd_exp','cm_exp');
    alpha_el = alpha_exp_cl;
    CN_e = cl_exp;
    alpha_ed = alpha_exp_cd;
    CD_e = cd_exp;
    alpha_em = alpha_exp_cm;
    CM_e = cm_exp;
elseif strcmp(airfoil,'NACA0012') && k == 0.099 && alphabase == deg2rad(10.92) && A_alpha == deg2rad(4.89) && H==0 && phi==0 && M>=0.3 && M<0.45
    load('experimental_data/NASA Data/frame_7119.mat',...
        'alpha_exp_cl','alpha_exp_cd','alpha_exp_cm','cl_exp','cd_exp','cm_exp');
    alpha_el = alpha_exp_cl;
    CN_e = cl_exp;
    alpha_ed = alpha_exp_cd;
    CD_e = cd_exp;
    alpha_em = alpha_exp_cm;
    CM_e = cm_exp;
elseif strcmp(airfoil,'NACA0012') && k == 0.149 && alphabase == deg2rad(10.96) && A_alpha == deg2rad(4.89) && H==0 && phi==0 && M>=0.3 && M<0.45
    load('experimental_data/NASA Data/frame_7120.mat',...
        'alpha_exp_cl','alpha_exp_cd','alpha_exp_cm','cl_exp','cd_exp','cm_exp');
    alpha_el = alpha_exp_cl;
    CN_e = cl_exp;
    alpha_ed = alpha_exp_cd;
    CD_e = cd_exp;
    alpha_em = alpha_exp_cm;
    CM_e = cm_exp;
elseif strcmp(airfoil,'NACA0012') && k == 0.198 && alphabase == deg2rad(10.88) && A_alpha == deg2rad(4.90) && H==0 && phi==0 && M>=0.3 && M<0.45
    load('experimental_data/NASA Data/frame_7121.mat',...
        'alpha_exp_cl','alpha_exp_cd','alpha_exp_cm','cl_exp','cd_exp','cm_exp');
    alpha_el = alpha_exp_cl;
    CN_e = cl_exp;
    alpha_ed = alpha_exp_cd;
    CD_e = cd_exp;
    alpha_em = alpha_exp_cm;
    CM_e = cm_exp;
elseif strcmp(airfoil,'NACA0012') && k == 0.025 && alphabase == deg2rad(11.92) && A_alpha == deg2rad(4.88) && H==0 && phi==0 && M>=0.3 && M<0.45
    load('experimental_data/NASA Data/frame_7200.mat',...
        'alpha_exp_cl','alpha_exp_cd','alpha_exp_cm','cl_exp','cd_exp','cm_exp');
    alpha_el = alpha_exp_cl;
    CN_e = cl_exp;
    alpha_ed = alpha_exp_cd;
    CD_e = cd_exp;
    alpha_em = alpha_exp_cm;
    CM_e = cm_exp;
elseif strcmp(airfoil,'NACA0012') && k == 0.050 && alphabase == deg2rad(11.93) && A_alpha == deg2rad(4.89) && H==0 && phi==0 && M>=0.3 && M<0.42
    load('experimental_data/NASA Data/frame_7202.mat',...
        'alpha_exp_cl','alpha_exp_cd','alpha_exp_cm','cl_exp','cd_exp','cm_exp');
    alpha_el = alpha_exp_cl;
    CN_e = cl_exp;
    alpha_ed = alpha_exp_cd;
    CD_e = cd_exp;
    alpha_em = alpha_exp_cm;
    CM_e = cm_exp;
elseif strcmp(airfoil,'NACA0012') && k == 0.099 && alphabase == deg2rad(11.92) && A_alpha == deg2rad(4.90) && H==0 && phi==0 && M>=0.3 && M<0.42
    load('experimental_data/NASA Data/frame_7205.mat',...
        'alpha_exp_cl','alpha_exp_cd','alpha_exp_cm','cl_exp','cd_exp','cm_exp');
    alpha_el = alpha_exp_cl;
    CN_e = cl_exp;
    alpha_ed = alpha_exp_cd;
    CD_e = cd_exp;
    alpha_em = alpha_exp_cm;
    CM_e = cm_exp;
elseif strcmp(airfoil,'NACA0012') && k == 0.198 && alphabase == deg2rad(11.88) && A_alpha == deg2rad(4.90) && H==0 && phi==0 && M>=0.3 && M<0.42
    load('experimental_data/NASA Data/frame_7207.mat',...
        'alpha_exp_cl','alpha_exp_cd','alpha_exp_cm','cl_exp','cd_exp','cm_exp');
    alpha_el = alpha_exp_cl;
    CN_e = cl_exp;
    alpha_ed = alpha_exp_cd;
    CD_e = cd_exp;
    alpha_em = alpha_exp_cm;
    CM_e = cm_exp;
elseif strcmp(airfoil,'NACA0012') && k == 0.149 && alphabase == deg2rad(8.71) && A_alpha == deg2rad(4.90) && H==0 && phi==0 && M>=0.3 && M<0.42
    load('experimental_data/NASA Data/frame_7212.mat',...
        'alpha_exp_cl','alpha_exp_cd','alpha_exp_cm','cl_exp','cd_exp','cm_exp');
    alpha_el = alpha_exp_cl;
    CN_e = cl_exp;
    alpha_ed = alpha_exp_cd;
    CD_e = cd_exp;
    alpha_em = alpha_exp_cm;
    CM_e = cm_exp;
elseif strcmp(airfoil,'NACA0012') && k == 0.099 && alphabase == deg2rad(8.71) && A_alpha == deg2rad(4.90) && H==0 && phi==0 && M>=0.3 && M<0.42
    load('experimental_data/NASA Data/frame_7214.mat',...
        'alpha_exp_cl','alpha_exp_cd','alpha_exp_cm','cl_exp','cd_exp','cm_exp');
    alpha_el = alpha_exp_cl;
    CN_e = cl_exp;
    alpha_ed = alpha_exp_cd;
    CD_e = cd_exp;
    alpha_em = alpha_exp_cm;
    CM_e = cm_exp;
elseif strcmp(airfoil,'NACA0012') && k == 0.050 && alphabase == deg2rad(8.71) && A_alpha == deg2rad(4.89) && H==0 && phi==0 && M>=0.3 && M<0.42
    load('experimental_data/NASA Data/frame_7216.mat',...
        'alpha_exp_cl','alpha_exp_cd','alpha_exp_cm','cl_exp','cd_exp','cm_exp');
    alpha_el = alpha_exp_cl;
    CN_e = cl_exp;
    alpha_ed = alpha_exp_cd;
    CD_e = cd_exp;
    alpha_em = alpha_exp_cm;
    CM_e = cm_exp;
elseif strcmp(airfoil,'NACA0012') && k == 0.051 && alphabase == deg2rad(9.76) && A_alpha == deg2rad(4.95) && H==0 && phi==0 && M>=0.3 && M<0.42
    load('experimental_data/NASA Data/frame_7222.mat',...
        'alpha_exp_cl','alpha_exp_cd','alpha_exp_cm','cl_exp','cd_exp','cm_exp');
    alpha_el = alpha_exp_cl;
    CN_e = cl_exp;
    alpha_ed = alpha_exp_cd;
    CD_e = cd_exp;
    alpha_em = alpha_exp_cm;
    CM_e = cm_exp;
elseif strcmp(airfoil,'NACA0012') && k == 0.151 && alphabase == deg2rad(9.93) && A_alpha == deg2rad(4.90) && H==0 && phi==0 && M<=0.3
    load('experimental_data/NASA Data/frame_7300.mat',...
        'alpha_exp_cl','alpha_exp_cd','alpha_exp_cm','cl_exp','cd_exp','cm_exp');
    alpha_el = alpha_exp_cl;
    CN_e = cl_exp;
    alpha_ed = alpha_exp_cd;
    CD_e = cd_exp;
    alpha_em = alpha_exp_cm;
    CM_e = cm_exp;
else
    alpha_el = [];
    CN_e = [];
    alpha_ed = [];
    CD_e = [];
    alpha_em = [];
    CM_e = [];
end

end