function [alpha_el,CN_e,alpha_ed,CD_e,alpha_em,CM_e] = experimental_results(airfoil,k,alphabase,A_alpha)

if strcmp(airfoil,'NACA0012') && k == 0.1 && alphabase == deg2rad(5) && A_alpha == deg2rad(10)
    Exp = 1;
elseif strcmp(airfoil,'HH02') && k == 0.1 && alphabase == deg2rad(5) && A_alpha == deg2rad(10)
    Exp = 2;
elseif strcmp(airfoil,'SC1095') && k == 0.1 && alphabase == deg2rad(5) && A_alpha == deg2rad(10)
    Exp = 3;
elseif strcmp(airfoil,'NACA0012') && k == 0.1 && alphabase == deg2rad(10) && A_alpha == deg2rad(10)
    Exp = 4;
elseif strcmp(airfoil,'HH02') && k == 0.1 && alphabase == deg2rad(10) && A_alpha == deg2rad(10)
    Exp = 5;
elseif strcmp(airfoil,'SC1095') && k == 0.1 && alphabase == deg2rad(10) && A_alpha == deg2rad(10)
    Exp = 6;
elseif strcmp(airfoil,'NACA0012') && k == 0.1 && alphabase == deg2rad(15) && A_alpha == deg2rad(10)
    Exp = 7;
elseif strcmp(airfoil,'HH02') && k == 0.1 && alphabase == deg2rad(15) && A_alpha == deg2rad(10)
    Exp = 8;
elseif strcmp(airfoil,'SC1095') && k == 0.1 && alphabase == deg2rad(15) && A_alpha == deg2rad(10)
    Exp = 9;
elseif strcmp(airfoil,'NACA0012') && k == 0.074 && alphabase == deg2rad(2.1) && A_alpha == deg2rad(8.2)
    Exp = 10;
elseif strcmp(airfoil,'NACA0012') && k == 0.075 && alphabase == deg2rad(10.3) && A_alpha == deg2rad(8.1)
    Exp = 11;
elseif strcmp(airfoil,'NACA0012') && k == 0.1 && alphabase == deg2rad(10) && A_alpha == deg2rad(15)
    Exp = 12;
elseif strcmp(airfoil,'NACA0012') && k == 0.1 && alphabase == deg2rad(0) && A_alpha == deg2rad(10)
    Exp = 13;
elseif strcmp(airfoil,'NACA0012') && k == 0.1 && alphabase == deg2rad(0) && A_alpha == deg2rad(15)
    Exp = 14;
else
    Exp = 0;
end

if Exp == 1
    load('experimental_data/CN_NACA0012.mat');
    alpha_el = CN_NACA0012(:,1);
    CN_e = CN_NACA0012(:,2);
    load('experimental_data/CD_NACA0012.mat');
    alpha_ed = CD_NACA0012(:,1);
    CD_e = CD_NACA0012(:,2);
    load('experimental_data/CM_NACA0012_510.mat');
    alpha_em = CM_NACA0012_510(:,1);
    CM_e = CM_NACA0012_510(:,2);
elseif Exp == 2
    load('experimental_data/CN_HH02.mat');
    alpha_el = CN_HH02(:,1);
    CN_e = CN_HH02(:,2);
    load('experimental_data/CD_HH02.mat');
    alpha_ed = CD_HH02(:,1);
    CD_e = CD_HH02(:,2);
    alpha_em = [];
    CM_e = [];
elseif Exp == 3
    load('experimental_data/CN_SC1095.mat');
    alpha_el = CN_SC1095(:,1);
    CN_e = CN_SC1095(:,2);
    load('experimental_data/CD_SC1095.mat');
    alpha_ed = CD_SC1095(:,1);
    CD_e = CD_SC1095(:,2);
    alpha_em = [];
    CM_e = [];
elseif Exp == 4
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
elseif Exp == 5
    load('experimental_data/CN_HH02_1010.mat');
    alpha_el = CN_HH02_1010(:,1);
    CN_e = CN_HH02_1010(:,2);
    load('experimental_data/CD_HH02_1010.mat');
    alpha_ed = CD_HH02_1010(:,1);
    CD_e = CD_HH02_1010(:,2);
    alpha_em = [];
    CM_e = [];
elseif Exp == 6
    load('experimental_data/CN_SC1095_1010.mat');
    alpha_el = CN_SC1095_1010(:,1);
    CN_e = CN_SC1095_1010(:,2);
    load('experimental_data/CD_SC1095_1010.mat');
    alpha_ed = CD_SC1095_1010(:,1);
    CD_e = CD_SC1095_1010(:,2);
    alpha_em = [];
    CM_e = [];
elseif Exp == 7
    load('experimental_data/CN_NACA0012_1510.mat');
    alpha_el = CN_NACA0012_1510(:,1);
    CN_e = CN_NACA0012_1510(:,2);
    load('experimental_data/CD_NACA0012_1510.mat');
    alpha_ed = CD_NACA0012_1510(:,1);
    CD_e = CD_NACA0012_1510(:,2);
    load('experimental_data/CM_NACA0012_1510.mat');
    alpha_em = CM_NACA0012_1510(:,1);
    CM_e = CM_NACA0012_1510(:,2);
elseif Exp == 8
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
elseif Exp == 9
    load('experimental_data/CN_SC1095_1510.mat');
    alpha_el = CN_SC1095_1510(:,1);
    CN_e = CN_SC1095_1510(:,2);
    load('experimental_data/CD_SC1095_1510.mat');
    alpha_ed = CD_SC1095_1510(:,1);
    CD_e = CD_SC1095_1510(:,2);
    alpha_em = [];
    CM_e = [];
elseif Exp == 10
    load('experimental_data/CN_NACA0012_28.mat');
    alpha_el = CN_NACA0012_28(:,1);
    CN_e = CN_NACA0012_28(:,2);
    alpha_ed = [];
    CD_e = [];
    alpha_em = [];
    CM_e = [];
elseif Exp == 11
    load('experimental_data/CN_NACA0012_108.mat');
    alpha_el = CN_NACA0012_108(:,1);
    CN_e = CN_NACA0012_108(:,2);
    alpha_ed = [];
    CD_e = [];
    alpha_em = [];
    CM_e = [];
elseif Exp == 12
    load('experimental_data/Cl_3D_DES_Re_135_10_15sinwt_shift_4_1cycle.mat');
    alpha_el = AOA_ds_Cl_3D_DES_Re_135_10_15sinwt_shift_4_1cycle;
    CN_e = Cl_3D_DES_Re_135_10_15sinwt_shift_4_1cycle;
    alpha_ed = [];
    CD_e = [];
    alpha_em = [];
    CM_e = [];
elseif Exp == 13
    load('experimental_data/Cl_3D_DES_Re_135_0_10sinwt_shift_2_1cycle.mat');
    alpha_el = AOA_ds_Cl_3D_DES_Re_135_0_10sinwt_shift_2_1cycle;
    CN_e = Cl_3D_DES_Re_135_0_10sinwt_shift_2_1cycle;
    alpha_ed = [];
    CD_e = [];
    alpha_em = [];
    CM_e = [];
elseif Exp == 14
    load('experimental_data/Cl_3D_DES_Re_135_0_15sinwt_shift_2_1cycle.mat');
    alpha_el = AOA_ds_Cl_3D_DES_Re_135_0_15sinwt_shift_2_1cycle;
    CN_e = Cl_3D_DES_Re_135_0_15sinwt_shift_2_1cycle;
    alpha_ed = [];
    CD_e = [];
    alpha_em = [];
    CM_e = [];
else
    alpha_el = [];
    CN_e = [];
    alpha_ed = [];
    CD_e = [];
    alpha_em = [];
    CM_e = [];
end

end