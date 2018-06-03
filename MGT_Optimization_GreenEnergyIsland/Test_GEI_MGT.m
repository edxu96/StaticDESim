% Title: Thermo-economic Optimization of Distributed Energy System in Green Energy Island.
% Based on the theory of nolinear equality and inequality constraints.
% Method: Genetic Algorithm within MATLAB Global Optimization Toolbox.
% Version: 1.6, 2018.6.3, Jie Xu.
% SubTitle: Test Calculation of MGT, 1.2.
% 1. Micro Gas Turbine (MGT)
% 2. Aqueous Lithium-Bromine Single-Effect Absorption Chiller (AC_ALB)
% 3. R123 Organic Recycle Cycle (ORC_R123)
% 4. R410a Heat Pump (HP_R410a)
% 5. R134a Vapor Compression Chiller (VCC_R134a)
clear; clc;
%% 1. Decision variable of MGT. ----------------------------------------------------------------------------------------------------------------------------------------
x(1) = 390.7*1000;       % p_2,       Pa,   Outlet Pressure of Air Compressor
x(2) = 0.813;            % Eta_MGTac,       Isentropic Efficiency of Air Compressor
x(3) = 0.910;            % Eta_MGTt,        Isentropic Efficiency of Gas Turbine
x(4) = 839.38;           % T_3,       K,    Outlet Temp of Air from Regenerator
x(5) = 1205;             % T_4,       K,    Inlet Temp of Gas Turbine
x(6) = 0.03;             % Ms_MGT,    kg/s, Fluid Rate of Saturated Steam from HRSG
x(7) = 60 * 1000;        % Wp_MGT,    W,    Net Power from Micro Gas Turbine
%% 2. General Constants ------------------------------------------------------------------------------------------------------------------------------------------------
R = 8.314472;                                     % J/(mol*K), Universial Gas Constant
p_0 = 101.325 * 1000;                             % Pa, Pressure of atmosphere.
T_0 = 25 + 273.15;                                % K, Temperature of atmosphere.
T_0H = 35 + 273.15;                               % K, Acceptable Highest Temp of atmosphere.
T_0L = 15 + 273.15;                               % K, Acceptable Lowest Temp of atmosphere.
% 2.1 Constants for MGT
LHV = 500000 * 1000;                              % 纯甲烷燃料的低热值
DELTA_p_cc = 0.05; DELTA_p_aRE = 0.05; ...
DELTA_p_gRE = 0.03; DELTA_p_HRSG = 0.05;          % 压力损失
p_1 = 101.3*1000; T_1 = 288.15;                   % 环境状态
p_8 = 500*1000; p_8p = p_8; p_9 = p_8; p_7 = p_1; % 余热锅炉水蒸气侧的热力学参数
T_8 = 343.15; T_8p = 410; T_9 = 425;              % 余热锅炉水蒸气侧的热力学参数
Z_f = 4 * 10^(-9);                                % 单位能量的燃料价格
A_MGT = 837.68; B_MGT = -143.22; C_MGT = 491.79;  % 回归分析统计的价格方程系数
AlphaC_MGT = 0.73;                                % 费用估价系数
Ka_MGT = 1.4;                                     % 空气的等熵系数
Kg_MGT = 1.33;                                    % 燃气的等熵系数
CRF = 0.182;                                      % 投资回收系数
PhiM = 1.06;                                      % 维护因子
c_pa = 1004;                                      % J/kg/K, air specific heat capacity at constant pressure
c_pg = 1170;                                      % J/kg/K, 实际燃气 specific heat capacity at constant pressure
c_pw = 4200;                                      % J/kg/K, water specific heat capacity at constant pressure
ETA_CC = 0.98877;
EtaG_MGTt = 0.98;
% 2.1 Decision Variables in MGT.
      p_2 = x(1);               % Pa,   Outlet Pressure of Air Compressor
Eta_MGTac = x(2);               %       Isentropic Efficiency of Air Compressor
 Eta_MGTt = x(3);               %       Isentropic Efficiency of Gas Turbine
      T_3 = x(4);               % K,    Outlet Temp of Air from Regenerator
      T_4 = x(5);               % K,    Inlet Temp of Gas Turbine
   Ms_MGT = x(6);               % kg/s, Fluid Rate of Saturated Steam from HRSG
   Wp_MGT = x(7);               % W,    Net Power from Micro Gas Turbine
%% 3. Pre-defined Condition. -------------------------------------------------------------------------------------------------------------------------------------------
 p_s = 500 * 1000;                                 % Pa, Pressure of supplying steam.
T_cw = 4 + 273.15; p_cw = p_0;                    % Supplying cooling water.
  Ms = 0.05;
 Mcw = 0.05;
  We = 100000;
%% 4. Mathematical Model Preparation  ----------------------------------------------------------------------------------------------------------------------------------
h_w0 = CoolProp.PropsSI('H', 'T', T_0, 'P', p_0, 'Water');
s_w0 = CoolProp.PropsSI('S', 'T', T_0, 'P', p_0, 'Water');
ex_ph0 = 0;
T_s = CoolProp.PropsSI('T', 'P', p_s, 'Q', 1, 'Water');
h_s = CoolProp.PropsSI('H', 'P', p_s, 'Q', 1, 'Water');
s_s = CoolProp.PropsSI('S', 'P', p_s, 'Q', 1, 'Water');
E_s = (h_s - h_w0) - T_0 * (s_s - s_w0);
h_cw = CoolProp.PropsSI('H', 'T', T_cw, 'P', p_cw, 'Water');
s_cw = CoolProp.PropsSI('S', 'T', T_cw, 'P', p_cw, 'Water');
E_cw = (h_cw - h_w0) + T_0 * (s_cw - s_w0);
h_w0H = CoolProp.PropsSI('H', 'T', T_0H, 'P', p_0, 'Water');
s_w0H = CoolProp.PropsSI('S', 'T', T_0H, 'P', p_0, 'Water');
E_w0H = (h_w0H - h_w0) + T_0 * (s_w0H - s_w0);
h_w0L = CoolProp.PropsSI('H', 'T', T_0L, 'P', p_0, 'Water');
s_w0L = CoolProp.PropsSI('S', 'T', T_0L, 'P', p_0, 'Water');
E_w0L = (h_w0L - h_w0) + T_0 * (s_w0L - s_w0);
h_8 = CoolProp.PropsSI('H', 'T', T_8, 'P', p_8, 'Water');
h_8p = CoolProp.PropsSI('H', 'T', T_8p, 'P', p_8p, 'Water');
% Calculate the h, s, E_ph of air, fuel and water in status 0.
h_a0 = CoolProp.PropsSI('H', 'T', T_0, 'P', p_0, 'Air');
s_a0 = CoolProp.PropsSI('S', 'T', T_0, 'P', p_0, 'Air');
h_f0 = CoolProp.PropsSI('H', 'T', T_0, 'P', p_0, 'Methane');
s_f0 = CoolProp.PropsSI('S', 'T', T_0, 'P', p_0, 'Methane');
h_f1 = h_f0; s_f1 = s_f0;
ex_phf1 = (h_f1 - h_f0) - T_0 * (s_f1 - h_f0);
%% 5 Mathematical Model of Micro Gas Turbine (MGT) ---------------------------------------------------------------------------------------------------------------------
% 5.1 Mathematical Model of Brayton Cycle.
T_2 = T_1 * (1 + 1/Eta_MGTac * ...
      ((p_2/p_1)^((Ka_MGT-1)/Ka_MGT) - 1));                        % (1)  AC
p_3 = p_2 * (1 - DELTA_p_aRE);                                     % (7)  RE
p_4 = p_3 * (1 - DELTA_p_cc);                                      % (5)  CC
p_6 = p_7 / (1 - DELTA_p_HRSG);                                    % (14) HRSG
p_5 = p_6 / (1 - DELTA_p_gRE);                                     % (8)  RE
T_5 = T_4 * (1 - Eta_MGTt * (1 - (p_4./p_5)^((1-Kg_MGT)/Kg_MGT))); % (9)  GT
H = (c_pg * T_4 - ETA_CC * LHV) / (c_pa * T_3 - ETA_CC * LHV);     % (3)(4) CC, H = m_a / m_g;
m_g = Wp_MGT / (c_pg*(T_4-T_5) - c_pa*(T_2-T_1)*H);                % (2) (10)
m_a = H * m_g;                                                     % H = m_a / m_g;
m_f = m_g - m_a;                                                   % (3)  CC
T_6 = T_5 - m_a * c_pa * (T_3 - T_2) / (m_g * c_pg);               % (6)  RE
h_9 = h_s;
T_7p = T_6 - Ms_MGT * (h_9 - h_8p) / (m_g * c_pg);                 % (12) HRSG
T_7 = T_6 - Ms_MGT * (h_9 - h_8) / (m_g * c_pg);                   % (13) HRSG
Wp_MGTac = m_a * c_pa * (T_2 - T_1);                               % (2)  AC
Wp_MGTt = Wp_MGT + Wp_MGTac;                                       % (11) GT
We_MGT = Wp_MGT * EtaG_MGTt;
%% 5.2 Calculate the h and s of mix gas before combustion.
h_m0 = m_a/m_g * h_a0 + m_f/m_g * h_f0;
s_m0 = m_a/m_g * s_a0 + m_f/m_g * s_f0;
% Calculate the h, s, E_ph of air in status 1.
h_a1 = CoolProp.PropsSI('H', 'T', T_1, 'P', p_1, 'Air');
s_a1 = CoolProp.PropsSI('S', 'T', T_1, 'P', p_1, 'Air');
ex_ph1 = (h_a1 - h_a0) - T_0 * (s_a1 - s_a0);
% Calculate the h, s, E_ph of air in status 2.
h_f2 = h_f0; s_f2 = s_f0;
h_a2 = CoolProp.PropsSI('H', 'T', T_2, 'P', p_2, 'Air');
s_a2 = CoolProp.PropsSI('S', 'T', T_2, 'P', p_2, 'Air');
h_m2 = m_a/m_g * h_a2 + m_f/m_g * h_f2;
s_m2 = m_a/m_g * s_a2 + m_f/m_g * s_f2;
ex_ph2 = (h_a2 - h_a0) - T_0 * (s_a2 - s_a0);
% Calculate the h, s, E_ph of gas in status 3.
h_a3 = CoolProp.PropsSI('H', 'T', T_3, 'P', p_3, 'Air');
s_a3 = CoolProp.PropsSI('S', 'T', T_3, 'P', p_3, 'Air');
h_m3 = m_a/m_g * h_a3 + m_f/m_g * h_f0;
s_m3 = m_a/m_g * s_a3 + m_f/m_g * s_f0;
ex_ph3 = (h_a3 - h_a0) - T_0 * (s_a3 - s_a0);
%% 5.3 Calculate the components of waste gas and water. -----------------------------------------------------------------------------------------------------
y_O = 0.21;                                    % Volume fraction of Oxygen in Air.
n_aC = 2 / y_O;                                % Needed air for 1 mol fuel combustion
n_gC = 1 + 2 + n_aC - 2;                       % Waste gas from combustion with just needed air
y_C = (1 + n_aC * 0.0004) / n_gC;
y_H = 2 / n_gC;
y_N = 1 - y_C - y_H;
M_gC = y_C * 44 + y_H * 18 + y_N * 28;
m_aC = n_aC * 29 / 16;
m_gC = n_gC * M_gC / 16;
r_water = (2257.2 * 1000);                             % J/kg, latent heat of vaporization
HHV = LHV + 1000/16 * n_gC * y_H * 18/1000 * r_water;  % J/kg
ex_f = HHV;                                            % J/kg
% Calculate the h, s, E_ph of gas in status 0.
m_ax = m_a - m_f/16*2/0.21*29;
h_H0 = CoolProp.PropsSI('H', 'T', T_0, 'P', p_0, 'Water');
s_H0 = CoolProp.PropsSI('S', 'T', T_0, 'P', p_0, 'Water');
h_C0 = CoolProp.PropsSI('H', 'T', T_0, 'P', p_0, 'CarbonDioxide');
s_C0 = CoolProp.PropsSI('S', 'T', T_0, 'P', p_0, 'CarbonDioxide');
h_N0 = CoolProp.PropsSI('H', 'T', T_0, 'P', p_0, 'Nitrogen');
s_N0 = CoolProp.PropsSI('S', 'T', T_0, 'P', p_0, 'Nitrogen');
h_g0 = m_ax/m_g * h_a0 + (m_a-m_ax)/m_g * y_H*18/M_gC * h_H0 + ...
       (m_a-m_ax)/m_g * y_C*44/M_gC * h_C0 + (m_a-m_ax)/m_g * y_N*28/M_gC * h_N0;
s_g0 = m_ax/m_g * s_a0 + (m_a-m_ax)/m_g * y_H*18/M_gC * s_H0 + ...
       (m_a-m_ax)/m_g * y_C*44/M_gC * s_C0 + (m_a-m_ax)/m_g * y_N*28/M_gC * s_N0;
% Calculate the h, s, E_ph of gas in status 4.
h_a4 = CoolProp.PropsSI('H', 'T', T_4, 'P', p_4, 'Air');
s_a4 = CoolProp.PropsSI('S', 'T', T_4, 'P', p_4, 'Air');
h_H4 = CoolProp.PropsSI('H', 'T', T_4, 'P', p_4, 'Water');
s_H4 = CoolProp.PropsSI('S', 'T', T_4, 'P', p_4, 'Water');
h_C4 = CoolProp.PropsSI('H', 'T', T_4, 'P', p_4, 'CarbonDioxide');
s_C4 = CoolProp.PropsSI('S', 'T', T_4, 'P', p_4, 'CarbonDioxide');
h_N4 = CoolProp.PropsSI('H', 'T', T_4, 'P', p_4, 'Nitrogen');
s_N4 = CoolProp.PropsSI('S', 'T', T_4, 'P', p_4, 'Nitrogen');
h_g4 = m_ax/m_g * h_a4 + (m_a-m_ax)/m_g * y_H*18/M_gC * h_H4 + ...
       (m_a-m_ax)/m_g * y_C*44/M_gC * h_C4 + (m_a-m_ax)/m_g * y_N*28/M_gC * h_N4;
s_g4 = m_ax/m_g * s_a4 + (m_a-m_ax)/m_g * y_H*18/M_gC * s_H4 + ...
       (m_a-m_ax)/m_g * y_C*44/M_gC * s_C4 + (m_a-m_ax)/m_g * y_N*28/M_gC * s_N4;
ex_ph4 = (h_g4 - h_g0) - T_0 * (s_g4 - s_g0);
% Calculate the h, s, E_ph of gas in status 5.
h_a5 = CoolProp.PropsSI('H', 'T', T_5, 'P', p_5, 'Air');
s_a5 = CoolProp.PropsSI('S', 'T', T_5, 'P', p_5, 'Air');
h_H5 = CoolProp.PropsSI('H', 'T', T_5, 'P', p_5, 'Water');
s_H5 = CoolProp.PropsSI('S', 'T', T_5, 'P', p_5, 'Water');
h_C5 = CoolProp.PropsSI('H', 'T', T_5, 'P', p_5, 'CarbonDioxide');
s_C5 = CoolProp.PropsSI('S', 'T', T_5, 'P', p_5, 'CarbonDioxide');
h_N5 = CoolProp.PropsSI('H', 'T', T_5, 'P', p_5, 'Nitrogen');
s_N5 = CoolProp.PropsSI('S', 'T', T_5, 'P', p_5, 'Nitrogen');
h_g5 = m_ax/m_g * h_a5 + (m_a-m_ax)/m_g * y_H*18/M_gC * h_H5 + ...
       (m_a-m_ax)/m_g * y_C*44/M_gC * h_C5 + (m_a-m_ax)/m_g * y_N*28/M_gC * h_N5;
s_g5 = m_ax/m_g * s_a5 + (m_a-m_ax)/m_g * y_H*18/M_gC * s_H5 + ...
       (m_a-m_ax)/m_g * y_C*44/M_gC * s_C5 + (m_a-m_ax)/m_g * y_N*28/M_gC * s_N5;
ex_ph5 = (h_g5 - h_g0) - T_0 * (s_g5 - s_g0);
% Calculate the h, s, E_ph of gas in status 6.
h_a6 = CoolProp.PropsSI('H', 'T', T_6, 'P', p_6, 'Air');
s_a6 = CoolProp.PropsSI('S', 'T', T_6, 'P', p_6, 'Air');
h_H6 = CoolProp.PropsSI('H', 'T', T_6, 'P', p_6, 'Water');
s_H6 = CoolProp.PropsSI('S', 'T', T_6, 'P', p_6, 'Water');
h_C6 = CoolProp.PropsSI('H', 'T', T_6, 'P', p_6, 'CarbonDioxide');
s_C6 = CoolProp.PropsSI('S', 'T', T_6, 'P', p_6, 'CarbonDioxide');
h_N6 = CoolProp.PropsSI('H', 'T', T_6, 'P', p_6, 'Nitrogen');
s_N6 = CoolProp.PropsSI('S', 'T', T_6, 'P', p_6, 'Nitrogen');
h_g6 = m_ax/m_g * h_a6 + (m_a-m_ax)/m_g * y_H*18/M_gC * h_H6 + ...
       (m_a-m_ax)/m_g * y_C*44/M_gC * h_C6 + (m_a-m_ax)/m_g * y_N*28/M_gC * h_N6;
s_g6 = m_ax/m_g * s_a6 + (m_a-m_ax)/m_g * y_H*18/M_gC * s_H6 + ...
       (m_a-m_ax)/m_g * y_C*44/M_gC * s_C6 + (m_a-m_ax)/m_g * y_N*28/M_gC * s_N6;
ex_ph6 = (h_g6 - h_g0) - T_0 * (s_g6 - s_g0);
% Calculate the h, s, E_ph of gas in status 7.
h_a7 = CoolProp.PropsSI('H', 'T', T_7, 'P', p_7, 'Air');
s_a7 = CoolProp.PropsSI('S', 'T', T_7, 'P', p_7, 'Air');
h_H7 = CoolProp.PropsSI('H', 'T', T_7, 'P', p_7, 'Water');
s_H7 = CoolProp.PropsSI('S', 'T', T_7, 'P', p_7, 'Water');
h_C7 = CoolProp.PropsSI('H', 'T', T_7, 'P', p_7, 'CarbonDioxide');
s_C7 = CoolProp.PropsSI('S', 'T', T_7, 'P', p_7, 'CarbonDioxide');
h_N7 = CoolProp.PropsSI('H', 'T', T_7, 'P', p_7, 'Nitrogen');
s_N7 = CoolProp.PropsSI('S', 'T', T_7, 'P', p_7, 'Nitrogen');
h_g7 = m_ax/m_g * h_a7 + (m_a-m_ax)/m_g * y_H*18/M_gC * h_H7 + ...
       (m_a-m_ax)/m_g * y_C*44/M_gC * h_C7 + (m_a-m_ax)/m_g * y_N*28/M_gC * h_N7;
s_g7 = m_ax/m_g * s_a7 + (m_a-m_ax)/m_g * y_H*18/M_gC * s_H7 + ...
       (m_a-m_ax)/m_g * y_C*44/M_gC * s_C7 + (m_a-m_ax)/m_g * y_N*28/M_gC * s_N7;
ex_ph7 = (h_g7 - h_g0) - T_0 * (s_g7 - s_g0);
% Set the temp of waste gas emitted to atmosphere. ??????
% T_wg =
% Calculate the h, s, E_ph of supplying steam in status 8.
s_8 = CoolProp.PropsSI('S', 'T', T_8, 'P', p_8, 'Water');
ex_ph8 = (h_8 - h_w0) - T_0 * (s_8 - s_w0);
% Calculate the h, s, E_ph of supplying steam in status 9.
s_9 = s_s;
ex_ph9 = E_s;
%% 5.1.4 Exergy Analysis of MGT.
Ed_MGTac = Wp_MGTac + m_a * ex_ph1 - m_a * ex_ph2;
Ed_MGTre = (ex_ph2 * m_a + ex_ph5 * m_g) - (ex_ph3 * m_a + ex_ph6 * m_g);
Ed_MGTcc = (ex_ph3 * m_a + ex_f * m_f) - ex_ph4 * m_g;
Ed_MGTt = ex_ph4 * m_g - (Wp_MGTt + ex_ph5 * m_g);
Ed_MGThrsg = (ex_ph8 * Ms_MGT + ex_ph6 * m_g) - (ex_ph9 * Ms_MGT + ex_ph7 * m_g);
% Exergy efficiency of components of MGT.
PhiE_MGTac = (ex_ph2 - ex_ph1) * m_a / Wp_MGTac;
PhiE_MGTre = (ex_ph3 * m_a + ex_ph6 * m_g) / (ex_ph2 * m_a + ex_ph5 * m_g);
PhiE_MGTcc = ex_ph4 * m_g / (ex_ph3 * m_a + ex_f * m_f);
PhiE_MGTt = Wp_MGT / ((ex_ph4 - ex_ph5) * m_g);
PhiE_MGThrsg = (ex_ph9 * Ms_MGT + ex_ph7 * m_g) / (ex_ph8 * Ms_MGT + ex_ph6 * m_g);
%% Plot the Diagram of MGT. -----------------------------------------------------------------------------------------------------------------------
% 6.1 Plot the T-s Diagram of Brayton Closed Cycle.
figure(1);
CycleT_MGT = [T_1 T_2 T_3 T_4 T_5 T_6 T_1]; CycleS_MGT = [s_m0 s_m2 s_m3 s_g4 s_g5 s_g6 s_m0];
plot(CycleS_MGT,CycleT_MGT);
ylabel('温度 T (K)'); xlabel('熵 s (J/kg/K)');
title('微型燃气轮机闭合循环温熵图 T-s Diagram of MGT Closed Cycle', ...
      'FontSize',20,'FontWeight','bold');
text(s_m0, T_1, '1', 'FontSize', 15);
text(s_m2, T_2, '2', 'FontSize', 15);
text(s_m3, T_3, '3', 'FontSize', 15);
text(s_g4, T_4, '4', 'FontSize', 15);
text(s_g5, T_5, '5', 'FontSize', 15);
text(s_g6, T_6, '6', 'FontSize', 15);
% 6.2 Plot the Temp Change Diagram of DER System in GEI 
figure(2);
ChgT_MGTy = [T_1, T_2, T_3, T_4, T_5, T_6, T_7];
ChgT_MGTx = [1,2,3,4,5,6,7];
ChgT_MGT = plot(ChgT_MGTx, ChgT_MGTy, '-*');
xlim([0, 8]);
title('微型燃气轮机温度变化图 Temp Change Diagram of MGT', ...
      'FontSize',20,'FontWeight','bold');
ylabel('温度 T (K)');
% 6.3 Plot the Pressure Change Diagram of DER System in GEI
figure(3);
ChgP_MGTy = [p_1, p_2, p_3, p_4, p_5, p_6, p_7];
ChgP_MGTx = [1,2,3,4,5,6,7];
ChgP_MGT = plot(ChgP_MGTx, ChgP_MGTy, '-*');
xlim([0, 8]);
title('微型燃气轮机压力变化图 Pressure Change Diagram of MGT', ...
      'FontSize',20,'FontWeight','bold');
ylabel('压力 p (Pa)');
% 6.4 Plot the bar chart of exergy efficiency of equipments in MGT.
figure(4);
BarEe_MGTx = categorical({'PhiE_M_G_T_a_c','PhiE_M_G_T_r_e','PhiE_M_G_T_c_c', ...
                         'PhiE_M_G_T_t','PhiE_M_G_T_h_r_s_g'});
BarEe_MGTy = [PhiE_MGTac, PhiE_MGTre, PhiE_MGTcc, PhiE_MGTt, PhiE_MGThrsg];
BarEe_MGT = bar(BarEe_MGTx, BarEe_MGTy, 0.5);
ylim([0 1.1]);
title('微型燃气轮机组件㶲效率条形图 Bar Chart of Exergy Efficiency of MGT', ...
      'FontSize',20,'FontWeight','bold');
% 6.5 Plot the bar chart of exergy damage of equipments in MGT.
figure(5);
BarEd_MGTx = categorical({'Ed_M_G_T_a_c','Ed_M_G_T_r_e','Ed_M_G_T_c_c', ...
                          'Ed_M_G_T_t','Ed_M_G_T_h_r_s_g'});
BarEd_MGTy = [Ed_MGTac, Ed_MGTre, Ed_MGTcc, Ed_MGTt, Ed_MGThrsg];
BarEd_MGT = bar(BarEd_MGTx, BarEd_MGTy, 0.5);
title('微型燃气轮机组件㶲损失条形图 Bar Chart of Exergy Damage of MGT', ...
      'FontSize',20,'FontWeight','bold');
% 6.6 Exergy amount of components of MGT. ---------------------------------------------------------------------------------------------------------------------------
figure(6);
BAR_MGTyE = [0,                                             ex_ph0*m_a, ex_f*m_f, Wp_MGTac+(ex_ph5-ex_ph6)*m_g+ex_ph8*Ms_MGT, 0; ...
             Ed_MGTac,                                      ex_ph2*m_a, ex_f*m_f, (ex_ph5-ex_ph6)*m_g+ex_ph8*Ms_MGT,          0; ...
             Ed_MGTac+Ed_MGTre,                             ex_ph3*m_a, ex_f*m_f, ex_ph8*Ms_MGT,                              0; ...
             Ed_MGTac+Ed_MGTre+Ed_MGTcc,                    ex_ph4*m_g, 0,        ex_ph8*Ms_MGT,                              0; ...
             Ed_MGTac+Ed_MGTre+Ed_MGTcc+Ed_MGTt,            ex_ph5*m_g, 0,        ex_ph8*Ms_MGT,                              Wp_MGTt; ...
             Ed_MGTac+Ed_MGTre+Ed_MGTcc+Ed_MGTt,            ex_ph6*m_g, 0,        ex_ph8*Ms_MGT,                              Wp_MGTt+(ex_ph5-ex_ph6)*m_g; ...
             Ed_MGTac+Ed_MGTre+Ed_MGTcc+Ed_MGTt+Ed_MGThrsg, ex_ph7*m_g, 0,        0,                                          Wp_MGTt+(ex_ph5-ex_ph6)*m_g+E_s*Ms_MGT];
bar(BAR_MGTyE,'stacked')
title('微型燃气轮机㶲流图 Exergy Flow of MGT', ...
      'FontSize',20,'FontWeight','bold');
legend('㶲损 Ed', '物理㶲 Eph', '化学㶲 Ech', '输入㶲 Ein', '输出㶲 Eout');
