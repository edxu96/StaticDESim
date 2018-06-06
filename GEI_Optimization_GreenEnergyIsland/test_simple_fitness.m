% Title: Thermo-economic Optimization of Distributed Energy System in Green Energy Island.
% Based on the theory of nolinear equality and inequality constraints.
% Method: Genetic Algorithm within MATLAB Global Optimization Toolbox.
% Version: 2.0, 2018.6.5, Jie Xu.
% Test
% SubTitle: Fitness Function for Exergy Damage Minimization
% 1. Micro Gas Turbine (MGT) with a fixed T_7 = 400K in Northern China, 1.0.
% 2. Aqueous Lithium-Bromine Single-Effect Absorption Chiller (AC_ALB)
% 3. R123 Organic Recycle Cycle (ORC_R123)
% 4. R410a Heat Pump (HP_R410a)
% 5. R134a Vapor Compression Chiller (VCC_R134a)
%% Input decision variable.
% Decision variable of MGT.
x(1) = 633186.25;
x(2) = 0.837;
x(3) = 0.917;
x(4) = 769.395751953125;
x(5) = 1300;
% Decision variable of AC_ALB.
x(6) = 4 + 273.15;
x(7) = 30 + 273.15;
x(8) = 0.5322;
x(11) = 0.6711;
x(12) = 0.05;
x(13) = 0.8;
% Decision variable of ORC_R123.
x(14) = 0.05;
x(15) = 200 + 273.15;
x(16) = 5 * 101.325 * 1000;
x(17) = 80 + 273.15;
% Decision variable of HP_R410a.
x(18) = 30 + 273.15;
x(19) = 70 + 273.15;
% Decision variable of VCC_R134a.
x(20) = 30 + 273.15;
x(21) = 70 + 273.15;
%% 3. Decision Variables. -------------------------------------------------------------------------------
% 3.1 Decision Variables in MGT.
      p_2 = x(1);               % Pa,   Outlet Pressure of Air Compressor
Eta_MGTac = x(2);               %       Isentropic Efficiency of Air Compressor
 Eta_MGTt = x(3);               %       Isentropic Efficiency of Gas Turbine
      T_3 = x(4);               % K,    Outlet Temp of Air from Regenerator
      T_4 = x(5);               % K,    Inlet Temp of Gas Turbine
   Ms_MGT = x(6);               % kg/s, Fluid Rate of Saturated Steam from HRSG
   Wp_MGT = x(7);               % W,    Net Power from Micro Gas Turbine
% 3.2 Decision Variables in AC_ALB.
   T_AC10 = x(8);               % K,    Outlet temperature of Evaporator.
    T_AC8 = x(9);               % K,    Outlet temperature of Condenser.
      y_1 = x(10);              %       Mass Fraction of Outlet Solution from Absorber.
      y_4 = x(11);              %       Mass Fraction of Outlet Solution from Desorber.
    M_AC1 = x(12);              % kg/s, Fluid Rate Outlet Solution from Absorber.
  Eta_shx = x(13);              %       Solution Heat Exchanger Ratio.
% 3.3 Decision Variables in ORC_R123.
    M_ORC = x(14);              % kg/s, Fluid Rate.
   T_ORC3 = x(15);              % K,    Inlet temperature of turbine
DeltaP_ORCp = x(16);              % K,    Outlet pressure of pump / inlet pressure of turbine
   T_ORC1 = x(17);              % K,    Outlet temperature of condenser.
% 3.4 Decision Variables in HP_R410a.
   T_HP1 = x(18);               % K,    Outlet temperature of evaporator.
   T_HP3 = x(19);               % K,    Outlet temperature of condenser.
%  T_HP2 = x(20);               % K,    Outlet temperature of compressor.
%   M_HP = x(20);               % kg/s, Fluid rate of HP_R410a.
% 3.5 Decision Variables in VCC_R134a.
  T_VCC1 = x(20);               % K,    Outlet temperature of evaporator.
  T_VCC3 = x(21);               % K,    Outlet temperature of condenser.
% T_VCC2 = x(24);               % K,    Outlet temperature of compressor.
%  M_VCC = x(23);               % kg/s, Fluid Rate.
%% 2. Pre-defined Condition. ----------------------------------------------------------------------------
 p_s = 500 * 1000;                                 % Pa, Pressure of supplying steam.
T_cw = 4 + 273.15; p_cw = P_0;                    % Supplying cooling water.
  Ms = 0.05;
 Mcw = 0.05;
  We = 100000;
 T_7 = 400;
% 1. General Constants ---------------------------------------------------------------------------------
R = 8.314472;                                     % J/(mol*K), Universial Gas Constant
P_0 = 101.325 * 1000;                             % Pa, Pressure of atmosphere.
T_0 = 25 + 273.15;                                % K, Temperature of atmosphere.
T_0H = 35 + 273.15;                               % K, Acceptable Highest Temp of atmosphere.
T_0L = 15 + 273.15;                               % K, Acceptable Lowest Temp of atmosphere.
c_pw = 4200;                                      % J/kg/K, water specific heat capacity at constant pressure
Cp_a = 1004;                                      % J/kg/K, air specific heat capacity at constant pressure
Cp_g = 1170;                                      % J/kg/K, 实际燃气 specific heat capacity at constant pressure
K = 10000;                                        % W / m^2 / K, Thermal conductivity.
R_w = (2257.2 * 1000);                            % J/kg, latent heat of vaporization
% 1.1 Constants for MGT
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
ETA_CC = 0.98877;
% 1.2 Constants for AC_ALB
Z_A = 100;                                        % RMB / m^2, Cost rate of area of heat transfer.
Z_W = 3.25E-3;                                    % RMB / kg, Cost rate of cooling water.
Z_cw = 0.5;                                       % RMB / kW*h, Profit rate of supplying cooling load.
Z_E = 0.6;                                        % RMB / kW*h, Cost of supplying electricity for pump.
Q_wAC = 44.1 * 1000;                              % W,   ????
Z_wAC = 235550;                                   % RMB, ????
ALPHA_wAC = 0.73;                                 %      ????
EtaS_ACv = 0.99;
EtaS_ACp = 0.99;
% 1.3 Constants for ORC_R123.
Rg_R123 = 0.0544 * 1000;                          % J/(mol*K), Gas Constant.
M_R123 = 152.93 / 1000;                           % kg/mol, Molecular Weight of R123.
Tc_R123 = 456.83;                                 % K, Critical Temp of R123.
Pc_R123 = 3668.0 * 1000;                          % Pa, Critical Pressure of R123.
EtaS_ORCs = 0.99;                                 % Isentropic efficiency of turbine.
EtaS_ORCp = 0.99;                                 % Efficiency of pump.
DeltaP_ORCc = 100 * 1000;                         % Pa, Pressure drop in condenser.
DeltaP_ORCe = 100 * 1000;                         % Pa, Pressure drop in evaporator.
OMEGA_R123 = 0.281922497036;                      % Acentric Factor, R123.
KAPPA_R123 = 0.37464 + (1.54226 - ...             % Dependent on OMEGA(working substance),
        0.26992 * OMEGA_R123) * OMEGA_R123;       % Temp-independent parameter in PR-EOS
aTc_R123 = 0.457235529 * (Rg_R123 * ...
           Tc_R123)^2 ./ Pc_R123;                 % Critical Point Restriction "a(Tc_R123)"
b_R123 = 0.077796074 * Rg_R123 * ...              % m^3/mol, Critical Point Restriction "b",
         Tc_R123 ./ Pc_R123;                      % Temp-independent parameter in PR-EOS
Q_wORC = 23.6 * 1000;                             % W,   ????
Z_wORC = 385600;                                  % RMB, ????
ALPHA_wORC = 0.73;                                %      ????
EtaG_MGTt = 0.9; EtaG_ORCt = 0.9;                 % Efficiency of generator.
%% 1.4 Constants for HP_R410a.
R = 8.314472;                                 % J/(mol*K), Universial Gas Constant
m_R410a = 72.58 / 1000;                       % kg / mol , Molar Mass
Rg_R410a = 0.11455 * 1000;                    % J/(K*kg) , Gas Constant - R134a
Tc_R410a = 345.28;                            % K , temperature in Critical Point.
Pc_R410a = 4926.1 * 1000;                     % Pa, pressure in Critical Point.
Dc_R410a = 488.90;                            % kg/m^3, critical density.
Tb_R410a = -26.06 + 273.15;                   % K , Boiling point at one atmosphere.
EtaS_HPp = 0.99;                               % Isentropic efficiency of pump.
EtaS_HPt = 0.99;                               % Isentropic efficiency of turbine.
DeltaP_HPc = 100 * 1000;                      % Pa, Pressure drop in condenser.
DeltaP_HPe = 100 * 1000;                      % Pa, Pressure drop in evaporator.
Omega_R410a = 0.296;                          % Acentric Factor.
Kappa_R410a = 0.37464 + ...                   % Dependent on Omega_R410a(working substance),
        (1.54226 - 0.26992 * Omega_R410a) * Omega_R410a;  % Temperature-independent parameter in PR-EOS
aTc_R410a = 0.457235529 * (Rg_R410a * Tc_R410a)^2 ./ Pc_R410a;  % Critical Point Restriction "a(Tc_R410a)"
b_R410a = 0.077796074 * Rg_R410a * Tc_R410a ./ Pc_R410a;        % m^3/mol, Critical Point Restriction "b_R410a",
                                                          % Temperature-independent parameter in PR-EOS
%% 1.5 Constants for VCC_R134a.
R = 8.314472;                % J/(mol*K), Universial Gas Constant
m_R134a = 102.03 / 1000;           % kg / mol , Molar Mass
Rg_R134a = 0.0815 * 1000;         % J/(K*kg) , Gas Constant - R134a
Tc_R134a = 374.23;                % K  , temperature in Critical Point.
Pc_R134a = 4060.3 * 1000;         % Pa , pressure in Critical Point.
Tb_R134a = -26.06 + 273.15;     % K  , Boiling point at one atmosphere.
ETA_ps = 0.7;                                 % Isentropic efficiency of pump.
ETA_ts = 0.7;                                 % Isentropic efficiency of turbine.
ETA_v = 0.8;                                  % Efficiency of throttle valve.
ETA_p = 0.8;                                  % Efficiency of compressor(p).
ETA_e = 0.9;                                  % Efficiency of evaporator.
ETA_c = 0.9;                                  % Efficiency of condenser.
DeltaP_C = 100 * 1000;                        % Pa, Pressure drop in condenser.
DeltaP_E = 100 * 1000;                        % Pa, Pressure drop in evaporator.
Omega_R134a = 0.332;                                % Acentric Factor.
Kappa_R134a = 0.37464 + ...                         % Dependent on Omega_R134a(working substance),
        (1.54226 - 0.26992 * Omega_R134a) * Omega_R134a;  % Temperature-independent parameter in PR-EOS
aTc_R134a = 0.457235529 * (Rg_R134a * Tc_R134a)^2 ./ Pc_R134a;    % Critical Point Restriction "a(Tc_R134a)"
b_R134a = 0.077796074 * Rg_R134a * Tc_R134a ./ Pc_R134a;           % m^3/mol, Critical Point Restriction "b_R134a",
                                              % Temperature-independent parameter in PR-EOS
%% 4. Refrence State. ----------------------------------------------------------
H_w0 = CoolProp.PropsSI('H', 'T', T_0, 'P', P_0, 'Water');
S_w0 = CoolProp.PropsSI('S', 'T', T_0, 'P', P_0, 'Water');
E_w0 = 0;
T_s = CoolProp.PropsSI('T', 'P', p_s, 'Q', 1, 'Water');
H_s = CoolProp.PropsSI('H', 'P', p_s, 'Q', 1, 'Water');
S_s = CoolProp.PropsSI('S', 'P', p_s, 'Q', 1, 'Water');
E_s = (H_s - H_w0) + T_0 * (S_s - S_w0);
H_cw = CoolProp.PropsSI('H', 'T', T_cw, 'P', p_cw, 'Water');
S_Cw = CoolProp.PropsSI('S', 'T', T_cw, 'P', p_cw, 'Water');
E_cw = (H_cw - H_w0) + T_0 * (S_Cw - S_w0);
H_w0H = CoolProp.PropsSI('H', 'T', T_0H, 'P', P_0, 'Water');
S_w0H = CoolProp.PropsSI('S', 'T', T_0H, 'P', P_0, 'Water');
E_w0H = (H_w0H - H_w0) + T_0 * (S_w0H - S_w0);
H_w0L = CoolProp.PropsSI('H', 'T', T_0L, 'P', P_0, 'Water');
S_w0L = CoolProp.PropsSI('S', 'T', T_0L, 'P', P_0, 'Water');
E_w0L = (H_w0L - H_w0) + T_0 * (S_w0L - S_w0);
% Calculate the h, s, E_ph of air, fuel and water in status 0.
H_a0 = CoolProp.PropsSI('H', 'T', T_0, 'P', P_0, 'Air');
S_a0 = CoolProp.PropsSI('S', 'T', T_0, 'P', P_0, 'Air');
H_f0 = CoolProp.PropsSI('H', 'T', T_0, 'P', P_0, 'Methane');
s_f0 = CoolProp.PropsSI('S', 'T', T_0, 'P', P_0, 'Methane');
H_f1 = H_f0; s_f1 = s_f0;
ex_phf1 = (H_f1 - H_f0) - T_0 * (s_f1 - H_f0);
%% 5.1 Mathematical Model of Micro Gas Turbine (MGT) ----------------------------------------------------
% 5.1.1 Mathematical Model of Brayton Cycle.
T_2 = T_1 * (1 + 1/Eta_MGTac * ...
      ((p_2/p_1)^((Ka_MGT-1)/Ka_MGT) - 1));                        % (1)  AC
p_3 = p_2 * (1 - DELTA_p_aRE);                                     % (7)  RE
p_4 = p_3 * (1 - DELTA_p_cc);                                      % (5)  CC
p_6 = p_7 / (1 - DELTA_p_HRSG);                                    % (14) HRSG
p_5 = p_6 / (1 - DELTA_p_gRE);                                     % (8)  RE
T_5 = T_4 * (1 - Eta_MGTt * (1 - (p_4./p_5)^((1-Kg_MGT)/Kg_MGT))); % (9)  GT
H = (Cp_g * T_4 - ETA_CC * LHV) / (Cp_a * T_3 - ETA_CC * LHV);     % (3)(4) CC, H = Ma / Mg;
Mg = Wp_MGT / (Cp_g*(T_4-T_5) - Cp_a*(T_2-T_1)*H);                % (2) (10)
Ma = H * Mg;                                                     % H = Ma / Mg;
Mf_MGT = Mg - Ma;                                                   % (3)  CC
T_6 = T_5 - Ma * Cp_a * (T_3 - T_2) / (Mg * Cp_g);               % (6)  RE
H_MGT9 = H_s;
H_MGT8 = CoolProp.PropsSI('H', 'T', 343.15, 'P', 500*1000, 'Water');
H_8p = CoolProp.PropsSI('H', 'T', 410, 'P', 500*1000, 'Water');
Ms_MGT = (T_6 - T_7) * (Mg * Cp_g) / (H_MGT9 - H_MGT8);                 % (13) HRSG
T_7p = T_6 - Ms_MGT * (H_MGT9 - H_8p) / (Mg * Cp_g);                 % (12) HRSG
Wp_MGTac = Ma * Cp_a * (T_2 - T_1);                               % (2)  AC
Wp_MGTt = Wp_MGT + Wp_MGTac;                                       % (11) GT
We_MGT = Wp_MGT * EtaG_MGTt;
%% 5.1.2 Calculate the h and s of air and gas before combustion.
% Calculate the h, s of air&fuel, gas in status 0.
H_m0 = Ma/Mg * H_a0 + Mf_MGT/Mg * H_f0;
S_m0 = Ma/Mg * S_a0 + Mf_MGT/Mg * s_f0;
Ma_x = Ma - Mf_MGT/16*2/0.21*29;
H_H0 = CoolProp.PropsSI('H', 'T', T_0, 'P', P_0, 'Water');
S_H0 = CoolProp.PropsSI('S', 'T', T_0, 'P', P_0, 'Water');
H_C0 = CoolProp.PropsSI('H', 'T', T_0, 'P', P_0, 'CarbonDioxide');
S_C0 = CoolProp.PropsSI('S', 'T', T_0, 'P', P_0, 'CarbonDioxide');
H_N0 = CoolProp.PropsSI('H', 'T', T_0, 'P', P_0, 'Nitrogen');
S_N0 = CoolProp.PropsSI('S', 'T', T_0, 'P', P_0, 'Nitrogen');
H_g0 = Ma_x/Mg * H_a0 + (Ma-Ma_x)/Mg * y_H*18/MgC * H_H0 + ...
       (Ma-Ma_x)/Mg * y_C*44/MgC * H_C0 + (Ma-Ma_x)/Mg * y_N*28/MgC * H_N0;
S_g0 = Ma_x/Mg * S_a0 + (Ma-Ma_x)/Mg * y_H*18/MgC * S_H0 + ...
       (Ma-Ma_x)/Mg * y_C*44/MgC * S_C0 + (Ma-Ma_x)/Mg * y_N*28/MgC * S_N0;
% Calculate the h, s, E_ph of air in status 1.
H_a1 = CoolProp.PropsSI('H', 'T', T_1, 'P', p_1, 'Air');
S_a1 = CoolProp.PropsSI('S', 'T', T_1, 'P', p_1, 'Air');
E_MGT1 = (H_a1 - H_a0) - T_0 * (S_a1 - S_a0);
% Calculate the h, s, E_ph of air in status 2.
H_a2 = CoolProp.PropsSI('H', 'T', T_2, 'P', p_2, 'Air');
S_a2 = CoolProp.PropsSI('S', 'T', T_2, 'P', p_2, 'Air');
E_MGT2 = (H_a2 - H_a0) - T_0 * (S_a2 - S_a0);
% Calculate the h, s, E_ph of gas in status 3.
H_a3 = CoolProp.PropsSI('H', 'T', T_3, 'P', p_3, 'Air');
S_a3 = CoolProp.PropsSI('S', 'T', T_3, 'P', p_3, 'Air');
H_m3 = Ma/Mg * H_a3 + Mf_MGT/Mg * H_f0;
s_m3 = Ma/Mg * S_a3 + Mf_MGT/Mg * s_f0;
E_MGT3 = (H_m3 - H_m0) - T_0 * (s_m3 - S_m0);
%% 5.1.3 Calculate the components of waste gas.
y_O = 0.21;                                    % Volume fraction of Oxygen in Air.
n_aC = 2 / y_O;                                % Needed air for 1 mol fuel combustion
n_gC = 1 + 2 + n_aC - 2;                       % Waste gas from combustion with just needed air
y_C = (1 + n_aC * 0.0004) / n_gC;
y_H = 2 / n_gC;
y_N = 1 - y_C - y_H;
MgC = y_C * 44 + y_H * 18 + y_N * 28;
M_AC = n_aC * 29 / 16;
MgC = n_gC * MgC / 16;
HHV = LHV + 1000/16 * n_gC * y_H * 18/1000 * R_w; % J/kg
E_f = HHV;                                            % J/kg
% Calculate the h, s, E_ph of gas in status 4.
H_a4 = CoolProp.PropsSI('H', 'T', T_4, 'P', p_4, 'Air');
S_a4 = CoolProp.PropsSI('S', 'T', T_4, 'P', p_4, 'Air');
H_H4 = CoolProp.PropsSI('H', 'T', T_4, 'P', p_4, 'Water');
S_H4 = CoolProp.PropsSI('S', 'T', T_4, 'P', p_4, 'Water');
H_C4 = CoolProp.PropsSI('H', 'T', T_4, 'P', p_4, 'CarbonDioxide');
S_C4 = CoolProp.PropsSI('S', 'T', T_4, 'P', p_4, 'CarbonDioxide');
H_N4 = CoolProp.PropsSI('H', 'T', T_4, 'P', p_4, 'Nitrogen');
S_N4 = CoolProp.PropsSI('S', 'T', T_4, 'P', p_4, 'Nitrogen');
H_g4 = Ma_x/Mg * H_a4 + (Ma-Ma_x)/Mg * y_H*18/MgC * H_H4 + ...
       (Ma-Ma_x)/Mg * y_C*44/MgC * H_C4 + (Ma-Ma_x)/Mg * y_N*28/MgC * H_N4;
S_g4 = Ma_x/Mg * S_a4 + (Ma-Ma_x)/Mg * y_H*18/MgC * S_H4 + ...
       (Ma-Ma_x)/Mg * y_C*44/MgC * S_C4 + (Ma-Ma_x)/Mg * y_N*28/MgC * S_N4;
E_MGT4 = (H_g4 - H_g0) - T_0 * (S_g4 - S_g0);
% Calculate the h, s, E_ph of gas in status 5.
H_a5 = CoolProp.PropsSI('H', 'T', T_5, 'P', p_5, 'Air');
S_a5 = CoolProp.PropsSI('S', 'T', T_5, 'P', p_5, 'Air');
H_H5 = CoolProp.PropsSI('H', 'T', T_5, 'P', p_5, 'Water');
S_H5 = CoolProp.PropsSI('S', 'T', T_5, 'P', p_5, 'Water');
H_C5 = CoolProp.PropsSI('H', 'T', T_5, 'P', p_5, 'CarbonDioxide');
S_C5 = CoolProp.PropsSI('S', 'T', T_5, 'P', p_5, 'CarbonDioxide');
H_N5 = CoolProp.PropsSI('H', 'T', T_5, 'P', p_5, 'Nitrogen');
S_N5 = CoolProp.PropsSI('S', 'T', T_5, 'P', p_5, 'Nitrogen');
H_g5 = Ma_x/Mg * H_a5 + (Ma-Ma_x)/Mg * y_H*18/MgC * H_H5 + ...
       (Ma-Ma_x)/Mg * y_C*44/MgC * H_C5 + (Ma-Ma_x)/Mg * y_N*28/MgC * H_N5;
S_g5 = Ma_x/Mg * S_a5 + (Ma-Ma_x)/Mg * y_H*18/MgC * S_H5 + ...
       (Ma-Ma_x)/Mg * y_C*44/MgC * S_C5 + (Ma-Ma_x)/Mg * y_N*28/MgC * S_N5;
E_MGT5 = (H_g5 - H_g0) - T_0 * (S_g5 - S_g0);
% Calculate tHe h, s, E_ph of gas in status 6.
H_a6 = CoolProp.PropsSI('H', 'T', T_6, 'P', p_6, 'Air');
S_a6 = CoolProp.PropsSI('S', 'T', T_6, 'P', p_6, 'Air');
H_H6 = CoolProp.PropsSI('H', 'T', T_6, 'P', p_6, 'Water');
S_H6 = CoolProp.PropsSI('S', 'T', T_6, 'P', p_6, 'Water');
H_C6 = CoolProp.PropsSI('H', 'T', T_6, 'P', p_6, 'CarbonDioxide');
S_C6 = CoolProp.PropsSI('S', 'T', T_6, 'P', p_6, 'CarbonDioxide');
H_N6 = CoolProp.PropsSI('H', 'T', T_6, 'P', p_6, 'Nitrogen');
S_N6 = CoolProp.PropsSI('S', 'T', T_6, 'P', p_6, 'Nitrogen');
H_g7 = Ma_x/Mg * H_a6 + (Ma-Ma_x)/Mg * y_H*18/MgC * H_H6 + ...
       (Ma-Ma_x)/Mg * y_C*44/MgC * H_C6 + (Ma-Ma_x)/Mg * y_N*28/MgC * H_N6;
S_g6 = Ma_x/Mg * S_a6 + (Ma-Ma_x)/Mg * y_H*18/MgC * S_H6 + ...
       (Ma-Ma_x)/Mg * y_C*44/MgC * S_C6 + (Ma-Ma_x)/Mg * y_N*28/MgC * S_N6;
E_MGT6 = (H_g7 - H_g0) - T_0 * (S_g6 - S_g0);
% Calculate the h, s, E_ph of gas in status 7.
H_a7 = CoolProp.PropsSI('H', 'T', T_7, 'P', p_7, 'Air');
S_a7 = CoolProp.PropsSI('S', 'T', T_7, 'P', p_7, 'Air');
H_H7 = CoolProp.PropsSI('H', 'T', T_7, 'P', p_7, 'Water');
S_H7 = CoolProp.PropsSI('S', 'T', T_7, 'P', p_7, 'Water');
H_C7 = CoolProp.PropsSI('H', 'T', T_7, 'P', p_7, 'CarbonDioxide');
S_C7 = CoolProp.PropsSI('S', 'T', T_7, 'P', p_7, 'CarbonDioxide');
H_N7 = CoolProp.PropsSI('H', 'T', T_7, 'P', p_7, 'Nitrogen');
S_N7 = CoolProp.PropsSI('S', 'T', T_7, 'P', p_7, 'Nitrogen');
H_g7 = Ma_x/Mg * H_a7 + (Ma-Ma_x)/Mg * y_H*18/MgC * H_H7 + ...
       (Ma-Ma_x)/Mg * y_C*44/MgC * H_C7 + (Ma-Ma_x)/Mg * y_N*28/MgC * H_N7;
S_g7 = Ma_x/Mg * S_a7 + (Ma-Ma_x)/Mg * y_H*18/MgC * S_H7 + ...
       (Ma-Ma_x)/Mg * y_C*44/MgC * S_C7 + (Ma-Ma_x)/Mg * y_N*28/MgC * S_N7;
E_MGT7 = (H_g7 - H_g0) - T_0 * (S_g7 - S_g0);
%% 5.1.4 Calculate the s, E_ph of supplying steam.
% Calculate the H, s, E_ph of supplying steam in status 8.
S_MGT8 = CoolProp.PropsSI('S', 'T', T_8, 'P', p_8, 'Water');
E_MGT8 = (H_MGT8 - H_w0) - T_0 * (S_MGT8 - S_w0);
% Calculate the H, s, E_ph of supplying steam in status 9.
S_MGT9 = S_s;
E_MGT9 = (H_MGT9 - H_w0) - T_0 * (S_MGT9 - S_w0);
%% 5.1.5 Exergy Analysis of MGT.
Ed_MGTac = Ma * E_MGT1 - Ma * E_MGT2;
Ed_MGTre = (E_MGT2 * Ma + E_MGT5 * Mg) - (E_MGT3 * Ma + E_MGT6 * Mg);
Ed_MGTcc = (E_MGT3 * Ma + E_f * Mf_MGT) - E_MGT4 * Mg;
Ed_MGTt = E_MGT4 * Mg - (Wp_MGT + E_MGT5 * Mg);
Ed_MGTre = (E_MGT2 * Ma + E_MGT5 * Mg) - (E_MGT3 * Ma + E_MGT6 * Mg);
Ed_MGThrsg = (E_MGT8 * Ms_MGT + E_MGT6 * Mg) - (E_MGT9 * Ms_MGT + E_MGT7 * Mg);
% Exergy efficiency of components of MGT.
PhiE_MGTcc = E_MGT4 * Mg / (E_MGT3 * Ma + E_f * Mf_MGT);
PhiE_MGTt = Wp_MGT / ((E_MGT4 - E_MGT5) * Mg);
PhiE_MGTac = (E_MGT2 - E_MGT1) * Ma / Wp_MGTac;
PhiE_MGTre = (E_MGT3 * Ma + E_MGT6 * Mg) / (E_MGT2 * Ma + E_MGT5 * Mg);
PhiE_MGThrsg = (E_MGT9 * Ms_MGT + E_MGT7 * Mg) / (E_MGT8 * Ms_MGT + E_MGT6 * Mg);
%% 5.2 Mathematical Model of Aqueous Lithium-Bromide Absorption Chiller (AC_ALB)-------------------------
% 5.2.1 Status 10 of AC_ALB.
P_AC10 = CoolProp.PropsSI('P','T', T_AC10, 'Q', 1, 'Water');
P_L = P_AC10;
S_AC10 = CoolProp.PropsSI('S', 'T', T_AC10, 'Q', 1, 'Water');
H_AC10 = CoolProp.PropsSI('H', 'T', T_AC10, 'Q', 1, 'Water');
E_AC10 = (H_AC10 - H_w0) - T_0 * (S_AC10 - S_w0);
% 5.2.2 Status 8 of AC_ALB.
P_AC8 = CoolProp.PropsSI('P','T', T_AC8, 'Q', 0, 'Water');
P_H = P_AC8;
S_AC8 = CoolProp.PropsSI('S', 'T', T_AC8, 'Q', 0, 'Water');
H_AC8 = CoolProp.PropsSI('H', 'T', T_AC8, 'Q', 0, 'Water');
E_AC8 = (H_AC8 - H_w0) - T_0 * (S_AC8 - S_w0);
% 5.2.3 Status 9 of AC_ALB.
S_AC9 = S_AC8 / EtaS_ACv;
P_AC9 = P_AC10;
T_AC9 = T_AC10;
H_AC9 = CoolProp.PropsSI('H', 'S', S_AC9, 'T', T_AC9, 'Water');
E_AC9 = (H_AC9 - H_w0) - T_0 * (S_AC9 - S_w0);
% 5.2.4 Status 1 of AC_ALB.
P_AC1 = P_L;
T_1wC = CoolProp.PropsSI('T', 'P', P_AC1, 'Q', 0, 'Water') - 273.15;
a0 = -2.00755; a1 = 0.16976; a2 = -3.13336E-3; a3 = 1.97668E-5;
b0 = 124.937;  b1 = -7.7162; b2 = 0.152286;    b3 = -7.9509E-4;
T_AC1 = T_1wC * (a0 + a1 * (y_1*100) + a2 * (y_1*100)^2 + a3 * (y_1*100)^3) + ...
      b0 + b1 * (y_1*100) + b2 * (y_1*100)^2 + b3 * (y_1*100)^3 + 273.15;
y_1str = num2str(y_1);
S_AC1 = CoolProp.PropsSI('S', 'T', T_AC1, 'Q', 0, strcat('INCOMP::LiBr[',y_1str,']'));
H_AC1 = CoolProp.PropsSI('H', 'T', T_AC1, 'Q', 0, strcat('INCOMP::LiBr[',y_1str,']'));
RHO_1 = CoolProp.PropsSI('D', 'T', T_AC1, 'P', P_AC1, strcat('INCOMP::LiBr[',y_1str,']'));
v_1 = 1 / RHO_1;
S_sT0 = CoolProp.PropsSI('S', 'T', T_0, 'P', P_0, strcat('INCOMP::LiBr[',y_1str,']'));
H_sT0 = CoolProp.PropsSI('H', 'T', T_0, 'P', P_0, strcat('INCOMP::LiBr[',y_1str,']'));
E_AC1 = (H_AC1 - H_sT0) - T_0 * (S_AC1 - S_sT0);
% 5.2.5 Status 4 of AC_ALB.
P_AC4 = P_H;
M_AC2 = M_AC1; y_2 = y_1;
M_AC3 = M_AC2; y_3 = y_2;
M_AC4 = M_AC3 * y_3 / y_4;
T_4wC = CoolProp.PropsSI('T', 'P', P_AC4, 'Q', 0, 'Water') - 273.15;
T_AC4 = T_4wC * (a0 + a1 * (y_4*100) + a2 * (y_4*100)^2 + a3 * (y_4*100)^3) + ...
        b0 + b1 * (y_4*100) + b2 * (y_4*100)^2 + b3 * (y_4*100)^3 + 273.15;
y_4str = num2str(y_4);
S_AC4 = CoolProp.PropsSI('S', 'T', T_AC4, 'Q', 0, strcat('INCOMP::LiBr[',y_4str,']'));
H_AC4 = CoolProp.PropsSI('H', 'T', T_AC4, 'Q', 0, strcat('INCOMP::LiBr[',y_4str,']'));
S_sD0 = CoolProp.PropsSI('S', 'T', T_0, 'P', P_0, strcat('INCOMP::LiBr[',y_4str,']'));
H_sD0 = CoolProp.PropsSI('H', 'T', T_0, 'P', P_0, strcat('INCOMP::LiBr[',y_4str,']'));
E_AC4 = (H_AC4 - H_sD0) - T_0 * (S_AC4 - S_sD0);
% 5.2.6 Status 2 of AC_ALB.
S_AC2 = S_AC1;
P_AC2 = P_H;
H_AC2 = H_AC1 + v_1 * (P_AC2 - P_AC1);
We_ACp = M_AC2 * H_AC2 - M_AC1 * H_AC1;
y_2str = num2str(y_2);
T_AC2 = CoolProp.PropsSI('T', 'P', P_AC2, 'H', H_AC2, strcat('INCOMP::LiBr[',y_2str,']'));
E_AC2 = (H_AC2 - H_sT0) - T_0 * (S_AC2 - S_sT0);
% 5.2.7 Status 5 of AC_ALB.
M_AC5 = M_AC4;
y_5 = y_4;
P_AC5 = P_AC4;
T_AC5 = T_AC4 - (T_AC4 - T_AC2) * Eta_shx;
y_5str = num2str(y_5);
S_AC5 = CoolProp.PropsSI('S', 'T', T_AC5, 'P', P_AC5, strcat('INCOMP::LiBr[',y_5str,']'));
H_AC5 = CoolProp.PropsSI('H', 'T', T_AC5, 'P', P_AC5, strcat('INCOMP::LiBr[',y_5str,']'));
E_AC5 = (H_AC5 - H_sD0) - T_0 * (S_AC5 - S_sD0);
% 5.2.8 Status 3 of AC_ALB.
P_AC3 = P_H;
Q_shx = M_AC4 * H_AC4 - M_AC5 * H_AC5;
H_AC3 = (Q_shx + M_AC2 * H_AC2) / M_AC3;
y_3str = num2str(y_3);
T_AC3 = CoolProp.PropsSI('T', 'H', H_AC3-5000, 'P', P_AC3, strcat('INCOMP::LiBr[',y_3str,']')); % ????
S_AC3 = CoolProp.PropsSI('S', 'T', T_AC3, 'P', P_AC3, strcat('INCOMP::LiBr[',y_3str,']'));
E_AC3 = (H_AC3 - H_sT0) - T_0 * (S_AC3 - S_sT0);
% 5.2.9 Status 6 of AC_ALB.
y_6 = y_5;
T_AC6 = T_AC5;
H_AC6 = H_AC5;
P_AC6 = P_AC1;
y_6str = num2str(y_6);
S_AC6 = CoolProp.PropsSI('S', 'T', T_AC6, 'P', P_AC6, strcat('INCOMP::LiBr[',y_6str,']'));
M_AC6 = M_AC4;
E_AC5 = (H_AC5 - H_sD0) - T_0 * (S_AC5 - S_sD0);
% 5.2.10 Status 7 of AC_ALB.
M_AC7 = M_AC3 - M_AC4;
M_AC8 = M_AC7;
M_AC9 = M_AC8;
M_AC10 = M_AC9;
P_AC7 = P_H;
% Vapor leaving desorber is assumed in equilibrium with ...
% ... incoming solution stream concentration (state 3).
% This is a standard assumption that represents the best possible case.
T_3wC = CoolProp.PropsSI('T', 'P', P_AC3, 'Q', 0, 'Water') - 273.15;
T_3sat = T_3wC * (a0 + a1 * (y_3*100) + a2 * (y_3*100)^2 + a3 * (y_3*100)^3) + ...
         b0 + b1 * (y_3*100) + b2 * (y_3*100)^2 + b3 * (y_3*100)^3 + 273.15;
T_AC7 = T_3sat;
H_AC7 = CoolProp.PropsSI('H', 'P', P_AC7, 'T', T_AC7, 'Water');
S_AC7 = CoolProp.PropsSI('S', 'P', P_AC7, 'T', T_AC7, 'Water');
E_AC7 = (H_AC7 - H_w0) - T_0 * (S_AC7 - S_w0);
% --------------------------------------------------------------------------------------------------------
%% 5.3 Exergy Analysis of AC_ALB. -------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------
% 5.3.1 Exergy damage in evaporator of AC_ALB.
Q_ACe = M_AC9 * (H_AC10 - H_AC9);
T_AC17 = T_0; P_AC17 = P_0; P_AC18 = P_AC17;
H_AC17 = H_w0; S_AC17 = S_w0; E_AC17 = E_w0;
T_AC18 = T_cw; E_AC18 = E_cw;
DeltaT_ACe = ((T_AC17 - T_AC10) - (T_AC18 - T_AC9)) ...
              / log((T_AC17 - T_AC10) / (T_AC18 - T_AC9));
A_ACe = Q_ACe / DeltaT_ACe / K;
Mcw_ACe = Q_ACe / c_pw / (T_AC17 - T_AC18);
Ed_ACe = (E_AC9 * M_AC9 + E_AC17 * Mcw_ACe) - (E_AC10 * M_AC10 + E_AC18 * Mcw_ACe);
PhiE_ACe = (E_AC17 * Mcw_ACe - E_AC18 * Mcw_ACe) / (E_AC10 * M_AC10 - E_AC9 * M_AC9);
% 5.3.2 Exergy damage in desorber A_ACd. ???????????????????????????????????????????????????
T_AC11 = T_7; P_AC11 = p_7; P_AC12 = P_AC11;          % Temp of High Temp Smoke from MGT.
H_AC11 = H_g7; S_AC11 = H_g7; E_AC11 = E_MGT7;
Q_ACd = M_AC4 * H_AC4 + M_AC7 * H_AC7 - M_AC3 * H_AC3;       % W, HT Rate in desorber
syms T_AC12 A_ACd DeltaT_ACd
eq_ACd(1) = DeltaT_ACd == ((T_AC11 - T_AC4) - (T_AC12 - T_AC3)) ...
                        / log((T_AC11 - T_AC4) / (T_AC12 - T_AC3));
eq_ACd(2) = Q_ACd == DeltaT_ACd * K * A_ACd;
eq_ACd(3) = Q_ACd == Cp_g * Mg * (T_AC11 - T_AC12);
[ST_AC12, SA_ACd, SDeltaT_ACd] = solve(eq_ACd);
    T_AC12 = double(ST_AC12);
     A_ACd = double(SA_ACd);
DeltaT_ACd = double(SDeltaT_ACd);
H_a12 = CoolProp.PropsSI('H', 'T', T_AC12, 'P', P_AC12, 'Air');
S_a12 = CoolProp.PropsSI('S', 'T', T_AC12, 'P', P_AC12, 'Air');
H_H12 = CoolProp.PropsSI('H', 'T', T_AC12, 'P', P_AC12, 'Water');
S_H12 = CoolProp.PropsSI('S', 'T', T_AC12, 'P', P_AC12, 'Water');
H_C12 = CoolProp.PropsSI('H', 'T', T_AC12, 'P', P_AC12, 'CarbonDioxide');
S_C12 = CoolProp.PropsSI('S', 'T', T_AC12, 'P', P_AC12, 'CarbonDioxide');
H_N12 = CoolProp.PropsSI('H', 'T', T_AC12, 'P', P_AC12, 'Nitrogen');
S_N12 = CoolProp.PropsSI('S', 'T', T_AC12, 'P', P_AC12, 'Nitrogen');
H_g12 = Ma_x/Mg * H_a12 + (Ma-Ma_x)/Mg * y_H*18/MgC * H_H12 + ...
       (Ma-Ma_x)/Mg * y_C*44/MgC * H_C12 + (Ma-Ma_x)/Mg * y_N*28/MgC * H_N12;
S_g12 = Ma_x/Mg * S_a12 + (Ma-Ma_x)/Mg * y_H*18/MgC * S_H12 + ...
       (Ma-Ma_x)/Mg * y_C*44/MgC * S_C12 + (Ma-Ma_x)/Mg * y_N*28/MgC * S_N12;
E_AC12 = (H_g12 - H_g0) - T_0 * (S_g12 - S_g0);
Ed_ACd = (E_AC11 * Mg + E_AC3 * M_AC3) - (E_AC7 * M_AC7 + E_AC4 * M_AC4 + E_AC12 * Mg);
PhiE_ACd = (E_AC7 * M_AC7 + E_AC4 * M_AC4 - E_AC3 * M_AC3) / (E_AC11 * Mg - E_AC12 * Mg);
% 5.3.3 Exergy damage in condenser A_ACc.
Q_ACc = M_AC7 * (H_AC7 - H_AC8);
T_AC15 = T_0;
T_AC16 = T_0H;
DeltaT_ACc = ((T_AC7 - T_AC16) - (T_AC8 - T_AC15)) / log((T_AC7 - T_AC16) / (T_AC8 - T_AC15));
A_ACc = Q_ACc / DeltaT_ACc / K;
H_AC15 = CoolProp.PropsSI('H', 'T', T_AC15, 'P', P_0, 'Water');
H_AC16 = CoolProp.PropsSI('H', 'T', T_AC16, 'P', P_0, 'Water');
M_ACc = Q_ACc / (H_AC16 - H_AC15);
Ed_ACc = (E_AC7 * M_AC7 + E_w0 * M_ACc) - (E_AC8 * M_AC8 + E_w0H * M_ACc);
PhiE_ACc = (E_w0H * M_ACc - E_w0 * M_ACc) / (E_AC7 * M_AC7 - E_AC8 * M_AC8);
% 5.3.4 Exergy damage in solution heat exchanger A_ACs.
Ed_ACs = (E_AC2 * M_AC2 + E_AC4 * M_AC4) - (E_AC3 * M_AC3 + E_AC5 * M_AC5);
PhiE_ACs = (E_AC3 * M_AC3 - E_AC2 * M_AC2) / (E_AC4 * M_AC4 - E_AC5 * M_AC5);
% 5.3.5 Exergy damage in absorber of A_ACa.
Q_ACa = M_AC10 * H_AC10 + M_AC6 * H_AC6 - M_AC1 * H_AC1; % W, HT rate of heat exchanger in desorber
T_AC13 = T_0;
T_AC14 = T_0H;
DeltaT_ACa = ((T_AC6 - T_AC14) - (T_AC1 - T_AC13)) / ...
              log((T_AC6 - T_AC14) / (T_AC1 - T_AC13));
A_ACa = Q_ACa / DeltaT_ACa / K;
H_AC13 = CoolProp.PropsSI('H', 'T', T_AC13, 'P', P_0, 'Water');
H_AC14 = CoolProp.PropsSI('H', 'T', T_AC14, 'P', P_0, 'Water');
M_ACa = Q_ACa / (H_AC14 - H_AC13);
Ed_ACa = (E_AC6 * M_AC6 + E_AC10 * M_AC10 + E_w0 * M_ACa) - (E_w0H * M_ACa + E_AC1 * M_AC1);
PhiE_ACa = (E_w0H * M_ACa - E_w0 * M_ACa) / (E_AC6 * M_AC6 + E_AC10 * M_AC10 - E_AC1 * M_AC1);
% 5.3.6 Exergy damage in other components of AC_ALB.
Ed_ACp = E_AC1 * M_AC1 - E_AC2 * M_AC2;
PhiE_ACp = (E_AC2 * M_AC2 - E_AC1 * M_AC1) / We_ACp;
% --------------------------------------------------------------------------------------------------------
%% 5.4 Mathematical Model of R123 Organic Rankine Cycle (ORC_R123) -------------------------------------
%--------------------------------------------------------------------------------------------------------
% 5.4.0 Status 0 of R123.
S_ORC0 = CoolProp.PropsSI('S', 'T', T_0, 'P', P_0, 'R123');
H_ORC0 = CoolProp.PropsSI('H', 'T', T_0, 'P', P_0, 'R123');
E_ORC0 = 0;
% 5.4.1 Status 1, 1-2 Pump
Af = 4.643358E2; Bf =  1.625985E3; Cf = -1.333543E3;
Df = 1.986142E3; Ef = -7.172430E2;
T_rORC1 = T_ORC1 ./ Tc_R123;                                          % Reduced Temerature
D_ORC1 = Af + Bf * (1-T_rORC1).^(1/3) + Cf * (1-T_rORC1).^(2/3) + ...
         Df * (1-T_rORC1) + Ef - (1-T_rORC1).^(4/3);
V_ORC1 = 1 / D_ORC1;
P_ORC1 = CoolProp.PropsSI('P', 'T', T_ORC1, 'Q', 0, 'R123');
S_ORC1 = CoolProp.PropsSI('S', 'T', T_ORC1, 'Q', 0, 'R123');
H_ORC1 = CoolProp.PropsSI('H', 'T', T_ORC1, 'Q', 0, 'R123');
E_ORC1 = (H_ORC1 - H_ORC0) - T_0 * (S_ORC1 - S_ORC0);
% 5.4.2 Status 2, 2-3 Evaporator
H_ORC2 = H_ORC1 + V_ORC1 * (P_ORC2 - P_ORC1);
P_ORC2 = P_ORC1 + DeltaP_ORCp;
T_ORC2 = CoolProp.PropsSI('T', 'H', H_ORC2, 'P', P_ORC2, 'R123');
S_ORC2 = S_ORC1 / EtaS_ORCp;                                             % ???
E_ORC2 = (H_ORC2 - H_ORC0) - T_0 * (S_ORC2 - S_ORC0);
%% 5.4.3 Status 3, 3-4 Gas Turbine
syms V_3sym
P_ORC3 = P_ORC2;
Tr_ORC3 = T_ORC3 ./ Tc_R123;                             % Reduced Temerature
ALPHASqrt3 = 1 + KAPPA_R123 * (1 - sqrt(Tr_ORC3));
ALPHA3 = ALPHASqrt3^2;                                % Temp-dependent para in PR-EOS
aT_ORC3 = aTc_R123 * ALPHA3;                             % Temp-dependent para in PR-EOS
SV_3sym = solve(P_ORC3 == R*T_ORC3 / (V_3sym-b_R123) - aT_ORC3 / (V_3sym*(V_3sym+b_R123) + b_R123*(V_3sym-b_R123)));
V_ORC3 = double(SV_3sym);
V_ORC3 = V_ORC3(imag(V_ORC3)==0);
S_ORC3 = CoolProp.PropsSI('S', 'T', T_ORC3, 'P', P_ORC3, 'R123');
H_ORC3 = CoolProp.PropsSI('H', 'T', T_ORC3, 'P', P_ORC3, 'R123');
E_ORC3 = (H_ORC3 - H_ORC0) - T_0 * (S_ORC3 - S_ORC0);
%% 5.4.4 Status 4, 4-1 Condenser
P_ORC4 = P_ORC1;
T_ORC4 = CoolProp.PropsSI('T', 'P', P_ORC4, 'Q', 1, 'R123');
%{
syms T_4sym
A =  1.656333E3;  B = -2.480583E6; C = 1.792522E1;
D = -8.868380E-2; E =  4.617861E2; F = 1.666667E3;
ST_4sym = solve(log10(P_ORC4/1000) == A + B/T_4sym + C * log10(T_4sym) + D * T_4sym + ...
                                       E * ((F-T_4sym)/T_4sym) * log10(F-T_4sym));
T_ORC4 = double(ST_4sym);
T_ORC4 = real(T_ORC4);
% T_ORC4 = T_ORC4(imag(T_ORC4)==0);
syms V_4sym
Tr_ORC4 = T_ORC4 ./ Tc_R123;                             % Reduced Temerature
ALPHASqrt4 = 1 + KAPPA_R123 * (1 - sqrt(Tr_ORC4));
ALPHA4 = ALPHASqrt4^2;                                % Temp-dependent parameter in PR-EOS
aT_ORC4 = aTc_R123 * ALPHA4;                             % Temp-dependent parameter in PR-EOS, on equation
SV_4sym = solve(P_ORC4 == R*T_ORC4 / (V_4sym-b_R123) - aT_ORC4 / (V_4sym*(V_4sym+b_R123) + b_R123*(V_4sym-b_R123)));
V_ORC4 = double(SV_4sym);
V_ORC4 = V_ORC4(imag(V_ORC4)==0);
%}
S_ORC4 = S_ORC3;
H_ORC4 = CoolProp.PropsSI('H', 'P', P_ORC4, 'Q', 1, 'R123');
E_ORC4 = (H_ORC4 - H_ORC0) - T_0 * (S_ORC4 - S_ORC0);
%% Exergy Analysis of ORC_R123.
% Exergy damage of evaporator in ORC_R123.
T_ORC19 = T_AC12; P_ORC19 = P_AC12; P_ORC20 = P_ORC19;
H_ORC19 = H_AC12; S_ORC19 = S_AC12; E_ORC19 = ex_AC12;
Q_ORC1 = M_ORC * (H_ORC3 - H_ORC2);
syms T_ORC20 A_ORC1 DELTA_T_ORC1
eq_ORC1(1) = DELTA_T_ORC1 == ((T_ORC19 - T_ORC3) - (T_ORC20 - T_ORC2)) / ...
                          log((T_ORC19 - T_ORC3) / (T_ORC20 - T_ORC2));
eq_ORC1(2) = Q_ORC1 == DELTA_T_ORC1 * K * A_ORC1;
eq_ORC1(3) = Q_ORC1 == Mg * Cp_g * (T_ORC19 - T_ORC20);
[ST_ORC20,SA_ORC1,SDELTA_T_ORC1] = solve(eq_ORC1);
     T_ORC20 = double(ST_ORC20);
      A_ORC1 = double(SA_ORC1);
DELTA_T_ORC1 = double(SDELTA_T_ORC1);
% Calculate the H, s, E_ph of gas in status 20 in ORC_R123.
H_a20 = CoolProp.PropsSI('H', 'T', T_ORC20, 'P', P_ORC20, 'Air');
S_a20 = CoolProp.PropsSI('S', 'T', T_ORC20, 'P', P_ORC20, 'Air');
H_H20 = CoolProp.PropsSI('H', 'T', T_ORC20, 'P', P_ORC20, 'Water');
S_H20 = CoolProp.PropsSI('S', 'T', T_ORC20, 'P', P_ORC20, 'Water');
H_C20 = CoolProp.PropsSI('H', 'T', T_ORC20, 'P', P_ORC20, 'CarbonDioxide');
S_C20 = CoolProp.PropsSI('S', 'T', T_ORC20, 'P', P_ORC20, 'CarbonDioxide');
H_N20 = CoolProp.PropsSI('H', 'T', T_ORC20, 'P', P_ORC20, 'Nitrogen');
S_N20 = CoolProp.PropsSI('S', 'T', T_ORC20, 'P', P_ORC20, 'Nitrogen');
H_g20 = Ma_x/Mg * H_a20 + (Ma-Ma_x)/Mg * y_H*18/MgC * H_H20 + ...
       (Ma-Ma_x)/Mg * y_C*44/MgC * H_C20 + (Ma-Ma_x)/Mg * y_N*28/MgC * H_N20;
S_g20 = Ma_x/Mg * S_a20 + (Ma-Ma_x)/Mg * y_H*18/MgC * S_H20 + ...
       (Ma-Ma_x)/Mg * y_C*44/MgC * S_C20 + (Ma-Ma_x)/Mg * y_N*28/MgC * S_N20;
E_ORC20 = (H_g20 - H_g0) - T_0 * (S_g20 - S_g0);
PhiEx_ORC1 = (E_ORC4 * M_ORC + E_ORC20 * Mg) / (E_ORC19 * Mg + E_ORC3 * M_ORC);
Ed_ORCe = (E_ORC19 * Mg + E_ORC3 * M_ORC) - (E_ORC4 * M_ORC + E_ORC20 * Mg);
% Exergy damage of condenser in ORC_R123.
Q_ORC2 = M_ORC * (H_ORC4 - H_ORC1);
DELTA_T_ORC2 = ((T_ORC4 - T_0H) - (T_ORC1 - T_0)) ...
                / log((T_ORC4 - T_0H) / (T_ORC1 - T_0));
A_ORC2 = Q_ORC2 / DELTA_T_ORC2 / K;
H_ORC21 = CoolProp.PropsSI('H', 'T', T_0, 'P', P_0, 'Water');
H_ORC22 = CoolProp.PropsSI('H', 'T', T_0H, 'P', P_0, 'Water');
M_ORCc = Q_ORC2 / (H_ORC22 - H_ORC21);                   % kg/s, Amount of supplying cooling water.
Ed_ORCc = (E_ORC4 * M_ORC + E_w0 * M_ORCc) - (E_ORC1 * M_ORC + E_w0H * M_ORCc);
% Exergy damage of gae turbine in ORC_R123.
WP_ORCt = M_ORC * (H_ORC3 - H_ORC4);                       % W,    Output work of ORC_R123.
We_ORCt = WP_ORCt * EtaG_ORCt;
Ed_ORCt = E_ORC3 * M_ORC - (E_ORC4 * M_ORC + We_ORCt);
% Exergy damage of pump in ORC_R123.
We_ORCp = M_ORC * (H_ORC2 - H_ORC1);                      % W,    Required Work of Pump in ORC_R123.
Ed_ORCp = (E_ORC1 * M_ORC + We_ORCp) - E_ORC2 * M_ORC;
% --------------------------------------------------------------------------------------------------------
%% 5.5 Mathematical Model of HP_R410a. -----------------------------------------------------------
% --------------------------------------------------------------------------------------------------------
MS_HPc = Ms - Ms_MGT;
Q_HPc = MS_HPc * (H_s - H_w0);
% 5.5.1 Status 0 of HP_R410a.
S_HP0 = CoolProp.PropsSI('S', 'T', T_0, 'P', P_0, 'R410a');
H_HP0 = CoolProp.PropsSI('H', 'T', T_0, 'P', P_0, 'R410a');
%% 5.5.2 Status 1 of HP_R410a.
P_HP1 = CoolProp.PropsSI('P', 'T', T_HP1, 'Q', 1, 'R410a');
% Solve for V_HP1
syms V_HP1sym
  Tr_HP1 = T_HP1 ./ Tc_R410a;                             % Reduced Temerature
  ALPHASqrt1 = 1 + Kappa_R410a * (1 - sqrt(Tr_HP1));
  ALPHA4 = ALPHASqrt1^2;                                  % Temp-dependent para in PR-EOS
  a_T1 = aTc_R410a * ALPHA4;                              % Temp-dependent para in PR-EOS
  SV_HP1sym = solve(P_HP1 == R*T_HP1 / (V_HP1sym-b_R410a) - ...
                         a_T1 / (V_HP1sym*(V_HP1sym+b_R410a) + b_R410a*(V_HP1sym-b_R410a)));
V_HP1 = double(SV_HP1sym); V_HP1 = V_HP1(imag(V_HP1)==0);
S_HP1 = CoolProp.PropsSI('S', 'T', T_HP1, 'Q', 1, 'R410a');
H_HP1 = CoolProp.PropsSI('H', 'T', T_HP1, 'Q', 1, 'R410a');
E_HP1 = (H_HP1 - H_HP0) - T_0 * (S_HP1 - S_HP0);
%% 5.5.3 Status 4 of HP_R410a.
P_HP4 = P_HP1; T_HP4 = T_HP1;
% Solve for V_HP4.
syms V_HP4sym
  Tr_HP4 = T_HP4 ./ Tc_R410a;                             % Reduced Temerature
  ALPHASqrt4 = 1 + Kappa_R410a * (1 - sqrt(Tr_HP4));
  ALPHA4 = ALPHASqrt4^2;                                 % Temp-dependent para in PR-EOS
  aT_ORC4 = aTc_R410a * ALPHA4;                          % Temp-dependent para in PR-EOS
  SV_HP4sym = solve(P_HP4 == R*T_HP4 / (V_HP4sym-b_R410a) - ...
                  aT_ORC4 / (V_HP4sym*(V_HP4sym+b_R410a) + b_R410a*(V_HP4sym-b_R410a)));
V_HP4 = double(SV_HP4sym);
V_HP4 = V_HP4(imag(V_HP4)==0);
%% 5.5.4 Status 3 of HP_R410a.
%{
% <8> solve V_HP3 through Equation for Density of the Saturated Liquid.
Af = 1.000000; Bf =  1.984734;    Cf = -1.767593E-01;
Df = 1.819972; Ef = -7.171684E-1;
Tr_HP3 = T_HP3 ./ Tc_R410a;                             % Reduced Temerature
D_HP3 = Dc_R410a * (Af + Bf * (1-Tr_HP3).^(1/3) + Cf * (1-Tr_HP3).^(2/3) + ...
                  Df * (1-Tr_HP3) + Ef - (1-Tr_HP3).^(4/3));
V_HP3 = D_HP3 / 1;
% <9> Solve for P_HP3.
syms P_HP3sym
  Tr_HP3 = T_HP3 ./ Tc_R410a;                             % Reduced Temerature
  ALPHASqrt3 = 1 + Kappa_R410a * (1 - sqrt(Tr_HP3));
  ALPHA3 = ALPHASqrt3^2;                                  % Temp-dependent para in PR-EOS
  aT_ORC3 = aTc_R410a * ALPHA3;                           % Temp-dependent para in PR-EOS
  SP_HP3sym = solve(P_HP3sym == R*T_HP3 / (V_HP3-b_R410a) - ...
                  aT_ORC3 / (V_HP3*(V_HP3+b_R410a) + b_R410a*(V_HP3-b_R410a)));
  P_HP3 = double(SP_HP3sym);
  P_HP3 = P_HP3(imag(P_HP3)==0);
%}
P_HP3 = CoolProp.PropsSI('P', 'T', T_HP3, 'Q', 0, 'R410a');
S_HP3 = CoolProp.PropsSI('S', 'T', T_HP3, 'Q', 0, 'R410a');
H_HP3 = CoolProp.PropsSI('H', 'T', T_HP3, 'Q', 0, 'R410a');
E_HP3 = (H_HP3 - H_HP0) - T_0 * (S_HP3 - S_HP0);
%
S_HP4 = S_HP3 / EtaS_HPv;
H_HP4 = H_HP3;
E_HP4 = (H_HP4 - H_HP0) - T_0 * (S_HP4 - S_HP0);
%% 5.5.5 Status 2 of HP_R410a.
P_HP2 = P_HP3;
S_HP2 = S_HP1 / EtaS_HPp;
T_HP2 = CoolProp.PropsSI('T', 'P', P_HP2, 'S', S_HP2, 'R410a');
%
syms V_HP2sym
  Tr_HP2 = T_HP2 ./ Tc_R410a;                         % Reduced Temerature
  ALPHASqrt2 = 1 + Kappa_R410a * (1 - sqrt(Tr_HP2));
  ALPHA2 = ALPHASqrt2^2;                              % Temp-dependent para in PR-EOS
  a_T2 = aTc_R410a * ALPHA2;                          % Temp-dependent para in PR-EOS
  SV_HP2sym = solve(P_HP2 == R*T_HP2 / (V_HP2sym-b_R410a) - ...
                  a_T2 / (V_HP2sym*(V_HP2sym+b_R410a) + b_R410a*(V_HP2sym-b_R410a)));
  V_HP2 = double(SV_HP2sym);
  V_HP2 = V_HP2(imag(V_HP2)==0);
H_HP2 = CoolProp.PropsSI('H', 'T', T_HP2, 'P', P_HP2, 'R410a');
E_HP2 = (H_HP2 - H_HP0) - T_0 * (S_HP2 - S_HP0);
%% Exergy Analysis of HP_R410a. ------------------------------------------------------------------------
% 5.5.6 Exergy damage in condenser of HP_R410a.
M_HP = Q_HPc / (H_HP2 - H_HP3);
DeltaT_HPc = ((T_HP2 - T_s) - (T_HP3 - T_8)) ...
              / log((T_HP2 - T_s) / (T_HP3 - T_8));
A_HPc = Q_HPc / DeltaT_HPc / K;
Ed_HPc = (E_HP2 * M_HP + E_MGT8 * MS_HPc) - (E_HP3 * M_HP + E_s * MS_HPc);
PhiE_HPc = (E_s * Ms_HPc - E_MGT8 * Ms_HPc) / (E_HP2 * M_HP - E_HP3 * M_HP);
% 5.5.7 Exergy damage in pump of HP_R410a.
We_HPp = M_HP * (H_HP2 - H_HP1);
Ed_HPp = (E_HP1 * M_HP + We_HPp) - E_HP2 * M_HP;
PhiE_HPp = (E_HP2 * M_HP - E_HP1 * M_HP) / We_HPp;
% 5.5.8 Area of heat exchange in evaporator A_HPe.
Q_HPe = M_HP * (H_HP1 - H_HP4);
M_HPe = Q_HPe / (H_w0 - H_w0L);
DeltaT_HPe = ((T_0 - T_HP1) - (T_0L - T_HP4)) ...
              / log((T_0 - T_HP1) / (T_0L - T_HP4));
A_HPe = Q_HPe / DeltaT_HPe / K;
Ed_HPe = (E_HP4 * M_HP + E_w0 * M_HPe) - (E_HP1 * M_HP + E_w0L * M_HPe);
PhiE_HPe = (E_HP1 * M_HP - E_HP4 * M_HP) / (E_w0 * M_HPe - E_w0L * M_HPe);
% ----------------------------------------------------------------------------------------------------
%% 5.6 Mathematical Model of VCC_R134a. -----------------------------------------------------------
% ----------------------------------------------------------------------------------------------------
Mcw_VCCe = Mcw - Mcw_ACe;
Q_VCCe = Mcw_VCCe * (H_w0 - H_cw);
% 5.6.1 Status 0 of HP_R410a.
S_VCC0 = CoolProp.PropsSI('S', 'T', T_0, 'P', P_0, 'R134a');
H_VCC0 = CoolProp.PropsSI('H', 'T', T_0, 'P', P_0, 'R134a');
E_VCC0 = 0;
% 5.6.2 Status 1 of VCC_R134a.
% <1> Solve T_VCC4 through Equation for Saturated Vapor Pressure.
A =  4.069889E1;  B = -2.362540E3;  C = -1.306883E1;
D =  7.616005E-3; E =  2.342564E-1; F =  3.761111E2;
P_VCC1 = 10^(A + B./T_VCC1 + C * log10(T_VCC1) + D * T_VCC1 + ...
         E * ((F-T_VCC1)./T_VCC1) * log10(F-T_VCC1)) * 1000;
% <2> Solve for V_VCC1
syms V_VCC1sym
  Tr_VCC1 = T_VCC1 ./ Tc_R134a;                             % Reduced Temerature
  ALPHASqrt1 = 1 + Kappa_R134a * (1 - sqrt(Tr_VCC1));
  ALPHA4 = ALPHASqrt1^2;                                   % Temp-dependent para in PR-EOS
  aT_VCC1 = aTc_R134a * ALPHA4;                             % Temp-dependent para in PR-EOS
  SV_VCC1sym = solve(P_VCC1 == R*T_VCC1 / (V_VCC1sym-b_R134a) - ...
                    aT_VCC1 / (V_VCC1sym*(V_VCC1sym+b_R134a) + b_R134a*(V_VCC1sym-b_R134a)));
  V_VCC1 = double(SV_VCC1sym);
  V_VCC1 = V_VCC1(imag(V_VCC1)==0);
S_VCC1 = CoolProp.PropsSI('S', 'T', T_VCC1, 'Q', 1, 'R134a');
H_VCC1 = CoolProp.PropsSI('H', 'T', T_VCC1, 'Q', 1, 'R134a');
E_VCC1 = (H_VCC1 - H_VCC0) - T_0 * (S_VCC1 - S_VCC0);
% 5.6.3(1) Status 4 of VCC_R134a.
P_VCC4 = P_VCC1; T_VCC4 = T_VCC1;
syms V_VCC4sym
  Tr_VCC4 = T_VCC4 ./ Tc_R134a;                             % Reduced Temerature
  ALPHASqrt4 = 1 + Kappa_R134a * (1 - sqrt(Tr_VCC4));
  ALPHA4 = ALPHASqrt4^2;                                   % Temp-dependent para in PR-EOS
  aT_VCC4 = aTc_R134a * ALPHA4;                             % Temp-dependent para in PR-EOS
  SV_VCC4sym = solve(P_VCC4 == R*T_VCC4 / (V_VCC4sym-b_R134a) - ...
                     aT_VCC4 / (V_VCC4sym*(V_VCC4sym+b_R134a) + b_R134a*(V_VCC4sym-b_R134a)));
  V_VCC4 = double(SV_VCC4sym);
  V_VCC4 = V_VCC4(imag(V_VCC4)==0);
% 5.6.4 Status 3 of VCC_R134a.
% <8> solve V_VCC3 through Equation for Density of the Saturated Liquid.
Af =  5.281464E2; Bf =  7.551834E2; Cf = 1.028676E3;
Df = -9.491172E2; Ef = 5.935660E2;
Tr_VCC3 = T_VCC3 ./ Tc_R134a;                             % Reduced Temerature
RHO_3 = Af + Bf * (1-Tr_VCC3).^(1/3) + Cf * (1-Tr_VCC3).^(2/3) + ...
        Df * (1-Tr_VCC3) + Ef - (1-Tr_VCC3).^(4/3);
V_VCC3 = RHO_3 / 1;
%{
% <9> Solve for P_VCC3.
syms P_VCC3sym
  Tr_VCC3 = T_VCC3 ./ Tc_R134a;                             % Reduced Temerature
  ALPHASqrt3 = 1 + Kappa_R134a * (1 - sqrt(Tr_VCC3));
  ALPHA3 = ALPHASqrt3^2;                                 % Temp-dependent para in PR-EOS
  aT_VCC3 = aTc_R134a * ALPHA3;                          % Temp-dependent para in PR-EOS
  SP_VCC3sym = solve(P_VCC3sym == R*T_VCC3 / (V_VCC3-b_R134a) - ...
                     aT_VCC3 / (V_VCC3*(V_VCC3+b_R134a) + b_R134a*(V_VCC3-b_R134a)));
  P_VCC3 = double(SP_VCC3sym);
  P_VCC3 = P_VCC3(imag(P_VCC3)==0);
%}
P_VCC3 = CoolProp.PropsSI('P', 'T', T_VCC3, 'Q', 0, 'R134a');
S_VCC3 = CoolProp.PropsSI('S', 'T', T_VCC3, 'Q', 0, 'R134a');
H_VCC3 = CoolProp.PropsSI('H', 'T', T_VCC3, 'Q', 0, 'R134a');
E_VCC3 = (H_VCC3 - H_VCC0) - T_0 * (S_VCC3 - S_VCC0);
% 5.6.3(2) Status 4 of VCC_R134a.
S_VCC4 = S_VCC3 / EtaS_VCCv;
H_VCC4 = H_VCC3;
E_VCC4 = (H_VCC4 - H_VCC0) - T_0 * (S_VCC4 - S_VCC0);
% 5.6.5 Status 2 of VCC_R134a.
P_VCC2 = P_VCC3;
S_VCC2 = S_VCC1 / EtaS_VCCp;
T_VCC2 = CoolProp.PropsSI('T', 'P', P_VCC2, 'S', S_VCC2, 'R134a');
% <11> Solve for V_VCC2.
syms V_VCC2sym
  Tr_VCC2 = T_VCC2 ./ Tc_R134a;                             % Reduced Temerature
  ALPHASqrt2 = 1 + Kappa_R134a * (1 - sqrt(Tr_VCC2));
  ALPHA2 = ALPHASqrt2^2;                                    % Temp-dependent para in PR-EOS
  aT_VCC2 = aTc_R134a * ALPHA2;                             % Temp-dependent para in PR-EOS
  SV_VCC2sym = solve(P_VCC2 == R*T_VCC2 / (V_VCC2sym-b_R134a) - ...
                    aT_VCC2 / (V_VCC2sym*(V_VCC2sym+b_R134a) + b_R134a*(V_VCC2sym-b_R134a)));
  V_VCC2 = double(SV_VCC2sym);
  V_VCC2 = V_VCC2(imag(V_VCC2)==0);
H_VCC2 = CoolProp.PropsSI('H', 'P', P_VCC2, 'S', S_VCC2, 'R134a');
E_VCC2 = (H_VCC2 - H_VCC0) - T_0 * (S_VCC2 - S_VCC0);
%% Exergy Analysis of VCC_R410a.
% 5.6.6 Exergy damage in evaporator of VCC_R134a. Absorb heat through evaporator.
M_VCC = Q_VCCe / (H_VCC1 - H_VCC4);
DeltaT_VCCe = ((T_0 - T_VCC1) - (T_cw - T_VCC4)) ...
               / log((T_0 - T_VCC1) / (T_cw - T_VCC4));
A_VCCe = Q_VCCe / DeltaT_VCCe / K;
Ed_VCCe = (E_VCC4 * M_VCC + E_w0 * Mcw_VCCe) - (E_VCC1 * M_VCC + E_cw * Mcw_VCCe);
PhiE_VCCe = (E_w0 * Mcw_VCCe - E_cw * Mcw_VCCe) / (E_VCC1 * M_VCC - E_VCC4 * M_VCC);
% 5.6.7 Exergy damage in pump of VCC_R134a.
We_VCCp = M_VCC * (H_VCC2 - H_VCC1);
Ed_VCCp = (E_VCC1 * M_VCC + We_VCCp) - E_VCC2 * M_VCC;
PhiE_VCCp = (E_VCC2 - E_VCC1) * M_VCC / We_VCCp;
% 5.6.8 Exergy damage in condenser of VCC_R134a.
Q_VCCc = M_VCC * (H_VCC2 - H_VCC3);                        % W, HT Rate in desorber
M_VCCc = Q_VCCc / (h_w0H - h_w0);
DeltaT_VCCc = ((T_VCC2 - T_0H) - (T_VCC3 - T_0)) ...
               / log((T_VCC2 - T_0H) / (T_VCC3 - T_0));
A_VCCc = Q_VCCc / DeltaT_VCCc / K;
Ed_VCCc = (E_VCC2 * M_VCC + E_w0 * M_VCCc) - (E_VCC3 * M_VCC + E_w0H * M_VCCc);
PhiE_VCCc = (E_VCC2 * M_VCC - E_VCC3 * M_VCC) / (E_w0H * M_VCCc - E_w0 * M_VCCc);
