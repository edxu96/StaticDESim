% Title: Thermo-economic Optimization of Distributed Energy System in Green Energy Island.
% Based on the theory of nolinear equality and inequality constraints.
% Method: Genetic Algorithm within MATLAB Global Optimization Toolbox.
% Version: 1.1, 2018.6.2, Jie Xu.
% SubTitle: Calculate and the result of GEI optim, then exergy analysis.
% 1. Micro Gas Turbine (MGT)
% 2. Aqueous Lithium-Bromine Single-Effect Absorption Chiller (AC_ALB)
% 3. R123 Organic Recycle Cycle (ORC_R123)
% 4. R410a Heat Pump (HP_R410a)
% 5. R134a Vapor Compression Chiller (VCC_R134a)
function f = simple_fitness(x)
%% Input decision variable.
% Decision variable of MGT.
x(1) = 633186.25;
x(2) = 0.837;
x(3) = 0.917;
x(4) = 769.395751953125;
x(5) = 1300;
x(6) = 0.03;
x(7) = 60 * 1000;
% Decision variable of AC_ALB.
x(8) = 4 + 273.15;
x(9) = 30 + 273.15;
x(10) = 0.5322;
x(11) = 0.6711;
x(12) = 0.05;
x(13) = 0.8;
% Decision variable of ORC_R123.
x(14) = 0.05;
x(15) = 200 + 273.15;
x(16) = 6 * 101.325 * 1000;
x(17) = 80 + 273.15;
% Decision variable of HP_R410a.
x(18) = 30 + 273.15;
x(19) = 70 + 273.15;
x(20) = 10 + 273.15;
x(21) = 0.005;
% Decision variable of VCC_R134a.
x(22) = 30 + 273.15;
x(23) = 70 + 273.15;
x(24) = 10 + 273.15;
x(25) = 0.005;
%% 1. General Constants ---------------------------------------------------------------------------------
R = 8.314472;                                     % J/(mol*K), Universial Gas Constant
p_0 = 101.325 * 1000;                             % Pa, Pressure of atmosphere.
T_0 = 25 + 273.15;                                % K, Temperature of atmosphere.
T_0H = 35 + 273.15;                               % K, Acceptable Highest Temp of atmosphere.
T_0L = 15 + 273.15;                               % K, Acceptable Lowest Temp of atmosphere.
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
h_8p = 640039.2;                                  % 回水的焓值 (503.97 + (589.30-503.97)/20*16.85) * 1000
h_8 = 293316;                                     % 余热回收后的焓值 (251.56 + (335.29-251.56)/20*10) * 1000
c_pa = 1004;                                      % J/kg/K, air specific heat capacity at constant pressure
c_pg = 1170;                                      % J/kg/K, 实际燃气 specific heat capacity at constant pressure
c_pw = 4200;                                      % J/kg/K, water specific heat capacity at constant pressure
ETA_CC = 0.98877;
% 1.2 Constants for AC_ALB
K = 1000;                                         % W / m^2 / K, Thermal conductivity.
Z_A = 100;                                        % RMB / m^2, Cost rate of area of heat transfer.
Z_W = 3.25E-3;                                    % RMB / kg, Cost rate of cooling water.
Z_cw = 0.5;                                       % RMB / kW*h, Profit rate of supplying cooling load.
Z_E = 0.6;                                        % RMB / kW*h, Cost of supplying electricity for pump.
Q_wAC = 44.1 * 1000;                              % W,   ????
Z_wAC = 235550;                                   % RMB, ????
ALPHA_wAC = 0.73;                                 %      ????
% 1.3 Constants for ORC_R123.
Rg_R123 = 0.0544 * 1000;                          % J/(mol*K), Gas Constant.
M_R123 = 152.93 / 1000;                           % kg/mol, Molecular Weight of R123.
Tc_R123 = 456.83;                                 % K, Critical Temp of R123.
Pc_R123 = 3668.0 * 1000;                          % Pa, Critical Pressure of R123.
ETA_ps = 0.9;                                     % Isentropic efficiency of pump.
ETA_ts = 0.9;                                     % Isentropic efficiency of turbine.
ETA_p = 0.8;                                      % Efficiency of pump.
ETA_t = 0.8;                                      % Efficiency of turbine.
ETA_e = 0.9;                                      % Efficiency of evaporator.
ETA_c = 0.9;                                      % Efficiency of condenser.
DELTA_p_C = 100 * 1000;                           % Pa, Pressure drop in condenser.
DELTA_p_E = 100 * 1000;                           % Pa, Pressure drop in evaporator.
OMEGA_R123 = 0.281922497036;                      % Acentric Factor, R123.
KAPPA_R123 = 0.37464 + (1.54226 - ...             % Dependent on OMEGA(working substance),
        0.26992 * OMEGA_R123) * OMEGA_R123;       % Temp-independent parameter in PR-EOS
aTc_R123 = 0.457235529 * (Rg_R123 * ...
           Tc_R123)^2 ./ Pc_R123;                 % Critical Point Restriction "a(T_c)"
b_R123 = 0.077796074 * Rg_R123 * ...              % m^3/mol, Critical Point Restriction "b",
         Tc_R123 ./ Pc_R123;                      % Temp-independent parameter in PR-EOS
Q_wORC = 23.6 * 1000;                             % W,   ????
Z_wORC = 385600;                                  % RMB, ????
ALPHA_wORC = 0.73;                                %      ????
Eta_G = 0.9;                                      % Efficiency of generator.
%% 1.4 Constants for HP_R410a.
R = 8.314472;                                 % J/(mol*K), Universial Gas Constant
m_R410a = 72.58 / 1000;                       % kg / mol , Molar Mass
Rg_R410a = 0.11455 * 1000;                    % J/(K*kg) , Gas Constant - R134a
Tc_R410a = 345.28;                            % K , temperature in Critical Point.
Pc_R410a = 4926.1 * 1000;                     % Pa, pressure in Critical Point.
Dc_R410a = 488.90;                            % kg/m^3, critical density.
Tb_R410a = -26.06 + 273.15;                   % K , Boiling point at one atmosphere.
EtaS_HPp = 0.7;                               % Isentropic efficiency of pump.
EtaS_HPt = 0.7;                               % Isentropic efficiency of turbine.
Eta_HPv = 0.8;                                % Efficiency of throttle valve.
Eta_HPp = 0.8;                                % Efficiency of compressor(p).
Eta_HPe = 0.9;                                % Efficiency of evaporator.
Eta_HPc = 0.9;                                % Efficiency of condenser.
DeltaP_HPc = 100 * 1000;                      % Pa, Pressure drop in condenser.
DeltaP_HPe = 100 * 1000;                      % Pa, Pressure drop in evaporator.
Omega_R410a = 0.296;                          % Acentric Factor.
Kappa_R410a = 0.37464 + ...                   % Dependent on Omega_R410a(working substance),
        (1.54226 - 0.26992 * Omega_R410a) * Omega_R410a;  % Temperature-independent parameter in PR-EOS
aTc_R410a = 0.457235529 * (R_gR410 * Tc_R140a)^2 ./ p_c;  % Critical Point Restriction "a(Tc_R140a)"
b_R410a = 0.077796074 * R_gR410 * Tc_R140a ./ p_c;        % m^3/mol, Critical Point Restriction "b_R410a",
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
DeltaP_C = 100 * 1000;                       % Pa, Pressure drop in condenser.
DeltaP_E = 100 * 1000;                       % Pa, Pressure drop in evaporator.
Omega_R134a = 0.332;                                % Acentric Factor.
Kappa_R134a = 0.37464 + ...                         % Dependent on Omega_R134a(working substance),
        (1.54226 - 0.26992 * Omega_R134a) * Omega_R134a;  % Temperature-independent parameter in PR-EOS
aTc_R134a = 0.457235529 * (Rg_R134a * Tc_R134a)^2 ./ Pc_R134a;    % Critical Point Restriction "a(Tc_R134a)"
b_R134a = 0.077796074 * Rg_R134a * Tc_R134a ./ Pc_R134a;           % m^3/mol, Critical Point Restriction "b_R134a",
                                              % Temperature-independent parameter in PR-EOS
%% 3. Pre-defined Condition. ----------------------------------------------------------------------------
N = 8000;                                         % Operating Hours in Unit Years
N_y = 10;                                         % Unit Years
p_s = 500 * 1000;                                 % Pa, Pressure of supplying steam.
T_cw = 4 + 273.15; p_cw = p_0;                    % Supplying cooling water.
  Ms = 10000;
 Mcw = 10000;
   E = 10000;
%% 4. Decision Variables. -------------------------------------------------------------------------------
% 4.1 Decision Variables in MGT.
      p_2 = x(1);               % Pa,   Outlet Pressure of Air Compressor
Eta_MGTac = x(2);               %       Isentropic Efficiency of Air Compressor
 Eta_MGTt = x(3);               %       Isentropic Efficiency of Gas Turbine
      T_3 = x(4);               % K,    Outlet Temp of Air from Regenerator
      T_4 = x(5);               % K,    Inlet Temp of Gas Turbine
   Ms_MGT = x(6);               % kg/s, Fluid Rate of Saturated Steam from HRSG
   Wp_MGT = x(7);               % W,    Net Power from Micro Gas Turbine
% 4.2 Decision Variables in AC_ALB.
   T_AC10 = x(8);               % K,    Outlet temperature of Evaporator.
    T_AC8 = x(9);               % K,    Outlet temperature of Condenser.
      y_1 = x(10);              %       Mass Fraction of Outlet Solution from Absorber.
      y_4 = x(11);              %       Mass Fraction of Outlet Solution from Desorber.
    M_AC1 = x(12);              % kg/s, Fluid Rate Outlet Solution from Absorber.
  Eta_shx = x(13);              %       Solution Heat Exchanger Ratio.
% 4.3 Decision Variables in ORC_R123.
    M_ORC = x(14);              % kg/s, Fluid Rate.
   T_ORC3 = x(15);              % K,    Inlet temperature of turbine
   p_ORC2 = x(16);              % K,    Outlet pressure of pump / inlet pressure of turbine
   T_ORC1 = x(17);              % K,    Outlet temperature of condenser.
% 4.4 Decision Variables in HP_R410a.
   T_HP1 = x(18);               % K,    Outlet temperature of evaporator.
   T_HP3 = x(19);               % K,    Outlet temperature of condenser.
   T_HP2 = x(20);               % K,    Outlet temperature of compressor.
    M_HP = x(21);               % kg/s, Fluid rate of HP_R410a.
% 4.5 Decision Variables in VCC_R134a.
  T_VCC1 = x(22);               % K,    Outlet temperature of evaporator.
  T_VCC3 = x(23);               % K,    Outlet temperature of condenser.
  T_VCC2 = x(24);               % K,    Outlet temperature of compressor.
   M_VCC = x(25);               % kg/s, Fluid Rate.
%% 5. Mathematical Model of Green Energy Island ----------------------------------------------------------
Ex_s  = Ms; %
Ex_cw = 1000;
Ex_e  = 1000;
h_w0 = CoolProp.PropsSI('H', 'T', T_0, 'P', p_0, 'Water');
s_w0 = CoolProp.PropsSI('S', 'T', T_0, 'P', p_0, 'Water');
ex_ph0 = 0;
h_s = CoolProp.PropsSI('H', 'P', p_s, 'Q', 1, 'Water');
s_s = CoolProp.PropsSI('S', 'P', p_s, 'Q', 1, 'Water');
%% 5.1 Mathematical Model of Micro Gas Turbine (MGT) ----------------------------------------------------
% 5.1.1 Mathematical Model of Brayton Cycle.
T_2 = T_1 .* (1 + 1./Eta_MGTac * ...
      ((p_2./p_1)^((Ka_MGT-1)./Ka_MGT) - 1));                         % (1)  AC
p_3 = p_2 * (1 - DELTA_p_aRE);                                  % (7)  RE
p_4 = p_3 * (1 - DELTA_p_cc);                                   % (5)  CC
p_6 = p_7 ./ (1 - DELTA_p_HRSG);                                % (14) HRSG
p_5 = p_6 ./ (1 - DELTA_p_gRE);                                 % (8)  RE
T_5 = T_4 * (1 - Eta_MGTt * (1 - (p_4./p_5)^((1-Kg_MGT)./Kg_MGT)));   % (9)  GT
H = (c_pg * T_4 - ETA_CC * LHV) ./ (c_pa * T_3 - ETA_CC * LHV); % (4) CC, H = m_a / m_g;
m_g = Wp_MGT ./ (c_pg*(T_4-T_5) - c_pa*(T_2-T_1)*H);            % (2) (10)
m_a = H * m_g;                                                  % H = m_a / m_g;
m_f = m_g - m_a;                                                % (3)  CC
T_6 = T_5 - m_a * c_pa * (T_3 - T_2) ./ (m_g * c_pg);           % (6)  RE
h_9 = h_s;
T_7p = T_6 - Ms_MGT * (h_9 - h_8p) ./ (m_g * c_pg);             % (12) HRSG
T_7 = T_6 - Ms_MGT * (h_9 - h_8) ./ (m_g * c_pg);               % (13) HRSG
W_AC = m_a * c_pa * (T_2 - T_1);                                % (2)  AC
W_GT = Wp_MGT + W_AC;                                           % (11) GT
%% 5.1.2 Calculate the h and s of air and gas before combustion.
% Calculate the h, s, E_ph of air, fuel and water in status 0.
h_a0 = CoolProp.PropsSI('H', 'T', T_0, 'P', p_0, 'Air');
s_a0 = CoolProp.PropsSI('S', 'T', T_0, 'P', p_0, 'Air');
h_f0 = CoolProp.PropsSI('H', 'T', T_0, 'P', p_0, 'Methane');
s_f0 = CoolProp.PropsSI('S', 'T', T_0, 'P', p_0, 'Methane');
% Calculate the h, s, E_ph of air in status 1.
h_a1 = CoolProp.PropsSI('H', 'T', T_1, 'P', p_1, 'Air');
s_a1 = CoolProp.PropsSI('S', 'T', T_1, 'P', p_1, 'Air');
ex_ph1 = (h_a1 - h_a1) + T_1 * (s_a1 - s_a1);
% Calculate the h, s, E_ph of air in status 2.
h_a2 = CoolProp.PropsSI('H', 'T', T_2, 'P', p_2, 'Air');
s_a2 = CoolProp.PropsSI('S', 'T', T_2, 'P', p_2, 'Air');
ex_ph2 = (h_a2 - h_a1) - T_1 * (s_a2 - s_a1);
% Calculate the h, s, E_ph of gas in status 3.
h_a3 = CoolProp.PropsSI('H', 'T', T_3, 'P', p_3, 'Air');
s_a3 = CoolProp.PropsSI('S', 'T', T_3, 'P', p_3, 'Air');
h_g3 = m_a/m_g * h_a3 + m_f/m_g * h_f0;
s_g3 = m_a/m_g * s_a3 + m_f/m_g * s_f0;
ex_ph3 = ((h_a3 - h_a1) - T_1 * (s_a3 - s_a1)) * m_a/m_g;
%% 5.1.3 Calculate the components of waste gas.
y_O = 0.21;                                    % Volume fraction of Oxygen in Air.
n_aC = 2 / y_O;                                % Needed air for 1 mol fuel combustion
n_gC = 1 + 2 + n_aC - 2;                       % Waste gas from combustion with just needed air
y_C = (1 + n_aC * 0.0004) / n_gC;
y_H = 2 / n_gC;
y_N = 1 - y_C - y_H;
M_gC = y_C * 44 + y_H * 18 + y_N * 28;
m_aC = n_aC * 29 / 16;
m_gC = n_gC * M_gC / 16;
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
%% 5.1.4 Calculate the s, E_ph of supplying steam.
% Calculate the h, s, E_ph of supplying steam in status 8.
s_8 = CoolProp.PropsSI('S', 'T', T_8, 'P', p_8, 'Water');
ex_ph8 = (h_8 - h_w0) - T_0 * (s_8 - s_w0);
% Calculate the h, s, E_ph of supplying steam in status 9.
s_9 = s_s;
ex_ph9 = (h_9 - h_w0) - T_0 * (s_9 - s_w0);
%% 5.1.5 Exergy Analysis of MGT.
r_water = (2257.2 * 1000);                            % J/kg, latent heat of vaporization
HHV = LHV + 1000/16 * n_gC * y_H * 18/1000 * r_water; % J/kg
ex_f = HHV;                                            % J/kg
PhiEx_CC = ex_ph4 / (ex_ph3 + ex_f);
PhiEx_GT = W_GT / ((ex_ph4 - ex_ph5) * m_a);
PhiEx_AC = (ex_ph2 - ex_ph1) * m_a / W_AC;
PhiEx_RE = (ex_ph3 * m_a + ex_ph6 * m_g) / (ex_ph2 * m_a + ex_ph5 * m_g);
PhiEx_HRSG = (ex_ph9 * Ms_MGT + ex_ph7 * m_g) / (ex_ph8 * Ms_MGT + ex_ph6 * m_g);
%% 5.2 Mathematical Model of Aqueous Lithium-Bromide Absorption Chiller (AC_ALB)-------------------------
% Reference State of Lithium Bromide Solution and Water in AC_ALB.
s_w0 = CoolProp.PropsSI('S', 'T', T_0, 'P', p_0, 'Water');
h_w0 = CoolProp.PropsSI('H', 'T', T_0, 'P', p_0, 'Water');
E_w0 = 0;
% 5.2.1 Status 10 of AC_ALB.
p_AC10 = CoolProp.PropsSI('P','T', T_AC10, 'Q', 1, 'Water');
p_L = p_AC10;
s_AC10 = CoolProp.PropsSI('S', 'T', T_AC10, 'P', p_AC10, 'Water');
h_AC10 = CoolProp.PropsSI('H', 'T', T_AC10, 'P', p_AC10, 'Water');
ex_AC10 = (h_AC10 - h_w0) - T_0 * (s_AC10 - s_w0);
% 5.2.2 Status 8 of AC_ALB.
p_AC8 = CoolProp.PropsSI('P','T', T_AC8, 'Q', 0, 'Water');
p_H = p_AC8;
s_AC8 = CoolProp.PropsSI('S', 'T', T_AC8, 'P', p_AC8, 'Water');
h_AC8 = CoolProp.PropsSI('H', 'T', T_AC8, 'P', p_AC8, 'Water');
ex_AC8 = (h_AC8 - h_w0) - T_0 * (s_AC8 - s_w0);
% 5.2.3 Status 9 of AC_ALB.
s_AC9 = s_AC8;
p_AC9 = p_AC10;
T_AC9 = T_AC10;
h_AC9 = CoolProp.PropsSI('H', 'S', s_AC9, 'T', T_AC9, 'Water');
ex_AC9 = (h_AC9 - h_w0) - T_0 * (s_AC9 - s_w0);
% 5.2.4 Status 1 of AC_ALB.
p_AC1 = p_L;
T_1wC = CoolProp.PropsSI('T', 'P', p_AC1, 'Q', 0, 'Water') - 273.15;
a0 = -2.00755; a1 = 0.16976; a2 = -3.13336E-3; a3 = 1.97668E-5;
b0 = 124.937;  b1 = -7.7162; b2 = 0.152286;    b3 = -7.9509E-4;
T_AC1 = T_1wC * (a0 + a1 * (y_1*100) + a2 * (y_1*100)^2 + a3 * (y_1*100)^3) + ...
      b0 + b1 * (y_1*100) + b2 * (y_1*100)^2 + b3 * (y_1*100)^3 + 273.15;
y_1str = num2str(y_1);
s_AC1 = CoolProp.PropsSI('S', 'T', T_AC1, 'Q', 0, strcat('INCOMP::LiBr[',y_1str,']'));
h_AC1 = CoolProp.PropsSI('H', 'T', T_AC1, 'Q', 0, strcat('INCOMP::LiBr[',y_1str,']'));
RHO_1 = CoolProp.PropsSI('D', 'T', T_AC1, 'P', p_AC1, strcat('INCOMP::LiBr[',y_1str,']'));
v_1 = 1 / RHO_1;
s_sT0 = CoolProp.PropsSI('S', 'T', T_0, 'P', p_0, strcat('INCOMP::LiBr[',y_1str,']'));
h_sT0 = CoolProp.PropsSI('H', 'T', T_0, 'P', p_0, strcat('INCOMP::LiBr[',y_1str,']'));
ex_AC1 = (h_AC1 - h_sT0) - T_0 * (s_AC1 - s_sT0);
% 5.2.5 Status 4 of AC_ALB.
p_AC4 = p_H;
M_AC2 = M_AC1; y_2 = y_1;
M_AC3 = M_AC2; y_3 = y_2;
M_AC4 = M_AC3 * y_3 / y_4;
T_4wC = CoolProp.PropsSI('T', 'P', p_AC4, 'Q', 0, 'Water') - 273.15;
T_AC4 = T_4wC * (a0 + a1 * (y_4*100) + a2 * (y_4*100)^2 + a3 * (y_4*100)^3) + ...
        b0 + b1 * (y_4*100) + b2 * (y_4*100)^2 + b3 * (y_4*100)^3 + 273.15;
y_4str = num2str(y_4);
s_AC4 = CoolProp.PropsSI('S', 'T', T_AC4, 'Q', 0, strcat('INCOMP::LiBr[',y_4str,']'));
h_AC4 = CoolProp.PropsSI('H', 'T', T_AC4, 'Q', 0, strcat('INCOMP::LiBr[',y_4str,']'));
s_sD0 = CoolProp.PropsSI('S', 'T', T_0, 'P', p_0, strcat('INCOMP::LiBr[',y_4str,']'));
h_sD0 = CoolProp.PropsSI('H', 'T', T_0, 'P', p_0, strcat('INCOMP::LiBr[',y_4str,']'));
ex_AC4 = (h_AC4 - h_sD0) - T_0 * (s_AC4 - s_sD0);
% 5.2.6 Status 2 of AC_ALB.
s_AC2 = s_AC1;
p_AC2 = p_H;
h_AC2 = h_AC1 + v_1 * (p_AC2 - p_AC1);
w_p = M_AC2 * h_AC2 - M_AC1 * h_AC1;
T_AC2 = T_AC1;
ex_AC2 = (h_AC2 - h_sT0) - T_0 * (s_AC2 - s_sT0);
% 5.2.7 Status 5 of AC_ALB.
M_AC5 = M_AC4;
y_5 = y_4;
p_AC5 = p_AC4;
T_AC5 = T_AC4 - (T_AC4 - T_AC2) * Eta_shx;
y_5str = num2str(y_5);
s_AC5 = CoolProp.PropsSI('S', 'T', T_AC5, 'P', p_AC5, strcat('INCOMP::LiBr[',y_5str,']'));
h_AC5 = CoolProp.PropsSI('H', 'T', T_AC5, 'P', p_AC5, strcat('INCOMP::LiBr[',y_5str,']'));
ex_AC5 = (h_AC5 - h_sD0) - T_0 * (s_AC5 - s_sD0);
% 5.2.8 Status 3 of AC_ALB.
p_AC3 = p_H;
Q_shx = M_AC4 * h_AC4 - M_AC5 * h_AC5;
h_AC3 = (Q_shx + M_AC2 * h_AC2) / M_AC3;
y_3str = num2str(y_3);
T_AC3 = CoolProp.PropsSI('T', 'H', h_AC3-5000, 'P', p_AC3, strcat('INCOMP::LiBr[',y_3str,']')); % ????
s_AC3 = CoolProp.PropsSI('S', 'T', T_AC3, 'P', p_AC3, strcat('INCOMP::LiBr[',y_3str,']'));
ex_AC3 = (h_AC3 - h_sT0) - T_0 * (s_AC3 - s_sT0);
% 5.2.9 Status 6 of AC_ALB.
y_6 = y_5;
T_AC6 = T_AC5;
h_AC6 = h_AC5;
p_AC6 = p_AC1;
y_6str = num2str(y_6);
s_AC6 = CoolProp.PropsSI('S', 'T', T_AC6, 'P', p_AC6, strcat('INCOMP::LiBr[',y_6str,']'));
M_AC6 = M_AC4;
ex_AC5 = (h_AC5 - h_sD0) - T_0 * (s_AC5 - s_sD0);
% 5.2.10 Status 7 of AC_ALB.
M_AC7 = M_AC3 - M_AC4;
M_AC8 = M_AC7;
M_AC9 = M_AC8;
M_AC10 = M_AC9;
p_AC7 = p_H;
% Vapor leaving desorber is assumed in equilibrium with ...
% ... incoming solution stream concentration (state 3).
% This is a standard assumption that represents the best possible case.
T_3wC = CoolProp.PropsSI('T', 'P', p_AC3, 'Q', 0, 'Water') - 273.15;
T_3sat = T_3wC * (a0 + a1 * (y_3*100) + a2 * (y_3*100)^2 + a3 * (y_3*100)^3) + ...
         b0 + b1 * (y_3*100) + b2 * (y_3*100)^2 + b3 * (y_3*100)^3 + 273.15;
T_AC7 = T_3sat;
h_AC7 = CoolProp.PropsSI('H', 'P', p_AC7, 'T', T_AC7, 'Water');
s_AC7 = CoolProp.PropsSI('S', 'P', p_AC7, 'T', T_AC7, 'Water');
ex_AC7 = (h_AC7 - h_w0) - T_0 * (s_AC7 - s_w0);
%% 5.3 Area of heat exchanger in desorber A_d.             ???
T_AC11 = T_7; p_AC11 = p_7; p_AC12 = p_AC11;          % Temp of High Temp Smoke from MGT.
h_AC11 = h_g7; s_AC11 = h_g7; ex_AC11 = ex_ph7;
Q_ACd = M_AC4 * h_AC4 + M_AC7 * h_AC7 - M_AC3 * h_AC3;       % W, HT Rate in desorber
syms T_AC12 A_ACd DeltaT_ACd
eq_ACd(1) = DeltaT_ACd == ((T_AC11 - T_AC4) - (T_AC12 - T_AC3)) ...
                        / log((T_AC11 - T_AC4) / (T_AC12 - T_AC3));
eq_ACd(2) = Q_ACd == DeltaT_ACd * K * A_ACd;
eq_ACd(3) = Q_ACd == c_pg * m_g * (T_AC11 - T_AC12);
[ST_AC12, SA_ACd, SDeltaT_ACd] = solve(eq_ACd);
   T_AC12 = double(ST_AC12);
      A_ACd = double(SA_ACd);
DeltaT_ACd = double(SDeltaT_ACd);
C_ACd = Z_A * A_ACd;
% Area of heat exchanger in evaporator A_ACe in AC_ALB.
Q_ACe = M_AC9 * (h_AC10 - h_AC9);
T_AC17 = T_0; p_AC17 = p_0; p_AC18 = p_AC17;
h_AC17 = h_w0; s_AC17 = s_w0; ex_AC17 = E_w0;
T_AC18 = 5 + 273.15;
syms Mcw_AC A_ACe DeltaT_ACe
eq_ACe(1) = DeltaT_ACe == ((T_AC17 - T_AC10) - (T_AC18 - T_AC9)) ...
                        / log((T_AC17 - T_AC10) / (T_AC18 - T_AC9));
eq_ACe(2) = Q_ACe == DELTA_T_e * K * A_ACe;
eq_ACe(3) = Q_ACe == c_pw * Mcw_AC * (T_AC17 - T_AC18);
[SMcw_AC, SA_e, SDeltaT_ACe] = solve(eq_ACe);
  Mcw_AC = double(SMcw_AC);
   A_ACe = double(SA_ACe);
DeltaT_e = double(SDeltaT_ACe);
C_ACcw = Z_cw * Mcw_AC * N * N_y * 3600;
% Calculate the h, s, E_ph of gas in status 12.
h_a12 = CoolProp.PropsSI('H', 'T', T_AC12, 'P', p_AC12, 'Air');
s_a12 = CoolProp.PropsSI('S', 'T', T_AC12, 'P', p_AC12, 'Air');
h_H12 = CoolProp.PropsSI('H', 'T', T_AC12, 'P', p_AC12, 'Water');
s_H12 = CoolProp.PropsSI('S', 'T', T_AC12, 'P', p_AC12, 'Water');
h_C12 = CoolProp.PropsSI('H', 'T', T_AC12, 'P', p_AC12, 'CarbonDioxide');
s_C12 = CoolProp.PropsSI('S', 'T', T_AC12, 'P', p_AC12, 'CarbonDioxide');
h_N12 = CoolProp.PropsSI('H', 'T', T_AC12, 'P', p_AC12, 'Nitrogen');
s_N12 = CoolProp.PropsSI('S', 'T', T_AC12, 'P', p_AC12, 'Nitrogen');
h_g12 = m_ax/m_g * h_a12 + (m_a-m_ax)/m_g * y_H*18/M_gC * h_H12 + ...
       (m_a-m_ax)/m_g * y_C*44/M_gC * h_C12 + (m_a-m_ax)/m_g * y_N*28/M_gC * h_N12;
s_g12 = m_ax/m_g * s_a12 + (m_a-m_ax)/m_g * y_H*18/M_gC * s_H12 + ...
       (m_a-m_ax)/m_g * y_C*44/M_gC * s_C12 + (m_a-m_ax)/m_g * y_N*28/M_gC * s_N12;
ex_AC12 = (h_g12 - h_g0) - T_0 * (s_g12 - s_g0);
C_Ad = Z_A * A_d;
PhiEx_ACd = (ex_AC7 * M_AC7 + ex_AC12 * m_g + ex_AC4 * m_4) / ...
          (ex_AC11 * m_g + ex_AC3 * M_AC3);
%%  5.4 Mathematical Model of R123 Organic Rankine Cycle (ORC_R123) -------------------------------------
% 5.4.0 Status 0 of R123.
s_ORC0 = CoolProp.PropsSI('S', 'T', T_0, 'P', p_0, 'R123');
h_ORC0 = CoolProp.PropsSI('H', 'T', T_0, 'P', p_0, 'R123');
E_ORCph0 = 0;
% 5.4.1 Status 1, 1-2 Pump
Af = 4.643358E2; Bf =  1.625985E3; Cf = -1.333543E3;
Df = 1.986142E3; Ef = -7.172430E2;
T_rORC1 = T_ORC1 ./ Tc_R123;                                          % Reduced Temerature
D_ORC1 = Af + Bf * (1-T_rORC1).^(1/3) + Cf * (1-T_rORC1).^(2/3) + ...
         Df * (1-T_rORC1) + Ef - (1-T_rORC1).^(4/3);
v_ORC1 = D_ORC1 / 1;
% p_ORC1 = 10 * CoolProp.PropsSI('P', 'T', T_ORC1, 'Q', 0, 'R123');
p_ORC1 = CoolProp.PropsSI('P', 'T', T_ORC1, 'Q', 0, 'R123');
s_ORC1 = CoolProp.PropsSI('S', 'T', T_ORC1, 'Q', 0, 'R123');
h_ORC1 = CoolProp.PropsSI('H', 'T', T_ORC1, 'Q', 0, 'R123');
E_ORCph1 = (h_ORC1 - h_ORC0) - T_0 * (s_ORC1 - s_ORC0);
%% 5.4.2 Status 2, 2-3 Evaporator
h_ORC2 = h_ORC1 + v_ORC1 * (p_ORC2 - p_ORC1);
T_ORC2 = T_ORC1;                                    % ???
% T_ORC2 = CoolProp.PropsSI('T', 'H', h_ORC2, 'P', p_ORC2, 'R123');
s_ORC2 = s_ORC1 / ETA_ps;                                             % ???
E_ORCph2 = (h_ORC2 - h_ORC0) - T_0 * (s_ORC2 - s_ORC0);
%% 5.4.3 Status 3, 3-4 Gas Turbine
syms v_3sym
p_ORC3 = p_ORC2;
T_r3 = T_ORC3 ./ Tc_R123;                             % Reduced Temerature
ALPHASqrt3 = 1 + KAPPA_R123 * (1 - sqrt(T_r3));
ALPHA3 = ALPHASqrt3^2;                                % Temp-dependent parameter in PR-EOS
a_T3 = aTc_R123 * ALPHA3;                             % Temp-dependent parameter in PR-EOS, on equation
Sv_3sym = solve(p_ORC3 == R*T_ORC3 / (v_3sym-b_R123) - a_T3 / (v_3sym*(v_3sym+b_R123) + b_R123*(v_3sym-b_R123)));
v_ORC3 = double(Sv_3sym);
v_ORC3 = v_ORC3(imag(v_ORC3)==0);
s_ORC3 = CoolProp.PropsSI('S', 'T', T_ORC3, 'P', p_ORC3, 'R123');
h_ORC3 = CoolProp.PropsSI('H', 'T', T_ORC3, 'P', p_ORC3, 'R123');
E_ORCph3 = (h_ORC3 - h_ORC0) - T_0 * (s_ORC3 - s_ORC0);
%% 5.4.4 Status 4, 4-1 Condenser
p_ORC4 = p_ORC1;
T_ORC4 = CoolProp.PropsSI('T', 'P', p_ORC4, 'Q', 1, 'R123');
%{
syms T_4sym
A =  1.656333E3;  B = -2.480583E6; C = 1.792522E1;
D = -8.868380E-2; E =  4.617861E2; F = 1.666667E3;
ST_4sym = solve(log10(p_ORC4/1000) == A + B/T_4sym + C * log10(T_4sym) + D * T_4sym + ...
                                       E * ((F-T_4sym)/T_4sym) * log10(F-T_4sym));
T_ORC4 = double(ST_4sym);
T_ORC4 = real(T_ORC4);
% T_ORC4 = T_ORC4(imag(T_ORC4)==0);
syms v_4sym
T_r4 = T_ORC4 ./ Tc_R123;                             % Reduced Temerature
ALPHASqrt4 = 1 + KAPPA_R123 * (1 - sqrt(T_r4));
ALPHA4 = ALPHASqrt4^2;                                % Temp-dependent parameter in PR-EOS
a_T4 = aTc_R123 * ALPHA4;                             % Temp-dependent parameter in PR-EOS, on equation
Sv_4sym = solve(p_ORC4 == R*T_ORC4 / (v_4sym-b_R123) - a_T4 / (v_4sym*(v_4sym+b_R123) + b_R123*(v_4sym-b_R123)));
v_ORC4 = double(Sv_4sym);
v_ORC4 = v_ORC4(imag(v_ORC4)==0);
%}
s_ORC4 = s_ORC3;
h_ORC4 = CoolProp.PropsSI('H', 'P', p_ORC4, 'Q', 1, 'R123');
E_ORCph4 = (h_ORC4 - h_ORC0) - T_0 * (s_ORC4 - s_ORC0);
%% 5.6 Area of Heat Exchange in Evaporator within ORC_R123.
T_ORC19 = T_AC12; p_ORC19 = p_AC12; p_ORC20 = p_ORC19;
h_ORC19 = h_AC12; s_ORC19 = s_AC12; E_ORCph19 = ex_AC12;
Q_ORC1 = M_ORC * (h_ORC3 - h_ORC2);
syms T_ORC20 A_ORC1 DELTA_T_ORC1
eq_ORC1(1) = DELTA_T_ORC1 == ((T_ORC19 - T_ORC3) - (T_ORC20 - T_ORC2)) / ...
                          log((T_ORC19 - T_ORC3) / (T_ORC20 - T_ORC2));
eq_ORC1(2) = Q_ORC1 == DELTA_T_ORC1 * K * A_ORC1;
eq_ORC1(3) = Q_ORC1 == m_g * c_pg * (T_ORC19 - T_ORC20);
[ST_ORC20,SA_ORC1,SDELTA_T_ORC1] = solve(eq_ORC1);
     T_ORC20 = double(ST_ORC20);
      A_ORC1 = double(SA_ORC1);
DELTA_T_ORC1 = double(SDELTA_T_ORC1);
% Calculate the h, s, E_ph of gas in status 20 in ORC_R123.
h_a20 = CoolProp.PropsSI('H', 'T', T_ORC20, 'P', p_ORC20, 'Air');
s_a20 = CoolProp.PropsSI('S', 'T', T_ORC20, 'P', p_ORC20, 'Air');
h_H20 = CoolProp.PropsSI('H', 'T', T_ORC20, 'P', p_ORC20, 'Water');
s_H20 = CoolProp.PropsSI('S', 'T', T_ORC20, 'P', p_ORC20, 'Water');
h_C20 = CoolProp.PropsSI('H', 'T', T_ORC20, 'P', p_ORC20, 'CarbonDioxide');
s_C20 = CoolProp.PropsSI('S', 'T', T_ORC20, 'P', p_ORC20, 'CarbonDioxide');
h_N20 = CoolProp.PropsSI('H', 'T', T_ORC20, 'P', p_ORC20, 'Nitrogen');
s_N20 = CoolProp.PropsSI('S', 'T', T_ORC20, 'P', p_ORC20, 'Nitrogen');
h_g20 = m_ax/m_g * h_a20 + (m_a-m_ax)/m_g * y_H*18/M_gC * h_H20 + ...
       (m_a-m_ax)/m_g * y_C*44/M_gC * h_C20 + (m_a-m_ax)/m_g * y_N*28/M_gC * h_N20;
s_g20 = m_ax/m_g * s_a20 + (m_a-m_ax)/m_g * y_H*18/M_gC * s_H20 + ...
       (m_a-m_ax)/m_g * y_C*44/M_gC * s_C20 + (m_a-m_ax)/m_g * y_N*28/M_gC * s_N20;
E_ORCph20 = (h_g20 - h_g0) - T_0 * (s_g20 - s_g0);
PhiEx_ORC1 = (E_ORCph4 * M_ORC + E_ORCph20 * m_g) / (E_ORCph19 * m_g + E_ORCph3 * M_ORC);
%% 5.5 Mathematical Model of HP_R410a. -----------------------------------------------------------
%% 3. (1)
P_HP1 = CoolProp.PropsSI('P', 'T', T_HP1, 'Q', 1, 'R410a');
% Solve for V_HP1
syms V_HP1sym
Tr_HP1 = T_HP1 ./ Tc_R140a;                             % Reduced Temerature
ALPHASqrt1 = 1 + Kappa_R410a * (1 - sqrt(Tr_HP1));
ALPHA4 = ALPHASqrt1^2;                         % Temp-dependent para in PR-EOS
a_T1 = aTc_R410a * ALPHA4;                          % Temp-dependent para in PR-EOS
SV_HP1sym = solve(P_HP1 == R*T_HP1 / (V_HP1sym-b_R410a) - ...
                       a_T1 / (V_HP1sym*(V_HP1sym+b_R410a) + b_R410a*(V_HP1sym-b_R410a)));
V_HP1 = double(SV_HP1sym);
V_HP1 = V_HP1(imag(V_HP1)==0);
%
S_HP1 = CoolProp.PropsSI('S', 'T', T_HP1, 'Q', 1, 'R134a');
h_HP1 = CoolProp.PropsSI('H', 'T', T_HP1, 'Q', 1, 'R134a');
%% 4. (4)
P_HP4 = P_HP1; T_HP4 = T_HP1;
% Solve for V_HP4.
syms V_HP4sym
Tr_HP4 = T_HP4 ./ Tc_R140a;                             % Reduced Temerature
ALPHASqrt4 = 1 + Kappa_R410a * (1 - sqrt(Tr_HP4));
ALPHA4 = ALPHASqrt4^2;                         % Temp-dependent para in PR-EOS
a_T4 = aTc_R410a * ALPHA4;                          % Temp-dependent para in PR-EOS
SV_HP4sym = solve(P_HP4 == R*T_HP4 / (V_HP4sym-b_R410a) - ...
                a_T4 / (V_HP4sym*(V_HP4sym+b_R410a) + b_R410a*(V_HP4sym-b_R410a)));
V_HP4 = double(SV_HP4sym);
V_HP4 = V_HP4(imag(V_HP4)==0);
%
S_HP4 = CoolProp.PropsSI('S', 'T', T_HP4, 'P', P_HP4, 'R410a');
h_HP4 = CoolProp.PropsSI('H', 'T', T_HP4, 'P', P_HP4, 'R410a');
%% (3)
% <8> solve V_HP3 through Equation for Density of the Saturated Liquid.
Af = 1.000000; Bf =  1.984734;    Cf = -1.767593E-01;
Df = 1.819972; Ef = -7.171684E-1;
Tr_HP3 = T_HP3 ./ Tc_R140a;                             % Reduced Temerature
D_HP3 = Dc_R410a * (Af + Bf * (1-Tr_HP3).^(1/3) + Cf * (1-Tr_HP3).^(2/3) + ...
                  Df * (1-Tr_HP3) + Ef - (1-Tr_HP3).^(4/3));
V_HP3 = D_HP3 / 1;
% <9> Solve for P_HP3.
syms P_HP3sym
Tr_HP3 = T_HP3 ./ Tc_R140a;                             % Reduced Temerature
ALPHASqrt3 = 1 + Kappa_R410a * (1 - sqrt(Tr_HP3));
ALPHA3 = ALPHASqrt3^2;                         % Temp-dependent para in PR-EOS
a_T3 = aTc_R410a * ALPHA3;                          % Temp-dependent para in PR-EOS
SP_HP3sym = solve(P_HP3sym == R*T_HP3 / (V_HP3-b_R410a) - ...
                a_T3 / (V_HP3*(V_HP3+b_R410a) + b_R410a*(V_HP3-b_R410a)));
P_HP3 = double(SP_HP3sym);
P_HP3 = P_HP3(imag(P_HP3)==0);
%
S_HP3 = CoolProp.PropsSI('S', 'T', T_HP3, 'Q', 0, 'R410a');
h_HP3 = CoolProp.PropsSI('H', 'T', T_HP3, 'Q', 0, 'R410a');
%% (2)
P_HP2 = P_HP3; S_HP2 = S_HP1;
% <11> Solve for V_HP2.
syms V_HP2sym
Tr_HP2 = T_HP2 ./ Tc_R140a;                             % Reduced Temerature
ALPHASqrt2 = 1 + Kappa_R410a * (1 - sqrt(Tr_HP2));
ALPHA2 = ALPHASqrt2^2;                         % Temp-dependent para in PR-EOS
a_T2 = aTc_R410a * ALPHA2;                          % Temp-dependent para in PR-EOS
SV_HP2sym = solve(P_HP2 == R*T_HP2 / (V_HP2sym-b_R410a) - ...
                a_T2 / (V_HP2sym*(V_HP2sym+b_R410a) + b_R410a*(V_HP2sym-b_R410a)));
V_HP2 = double(SV_HP2sym);
V_HP2 = V_HP2(imag(V_HP2)==0);
%
S_HP2 = CoolProp.PropsSI('S', 'T', T_HP2, 'P', P_HP2, 'R410a');
h_HP2 = CoolProp.PropsSI('H', 'T', T_HP2, 'P', P_HP2, 'R410a');
% Input work of compressor Wi_HPc.
Wi_HPc = M_HP * (h_HP2 - h_HP1);
%% 5.5 Mathematical Model of HP_R410a. -----------------------------------------------------------
% 3. (1)
% <1> Solve T_VCC4 through Equation for Saturated Vapor Pressure.
A =  4.069889E1;  B = -2.362540E3;  C = -1.306883E1;
D =  7.616005E-3; E =  2.342564E-1; F =  3.761111E2;
P_HP1 = 10^(A + B./T_VCC1 + C * log10(T_VCC1) + D * T_VCC1 + ...
          E * ((F-T_VCC1)./T_VCC1) * log10(F-T_VCC1)) * 1000;
% <2> Solve for V_HP1
syms V_HP1sym
Tr_HP1 = T_VCC1 ./ Tc_R134a;                             % Reduced Temerature
ALPHASqrt1 = 1 + Kappa_R134a * (1 - sqrt(Tr_HP1));
ALPHA4 = ALPHASqrt1^2;                                   % Temperature-dependent parameter in PR-EOS
aT_HP1 = aTc_R134a * ALPHA4;                             % Temperature-dependent parameter in PR-EOS (dimensions depend on equation)
SV_HP1sym = solve(P_HP1 == R*T_VCC1 / (V_HP1sym-b_R134a) - ...
                  aT_HP1 / (V_HP1sym*(V_HP1sym+b_R134a) + b_R134a*(V_HP1sym-b_R134a)));
V_HP1 = double(SV_HP1sym);
V_HP1 = V_HP1(imag(V_HP1)==0);
%
S_VCC1 = CoolProp.PropsSI('S', 'T', T_VCC1, 'Q', 1, 'R123');
H_VCC1 = CoolProp.PropsSI('H', 'T', T_VCC1, 'Q', 1, 'R123');
% 4. (4)
P_HP4 = P_HP1; T_VCC4 = T_VCC1;
% <5> Solve for V_HP4.
syms V_HP4sym
Tr_HP4 = T_VCC4 ./ Tc_R134a;                             % Reduced Temerature
ALPHASqrt4 = 1 + Kappa_R134a * (1 - sqrt(Tr_HP4));
ALPHA4 = ALPHASqrt4^2;                                   % Temperature-dependent parameter in PR-EOS
aT_HP4 = aTc_R134a * ALPHA4;                             % Temperature-dependent parameter in PR-EOS (dimensions depend on equation)
SV_HP4sym = solve(P_HP4 == R*T_VCC4 / (V_HP4sym-b_R134a) - ...
                  aT_HP4 / (V_HP4sym*(V_HP4sym+b_R134a) + b_R134a*(V_HP4sym-b_R134a)));
V_HP4 = double(SV_HP4sym);
V_HP4 = V_HP4(imag(V_HP4)==0);
%
S_VCC4 = CoolProp.PropsSI('S', 'T', T_VCC4, 'P', P_HP4, 'R123');
H_VCC4 = CoolProp.PropsSI('H', 'T', T_VCC4, 'P', P_HP4, 'R123');
% (3)
% <8> solve V_HP3 through Equation for Density of the Saturated Liquid.
Af =  5.281464E2; Bf =  7.551834E2; Cf = 1.028676E3;
Df = -9.491172E2; Ef = 5.935660E2;
Tr_HP3 = T_VCC3 ./ Tc_R134a;                             % Reduced Temerature
RHO_3 = Af + Bf * (1-Tr_HP3).^(1/3) + Cf * (1-Tr_HP3).^(2/3) + ...
        Df * (1-Tr_HP3) + Ef - (1-Tr_HP3).^(4/3);
V_HP3 = RHO_3 / 1;
% <9> Solve for P_HP3.
syms P_HP3sym
Tr_HP3 = T_VCC3 ./ Tc_R134a;                             % Reduced Temerature
ALPHASqrt3 = 1 + Kappa_R134a * (1 - sqrt(Tr_HP3));
ALPHA3 = ALPHASqrt3^2;                         % Temperature-dependent parameter in PR-EOS
aT_HP3 = aTc_R134a * ALPHA3;                          % Temperature-dependent parameter in PR-EOS (dimensions depend on equation)
SP_HP3sym = solve(P_HP3sym == R*T_VCC3 / (V_HP3-b_R134a) - ...
                aT_HP3 / (V_HP3*(V_HP3+b_R134a) + b_R134a*(V_HP3-b_R134a)));
P_HP3 = double(SP_HP3sym);
P_HP3 = P_HP3(imag(P_HP3)==0);
%
S_VCC3 = CoolProp.PropsSI('S', 'T', T_VCC3, 'Q', 0, 'R134a');
H_VCC3 = CoolProp.PropsSI('H', 'T', T_VCC3, 'Q', 0, 'R134a');
% (2)
P_HP2 = P_HP3; S_VCC2 = S_VCC1;
% <11> Solve for V_HP2.
syms V_HP2sym
Tr_HP2 = T_VCC2 ./ Tc_R134a;                             % Reduced Temerature
ALPHASqrt2 = 1 + Kappa_R134a * (1 - sqrt(Tr_HP2));
ALPHA2 = ALPHASqrt2^2;                         % Temperature-dependent parameter in PR-EOS
aT_HP2 = aTc_R134a * ALPHA2;                          % Temperature-dependent parameter in PR-EOS (dimensions depend on equation)
SV_HP2sym = solve(P_HP2 == R*T_VCC2 / (V_HP2sym-b_R134a) - ...
                  aT_HP2 / (V_HP2sym*(V_HP2sym+b_R134a) + b_R134a*(V_HP2sym-b_R134a)));
V_HP2 = double(SV_HP2sym);
V_HP2 = V_HP2(imag(V_HP2)==0);
%
S_VCC2 = CoolProp.PropsSI('S', 'T', T_VCC2, 'P', P_HP2, 'R134a');
H_VCC2 = CoolProp.PropsSI('H', 'T', T_VCC2, 'P', P_HP2, 'R134a');
% Input work of compressor Wi_HPp.
Wi_HPp = M_VCC * (H_VCC2 - H_VCC1);
%% 6. Economic Model ------------------------------------------------------------------------------------
%% 6.1 Economic Model of MGT.
ETA_t = Wp_MGT ./ (m_f * LHV);                             % Thermo Efficiency of MGT
P = A_MGT + B_MGT * log(Wp_MGT/1000) + C_MGT * exp(ETA_t); % 单位功率价格, P = 0.9;
z_MT = P * Wp_MGT / 1000;                                  % 燃气轮机的购置费
Q_HRSG = (h_9 - h_8) * Ms_MGT;                                % 余热锅炉热负荷
% z_HRSG = z_w * (Q_HRSG ./ Q_w)^AlphaC_MGT;                    % 余热锅炉购置费 ???
z_HRSG = (Q_HRSG * 4.851)^AlphaC_MGT;                           % 余热锅炉购置费
C_E_MGT = Wp_MGT * ETA_G * Z_E * N * N_y * 1000;           % Profit from Electricity generated by MGT.

%% 6.2 Economic Model of AC_ALB. ---------------------------------------------------------------------------
% Area of heat exchanger in solution heat exchanger A_ACs in AC_ALB.
Q_ACs = M_AC2 * (h_AC3 - h_AC2);
DeltaT_ACs = ((T_AC4 - T_AC3) - (T_AC5 - T_AC2)) / log((T_AC4 - T_AC3) / (T_AC5 - T_AC2));
A_ACs = Q_ACs / (DeltaT_ACs * K);
Ca_ACs = Z_A * A_ACs;
% Area of heat exchanger in condenser A_ACc in AC_ALB.
Q_c = M_AC7 * (h_AC7 - h_AC8);                             % W, HT Rate of heat exchanger in desorber
T_AC15 = T_0;
T_AC16 = T_0H;
DeltaT_ACc = ((T_AC7 - T_AC16) - (T_AC8 - T_AC15)) / log((T_AC7 - T_AC16) / (T_AC8 - T_AC15));
A_c = Q_c / DeltaT_ACc / K;
h_AC15 = CoolProp.PropsSI('H', 'T', T_AC15, 'P', p_0, 'Water');
h_AC16 = CoolProp.PropsSI('H', 'T', T_AC16, 'P', p_0, 'Water');
m_c = Q_c / (h_AC16 - h_AC15);
C_Ac = Z_A * A_c;
C_Wc = Z_W * m_c * N * N_y * 3600;
% Calculate the h, s, E_ph of water in status 18.
h_AC18 = CoolProp.PropsSI('H', 'T', T_AC18, 'P', p_AC18, 'Water');
s_AC18 = CoolProp.PropsSI('S', 'T', T_AC18, 'P', p_AC18, 'Water');
ex_AC18 = (h_g18 - h_w0) - T_0 * (s_g18 - s_w0);
E_ACe = Q_e * (T_0/T_AC18 - 1);
PhiEx_ACe = E_ACe * Mcw_AC / (ex_AC9 * M_AC9 + ex_AC17 * Mcw_AC);                 % ???
C_Ae = Z_A * A_e;
% Area of heat exchanger in absorber A_ACa in AC_ALB.
Q_ACa = M_AC10 * h_AC10 + M_AC6 * h_AC6 - M_AC1 * h_AC1; % W, HT rate of heat exchanger in desorber
T_AC13 = T_0;
T_AC14 = T_0H;
DeltaT_ACa = ((T_AC6 - T_AC14) - (T_AC1 - T_AC13)) / log((T_AC6 - T_AC14) / (T_AC1 - T_AC13));
A_ACa = Q_ACa / DeltaT_ACa / K;
h_AC13 = CoolProp.PropsSI('H', 'T', T_AC13, 'P', p_0, 'Water');
h_AC14 = CoolProp.PropsSI('H', 'T', T_AC14, 'P', p_0, 'Water');
m_ACa = Q_ACa / (h_AC14 - h_AC13);
C_ACa = Z_A * A_ACa;
Cw_ACa = Z_W * m_ACa * N * N_y * 3600;
W_ACa = M_AC1 * (h_AC2 - h_AC1);                        % W,   Electricity required for pump.
Ce_ACa = Z_E * W_ACa * N * N_y * 1000;                    % RMB, Cost of electricity required for pump.
z_AC = Q_ACe * (Q_wAC / Z_wAC)^ALPHA_wAC;                 % RMB, Estimated purchased price of AC_ALB.
%% 6.3 Economic Model of ORC_R123. -----------------------------------------------------------------------
% Area of Heat Exchange in Condenser within ORC_R123.
T_ORC21 = T_0;
T_ORC22 = T_0H;
Q_ORC2 = M_ORC * (h_ORC4 - h_ORC1);
DELTA_T_ORC2 = ((T_ORC4 - T_ORC22) - (T_ORC1 - T_ORC21)) / log((T_ORC4 - T_ORC22) / (T_ORC1 - T_ORC21));
A_ORC2 = Q_ORC2 / DELTA_T_ORC2 / K;
h_ORC21 = CoolProp.PropsSI('H', 'T', T_ORC21, 'P', p_0, 'Water');
h_ORC22 = CoolProp.PropsSI('H', 'T', T_ORC22, 'P', p_0, 'Water');
m_ORC2 = Q_ORC2 / (h_ORC22 - h_ORC21);                   % kg/s, Amount of supplying cooling water.
W_ORC = M_ORC * (h_ORC3 - h_ORC4);                       % W,    Output work of ORC_R123.
C_E_ORC = W_ORC * ETA_G * Z_E * N * N_y * 1000;          % RMB,  Profit of output work of ORC_R123.
W_ORCp = M_ORC * (h_ORC2 - h_ORC1);                      % W,    Required Work of Pump in ORC_R123.
C_EpORC = Z_E * W_ORCp * N * N_y * 1000;                 % RMB,  Cost of required Work of Pump in ORC_R123.
C_ORC1 = A_ORC1 * Z_A;                                   % RMB,  Cost of area for heat exchange in evaporator.
C_ORC2 = A_ORC2 * Z_A;                                   % RMB,  Cost of area for heat exchange in evaporator.
C_ORCw = m_ORC2 * Z_W;                                   % RMB,  Cost of supplying cooling water.
z_ORC = Z_wORC * (W_ORC / Q_wORC)^ALPHA_wORC;            % RMB,  Non-energy cost of ORC_R123. ??????
PhiEx_ORC = W_ORC / ((ex_ph3 - ex_ph4) * M_ORC);
%% 6.4 Economic Model of HP_R410a.
% Area of heat exchange in condenser A_HPc.
Q_HPc = M_HP * (h_HP2 - h_HP3);  % W, HT Rate in desorber
Ms_HP = Q_HPc / (h_9 - h_0);
DeltaT_HPc == ((T_HP2 - T_9) - (T_HP3 - T_0)) ...
                / log((T_HP2 - T_9) / (T_HP3 - T_0));
A_HPc = Q_HPc / DeltaT_HPc / K;
% Area of heat exchange in evaporator A_HPe.
Q_HPe = M_HP * (h_HP1 - h_HP4);
DeltaT_HPe == ((T_0 - T_HP1) - (T_0L - T_HP4)) ...
                / log((T_0 - T_HP1) / (T_0L - T_HP4));
A_HPe = Q_HPe / DeltaT_HPe / K;
%% 4. 经济模型 ---------------------------------------------------------------
% Area of heat exchange in condenser Q_VCCc
Q_VCCc = M_VCC * (H_VCC2 - H_VCC3);                        % W, HT Rate in desorber
Mcw_VCC = Q_VCCc / (h_0H - h_0);
DeltaT_VCCc == ((T_VCC2 - T_0H) - (T_VCC3 - T_0)) ...
                / log((T_VCC2 - T_0H) / (T_VCC3 - T_0));
A_VCCc = Q_VCCc / DeltaT_VCCc / K;
% Area of heat exchange in evaporator Q_VCCe
Q_VCCe = M_VCC * (H_VCC1 - H_VCC4);
m_VCCcw = Q_VCCc / (h_0 - h_cw);
DeltaT_VCCc == ((T_0 - T_VCC1) - (T_cw - T_VCC4)) ...
                / log((T_0 - T_VCC1) / (T_cw - T_VCC4));
A_VCCe = Q_VCCe / DeltaT_VCCc / K;
%% 5. Define Objective Function -------------------------------------------------------------------------
f = 1000000 - Z_f * m_f * LHV - CRF * PhiM * (z_MT + z_HRSG) / (3600 * N) + C_E_MGT + ... % MGT Objective
    500000 - CRF * PhiM * z_AC / (3600 * N) + ...                                    % AC_ALB Objective 1
    500000 - (Ca_ACs + Ca_ACd + Cw_ACc + Ca_ACc + Cw_ACa + Ca_ACa) - Ce_ACa + C_Ce + ...          % AC_ALB Objective 2
    500000 - CRF * PhiM * z_ORC / (3600 * N) + ...                                 % ORC_R123 Objective 1
    500000 - (C_ORC1 + C_ORC2 + C_ORCw + C_EpORC) + C_E_ORC;                      % ORC_R123 Objective 2
%% Plot the result of exergy analysis.
% Plot the bar chart of exergy efficiency of equipments in MGT.
BAR_MGTxE = categorical({'PhiEx_C_C','PhiEx_G_T','PhiEx_A_C','PhiEx_R_E','PhiEx_H_R_S_G', ...
                        'PhiEx_ACd','PhiEx_ACe','PhiEx_ORC1','PhiEx_ORC'});
BAR_MGTyE = [PhiEx_CC, PhiEx_GT, PhiEx_AC, PhiEx_RE, PhiEx_HRSG, PhiEx_ACd, PhiEx_ACe, PhiEx_ORC1, PhiEx_ORC];
BAR_MGT = bar(BAR_MGTxE, BAR_MGTyE, 0.5);
ylim([0 1.1]);
title('绿色能源岛各组分㶲效率条形图 Bar Chart of Exergy Efficiency of GEI', ...
      'FontSize',20,'FontWeight','bold');
