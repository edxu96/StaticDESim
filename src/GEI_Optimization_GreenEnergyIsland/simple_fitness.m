% Title: Thermo-economic Optimization of Distributed Energy System in Green Energy Island.
% Based on the theory of nolinear equality and inequality constraints.
% Method: Genetic Algorithm within MATLAB Global Optimization Toolbox.
% Version: 1.4, 2018.6.1, Jie Xu.
% SubTitle: Define Objective Function.
% 1. Micro Gas Turbine (MGT)
% 2. Aqueous Lithium-Bromine Single-Effect Absorption Chiller (AC_ALB)
% 3. R123 Organic Recycle Cycle (ORC_R123)

clear; clc;

% function f = simple_fitness(x)

x(1) = 633186.25;
x(2) = 0.837;
x(3) = 0.917;
x(4) = 769.395751953125;
x(5) = 1300;
x(6) = 0.03;
x(7) = 60 * 1000;

x(8) = 4 + 273.15;
x(9) = 30 + 273.15;
x(10) = 0.5322;
x(11) = 0.6711;
x(12) = 0.05;           %
x(13) = 0.8;         %

x(14) = 0.05;
x(15) = 200 + 273.15;
x(16) = 6 * 101.325 * 1000;
x(17) = 80 + 273.15;

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
c_f = 4 * 10^(-9);                                % 单位能量的燃料价格
A_MGT = 837.68; B_MGT = -143.22; C_MGT = 491.79;  % 回归分析统计的价格方程系数
ALPHA = 0.73;                                     % 费用估价系数
K_a = 1.4;                                        % 空气的等熵系数
K_g = 1.33;                                       % 燃气的等熵系数
CRF = 0.182;                                      % 投资回收系数
phi = 1.06;                                       % 维护因子
h_8p = 640039.2;                                  % 回水的焓值 (503.97 + (589.30-503.97)/20*16.85) * 1000
h_8 = 293316;                                     % 余热回收后的焓值 (251.56 + (335.29-251.56)/20*10) * 1000
h_9 = 2748.6 * 1000;                              % 过热水蒸汽焓值
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
R_gR123 = 0.0544 * 1000;                          % J/(mol*K), Gas Constant.
M_R123 = 152.93 / 1000;                           % kg/mol, Molecular Weight of R123.
T_cR123 = 456.83;                                 % K, Critical Temp of R123.
p_cR123 = 3668.0 * 1000;                          % Pa, Critical Pressure of R123.
T_Cin = 25 + 273.15;                              % K, Inlet temperature of cold source.
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
a_TcR123 = 0.457235529 * (R_gR123 * ...
           T_cR123)^2 ./ p_cR123;                 % Critical Point Restriction "a(T_c)"
b_R123 = 0.077796074 * R_gR123 * ...              % m^3/mol, Critical Point Restriction "b",
         T_cR123 ./ p_cR123;                      % Temp-independent parameter in PR-EOS
Q_wORC = 23.6 * 1000;                             % W,   ????
Z_wORC = 385600;                                  % RMB, ????
ALPHA_wORC = 0.73;                                %      ????
ETA_G = 0.9;
%% 3. Pre-defined Condition. ----------------------------------------------------------------------------
N = 8000;                                         % Operating Hours in Unit Years
N_y = 10;                                         % Unit Years
T_s = T_9; p_s = p_9;                             % Supplying steam.
T_cw = 4 + 273.15; p_cw = p_0;                    % Supplying cooling water.
%% Constants after Calculation.
h_0 = CoolProp.PropsSI('H', 'T', T_0, 'P', p_0, 'Water');
h_s = CoolProp.PropsSI('H', 'T', T_s, 'P', p_s, 'Water');
h_cw = CoolProp.PropsSI('H', 'T', T_cw, 'P', p_cw, 'Water');
%% 4. Decision Variables. -------------------------------------------------------------------------------
% 4.1 Decision Variables in MGT.
      p_2 = x(1);               % Outlet Pressure of Air Compressor
Eta_MGTac = x(2);               % Isentropic Efficiency of Air Compressor
 Eta_MGTt = x(3);               % Isentropic Efficiency of Gas Turbine
      T_3 = x(4);               % Outlet Temp of Air from Regenerator
      T_4 = x(5);               % Inlet Temp of Gas Turbine
      m_MGTs = x(6);               % Fluid Rate of Saturated Steam from HRSG
   W_pMGT = x(7);               % Net Power from Micro Gas Turbine
% 4.2 Decision Variables in AC_ALB.
   T_AC10 = x(8);               % K,    Outlet temperature of Evaporator.
    T_AC8 = x(9);               % K,    Outlet temperature of Condenser.
      y_1 = x(10);              %       Mass Fraction of Outlet Solution from Absorber.
      y_4 = x(11);              %       Mass Fraction of Outlet Solution from Desorber.
    m_AC1 = x(12);              % kg/s, Fluid Rate Outlet Solution from Absorber.
  ETA_shx = x(13);              % Solution Heat Exchanger Ratio.
% 4.3 Decision Variables in ORC_R123.
    m_ORC = x(14);              % Fluid Rate.
   T_ORC3 = x(15);              % inlet temperature of turbine
   p_ORC2 = x(16);              % outlet pressure of pump / inlet pressure of turbine
   T_ORC1 = x(17);              % Outlet temperature of condenser.


%% 5.1 Mathematical Model of Micro Gas Turbine (MGT) ----------------------------------------------------
T_2 = T_1 .* (1 + 1./Eta_MGTac * ...
      ((p_2./p_1)^((K_a-1)./K_a) - 1));                         % (1)  AC
p_3 = p_2 * (1 - DELTA_p_aRE);                                  % (7)  RE
p_4 = p_3 * (1 - DELTA_p_cc);                                   % (5)  CC
p_6 = p_7 ./ (1 - DELTA_p_HRSG);                                % (14) HRSG
p_5 = p_6 ./ (1 - DELTA_p_gRE);                                 % (8)  RE
T_5 = T_4 * (1 - Eta_MGTt * (1 - (p_4./p_5)^((1-K_g)./K_g)));     % (9)  GT
H = (c_pg * T_4 - ETA_CC * LHV) ./ (c_pa * T_3 - ETA_CC * LHV); % (4) CC, H = m_a / m_g;
m_g = W_pMGT ./ (c_pg*(T_4-T_5) - c_pa*(T_2-T_1)*H);            % (2) (10)
m_a = H * m_g;                                                  % H = m_a / m_g;
m_f = m_g - m_a;                                                % (3)  CC
T_6 = T_5 - m_a * c_pa * (T_3 - T_2) ./ (m_g * c_pg);           % (6)  RE
T_7p = T_6 - m_MGTs * (h_9 - h_8p) ./ (m_g * c_pg);                % (12) HRSG
T_7 = T_6 - m_MGTs * (h_9 - h_8) ./ (m_g * c_pg);                  % (13) HRSG
W_AC = m_a * c_pa * (T_2 - T_1);                                % (2)  AC
W_GT = W_pMGT + W_AC;                                           % (11) GT

%% 5.2 Mathematical Model of Aqueous Lithium-Bromide Absorption Chiller (AC_ALB)-------------------------
% 5.2.1 Status 10 of AC_ALB.
p_AC10 = CoolProp.PropsSI('P','T', T_AC10, 'Q', 1, 'Water');
p_L = p_AC10;
s_AC10 = CoolProp.PropsSI('S', 'T', T_AC10, 'P', p_AC10, 'Water');
h_AC10 = CoolProp.PropsSI('H', 'T', T_AC10, 'P', p_AC10, 'Water');
% 5.2.2 Status 8 of AC_ALB.
p_AC8 = CoolProp.PropsSI('P','T', T_AC8, 'Q', 0, 'Water');
p_H = p_AC8;
s_AC8 = CoolProp.PropsSI('S', 'T', T_AC8, 'P', p_AC8, 'Water');
h_AC8 = CoolProp.PropsSI('H', 'T', T_AC8, 'P', p_AC8, 'Water');
% 5.2.3 Status 9 of AC_ALB.
s_AC9 = s_AC8;
p_AC9 = p_AC10;
T_AC9 = T_AC10;
h_AC9 = CoolProp.PropsSI('H', 'S', s_AC9, 'T', T_AC9, 'Water');
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
% 5.2.5 Status 4 of AC_ALB.
p_AC4 = p_H;
m_AC2 = m_AC1; y_2 = y_1;
m_AC3 = m_AC2; y_3 = y_2;
m_AC4 = m_AC3 * y_3 / y_4;
T_4wC = CoolProp.PropsSI('T', 'P', p_AC4, 'Q', 0, 'Water') - 273.15;
T_AC4 = T_4wC * (a0 + a1 * (y_4*100) + a2 * (y_4*100)^2 + a3 * (y_4*100)^3) + ...
        b0 + b1 * (y_4*100) + b2 * (y_4*100)^2 + b3 * (y_4*100)^3 + 273.15;
y_4str = num2str(y_4);
s_AC4 = CoolProp.PropsSI('S', 'T', T_AC4, 'Q', 0, strcat('INCOMP::LiBr[',y_4str,']'));
h_AC4 = CoolProp.PropsSI('H', 'T', T_AC4, 'Q', 0, strcat('INCOMP::LiBr[',y_4str,']'));
% 5.2.6 Status 2 of AC_ALB.
s_AC2 = s_AC1;
p_AC2 = p_H;
h_AC2 = h_AC1 + v_1 * (p_AC2 - p_AC1);
w_p = m_AC2 * h_AC2 - m_AC1 * h_AC1;
T_AC2 = T_AC1;
% 5.2.7 Status 5 of AC_ALB.
m_AC5 = m_AC4;
y_5 = y_4;
p_AC5 = p_AC4;
T_AC5 = T_AC4 - (T_AC4 - T_AC2) * ETA_shx;
y_5str = num2str(y_5);
s_AC5 = CoolProp.PropsSI('S', 'T', T_AC5, 'P', p_AC5, strcat('INCOMP::LiBr[',y_5str,']'));
h_AC5 = CoolProp.PropsSI('H', 'T', T_AC5, 'P', p_AC5, strcat('INCOMP::LiBr[',y_5str,']'));
% 5.2.8 Status 3 of AC_ALB.
p_AC3 = p_H;
Q_shx = m_AC4 * h_AC4 - m_AC5 * h_AC5;
h_AC3 = (Q_shx + m_AC2 * h_AC2) / m_AC3;
y_3str = num2str(y_3);
T_AC3 = CoolProp.PropsSI('T', 'H', h_AC3-5000, 'P', p_AC3, strcat('INCOMP::LiBr[',y_3str,']')); % ????
s_AC3 = CoolProp.PropsSI('S', 'T', T_AC3, 'P', p_AC3, strcat('INCOMP::LiBr[',y_3str,']'));
% 5.2.9 Status 6 of AC_ALB.
y_6 = y_5;
T_AC6 = T_AC5;
h_AC6 = h_AC5;
p_AC6 = p_AC1;
y_6str = num2str(y_6);
s_AC6 = CoolProp.PropsSI('S', 'T', T_AC6, 'P', p_AC6, strcat('INCOMP::LiBr[',y_6str,']'));
m_AC6 = m_AC4;
% 5.2.10 Status 7 of AC_ALB.
m_AC7 = m_AC3 - m_AC4;
m_AC8 = m_AC7;
m_AC9 = m_AC8;
m_AC10 = m_AC9;
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
%% Area of heat exchanger in desorber A_ACd.             ???
T_AC11 = T_7; p_AC11 = p_7; p_AC12 = p_AC11;               % Temp of High Temp Smoke from MGT.
Q_ACd = m_AC4 * h_AC4 + m_AC7 * h_AC7 - m_AC3 * h_AC3;       % W, HT Rate in desorber
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

%%  5.3 Mathematical Model of R123 Organic Rankine Cycle (ORC_R123) -------------------------------------
% 5.3.1 Status 1, 1-2 Pump
Af = 4.643358E2; Bf =  1.625985E3; Cf = -1.333543E3;
Df = 1.986142E3; Ef = -7.172430E2;
T_rORC1 = T_ORC1 ./ T_cR123;                                          % Reduced Temerature
D_ORC1 = Af + Bf * (1-T_rORC1).^(1/3) + Cf * (1-T_rORC1).^(2/3) + ...
         Df * (1-T_rORC1) + Ef - (1-T_rORC1).^(4/3);
v_ORC1 = D_ORC1 / 1;
% p_ORC1 = 10 * CoolProp.PropsSI('P', 'T', T_ORC1, 'Q', 0, 'R123');
p_ORC1 = CoolProp.PropsSI('P', 'T', T_ORC1, 'Q', 0, 'R123');
s_ORC1 = CoolProp.PropsSI('S', 'T', T_ORC1, 'Q', 0, 'R123');
h_ORC1 = CoolProp.PropsSI('H', 'T', T_ORC1, 'Q', 0, 'R123');
%% 5.3.2 Status 2, 2-3 Evaporator
h_ORC2 = h_ORC1 + v_ORC1 * (p_ORC2 - p_ORC1);
T_ORC2 = T_ORC1;                                    % ???
% T_ORC2 = CoolProp.PropsSI('T', 'H', h_ORC2, 'P', p_ORC2, 'R123');
s_ORC2 = s_ORC1 / ETA_ps;                                             % ???
%% 5.3.3 Status 3, 3-4 Gas Turbine
syms v_3sym
p_ORC3 = p_ORC2;
T_r3 = T_ORC3 ./ T_cR123;                             % Reduced Temerature
ALPHASqrt3 = 1 + KAPPA_R123 * (1 - sqrt(T_r3));
ALPHA3 = ALPHASqrt3^2;                                % Temp-dependent parameter in PR-EOS
a_T3 = a_TcR123 * ALPHA3;                             % Temp-dependent parameter in PR-EOS, on equation
Sv_3sym = solve(p_ORC3 == R*T_ORC3 / (v_3sym-b_R123) - a_T3 / (v_3sym*(v_3sym+b_R123) + b_R123*(v_3sym-b_R123)));
v_ORC3 = double(Sv_3sym);
v_ORC3 = v_ORC3(imag(v_ORC3)==0);
s_ORC3 = CoolProp.PropsSI('S', 'T', T_ORC3, 'P', p_ORC3, 'R123');
h_ORC3 = CoolProp.PropsSI('H', 'T', T_ORC3, 'P', p_ORC3, 'R123');
%% 5.3.4 Status 4, 4-1 Condenser
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
T_r4 = T_ORC4 ./ T_cR123;                             % Reduced Temerature
ALPHASqrt4 = 1 + KAPPA_R123 * (1 - sqrt(T_r4));
ALPHA4 = ALPHASqrt4^2;                                % Temp-dependent parameter in PR-EOS
a_T4 = a_TcR123 * ALPHA4;                             % Temp-dependent parameter in PR-EOS, on equation
Sv_4sym = solve(p_ORC4 == R*T_ORC4 / (v_4sym-b_R123) - a_T4 / (v_4sym*(v_4sym+b_R123) + b_R123*(v_4sym-b_R123)));
v_ORC4 = double(Sv_4sym);
v_ORC4 = v_ORC4(imag(v_ORC4)==0);
%}
s_ORC4 = s_ORC3;
h_ORC4 = CoolProp.PropsSI('H', 'P', p_ORC4, 'Q', 1, 'R123');
%% Area of Heat Exchange in Evaporator within ORC_R123.
T_ORC19 = T_AC12; p_ORC19 = p_AC12; p_ORC20 = p_ORC19;
Q_ORC1 = m_ORC * (h_ORC3 - h_ORC2);
syms T_ORC20 A_ORC1 DELTA_T_ORC1
eq_ORC1(1) = DELTA_T_ORC1 == ((T_ORC19 - T_ORC3) - (T_ORC20 - T_ORC2)) / ...
                          log((T_ORC19 - T_ORC3) / (T_ORC20 - T_ORC2));
eq_ORC1(2) = Q_ORC1 == DELTA_T_ORC1 * K * A_ORC1;
eq_ORC1(3) = Q_ORC1 == m_g * c_pg * (T_ORC19 - T_ORC20);
[ST_ORC20,SA_ORC1,SDELTA_T_ORC1] = solve(eq_ORC1);
     T_ORC20 = double(ST_ORC20);
      A_ORC1 = double(SA_ORC1);
DELTA_T_ORC1 = double(SDELTA_T_ORC1);

%% 6. Economic Model ------------------------------------------------------------------------------------
%% 6.1 Economic Model of MGT.
ETA_t = W_pMGT ./ (m_f * LHV);                             % Thermo Efficiency of MGT
P = A_MGT + B_MGT * log(W_pMGT/1000) + C_MGT * exp(ETA_t); % 单位功率价格, P = 0.9;
z_MT = P * W_pMGT / 1000;                                  % 燃气轮机的购置费
Q_HRSG = (h_9 - h_8) * m_MGTs;                                % 余热锅炉热负荷
% z_HRSG = z_w * (Q_HRSG ./ Q_w)^ALPHA;                    % 余热锅炉购置费 ???
z_HRSG = (Q_HRSG * 4.851)^ALPHA;                           % 余热锅炉购置费
C_E_MGT = W_pMGT * ETA_G * Z_E * N * N_y * 1000;           % Profit from Electricity generated by MGT.
%% 6.2 Economic Model of AC_ALB.
% Area of heat exchanger in solution heat exchanger A_ACs in AC_ALB.
Q_ACs = m_AC2 * (h_AC3 - h_AC2);
DeltaT_ACs = ((T_AC4 - T_AC3) - (T_AC5 - T_AC2)) / log((T_AC4 - T_AC3) / (T_AC5 - T_AC2));
A_ACs = Q_ACs / (DeltaT_ACs * K);
Ca_ACs = Z_A * A_ACs;
% Area of heat exchanger in condenser A_ACc in AC_ALB.
Q_c = m_AC7 * (h_AC7 - h_AC8);                             % W, HT Rate of heat exchanger in desorber
T_AC15 = T_0;
T_AC16 = T_0H;
DeltaT_ACc = ((T_AC7 - T_AC16) - (T_AC8 - T_AC15)) / log((T_AC7 - T_AC16) / (T_AC8 - T_AC15));
A_c = Q_c / DeltaT_ACc / K;
h_AC15 = CoolProp.PropsSI('H', 'T', T_AC15, 'P', p_0, 'Water');
h_AC16 = CoolProp.PropsSI('H', 'T', T_AC16, 'P', p_0, 'Water');
m_c = Q_c / (h_AC16 - h_AC15);
C_Ac = Z_A * A_c;
C_Wc = Z_W * m_c * N * N_y * 3600;
% Area of heat exchanger in evaporator A_ACe in AC_ALB.
Q_ACe = m_AC9 * (h_AC10 - h_AC9);                        % W, Power of supplying cooling energy.
T_AC17 = T_0; p_AC17 = p_0; p_AC18 = p_AC17;
h_AC17 = h_w0; s_AC17 = s_w0; E_ACph17 = E_w0;
T_AC18 = 5 + 273.15;
syms m_ACcw A_ACe DeltaT_ACe
eq_ACe(1) = DeltaT_ACe == ((T_AC17 - T_AC10) - (T_AC18 - T_AC9)) ...
                        / log((T_AC17 - T_AC10) / (T_AC18 - T_AC9));
eq_ACe(2) = Q_ACe == DELTA_T_e * K * A_ACe;
eq_ACe(3) = Q_ACe == c_pw * m_ACcw * (T_AC17 - T_AC18);
[Sm_ACcw, SA_e, SDeltaT_ACe] = solve(eq_ACe);
  m_ACcw = double(Sm_ACcw);
   A_ACe = double(SA_ACe);
DeltaT_e = double(SDeltaT_ACe);
C_ACcw = Z_cw * m_ACcw * N * N_y * 3600;
% Area of heat exchanger in absorber A_ACa in AC_ALB.
Q_ACa = m_AC10 * h_AC10 + m_AC6 * h_AC6 - m_AC1 * h_AC1; % W, HT rate of heat exchanger in desorber
T_AC13 = T_0;
T_AC14 = T_0H;
DeltaT_ACa = ((T_AC6 - T_AC14) - (T_AC1 - T_AC13)) / log((T_AC6 - T_AC14) / (T_AC1 - T_AC13));
A_ACa = Q_ACa / DeltaT_ACa / K;
h_AC13 = CoolProp.PropsSI('H', 'T', T_AC13, 'P', p_0, 'Water');
h_AC14 = CoolProp.PropsSI('H', 'T', T_AC14, 'P', p_0, 'Water');
m_ACa = Q_ACa / (h_AC14 - h_AC13);
C_ACa = Z_A * A_ACa;
Cw_ACa = Z_W * m_ACa * N * N_y * 3600;
W_ACa = m_AC1 * (h_AC2 - h_AC1);                        % W,   Electricity required for pump.
Ce_ACa = Z_E * W_ACa * N * N_y * 1000;                    % RMB, Cost of electricity required for pump.
z_AC = Q_ACe * (Q_wAC / Z_wAC)^ALPHA_wAC;                 % RMB, Estimated purchased price of AC_ALB.
%% 6.3 Economic Model of ORC_R123.
% Area of Heat Exchange in Condenser within ORC_R123.
T_ORC21 = T_0;
T_ORC22 = T_0H;
Q_ORC2 = m_ORC * (h_ORC4 - h_ORC1);
DELTA_T_ORC2 = ((T_ORC4 - T_ORC22) - (T_ORC1 - T_ORC21)) / log((T_ORC4 - T_ORC22) / (T_ORC1 - T_ORC21));
A_ORC2 = Q_ORC2 / DELTA_T_ORC2 / K;
h_ORC21 = CoolProp.PropsSI('H', 'T', T_ORC21, 'P', p_0, 'Water');
h_ORC22 = CoolProp.PropsSI('H', 'T', T_ORC22, 'P', p_0, 'Water');
m_ORC2 = Q_ORC2 / (h_ORC22 - h_ORC21);                   % kg/s, Amount of supplying cooling water.
W_ORC = m_ORC * (h_ORC3 - h_ORC4);                       % W,    Output work of ORC_R123.
C_E_ORC = W_ORC * ETA_G * Z_E * N * N_y * 1000;          % RMB,  Profit of output work of ORC_R123.
W_ORCp = m_ORC * (h_ORC2 - h_ORC1);                      % W,    Required Work of Pump in ORC_R123.
C_EpORC = Z_E * W_ORCp * N * N_y * 1000;                 % RMB,  Cost of required Work of Pump in ORC_R123.
C_ORC1 = A_ORC1 * Z_A;                                   % RMB,  Cost of area for heat exchange in evaporator.
C_ORC2 = A_ORC2 * Z_A;                                   % RMB,  Cost of area for heat exchange in evaporator.
C_ORCw = m_ORC2 * Z_W;                                   % RMB,  Cost of supplying cooling water.
z_ORC = Z_wORC * (W_ORC / Q_wORC)^ALPHA_wORC;            % RMB,  Non-energy cost of ORC_R123. ??????

%% 5. Define Objective Function -------------------------------------------------------------------------
f = 1000000 - c_f * m_f * LHV - CRF * phi * (z_MT + z_HRSG) / (3600 * N) + C_E_MGT + ... % MGT Objective
    500000 - CRF * phi * z_AC / (3600 * N) + ...                                    % AC_ALB Objective 1
    500000 - (Ca_ACs + Ca_ACd + Cw_ACc + Ca_ACc + Cw_ACa + Ca_ACa) - Ce_ACa + C_Ce + ...          % AC_ALB Objective 2
    500000 - CRF * phi * z_ORC / (3600 * N) + ...                                 % ORC_R123 Objective 1
    500000 - (C_ORC1 + C_ORC2 + C_ORCw + C_EpORC) + C_E_ORC;                      % ORC_R123 Objective 2

% end
