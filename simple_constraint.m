% Title: Thermo-economic Optimization of Distributed Energy System in Green Energy Island.
% Based on the theory of nolinear equality and inequality constraints.
% Method: Genetic Algorithm within MATLAB Global Optimization Toolbox.
% Version: 1.0, 2018.6.2, Jie Xu.
% SubTitle: Non-linear Equality and Inequality Constants
% 1. Micro Gas Turbine (MGT)
% 2. Aqueous Lithium-Bromine Single-Effect Absorption Chiller (AC_ALB)
% 3. R123 Organic Recycle Cycle (ORC_R123)
% 4. R410a Heat Pump (HP_R410a)
% 5. R134a Vapor Compression Chiller (VCC_R134a)
function [c,ceq] = simple_constraint(x)
%% 1. Decision Variables. ------------------
% 1.1 Decision Variables in MGT.
      p_2 = x(1);               % Pa,   Outlet Pressure of Air Compressor
Eta_MGTac = x(2);               %       Isentropic Efficiency of Air Compressor
 Eta_MGTt = x(3);               %       Isentropic Efficiency of Gas Turbine
      T_3 = x(4);               % K,    Outlet Temp of Air from Regenerator
      T_4 = x(5);               % K,    Inlet Temp of Gas Turbine
   Ms_MGT = x(6);               % kg/s, Fluid Rate of Saturated Steam from HRSG
   Wp_MGT = x(7);               % W,    Net Power from Micro Gas Turbine
% 1.2 Decision Variables in AC_ALB.
   T_AC10 = x(8);               % K,    Outlet temperature of Evaporator.
    T_AC8 = x(9);               % K,    Outlet temperature of Condenser.
      y_1 = x(10);              %       Mass Fraction of Outlet Solution from Absorber.
      y_4 = x(11);              %       Mass Fraction of Outlet Solution from Desorber.
    M_AC1 = x(12);              % kg/s, Fluid Rate Outlet Solution from Absorber.
  Eta_shx = x(13);              %       Solution Heat Exchanger Ratio.
% 1.3 Decision Variables in ORC_R123.
    M_ORC = x(14);              % kg/s, Fluid Rate.
   T_ORC3 = x(15);              % K,    Inlet temperature of turbine
   p_ORC2 = x(16);              % K,    Outlet pressure of pump / inlet pressure of turbine
   T_ORC1 = x(17);              % K,    Outlet temperature of condenser.
% 1.4 Decision Variables in HP_R410a.
   T_HP1 = x(18);               % K,    Outlet temperature of evaporator.
   T_HP3 = x(19);               % K,    Outlet temperature of condenser.
   T_HP2 = x(20);               % K,    Outlet temperature of compressor.
    M_HP = x(21);               % kg/s, Fluid rate of HP_R410a.
% 1.5 Decision Variables in VCC_R134a.
  T_VCC1 = x(22);               % K,    Outlet temperature of evaporator.
  T_VCC3 = x(23);               % K,    Outlet temperature of condenser.
  T_VCC2 = x(24);               % K,    Outlet temperature of compressor.
   M_VCC = x(25);               % kg/s, Fluid Rate.
%% 2.1 Constants for MGT
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
% 2.2 Constants calculated for MGT.
p_s = 500 * 1000;                                 % Pa, Pressure of supplying steam.
T_s = CoolProp.PropsSI('T', 'P', p_s, 'Q', 1, 'Water');
h_s = CoolProp.PropsSI('H', 'P', p_s, 'Q', 1, 'Water');
% 3. Mathematical Model of Brayton Cycle.
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
%% 4. Nonlinear inequality constraints <= 0, A - B, A <= B. -----------------------------
% 4.1 Nonlinear inequality constraints from MGT.
c(1) = T_3 - T_5;
c(2) = T_2 - T_6;
Epsilon_MGTre = (T_3 - T_2) ./ (T_5 - T_2);
c(3) = Epsilon_MGTre - 1;                       % 回热度
DeltaT_MGTp = T_7p - T_9;
c(4) = 15 - DeltaT_MGTp;                        % 窄点温差
Eta_MGThrsg = (T_6 - T_7) ./ (T_6 - T_8);
c(5) = Eta_MGThrsg - 1;                         % 余热锅炉效率
c(6) = T_3 - T_4;                               % 机组净功率
c(7) = 400 - T_7;                               % To avoid formation of sulfuric acid in exhaust gases
% 4.2 Nonlinear inequality constraints from AC_ALB.
c(8) = y_1 - y_4;
% 4.3 Nonlinear inequality constraints from ORC_R123.
c(9) = T_ORC1 - T_ORC3;
p_ORC1 = CoolProp.PropsSI('P', 'T', T_ORC1, 'Q', 0, 'R123');
c(10) = p_ORC1 - p_ORC2;
% 4.4 Nonlinear inequality constraints from HP_R410a.
c(11) =  T_HP3 - T_HP2;
c(12) =  T_HP1 - T_HP3;
% 4.5 Nonlinear inequality constraints from VCC_R134a.
c(13) =  T_VCC3 - T_VCC2;
c(14) =  T_VCC1 - T_VCC3;
%% 5. Larger than zero, Nonlinear inequality constraints <= 0 ---------------------
% 5. Larger than zero from MGT.
c(15) = - W_AC;
c(16) = - W_GT;
c(17) = - m_g;
c(18) = - m_f;
%% 6. Nonlinear equality constraints -----------------------------------
ceq = [];
end
