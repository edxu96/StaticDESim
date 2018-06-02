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
  % 4.1 Decision Variables in MGT.
        p_2 = x(1);               % Pa,   Outlet Pressure of Air Compressor
  Eta_MGTac = x(2);               %       Isentropic Efficiency of Air Compressor
   Eta_MGTt = x(3);               %       Isentropic Efficiency of Gas Turbine
        T_3 = x(4);               % K,    Outlet Temp of Air from Regenerator
        T_4 = x(5);               % K,    Inlet Temp of Gas Turbine
     Ms_MGT = x(6);               % kg/s, Fluid Rate of Saturated Steam from HRSG
     Wp_MGT = x(7);               % W,    Net Power from Micro Gas Turbine
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
%% 1. Nonlinear inequality constraints <= 0, A - B, A <= B. -----------------------------
% 1.1 Nonlinear inequality constraints from MGT.
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
% 1.2 Nonlinear inequality constraints from AC_ALB.
c(8) = y_1 - y_4;
% 1.3 Nonlinear inequality constraints from ORC_R123.
c(9) = T_ORC1 - T_ORC3;
c(10) = p_ORC1 - p_ORC2;
% 1.4 Nonlinear inequality constraints from HP_R410a.
c(11) =  T_HP3 - T_HP2;
c(12) =  T_HP1 - T_HP3;
% 1.5 Nonlinear inequality constraints from VCC_R134a.
c(13) =  T_VCC3 - T_VCC2;
c(14) =  T_VCC1 - T_VCC3;
%% 2. Larger than zero, Nonlinear inequality constraints <= 0 ---------------------
% 2. Larger than zero from MGT.
c(101) = - W_AC;
c(102) = - W_GT;
c(103) = - m_g;
c(104) = - m_f;
%% Nonlinear equality constraints -----------------------------------
end
