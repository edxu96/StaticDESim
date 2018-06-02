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
%% Nonlinear inequality constraints <= 0 -----------------------------------
c()

c() = -

c(1) = x(4) - T_5;
c(2) = T_2 - T_6;
Epsilon_MGTre = (x(4) - T_2) ./ (T_5 - T_2);
c(3) = Epsilon_MGTre - 1;                       % 回热度
DELTA_T_p = T_7p - T_9;
c(4) = 15 - DELTA_T_p;                       % 窄点温差
Eta_MGThrsg = (T_6 - T_7) ./ (T_6 - T_8);
c(5) = Eta_MGThrsg - 1;                         % 余热锅炉效率
c(6) = x(4) - x(5);                          % 机组净功率
c(7) = 400 - T_7;                            % To avoid formation of sulfuric acid in exhaust gases
%% Larger than zero, Nonlinear inequality constraints <= 0 ---------------------
c(101) = - W_AC;
c(102) = - W_GT;
c(103) = - m_g;
c(104) = - m_f;

%% Nonlinear equality constraints -----------------------------------
ceq() =
end
