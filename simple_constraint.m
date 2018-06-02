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
%% Larger than zero, Nonlinear inequality constraints <= 0 ---------------------
c() = -
%% Nonlinear equality constraints -----------------------------------
ceq() =
end
