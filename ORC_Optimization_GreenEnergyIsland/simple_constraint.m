% Title: Ideal Organic Rankine Cycle Modeling and Optimization.
% Version: 1.0, Edward Xu, 2018.4.25.
% Subtitle: Constraints.
% 2-3  定压吸热 Evaporator 蒸发器
% 3-4  膨胀做功 Turbine    汽轮机
% 4-1  定压放热 Condenser  冷凝器
% 1-2  定熵加压 Pump       泵

function [c,ceq] = simple_constraint(x)

  q = x(1); % Fluid Rate.
T_3 = x(2); % Inlet temperature of turbine.     T_e
p_2 = x(3); % Outlet Pressure of pump.
T_1 = x(4); % Outlet temperature of Condenser.  
T_2 = x(5); % Outlet Temperature of Pump.

%% Nonlinear inequality constraints <= 0 -----------------------------------
c(1) = T_3 - T_Hin;
c(2) = T_2 - T_Hout;

%% Nonlinear equality constraints -----------------------------------
ceq(1) = s_4 == s_3;    % 
ceq(2) = s_2 == s_1;    % 
ceq(3) = W1 - W2 == 50; % Output work of organic rankine cycle ???

end
