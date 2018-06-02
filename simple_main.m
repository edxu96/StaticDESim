% Title: Thermo-economic Optimization of Distributed Energy System in Green Energy Island.
% Based on the theory of nolinear equality and inequality constraints.
% Method: Genetic Algorithm within MATLAB Global Optimization Toolbox.
% Version: 1.0, 2018.6.2, Jie Xu.
% SubTitle: Main Function.
% 1. Micro Gas Turbine (MGT)
% 2. Aqueous Lithium-Bromine Single-Effect Absorption Chiller (AC_ALB)
% 3. R123 Organic Recycle Cycle (ORC_R123)
% 4. R410a Heat Pump (HP_R410a)
% 5. R134a Vapor Compression Chiller (VCC_R134a)
clear;clc
%% Set Optimization Parameters
ObjectiveFunction = @simple_fitness;
nvars = 25;                                           % Number of variables
lb = [250*1000; 0.8; 0.8;  500; 500;  0.0001; 25000;  ...
      10+273.15;  10+273.15;  0.2; 0.2; 0.0001; 0.2; ...
      0.0001; 10+273.15;  101*1000; 10+273.15;  ...
      10+273.15;  10+273.15;  10+273.15;  0.0001];                % Lower bound
ub = [650*1000; 0.9; 0.93; 900; 1300; 0.01;   100000; ...
      300+273.15; 300+273.15; 0.8; 0.8; 0.01;  0.8; ...
      0.01;   300+273.15; 606*1000; 300+273.15; ...
      300+273.15; 300+273.15; 300+273.15; 0.01];               % Upper bound
ConstraintFunction = @simple_constraint;
rng(1,'twister')                                     % for reproducibility
X0(:,1) = round(rand(1,100)*400*1000) + 250000;
X0(:,2) = (round(rand(1,100)*100) + 800) / 1000;
X0(:,3) = (round(rand(1,100)*130) + 800) / 1000;
X0(:,4) = round(rand(1,100)*400) + 500;
X0(:,5) = round(rand(1,100)*800) + 500;
X0(:,6) = (round(rand(1,100)*990) + 10) / 100000;
X0(:,7) = round(rand(1,100)*(100000-2500)) + 2500;
X0(:,8) = round(rand(1,100)*290) + (10+273.15);
X0(:,9) = round(rand(1,100)*290) + (10+273.15);
X0(:,10) = (round(rand(1,100)*600) + 200) / 1000;
X0(:,11) = (round(rand(1,100)*600) + 200) / 1000;
X0(:,12) = (round(rand(1,100)*990) + 10) / 100000;
X0(:,13) = (round(rand(1,100)*600) + 200) / 1000;
X0(:,14) = (round(rand(1,100)*990) + 10) / 100000;
X0(:,15) = round(rand(1,100)*290) + (10+273.15);
X0(:,16) = round(rand(1,100)*505*1000) + 101*1000;
X0(:,17) = round(rand(1,100)*290) + (10+273.15);
X0(:,18) = round(rand(1,100)*290) + (10+273.15);
X0(:,19) = round(rand(1,100)*290) + (10+273.15);
X0(:,20) = round(rand(1,100)*290) + (10+273.15);
X0(:,21) = (round(rand(1,100)*990) + 10) / 100000;
X0(:,22) = round(rand(1,100)*290) + (10+273.15);
X0(:,23) = round(rand(1,100)*290) + (10+273.15);
X0(:,24) = round(rand(1,100)*290) + (10+273.15);
X0(:,25) = (round(rand(1,100)*990) + 10) / 100000;

% options.InitialPopulationMatrix = rand(20,7);
% X0 = [3,0.8,0.9,821,1164,190,2];                   % Start point (row vector)
% 'InitialPopulationRange', [101 0.5 0.5 0 0 100 1; 1616 0.9 0.93 1000 1600 200 10], ...
options = optimoptions(@ga, ...
          'PopulationSize',100,'InitialPopulationMatrix',X0,'PlotFcn', ...
          {@gaplotbestf,@gaplotmaxconstr,@gaplotrange},'Display','iter');
[x,fval,exitFlag,Output] = ga(ObjectiveFunction,nvars, ...
           [],[],[],[],lb,ub,ConstraintFunction,options);
%% Print the Result.
fprintf('The number of generations was : %d\n', Output.generations);
fprintf('The number of function evaluations was : %d\n', Output.funccount);
fprintf('The best function value found was : %g\n', fval);
