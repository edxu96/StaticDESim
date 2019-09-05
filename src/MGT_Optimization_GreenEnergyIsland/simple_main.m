% Title: Thermo-economic Optimization of Micro Gas Turbine.
% Based on the theory of nolinear equality and inequality constraints.
% Method: Genetic Algorithm within MATLAB Global Optimization Toolbox.
% Version: 4.5, 2018.4.17, Jie Xu.
% SubTitle: Main Function.
% p_2     -> x(1)
% ETA_AC  -> x(2)
% ETA_GT  -> x(3)
% T_3     -> x(4)
% T_4     -> x(5)
clear;clc
%% Set Optimization Parameters
ObjectiveFunction = @simple_fitness;
nvars = 5;                                           % Number of variables
lb = [250*1000; 0.8; 0.8;  500; 500];                % Lower bound
ub = [650*1000; 0.9; 0.93; 900; 1300];               % Upper bound
ConstraintFunction = @simple_constraint;
rng(1,'twister')                                     % for reproducibility
X0(:,1) = round(rand(1,100)*400*1000) + 250000;
% X0(:,1) = [120 120 120 120 120 120 120 120 120 120 ...
%            120 120 120 120 120 120 120 120 120 120];
X0(:,2) = round(rand(1,100)*100) + 800;
X0(:,2) = X0(:,2) / 1000;
X0(:,3) = round(rand(1,100)*130) + 800;
X0(:,3) = X0(:,3) / 1000;
X0(:,4) = round(rand(1,100)*400) + 500;
X0(:,5) = round(rand(1,100)*800) + 500;
% options.InitialPopulationMatrix = rand(20,7);
% X0 = [3,0.8,0.9,821,1164,190,2];                   % Start point (row vector)
% 'InitialPopulationRange', [101 0.5 0.5 0 0 100 1; 1616 0.9 0.93 1000 1600 200 10], ...
options = optimoptions(@ga, ...
          'PopulationSize',100,'InitialPopulationMatrix',X0,'PlotFcn', ...
          {@gaplotbestf,@gaplotmaxconstr,@gaplotrange},'Display','iter');
[x,fval,exitFlag,Output] = ga(ObjectiveFunction,nvars, ...
           [],[],[],[],lb,ub,ConstraintFunction,options)
%% Print the Result.
fprintf('The number of generations was : %d\n', Output.generations);
fprintf('The number of function evaluations was : %d\n', Output.funccount);
fprintf('The best function value found was : %g\n', fval);
