% Title: Vapor Compression Refrigeration Modeling and Optimization.
% Version: 1.0, Edward Xu, 2018.4.26.
% Subtitle: Define the Fitness Function.
% 2-3  定压放热 Condenser      冷凝器
% 3-4  绝热膨胀 Throttle Valve 节流阀
% 4-1  定压吸热 Evaporator     蒸发器
% 1-2  定熵加压 Compressor(p)  压缩机

function f = simple_fitness(x)
%% 1. Constant
global p_0 T_0 T_Cin ETA_ps ETA_ts ETA_v ETA_p ETA_e ETA_c DELTA_p_C DELTA_p_E ...
       OMEGA KAPPA a_Tc b

R = 8.314472;                % J/(mol*K), Universial Gas Constant
M = 102.03 / 1000;           % kg / mol , Molar Mass
R_G = 0.0815 * 1000;         % J/(K*kg) , Gas Constant - R134a
AAA =  39.49E-3;             % C_p1, heat capacity calculation parameter of R123, from《低品位热能有机朗肯循环特性与效率优化的研究》
BBB = -2.743E-5;             % C_p2, heat capacity calculation parameter of R123, from《低品位热能有机朗肯循环特性与效率优化的研究》
CCC = -0.122E-8;             % C_p3, heat capacity calculation parameter of R123, from《低品位热能有机朗肯循环特性与效率优化的研究》
DDD =  0.572E-11;            % C_p4, heat capacity calculation parameter of R123, from《低品位热能有机朗肯循环特性与效率优化的研究》
T_c = 374.23;                % K  , temperature in Critical Point.
p_c = 4060.3 * 1000;         % Pa , pressure in Critical Point.
Tboil = -26.06 + 273.15;     % K  , Boiling point at one atmosphere.
T_ref = 0 + 273.15;          % K  , temperature in reference state
p_ref = 101.325 * 1000;      % Pa , pressure in reference state

p_0 = 101.325;                                % kP, Pressure of atmosphere.
T_0 = 15 + 273.15;                            % K, Temperature of atmosphere.
% T_Hin = 275 + 273.15;                       % K, Inlet temperature of heat source.
T_Cin = 15 + 273.15;                          % K, Inlet temperature of cold source.
ETA_ps = 0.7;                                 % Isentropic efficiency of pump.
ETA_ts = 0.7;                                 % Isentropic efficiency of turbine.
ETA_v = 0.8;                                  % Efficiency of throttle valve.
ETA_p = 0.8;                                  % Efficiency of compressor(p).
ETA_e = 0.9;                                  % Efficiency of evaporator.
ETA_c = 0.9;                                  % Efficiency of condenser.
DELTA_p_C = 100 * 1000;                       % Pa, Pressure drop in condenser.
DELTA_p_E = 100 * 1000;                       % Pa, Pressure drop in evaporator.
OMEGA = 0.332;                                % Acentric Factor.
KAPPA = 0.37464 + ...                         % Dependent on OMEGA(working substance), 
        (1.54226 - 0.26992 * OMEGA) * OMEGA;  % Temperature-independent parameter in PR-EOS
a_Tc = 0.457235529 * (R_G * T_c)^2 ./ p_c;    % Critical Point Restriction "a(T_c)"
b = 0.077796074 * R_G * T_c ./ p_c;           % m^3/mol, Critical Point Restriction "b", 
                                              % Temperature-independent parameter in PR-EOS

%% 2. Data from ORC.

% T_Hin2 = T_Hout; % Waste smoke from ORC.

x(1) = 30 + 273.15;
x(2) = 70 + 273.15;
x(3) = 10 + 273.15;
x(4) = 0.005;

T_1 = x(1); % Outlet temperature of evaporator.
T_3 = x(2); % Outlet temperature of condenser.
T_2 = x(3); % Outlet temperature of compressor.
  q2 = x(4); % Fluid Rate.

%% 3. (1) 

% <1> Solve T_4 through Equation for Saturated Vapor Pressure.
A =  4.069889E1;  B = -2.362540E3;  C = -1.306883E1; 
D =  7.616005E-3; E =  2.342564E-1; F =  3.761111E2;
p_1 = 10^(A + B./T_1 + C * log10(T_1) + D * T_1 + ...
          E * ((F-T_1)./T_1) * log10(F-T_1)) * 1000;

% <2> Solve for v_1
syms v_1sym
T_r1 = T_1 ./ T_c;                             % Reduced Temerature
ALPHASqrt1 = 1 + KAPPA * (1 - sqrt(T_r1));
ALPHA4 = ALPHASqrt1^2;                         % Temperature-dependent parameter in PR-EOS
a_T1 = a_Tc * ALPHA4;                          % Temperature-dependent parameter in PR-EOS (dimensions depend on equation)
Sv_1sym = solve(p_1 == R*T_1 / (v_1sym-b) - a_T1 / (v_1sym*(v_1sym+b) + b*(v_1sym-b)));
v_1 = double(Sv_1sym);
v_1 = v_1(imag(v_1)==0);

% <3> <4> Calculate s_1 h_1 through PR EOS.
T = T_1; p = p_1;
[H,S] = ThermoProp_R123_EdXu_3(T,p); s_1 = S(1); h_1 = H(1);
clear T p H S; 

%% 4. (4)

p_4 = p_1;
T_4 = T_1;

% <5> Solve for v_4.
syms v_4sym
T_r4 = T_4 ./ T_c;                             % Reduced Temerature
ALPHASqrt4 = 1 + KAPPA * (1 - sqrt(T_r4));
ALPHA4 = ALPHASqrt4^2;                         % Temperature-dependent parameter in PR-EOS
a_T4 = a_Tc * ALPHA4;                          % Temperature-dependent parameter in PR-EOS (dimensions depend on equation)
Sv_4sym = solve(p_4 == R*T_4 / (v_4sym-b) - ...
                a_T4 / (v_4sym*(v_4sym+b) + b*(v_4sym-b)));
v_4 = double(Sv_4sym);
v_4 = v_4(imag(v_4)==0);

% <6> <7> Calculate s_4 h_4 through PR EOS.
T = T_4; p = p_4;
[H,S] = ThermoProp_R123_EdXu_3(T,p); s_4 = S(1); h_4 = H(1);
clear T p H S; 

%% (3)

% <8> solve v_3 through Equation for Density of the Saturated Liquid.
Af =  5.281464E2; Bf =  7.551834E2; Cf = 1.028676E3;
Df = -9.491172E2; Ef = 5.935660E2; 
T_r3 = T_3 ./ T_c;                             % Reduced Temerature
RHO_3 = Af + Bf * (1-T_r3).^(1/3) + Cf * (1-T_r3).^(2/3) + ...
        Df * (1-T_r3) + Ef - (1-T_r3).^(4/3);
v_3 = RHO_3 / 1;

% <9> Solve for p_3.
syms p_3sym
T_r3 = T_3 ./ T_c;                             % Reduced Temerature
ALPHASqrt3 = 1 + KAPPA * (1 - sqrt(T_r3));
ALPHA3 = ALPHASqrt3^2;                         % Temperature-dependent parameter in PR-EOS
a_T3 = a_Tc * ALPHA3;                          % Temperature-dependent parameter in PR-EOS (dimensions depend on equation)
Sp_3sym = solve(p_3sym == R*T_3 / (v_3-b) - ...
                a_T3 / (v_3*(v_3+b) + b*(v_3-b)));
p_3 = double(Sp_3sym);
p_3 = p_3(imag(p_3)==0);

% <10> Calculate s_3 through PR EOS.
T = T_3; p = p_3;
[H,S] = ThermoProp_R123_EdXu_3(T,p); s_3 = S(1);
% h_3 = H(1); ???
clear T p H S;

%% (2)

p_2 = p_3;
s_2 = s_1;

% <11> Solve for v_2.
syms v_2sym
T_r2 = T_2 ./ T_c;                             % Reduced Temerature
ALPHASqrt2 = 1 + KAPPA * (1 - sqrt(T_r2));
ALPHA2 = ALPHASqrt2^2;                         % Temperature-dependent parameter in PR-EOS
a_T2 = a_Tc * ALPHA2;                          % Temperature-dependent parameter in PR-EOS (dimensions depend on equation)
Sv_2sym = solve(p_2 == R*T_2 / (v_2sym-b) - ...
                a_T2 / (v_2sym*(v_2sym+b) + b*(v_2sym-b)));
v_2 = double(Sv_2sym);
v_2 = v_2(imag(v_2)==0);

% <12> Calculate h_2 through PR EOS.
T = T_2; p = p_2;
[H,S] = ThermoProp_R123_EdXu_3(T,p); h_2 = H(1);
% s_2 = S(1); ???
clear T p H S;

% Input work of compressor W_1.
W_1 = q2 * (h_2 - h_1);

%{
%% 4. 经济模型 ---------------------------------------------------------------

% Area of heat exchange in condenser Q_2
Q_2 = q2 * (h_2 - h_3);
q_C = 0.002;                                 % Assume the rate of cooling flow.
h_Cin = T_Cin * 1000;                        % Calculate the inlet enthalpy of air.
h_Cout = h_Cin + Q_2 / q_C / ETA_c;          % Exchange heat of evaporator.
T_Cout = h_Cout / 1000;                      % Calculate the outlet temperature of air.
DELTA_T_m2 = ((T_4 - T_Cout) - (T_1 - T_Cin)) / ...
             log((T_4 - T_Cout) / (T_1 - T_Cin));
K_2 = 1000;                                  % Define the heat transfer efficiency of evaporator.
A_2 = Q_2 / K_2 / DELTA_T_m2;

C_1 = A_1 * 1000; % Cost of area for heat exchange in evaporator.
C_2 = A_2 * 1000; % Cost of area for heat exchange in evaporator.

%% 5. 目标函数定义 -----------------------------------------------------------
f = C_1 + C_2;
%}

%% Plot the T-s Diagram.
T = [T_1 T_2 T_3 T_4 T_1]; s = [s_1 s_2 s_3 s_4 s_1];
plot(real(s),T);
axis([-40 100 273.15 350]);
ylabel('温度 T (K)'); xlabel('熵 s (J/kg/K)');
title('有机朗肯循环温熵图 T-s Diagram of ORC','FontSize',10,'FontWeight','normal');
text(s_1,T_1,'1');
text(s_2,T_2,'2');
text(s_3,T_3,'3');
text(s_4,T_4,'4');

end