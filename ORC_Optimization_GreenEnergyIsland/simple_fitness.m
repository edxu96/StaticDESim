% Title: Ideal Organic Rankine Cycle Modeling and Optimization.
% Version: 1.0, Edward Xu, 2018.4.25.
% Subtitle: Fitness Function.
% 2-3  定压吸热 Evaporator 蒸发器
% 3-4  膨胀做功 Turbine    汽轮机
% 4-1  定压放热 Condenser  冷凝器
% 1-2  定熵加压 Pump       泵

function f = simple_fitness(x)
%% 1. Constant
global p_0 T_0 T_Cin ETA_ps ETA_ts ETA_p ETA_t ETA_e ETA_c DELTA_p_C DELTA_p_E ...
       OMEGA KAPPA a_Tc b K_1 K_2

R = 8.314472;                % J/(mol*K), Universial Gas Constant
M = 152.93 / 1000;           % kg / mol , Molar Mass
R_G = 0.0544 * 1000;         % J/(K*kg) , Gas Constant - R123
AAA = 39.49E-3;              % C_p1, heat capacity calculation parameter of R123
BBB = -2.743E-5;             % C_p2, heat capacity calculation parameter of R123
CCC = -0.122E-8;             % C_p3, heat capacity calculation parameter of R123
DDD = 0.572E-11;             % C_p4, heat capacity calculation parameter of R123
T_c = 456.83;                % K  , temperature in Critical Point.
p_c = 3668.0 * 1000;         % Pa , pressure in Critical Point.
Tboil = 273.15 + 27.85;      % K  , Boiling point at one atmosphere.
T_ref = 273.15;              % K  , temperature in reference state
p_ref = 101.325 * 1000;      % Pa , pressure in reference state
   
p_0 = 101.325;                                % kP, Pressure of atmosphere.
T_0 = 15 + 273.15;                            % K, Temperature of atmosphere.
% T_Hin = 275 + 273.15;                       % K, Inlet temperature of heat source.
T_Cin = 15 + 273.15;                          % K, Inlet temperature of cold source.
ETA_ps = 0.7;                                 % Isentropic efficiency of pump.
ETA_ts = 0.7;                                 % Isentropic efficiency of turbine.
ETA_p = 0.8;                                  % Efficiency of pump.
ETA_t = 0.8;                                  % Efficiency of turbine.
ETA_e = 0.9;                                  % Efficiency of evaporator.
ETA_c = 0.9;                                  % Efficiency of condenser.
DELTA_p_C = 100 * 1000;                       % Pa, Pressure drop in condenser.
DELTA_p_E = 100 * 1000;                       % Pa, Pressure drop in evaporator.
OMEGA = 0.281922497036;                       % Acentric Factor.
KAPPA = 0.37464 + ...                         % Dependent on OMEGA(working substance), 
        (1.54226 - 0.26992 * OMEGA) * OMEGA;  % Temperature-independent parameter in PR-EOS
a_Tc = 0.457235529 * (R_G * T_c)^2 ./ p_c;    % Critical Point Restriction "a(T_c)"
b = 0.077796074 * R_G * T_c ./ p_c;           % m^3/mol, Critical Point Restriction "b", 
                                              % Temperature-independent parameter in PR-EOS

%% 2. Data from MGT
  q_H = m_g; % 0.0043 + 0.38477 kg/s
T_Hin = T_7; % 418.51 K

  q = x(1); % Fluid Rate.
T_3 = x(2); % inlet temperature of turbine
p_2 = x(3); % outlet pressure of pump / inlet pressure of turbine
T_1 = x(4); % Outlet temperature of Condenser.  
T_2 = x(5); % Outlet Temperature of Pump.

%% 3. 2-3 Evaporator 蒸发器 定压吸热 (2)
%{
DELTA_tm = Q_1 ./ (K * A);
syms T_2
T_2 = solve(DELTA_tm = ((T_Hin-T_3)-(T_Hout-T_2))./ ...
      log((T_Hin-T_3)./(T_Hout-T_2)));
%}

% Solve for v_2
syms v_2sym
T_r2 = T_2 ./ T_c;                             % Reduced Temerature
ALPHASqrt2 = 1 + KAPPA * (1 - sqrt(T_r2));
ALPHA2 = ALPHASqrt2^2;                         % Temperature-dependent parameter in PR-EOS
a_T2 = a_Tc * ALPHA2;                          % Temperature-dependent parameter in PR-EOS (dimensions depend on equation)
Sv_2sym = solve(p_2 == R*T_2 / (v_2sym-b) - a_T2 / (v_2sym*(v_2sym+b) + b*(v_2sym-b)));
v_2 = double(Sv_2sym);

% ??? check the s_2 = s_1
l = 2; T(l) = T_2;     
m = 2; p(m) = p_2;
ThermoProp_R123_EdXu; s_2 = S(1); h_2 = H(1);
clear T p H S;

%% 3. 3-4 Turbine 汽轮机（定熵）膨胀做功 (3)
 
% Solve for v_3
syms v_3sym
p_3 = p_2;
T_r3 = T_3 ./ T_c;                            % Reduced Temerature
ALPHASqrt3 = 1 + KAPPA * (1 - sqrt(T_r3));
ALPHA3 = ALPHASqrt3^2;                         % Temperature-dependent parameter in PR-EOS
a_T3 = a_Tc * ALPHA3;                          % Temperature-dependent parameter in PR-EOS (dimensions depend on equation)
Sv_3sym = solve(p_3 == R*T_3 / (v_3sym-b) - a_T3 / (v_3sym*(v_3sym+b) + b*(v_3sym-b)));
v_3 = double(Sv_3sym);

l = 3; T(l) = T_3;     
m = 3; p(m) = p_3;
ThermoProp_R123_EdXu; s_3 = S(1); h_3 = H(1);
clear T p H S;

% Area of heat exchange in evaporator
Q_1 = q * (h_3 - h_2);
h_Hin = T_Hin * 1000;                   % Calculate the inlet enthalpy of smoke.
h_Hout = h_Hin - Q_1 / q_H / ETA_e;     % Exchange heat of evaporator.
T_Hout = h_Hout / 1000;                 % Calculate the outlet temperature of smoke.
DELTA_T_m1 = ((T_Hin - T_3) - (T_Hout - T_2)) / ...
             log((T_Hin - T_3)/(T_Hout - T_2));
K_1 = 1000;                             % Define the heat transfer efficiency of evaporator.
A_1 = Q_1 / K_1 / DELTA_T_m1;


%% 3. 1-2 Pump 泵（定熵）加压 (1)

% solve v_1 through Equation for Density of the Saturated Liquid. -> 2
Af = 4.643358E2; Bf =  1.625985E3; Cf = -1.333543E3;
Df = 1.986142E3; Ef = -7.172430E2; 
T_r1 = T_1 ./ T_c;                             % Reduced Temerature
RHO_1 = Af + Bf * (1-T_r1).^(1/3) + Cf * (1-T_r1).^(2/3) + ...
        Df * (1-T_r1) + Ef - (1-T_r1).^(4/3);
v_1 = RHO_1 / 1;

% Solve p_1 through PR EOS.
ALPHASqrt1 = 1 + KAPPA * (1 - sqrt(T_r1));
ALPHA1 = ALPHASqrt1^2;                         % Temperature-dependent parameter in PR-EOS
a_T1 = a_Tc * ALPHA1;                          % Temperature-dependent parameter in PR-EOS (dimensions depend on equation)
p_1 = R*T_1 / (v-b) - a_T1 ./ (v*(v+b) + b*(v-b));

l = 1; T(l) = T_1;     
m = 1; p(m) = p_1;
ThermoProp_R123_EdXu; s_1 = S(2); h_1 = H(2);
clear T p H S;

p_4 = p_1;

% Used work of pump.
W_2 = q * (h_2 - h_1);

%% 3. 4-1 Condenser 冷凝器（定压定温）冷凝放热 (4)

% Solve T_4 through Equation for Saturated Vapor Pressure.
syms T_4sym
A =  1.656333E3; B = -2.480583E6; C = 1.792522E1; 
D = -8.868380E2; E =  4.617861E2; F = 1.666667E3;
ST_4sym =solve(log10(p_4./1000) == A + B./T_4sym + C * log10(T_4sym) + D * T_4sym + ...
                                   E * ((F-T_4sym)./T_4sym) * log10(F-T_4sym));
T_4 = double(ST_4sym);
T_4 = T_4(imag(T_4)==0);
                            
% Solve for v_4
syms v_4sym
T_r4 = T_4 ./ T_c;                             % Reduced Temerature
ALPHASqrt4 = 1 + KAPPA * (1 - sqrt(T_r4));
ALPHA4 = ALPHASqrt4^2;                         % Temperature-dependent parameter in PR-EOS
a_T4 = a_Tc * ALPHA4;                          % Temperature-dependent parameter in PR-EOS (dimensions depend on equation)
Sv_4sym = solve(p_4 == R*T_4 / (v_4sym-b) - a_T4 / (v_4sym*(v_4sym+b) + b*(v_4sym-b)));
v_4 = double(Sv_4sym);
v_4 = v_4(imag(v_4)==0);

% ??? check the s_4 = s_3
l = 4; T(l) = T_4;     
m = 4; p(m) = p_4;
ThermoProp_R123_EdXu; s_4 = S(1); h_4 = H(1);
clear T p H S;
               
% Area of heat exchange in condenser
Q_2 = q * (h_4 - h_1);
q_C = 0.002;                                 % Assume the rate of cooling flow.
h_Cin = T_Cin * 1000;                        % Calculate the inlet enthalpy of air.
h_Cout = h_Cin + Q_2 / q_C / ETA_c;          % Exchange heat of evaporator.
T_Cout = h_Cout / 1000;                      % Calculate the outlet temperature of air.
DELTA_T_m2 = ((T_4 - T_Cout) - (T_1 - T_Cin)) / ...
             log((T_4 - T_Cout) / (T_1 - T_Cin));
K_2 = 1000;                                  % Define the heat transfer efficiency of evaporator.
A_2 = Q_2 / K_2 / DELTA_T_m2;

% Output work of turbine.
W_1 = q * (h_3 - h_4);

%% 4. 经济模型 ---------------------------------------------------------------

C_1 = A_1 * 1000; % Cost of area for heat exchange in evaporator.
C_2 = A_2 * 1000; % Cost of area for heat exchange in evaporator.

%% 5. 目标函数定义 -----------------------------------------------------------
f = C_1 + C_2;

end