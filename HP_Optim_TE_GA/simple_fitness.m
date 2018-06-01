% Title: R410a Heat Pump (HP_R410a) Modeling.
% Version: 1.0, Edward Xu, 2018.6.1.
% Subtitle: Define the Fitness Function.
% 2-3  定压放热 Condenser      冷凝器
% 3-4  绝热膨胀 Throttle Valve 节流阀
% 4-1  定压吸热 Evaporator     蒸发器
% 1-2  定熵加压 Compressor(p)  压缩机

% function f = simple_fitness(x)
%% 1. Constant
R = 8.314472;                                 % J/(mol*K), Universial Gas Constant
M = 72.58 / 1000;                             % kg / mol , Molar Mass
R_gR410R410 = 0.11455 * 1000;                 % J/(K*kg) , Gas Constant - R134a
T_cR410a = 345.28;                            % K , temperature in Critical Point.
p_cR410a = 4926.1 * 1000;                     % Pa, pressure in Critical Point.
d_cR410a = 488.90;                            % kg/m^3, critical density.
T_bR410 = -26.06 + 273.15;                    % K , Boiling point at one atmosphere.
T_ref = 25 + 273.15;                          % K , temperature in reference state
p_ref = 101.325 * 1000;                       % Pa, pressure in reference state
ETA_ps = 0.7;                                 % Isentropic efficiency of pump.
ETA_ts = 0.7;                                 % Isentropic efficiency of turbine.
ETA_v = 0.8;                                  % Efficiency of throttle valve.
ETA_p = 0.8;                                  % Efficiency of compressor(p).
ETA_e = 0.9;                                  % Efficiency of evaporator.
ETA_c = 0.9;                                  % Efficiency of condenser.
DELTA_p_C = 100 * 1000;                       % Pa, Pressure drop in condenser.
DELTA_p_E = 100 * 1000;                       % Pa, Pressure drop in evaporator.
OMEGA = 0.296;                                % Acentric Factor.
KAPPA = 0.37464 + ...                         % Dependent on OMEGA(working substance),
        (1.54226 - 0.26992 * OMEGA) * OMEGA;  % Temperature-independent parameter in PR-EOS
a_Tc = 0.457235529 * (R_gR410 * T_c)^2 ./ p_c; % Critical Point Restriction "a(T_c)"
b = 0.077796074 * R_gR410 * T_c ./ p_c;        % m^3/mol, Critical Point Restriction "b",
                                              % Temperature-independent parameter in PR-EOS
%% 2. Input
x(1) = 30 + 273.15;
x(2) = 70 + 273.15;
x(3) = 10 + 273.15;
x(4) = 0.005;
 T_1 = x(1); % Outlet temperature of evaporator.
 T_3 = x(2); % Outlet temperature of condenser.
 T_2 = x(3); % Outlet temperature of compressor.
q_HP = x(4); % Fluid Rate.
%% 3. (1)
p_1 = CoolProp.PropsSI('P', 'T', T_1, 'Q', 1, 'R410a');
% Solve for v_1
syms v_1sym
T_r1 = T_1 ./ T_c;                             % Reduced Temerature
ALPHASqrt1 = 1 + KAPPA * (1 - sqrt(T_r1));
ALPHA4 = ALPHASqrt1^2;                         % Temp-dependent para in PR-EOS
a_T1 = a_Tc * ALPHA4;                          % Temp-dependent para in PR-EOS
Sv_1sym = solve(p_1 == R*T_1 / (v_1sym-b) - ...
                       a_T1 / (v_1sym*(v_1sym+b) + b*(v_1sym-b)));
v_1 = double(Sv_1sym);
v_1 = v_1(imag(v_1)==0);
%
s_1 = CoolProp.PropsSI('S', 'T', T_1, 'Q', 1, 'R134a');
h_1 = CoolProp.PropsSI('H', 'T', T_1, 'Q', 1, 'R134a');
%% 4. (4)
p_4 = p_1; T_4 = T_1;
% Solve for v_4.
syms v_4sym
T_r4 = T_4 ./ T_c;                             % Reduced Temerature
ALPHASqrt4 = 1 + KAPPA * (1 - sqrt(T_r4));
ALPHA4 = ALPHASqrt4^2;                         % Temp-dependent para in PR-EOS
a_T4 = a_Tc * ALPHA4;                          % Temp-dependent para in PR-EOS
Sv_4sym = solve(p_4 == R*T_4 / (v_4sym-b) - ...
                a_T4 / (v_4sym*(v_4sym+b) + b*(v_4sym-b)));
v_4 = double(Sv_4sym);
v_4 = v_4(imag(v_4)==0);
%
s_4 = CoolProp.PropsSI('S', 'T', T_4, 'P', p_4, 'R410a');
h_4 = CoolProp.PropsSI('H', 'T', T_4, 'P', p_4, 'R410a');
%% (3)
% <8> solve v_3 through Equation for Density of the Saturated Liquid.
Af = 1.000000; Bf =  1.984734;    Cf = -1.767593E-01;
Df = 1.819972; Ef = -7.171684E-1;
T_r3 = T_3 ./ T_c;                             % Reduced Temerature
d_3 = d_cR410a * (Af + Bf * (1-T_r3).^(1/3) + Cf * (1-T_r3).^(2/3) + ...
                  Df * (1-T_r3) + Ef - (1-T_r3).^(4/3));
v_3 = d_3 / 1;
% <9> Solve for p_3.
syms p_3sym
T_r3 = T_3 ./ T_c;                             % Reduced Temerature
ALPHASqrt3 = 1 + KAPPA * (1 - sqrt(T_r3));
ALPHA3 = ALPHASqrt3^2;                         % Temp-dependent para in PR-EOS
a_T3 = a_Tc * ALPHA3;                          % Temp-dependent para in PR-EOS
Sp_3sym = solve(p_3sym == R*T_3 / (v_3-b) - ...
                a_T3 / (v_3*(v_3+b) + b*(v_3-b)));
p_3 = double(Sp_3sym);
p_3 = p_3(imag(p_3)==0);
%
s_3 = CoolProp.PropsSI('S', 'T', T_3, 'Q', 0, 'R410a');
h_3 = CoolProp.PropsSI('H', 'T', T_3, 'Q', 0, 'R410a');
%% (2)
p_2 = p_3; s_2 = s_1;
% <11> Solve for v_2.
syms v_2sym
T_r2 = T_2 ./ T_c;                             % Reduced Temerature
ALPHASqrt2 = 1 + KAPPA * (1 - sqrt(T_r2));
ALPHA2 = ALPHASqrt2^2;                         % Temp-dependent para in PR-EOS
a_T2 = a_Tc * ALPHA2;                          % Temp-dependent para in PR-EOS
Sv_2sym = solve(p_2 == R*T_2 / (v_2sym-b) - ...
                a_T2 / (v_2sym*(v_2sym+b) + b*(v_2sym-b)));
v_2 = double(Sv_2sym);
v_2 = v_2(imag(v_2)==0);
%
s_2 = CoolProp.PropsSI('S', 'T', T_2, 'P', p_2, 'R410a');
h_2 = CoolProp.PropsSI('H', 'T', T_2, 'P', p_2, 'R410a');
% Input work of compressor W_1.
W_1 = q_HP * (h_2 - h_1);
%% 4. 经济模型 ---------------------------------------------------------------
% Area of heat exchange in condenser Q_2
Q_HPc = q_HP * (h_2 - h_3);  % W, HT Rate in desorber
m_s2 = Q_HPc / (h_9 - h_0);
DELTA_T_HPc == ((T_HP2 - T_9) - (T_HP3 - T_0)) ...
                / log((T_HP2 - T_9) / (T_HP3 - T_0));
A_HPc = Q_HPc / DELTA_T_HPc / K;
% Area of heat exchange in condenser Q_2
Q_HPe = q_HP * (h_1 - h_4);
DELTA_T_HPe == ((T_0 - T_1) - (T_0L - T_4)) ...
                / log((T_0 - T_1) / (T_0L - T_4));
A_HPe = Q_HPe / DELTA_T_HPe / K;
%% Economic Analysis of HP_R410a.
C_1 = A_1 * 1000; % Cost of area for heat exchange in evaporator.
C_2 = A_2 * 1000; % Cost of area for heat exchange in evaporator.
%% 5. 目标函数定义 -----------------------------------------------------------
f = C_1 + C_2;
%% Plot the T-s Diagram.
T = [T_1 T_2 T_3 T_4 T_1]; s = [s_1 s_2 s_3 s_4 s_1];
plot(real(s),T);
axis([950 1100 273.15 350]);
ylabel('温度 T (K)'); xlabel('熵 s (J/kg/K)');
title('蒸汽压缩制冷循环温熵图 T-s Diagram of RC','FontSize',10,'FontWeight','normal');
text(s_1,T_1,'1');
text(s_2,T_2,'2');
text(s_3,T_3,'3');
text(s_4,T_4,'4');

% end
