% Title: R134a Vapor Compression Chiller (VCC_R134a) Modeling and Optimization.
% Version: 3.0, Edward Xu, 2018.6.1.
% Subtitle: Define the Fitness Function.
% 2-3  定压放热 Condenser      冷凝器
% 3-4  绝热膨胀 Throttle Valve 节流阀
% 4-1  定压吸热 Evaporator     蒸发器
% 1-2  定熵加压 Compressor(p)  压缩机
% function f = simple_fitness(x)
%% 1. Constant
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
p_0 = 101.325;                                % kP, Pressure of atmosphere.
T_0 = 25 + 273.15;                            % K, Temperature of atmosphere.
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
%% 2. Input Data.
x(1) = 30 + 273.15;
x(2) = 70 + 273.15;
x(3) = 10 + 273.15;
x(4) = 0.005;
T_VCC1 = x(1); % Outlet temperature of evaporator.
T_VCC3 = x(2); % Outlet temperature of condenser.
T_VCC2 = x(3); % Outlet temperature of compressor.
q_2 = x(4); % Fluid Rate.
%% 3. (1)
% <1> Solve T_VCC4 through Equation for Saturated Vapor Pressure.
A =  4.069889E1;  B = -2.362540E3;  C = -1.306883E1;
D =  7.616005E-3; E =  2.342564E-1; F =  3.761111E2;
p_1 = 10^(A + B./T_VCC1 + C * log10(T_VCC1) + D * T_VCC1 + ...
          E * ((F-T_VCC1)./T_VCC1) * log10(F-T_VCC1)) * 1000;
% <2> Solve for v_1
syms v_1sym
T_r1 = T_VCC1 ./ T_c;                             % Reduced Temerature
ALPHASqrt1 = 1 + KAPPA * (1 - sqrt(T_r1));
ALPHA4 = ALPHASqrt1^2;                         % Temperature-dependent parameter in PR-EOS
a_T1 = a_Tc * ALPHA4;                          % Temperature-dependent parameter in PR-EOS (dimensions depend on equation)
Sv_1sym = solve(p_1 == R*T_VCC1 / (v_1sym-b) - a_T1 / (v_1sym*(v_1sym+b) + b*(v_1sym-b)));
v_1 = double(Sv_1sym);
v_1 = v_1(imag(v_1)==0);
%
s_VCC1 = CoolProp.PropsSI('S', 'T', T_VCC1, 'Q', 1, 'R123');
h_VCC1 = CoolProp.PropsSI('H', 'T', T_VCC1, 'Q', 1, 'R123');
%% 4. (4)
p_4 = p_1; T_VCC4 = T_VCC1;
% <5> Solve for v_4.
syms v_4sym
T_r4 = T_VCC4 ./ T_c;                             % Reduced Temerature
ALPHASqrt4 = 1 + KAPPA * (1 - sqrt(T_r4));
ALPHA4 = ALPHASqrt4^2;                         % Temperature-dependent parameter in PR-EOS
a_T4 = a_Tc * ALPHA4;                          % Temperature-dependent parameter in PR-EOS (dimensions depend on equation)
Sv_4sym = solve(p_4 == R*T_VCC4 / (v_4sym-b) - ...
                a_T4 / (v_4sym*(v_4sym+b) + b*(v_4sym-b)));
v_4 = double(Sv_4sym);
v_4 = v_4(imag(v_4)==0);
%
s_VCC4 = CoolProp.PropsSI('S', 'T', T_VCC4, 'P', p_4, 'R123');
h_VCC4 = CoolProp.PropsSI('H', 'T', T_VCC4, 'P', p_4, 'R123');
%% (3)
% <8> solve v_3 through Equation for Density of the Saturated Liquid.
Af =  5.281464E2; Bf =  7.551834E2; Cf = 1.028676E3;
Df = -9.491172E2; Ef = 5.935660E2;
T_r3 = T_VCC3 ./ T_c;                             % Reduced Temerature
RHO_3 = Af + Bf * (1-T_r3).^(1/3) + Cf * (1-T_r3).^(2/3) + ...
        Df * (1-T_r3) + Ef - (1-T_r3).^(4/3);
v_3 = RHO_3 / 1;
% <9> Solve for p_3.
syms p_3sym
T_r3 = T_VCC3 ./ T_c;                             % Reduced Temerature
ALPHASqrt3 = 1 + KAPPA * (1 - sqrt(T_r3));
ALPHA3 = ALPHASqrt3^2;                         % Temperature-dependent parameter in PR-EOS
a_T3 = a_Tc * ALPHA3;                          % Temperature-dependent parameter in PR-EOS (dimensions depend on equation)
Sp_3sym = solve(p_3sym == R*T_VCC3 / (v_3-b) - ...
                a_T3 / (v_3*(v_3+b) + b*(v_3-b)));
p_3 = double(Sp_3sym);
p_3 = p_3(imag(p_3)==0);
%
s_VCC3 = CoolProp.PropsSI('S', 'T', T_VCC3, 'Q', 0, 'R134a');
h_VCC3 = CoolProp.PropsSI('H', 'T', T_VCC3, 'Q', 0, 'R134a');
%% (2)
p_2 = p_3; s_VCC2 = s_VCC1;
% <11> Solve for v_2.
syms v_2sym
T_r2 = T_VCC2 ./ T_c;                             % Reduced Temerature
ALPHASqrt2 = 1 + KAPPA * (1 - sqrt(T_r2));
ALPHA2 = ALPHASqrt2^2;                         % Temperature-dependent parameter in PR-EOS
a_T2 = a_Tc * ALPHA2;                          % Temperature-dependent parameter in PR-EOS (dimensions depend on equation)
Sv_2sym = solve(p_2 == R*T_VCC2 / (v_2sym-b) - ...
                a_T2 / (v_2sym*(v_2sym+b) + b*(v_2sym-b)));
v_2 = double(Sv_2sym);
v_2 = v_2(imag(v_2)==0);
%
s_VCC2 = CoolProp.PropsSI('S', 'T', T_VCC2, 'P', p_2, 'R134a');
h_VCC2 = CoolProp.PropsSI('H', 'T', T_VCC2, 'P', p_2, 'R134a');
% Input work of compressor W_1.
W_1 = q2 * (h_VCC2 - h_VCC1);
%% 4. 经济模型 ---------------------------------------------------------------
% Area of heat exchange in condenser Q_VCCc
Q_VCCc = q_VCC * (h_VCC2 - h_VCC3);                        % W, HT Rate in desorber
m_VCCc = Q_VCCc / (h_0H - h_0);
DeltaT_VCCc == ((T_VCC2 - T_0H) - (T_VCC3 - T_0)) ...
                / log((T_VCC2 - T_0H) / (T_VCC3 - T_0));
A_VCCc = Q_VCCc / DeltaT_VCCc / K;
% Area of heat exchange in evaporator Q_VCCe
Q_VCCe = q_VCC * (h_VCC1 - h_VCC4);
m_VCCcw = Q_VCCc / (h_0 - h_cw);
DeltaT_VCCc == ((T_0 - T_VCC1) - (T_cw - T_VCC4)) ...
                / log((T_0 - T_VCC1) / (T_cw - T_VCC4));
A_VCCe = Q_VCCe / DeltaT_VCCc / K;
C_1 = A_VCCc * Z_A; % Cost of area for heat exchange in evaporator.
C_2 = A_VCCe * Z_A; % Cost of area for heat exchange in evaporator.
%% 5. 目标函数定义 -----------------------------------------------------------
f = C_1 + C_2;
%% Plot the T-s Diagram.
T = [T_VCC1 T_VCC2 T_VCC3 T_VCC4 T_VCC1]; s = [s_VCC1 s_VCC2 s_VCC3 s_VCC4 s_VCC1];
plot(real(s),T);
axis([950 1100 273.15 350]);
ylabel('温度 T (K)'); xlabel('熵 s (J/kg/K)');
title('蒸汽压缩制冷循环温熵图 T-s Diagram of RC','FontSize',10,'FontWeight','normal');
text(s_VCC1,T_VCC1,'1');
text(s_VCC2,T_VCC2,'2');
text(s_VCC3,T_VCC3,'3');
text(s_VCC4,T_VCC4,'4');

% end
