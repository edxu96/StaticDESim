% Title: Aqueous Lithium-Bromide Absorption Refrigeration Modeling and Optimization.
% Version: 3.1, Edward Xu, 2018.5.26.
% Subtitle: Define the Fitness Function.

function f = simple_fitness(x)

%% 1. Constant ------------------------------------------------------------------------------------------

R = 8.314472;                  % J/(mol*K), Universial Gas Constant
T_ref = 0 + 273.15;            % K  , temperature in reference state
p_ref = 101.325 * 1000;        % Pa , pressure in reference state

p_0 = 101.325;                 % kP, Pressure of atmosphere.
T_0 = 15 + 273.15;             % K, Temperature of atmosphere.
%{
ETA_ps = 0.7;                  % Isentropic efficiency of pump.
ETA_ts = 0.7;                  % Isentropic efficiency of turbine.
ETA_v = 0.8;                   % Efficiency of throttle valve.
ETA_p = 0.8;                   % Efficiency of compressor(p).
ETA_e = 0.9;                   % Efficiency of evaporator.
ETA_c = 0.9;                   % Efficiency of condenser.
DELTA_p_C = 100 * 1000;        % Pa, Pressure drop in condenser.
DELTA_p_E = 100 * 1000;        % Pa, Pressure drop in evaporator.
%}

K = ???;                       % W / m2 / K, Thermal conductivity.
Z_A = ???;                     % RMB / m2, Cost rate of area of heat transfer.
Z_W = ???;                     % RMB / kg, Cost rate of cooling water.
Z_C = ???;                     % RMB / J / s, Profit rate of supplying cooling load.
Z_E = ???;                     % RMB / J / s, Cost of supplying electricity for pump.

N = 10;                        % Years, operating time of system.

%% 2. Data from ORC. ------------------------------------------------------------------------------------

x(1) = 4 + 273.15;
x(2) = 30 + 273.15;
x(3) = 0.5322;
x(4) = 0.6711;
x(5) = 1;           %
x(6) = 0.8;         %

   T_10 = x(1);     % K,    Outlet temperature of Evaporator.
    T_8 = x(2);     % K,    Outlet temperature of Condenser.
    y_1 = x(3);     %       Mass Fraction of Outlet Solution from Absorber.
    y_4 = x(4);     %       Mass Fraction of Outlet Solution from Desorber.
    m_1 = x(5);     % kg/s, Fluid Rate Outlet Solution from Absorber.
ETA_shx = x(6);

%% 3. (10) ----------------------------------------------------------------------------------------------

% Pa. Calculate p_10 through equation of saturated ammonia vapor pressure.
p_10 = CoolProp.PropsSI('P','T', T_10, 'Q', 1, 'water');

p_L = p_10;

% J. Calculate s_10 & h_10 through departure equation.
s_10 = CoolProp.PropsSI('S', 'T', T_10, 'P', p_10, 'water');
h_10 = CoolProp.PropsSI('H', 'T', T_10, 'P', p_10, 'water');

%% 3. (8) -----------------------------------------------------------------------------------------------

% Pa. Calculate p_8 through equation of saturated liquid ammonia density.
p_8 = CoolProp.PropsSI('P','T', T_8, 'Q', 0, 'water');

p_H = p_8;

% J. Calculate s_8 & h_8 through departure equation.
s_8 = CoolProp.PropsSI('S', 'T', T_8, 'P', p_8, 'water');
h_8 = CoolProp.PropsSI('H', 'T', T_8, 'P', p_8, 'water');

%% 3. (9) -----------------------------------------------------------------------------------------------

s_9 = s_8;
p_9 = p_10;
T_9 = T_10;

% J. Calculate h_9
h_9 = CoolProp.PropsSI('H', 'S', s_9, 'T', T_9, 'water');

%% 3. (1) -----------------------------------------------------------------------------------------------

p_1 = p_L;

% Pa. Calculate T_1 through equation of saturated liquid solution density.
T_1wC = CoolProp.PropsSI('T', 'P', p_1, 'Q', 0, 'water') - 273.15;
a0 = -2.00755; a1 = 0.16976; a2 = -3.13336E-3; a3 = 1.97668E-5;
b0 = 124.937;  b1 = -7.7162; b2 = 0.152286;    b3 = -7.9509E-4;
T_1 = T_1wC * (a0 + a1 * (y_1*100) + a2 * (y_1*100)^2 + a3 * (y_1*100)^3) + ...
      b0 + b1 * (y_1*100) + b2 * (y_1*100)^2 + b3 * (y_1*100)^3 + 273.15;

% J. Calculate s_1 & h_1 through departure equation.
y_1str = num2str(y_1);
s_1 = CoolProp.PropsSI('S', 'T', T_1, 'Q', 0, strcat('INCOMP::LiBr[',y_1str,']'));
h_1 = CoolProp.PropsSI('H', 'T', T_1, 'Q', 0, strcat('INCOMP::LiBr[',y_1str,']'));

RHO_1 = CoolProp.PropsSI('D', 'T', T_1, 'P', p_1, strcat('INCOMP::LiBr[',y_1str,']'));
v_1 = 1 / RHO_1;

%% 3. (4) -----------------------------------------------------------------------------------------------

p_4 = p_H;

m_2 = m_1;
m_3 = m_2;
y_2 = y_1;
y_3 = y_2;
m_4 = m_3 * y_3 / y_4;

% Pa. Calculate T_4 through equation of saturated liquid solution density.
T_4wC = CoolProp.PropsSI('T', 'P', p_4, 'Q', 0, 'water') - 273.15;
T_4 = T_4wC * (a0 + a1 * (y_4*100) + a2 * (y_4*100)^2 + a3 * (y_4*100)^3) + ...
      b0 + b1 * (y_4*100) + b2 * (y_4*100)^2 + b3 * (y_4*100)^3 + 273.15;

% J. Calculate s_4 & h_4 through departure equation.
y_4str = num2str(y_4);
s_4 = CoolProp.PropsSI('S', 'T', T_4, 'Q', 0, strcat('INCOMP::LiBr[',y_4str,']'));
h_4 = CoolProp.PropsSI('H', 'T', T_4, 'Q', 0, strcat('INCOMP::LiBr[',y_4str,']'));

%% 3. (2) -----------------------------------------------------------------------------------------------

s_2 = s_1;
p_2 = p_H;
h_2 = h_1 + v_1 * (p_2 - p_1);
w_p = m_2 * h_2 - m_1 * h_1;

T_2 = T_1;
% T_2 = CoolProp.PropsSI('T', 'P', p_2, 'H', h_2, 'Q', 0, 'INCOMP::LiBr[0.5322]');

%% 3. (5) -----------------------------------------------------------------------------------------------

m_5 = m_4;
y_5 = y_4;
p_5 = p_4;

T_5 = T_4 - (T_4 - T_2) * ETA_shx;

% J. Calculate s_5 & h_5 through departure equation.
y_5str = num2str(y_5);
s_5 = CoolProp.PropsSI('S', 'T', T_5, 'P', p_5, strcat('INCOMP::LiBr[',y_5str,']'));
h_5 = CoolProp.PropsSI('H', 'T', T_5, 'P', p_5, strcat('INCOMP::LiBr[',y_5str,']'));

%% 3. (3) -----------------------------------------------------------------------------------------------

p_3 = p_H;

Q_shx = m_4 * h_4 - m_5 * h_5;
h_3 = (Q_shx + m_2 * h_2) / m_3;

y_3str = num2str(y_3);

T_3 = CoolProp.PropsSI('T', 'H', h_3-5000, 'P', p_3, strcat('INCOMP::LiBr[',y_3str,']')); % ????

% J. Calculate s_3 through departure equation and h.
s_3 = CoolProp.PropsSI('S', 'T', T_3, 'P', p_3, strcat('INCOMP::LiBr[',y_3str,']'));

%% 3. (6) -----------------------------------------------------------------------------------------------

y_6 = y_5;
T_6 = T_5;
h_6 = h_5;
p_6 = p_1;

% J. Calculate s_5 through departure equation and h.
y_6str = num2str(y_6);
s_6 = CoolProp.PropsSI('S', 'T', T_6, 'P', p_6, strcat('INCOMP::LiBr[',y_6str,']'));

m_6 = m_4;

%% 3. (7) -----------------------------------------------------------------------------------------------

m_7 = m_3 - m_4;
m_8 = m_7;
m_9 = m_8;
m_10 = m_9;

p_7 = p_H;

% Vapor leaving desorber is assumed in equilibrium with incoming solution stream concentration (state 3).
% This is a standard assumption that represents the best possible case.
T_3wC = CoolProp.PropsSI('T', 'P', p_3, 'Q', 0, 'water') - 273.15;
T_3sat = T_3wC * (a0 + a1 * (y_3*100) + a2 * (y_3*100)^2 + a3 * (y_3*100)^3) + ...
      b0 + b1 * (y_3*100) + b2 * (y_3*100)^2 + b3 * (y_3*100)^3 + 273.15;
T_7 = T_3sat;

% J. Calculate s_7 & h_7 through departure equation.
h_7 = CoolProp.PropsSI('H', 'P', p_7, 'T', T_7, 'water');
s_7 = CoolProp.PropsSI('S', 'P', p_7, 'T', T_7, 'water');

Q_d = m_4 * h_4 + m_7 * h_7 + m_3 * h_3;

Q_a = m_1 * h_1 + m_10 * h_10 + m_6 * h_6;

%% 4. Economic Model ------------------------------------------------------------------------------------

% Area of heat exchanger in solution heat exchanger A_s.
Q_s = m_2 * (h_3 - h_2);
DELTA_T_s = ((T_4 - T_3) - (T_5 - T_2)) / log((T_4 - T_3) / (T_5 - T_2));
A_s = Q_s / (DELTA_T_s * K);
C_As = Z_A * A_s;

% Area of heat exchanger in desorber A_d.
Q_d = m_4 * h_4 + m_7 * h_7 - m_3 * h_3;        % W, Heat transfer rate of heat exchanger in desorber
m_d = 0.1;                                      % Rate of waste smoke flow.
syms T_12 A_d
eq_d(1) = DELTA_T_d == ((T_11 - T_4) - (T_12 - T_3)) / log((T_11 - T_4) / (T_12 - T_3));
eq_d(2) = Q_d == DELTA_T_d * K * A_d;
[ST_12 SA_d] = solve(eq_d);
T_12 = double(ST_12);
A_d = double(SA_d);
C_Ad = Z_A * A_d;

% Area of heat exchanger in condenser A_c.
Q_c = m_7 * (h_7 - h_8);                        % W, Heat transfer rate of heat exchanger in desorber
T_15 = 25 + 273.15;
T_16 = 30 + 273.15;                             % Outlet cooling temp in condenser, required in 30'C
DELTA_T_c = ((T_7 - T_16) - (T_8 - T_15)) / log((T_7 - T_16) / (T_8 - T_15));
A_c = Q_c / DELTA_T_c / K;
h_15 = CoolProp.PropsSI('H', 'T', T_15, 'P', p_0, 'water');
h_16 = CoolProp.PropsSI('H', 'T', T_16, 'P', p_0, 'water');
m_c = Q_c / (h_16 - h_15);
C_Ac = Z_A * A_c;
C_Wc = Z_W * m_c * 60 * 60 * 24 * 365 * N;

% Area of heat exchanger in evaporator A_e.
Q_e = m_9 * (h_10 - h_9);
C_Ce = Z_c * Q_e * 60 * 60 * 24 * 365 * N;      % Profit of supplying cooling energy.

% Area of heat exchanger in absorber A_a.
Q_a = m_10 * h_10 + m_6 * h_6 - m_1 * h_1;      % W, Heat transfer rate of heat exchanger in desorber
T_13 = 25 + 273.15;
T_14 = 30 + 273.15;                             % Outlet cooling temp in absorber, required in 30'C
DELTA_T_a = ((T_6 - T_14) - (T_1 - T_13)) / log((T_6 - T_14) / (T_1 - T_13));
A_a = Q_a / DELTA_T_a / K;
h_13 = CoolProp.PropsSI('H', 'T', T_13, 'P', p_0, 'water');
h_14 = CoolProp.PropsSI('H', 'T', T_14, 'P', p_0, 'water');
m_a = Q_a / (h_14 - h_13);
C_Aa = Z_A * A_a;
C_Wa = Z_W * m_a * 60 * 60 * 24 * 365 * N;

% Electricity required for pump
W_p = m_1 * (h_2 - h_1);
C_Ep = Z_E * W_p * 60 * 60 * 24 * 365 * N;

%% 5. Define Fitness Function ---------------------------------------------------------------------------
f = 1000000 - (C_As + C_Ad + C_Wc + C_Ac + C_Wa + C_Aa) - C_Ep + C_Ce;
% Rest of the total budget, would like the maximum value

%{
%% Plot the T-s Diagram.
T = [T_10 T_7 T_8 T_9 T_10]; s = [s_10 s_7 s_8 s_9 s_10];
plot(real(s),T);
% axis([950 1100 273.15 350]);
ylabel('温度 T (K)'); xlabel('熵 s (J/kg/K)');
title('溴化锂吸收式制冷循环温熵图 T-s Diagram of Aqueous Lithium Bromide Absorption Chiller', ...
      'FontSize',20,'FontWeight','bold');
text(s_10,T_10,'10','FontSize',15);
text(s_7,T_7,'7','FontSize',15);
text(s_8,T_8,'8','FontSize',15);
text(s_9,T_9,'9','FontSize',15);
%}

end
