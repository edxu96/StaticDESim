% Title: Calculate the other three variables based on two knowing
%        thermodynamic variables. And plot the T-s Diagram of R123.
% Based on: MATLAB program from << Chemical, Biochemical, and Engineering
%                                  Thermodynamics >>
%          《用PR状态方程计算制冷工质的热力学性质及循环性能》
% Method 1, Version: 1.0, Edward Xu, 18.5.11

% Input from x[] ----------------------------------------------------------
%{
T = x(1);     %
p = x(2);     %
%}

function [H,S] = ThermoProp_R123_EdXu_1(T,p)
%% Constants --------------------------------------------------------------
R = 8.314472;                % J/(mol*K), Universial Gas Constant
M = 18 / 1000;               % kg / mol , Molar Mass
R_G = 0.0544 * 1000;         % J/(K*kg) , Gas Constant - R123
AAA = 39.49E-3;              % C_p1, heat capacity calculation parameter of R123
BBB = -2.743E-5;             % C_p2, heat capacity calculation parameter of R123
CCC = -0.122E-8;             % C_p3, heat capacity calculation parameter of R123
DDD = 0.572E-11;             % C_p4, heat capacity calculation parameter of R123
T_c = 456.83;                % K  , temperature in Critical Point.
p_c = 3668.0 * 1000;         % Pa , pressure in Critical Point.
OMEGA = 0.281922497036;      % Acentric Factor.
Tboil = 273.15 + 27.85;      % K  , Boiling point at one atmosphere.
T_ref = 273.15;              % K  , temperature in reference state
p_ref = 101.325 * 1000;      % Pa , pressure in reference state

%% Part1: Peng-Robinson Constant Calculation. -----------------------------
% Reference: << Chemical, Biochemical, and Engineering Thermodynamics, 5e >>
%            E6.4-2, P221, Ch6
a_Tc = 0.457235529 * (R_G * T_c)^2 ./ p_c;              % Critical Point Restriction "a(T_c)"
KAPPA = 0.37464 + (1.54226 - 0.26992 * OMEGA) * OMEGA;  % Dependent on OMEGA(working substance), Temperature-independent parameter in PR-EOS
T_r = T ./ T_c;                                         % Reduced Temerature
ALPHASqrt = 1 + KAPPA * (1 - sqrt(T_r));
ALPHA = ALPHASqrt^2;                                    % Temperature-dependent parameter in PR-EOS
a_T = a_Tc * ALPHA;                                     % Temperature-dependent parameter in PR-EOS (dimensions depend on equation)
b = 0.077796074 * R_G * T_c ./ p_c;                     % m^3/mol, Critical Point Restriction "b", Temperature-independent parameter in PR-EOS
DADT = - a_Tc * KAPPA * ALPHASqrt ./ sqrt(T_c * T);     % partial(a) / partial(T)
A = a_T * p ./ ((R_G * T)^2);                           % Parameters for Cubic Form of Equation of State, A
B = b * p ./ (R_G * T);                                 % Parameters for Cubic Form of Equation of State, B
% P-R EOS: Z^3 + (-1 + B) * Z^2 + (A - 3*B^2 - 2*B) * Z + (-A*B + B^2 + B^3) = 0
% And it can be transformed to : Z^3 + r * Z + q = 0
syms Z
SZ = solve(Z^3 + (-1 + B) * Z^2 + (A - 3*B^2 - 2*B) * Z + (-A*B + B^2 + B^3) == 0);
SZ

%% Part2: Solve Peng-Robinson EOS to get compressibility factor. ----------
% Reference: << Chemical, Biochemical, and Engineering Thermodynamics, 5e >>
%            E6.4-4, P222-223, Ch6
% Root = SolveCubic(Para_TcF(1),Para_TcF(2),Para_TcF(3));
Z(1) = Z_g;
Z(2) = Z_l;

%% Part3: Solve for Peng-Robinson compressibility factor. -----------------
TT = T;
pp = p;
Fugacity = SolveFugacity(A,B,Z,pp);
DepartHS = SolveDepartHS(a_T,b,B,Z,DADT,TT,pp,R,R_G);
DH(1) = DepartHS(1);
DH(2) = DepartHS(2);
DS(1) = DepartHS(3);
DS(2) = DepartHS(4);

%% Part4: Enthalpy and Entropy Change
H_IG = AAA * (T-T_ref) + BBB * (T^2-T_ref^2)./2 + ...
       CCC * (T^3-T_ref^3)./3 + DDD * (T^4-T_ref^4)./4;     % Ideal Gas
S_IG = AAA * log(T./T_ref) + BBB * (T-T_ref) + ...
       CCC * (T^2-T_ref^2)./2 + DDD * (T^3-T_ref^3)./3;     % Ideal Gas
S_IG = S_IG - R * log(p/p_ref);                                % Ideal Gas
H = DH + H_IG;
S = DS + S_IG;

%% Define SubFunction area ------------------------------------------------
% SubFunction1 SolveFugacity:
function Fugacity = SolveFugacity(A,B,Z,pp)
ParaFuga2 = A/B/sqrt(8);
for i=1:2
    ParaFuga1(i) = log((Z(i)+(1+sqrt(2))*B)/(Z(i)+(1-sqrt(2))*B));
    LogFugacityCoeff(i) = (Z(i)-1)-log(Z(i)-B)-ParaFuga1(i)*ParaFuga2;
    FugacityCoeff(i) = exp(LogFugacityCoeff(i));
    Fugacity(i) = FugacityCoeff(i) * pp;
end
clear i;
end

% SubFunction2 SolveDepartHS:
function DepartHS = SolveDepartHS(a_T,b,B,Z,DADT,TT,pp,R,R_G)
for i=1:2
    ParaFuga1(i) = log((Z(i) + (1+sqrt(2))*B) / (Z(i) + (1-sqrt(2))*B));
    DH(i) = (TT * DADT - a_T) * ParaFuga1(i) / b / sqrt(8);
    DH(i) = R * TT * (Z(i)-1) + DH(i) * R / R_G;
    % Gas Enthalpy Departure , EQN 6.4-29
    DS(i) = DADT * ParaFuga1(i) / b / sqrt(8);
    DS(i) = R * log((Z(i)-B)) + DS(i) * R / R_G;
    % Gas Entropy Departure  , EQN 6.4-30
end
clear i;
DepartHS = [DH DS];
end

end