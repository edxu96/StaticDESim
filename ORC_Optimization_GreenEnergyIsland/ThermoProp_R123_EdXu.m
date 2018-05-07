% Title: Calculate the other three variables based on two knowing
%        thermodynamic variables. And plot the T(l)-s Diagram of R123.
% Based on: MATLAB program from << Chemical, Biochemical, and Engineering
%                                  Thermodynamics >>
%          《用PR状态方程计算制冷工质的热力学性质及循环性能》
% Version: 3.1, Edward Xu, 18.4.26

% Input from x[] ----------------------------------------------------------
%{
T(l) = x(1);     %
p(m) = x(2);     %
%}

%{
clear;
T(l) = 273.15 + (370);        % K
p(m) = 200 * 1000; % Pa
%}

%% Constants --------------------------------------------------------------
R = 8.314472;                % J/(mol*K), Universial Gas Constant
M = 18 / 1000;               % kg / mol , Molar Mass
R_G = 0.0544 * 1000;         % J/(K*kg) , Gas Constant - R123
AAA =  39.49E-3;             % C_p1, heat capacity calculation parameter of R123, from《低品位热能有机朗肯循环特性与效率优化的研究》
BBB = -2.743E-5;             % C_p2, heat capacity calculation parameter of R123, from《低品位热能有机朗肯循环特性与效率优化的研究》
CCC = -0.122E-8;             % C_p3, heat capacity calculation parameter of R123, from《低品位热能有机朗肯循环特性与效率优化的研究》
DDD =  0.572E-11;            % C_p4, heat capacity calculation parameter of R123, from《低品位热能有机朗肯循环特性与效率优化的研究》
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
T_r = T(l) ./ T_c;                                      % Reduced Temerature
ALPHASqrt = 1 + KAPPA * (1 - sqrt(T_r));
ALPHA = ALPHASqrt^2;                                    % Temperature-dependent parameter in PR-EOS
a_T = a_Tc * ALPHA;                                     % Temperature-dependent parameter in PR-EOS (dimensions depend on equation)
b = 0.077796074 * R_G * T_c ./ p_c;                     % m^3/mol, Critical Point Restriction "b", Temperature-independent parameter in PR-EOS
DADT = - a_Tc * KAPPA * ALPHASqrt ./ sqrt(T_c * T(l));  % partial(a) / partial(T)
A = a_T * p(m) ./ ((R_G * T(l))^2);                     % Parameters for Cubic Form of Equation of State, A
B = b * p(m) ./ (R_G * T(l));                           % Parameters for Cubic Form of Equation of State, B
% ZZ = ZZroot(A,B);
[Z_g,Z_l] = ZZroot2(A,B);
Z_g
Z_l

%% Peng-Robinson EOS
%{
syms p T v
OMEGA = 0.281922497036;      % Acentric Factor.
KAPPA = 0.37464 + (1.54226 - 0.26992 * OMEGA) * OMEGA;  % Dependent on OMEGA(working substance), Temperature-independent parameter in PR-EOS
ALPHASqrt = 1 + KAPPA * (1 - sqrt(T_r));
ALPHA = ALPHASqrt^2;                                    % Temperature-dependent parameter in PR-EOS
a_Tc = 0.457235529 * (R_G * T_c)^2 ./ p_c;              % Critical Point Restriction "a(T_c)"
a_T = a_Tc * ALPHA;                                     % Temperature-dependent parameter in PR-EOS (dimensions depend on equation)
b = 0.077796074 * R_G * T_c ./ p_c;                     % m^3/mol, Critical Point Restriction "b", Temperature-independent parameter in PR-EOS
T_r = T(l) ./ T_c;                                      % Reduced Temerature
Sv = solve(p == R*T / (v-b) - a_T / (v*(v+b) + b*(v-b)))
%}

%% Part2: Solve Peng-Robinson EOS to get compressibility factor. ----------
% Reference: << Chemical, Biochemical, and Engineering Thermodynamics, 5e >>
%            E6.4-4, P222-223, Ch6
% Root = SolveCubic(Para_TcF(1),Para_TcF(2),Para_TcF(3));
% Root,
% Z(1) = max(ZZ); % Vapor Phase, most compressible;
% Z(2) = min(ZZ); % Liquid Phase, least compressible;
Z(1) = Z_g;
Z(2) = Z_l;

%% Part3: Solve for Peng-Robinson compressibility factor. -----------------
TT = T(l);
pp = p(m);
Fugacity = SolveFugacity(A,B,Z,pp);
DepartHS = SolveDepartHS(a_T,b,B,Z,DADT,TT,pp,R,R_G);
DH(1) = DepartHS(1);
DH(2) = DepartHS(2);
DS(1) = DepartHS(3);
DS(2) = DepartHS(4);

%% Part4: Enthalpy and Entropy Change
H_IG = AAA * (T(l)-T_ref) + BBB * (T(l)^2-T_ref^2)./2 + ...
       CCC * (T(l)^3-T_ref^3)./3 + DDD * (T(l)^4-T_ref^4)./4;     % Ideal Gas
S_IG = AAA * log(T(l)./T_ref) + BBB * (T(l)-T_ref) + ...
       CCC * (T(l)^2-T_ref^2)./2 + DDD * (T(l)^3-T_ref^3)./3;     % Ideal Gas
S_IG = S_IG - R * log(p(m)/p_ref);                                % Ideal Gas
H = DH + H_IG;
S = DS + S_IG;

%{

%% 2. Equation for Saturated Vapor Pressure. ------------------------------
syms p_sat
A =  1.656333E3; B = -2.480583E6; C = 1.792522E1; 
D = -8.868380E2; E =  4.617861E2; F = 1.666667E3;
eqn(2) = log10(p_sat./1000) == A + B./T(l) + C * log10(T) + D * T(l) + ...
                               E * ((F-T(l))./T(l)) * log10(F-T(l));

%% 3. Equation for Density of the Saturated Liquid. -----------------------
syms RHO_f
Af = 4.643358E2; Bf =  1.625985E3; Cf = -1.333543E3;
Df = 1.986142E3; Ef = -7.172430E2; 
eqn(3) = RHO_f == Af + Bf * (1-T_r).^(1/3) + Cf * (1-T_r).^(2/3) + ...
                  Df * (1-T_r) + Ef - (1-T_r).^(4/3);

%% 4. Equation for Ideal Gas Heat Capacity at constant pressure. ----------
syms c_p0
cp1 =  2.89811E1;
cp3 = -1.95477E-1;
cp2 =  3.04711E-1;
eqn(4) = c_p0 == cp1 + cp2 * T(l) + cp3 * T(l).^2;           % (J * mole.^-1 * K.^-1)
% Real Gas Heat Capacity Equation (at constant pressure).
% c_p  = c_p0 - T * int(diff(exp_v_1, T_1, 2), p, p_0, p_1);      

%% 5. Equation for Ideal Gas Heat Capacity at constant vapor. -------------
syms c_v0
a = -5.397695;
b =  3.275570E-2;
d =  4.340077E-8;
c = -6.358090E-5;
f =  6.622145E4; 
eqn(5) = c_v0 * 1000 == (a + b * T(l) + c * T(l).^2 + ...
                        d * T(l).^3 + f ./ T(l).^2) ;        % J * kg.^-1 * K.^-1

%}

%% Optput the result. -----------------------------------------------------
Fugacity;                                                                                 % bar
Compressibility = Z;
fprintf('                      Temperature in this condition is %4.1f K.\n', T(l));       % K
fprintf('                         Pressure in this condition is %4.1f Pa.\n', p(m));      % kPa
fprintf('      Enthalpy of saturated vapor in this condition is %f J/mol.\n', H(1));      % J/mol
fprintf('     Enthalpy of saturated liquid in this condition is %f J/mol.\n', H(2));      % J/mol
fprintf('       Entropy of saturated vapor in this condition is %f J/(mol*K).\n', S(1));  % J/(mol*K)
fprintf('      Entropy of saturated liquid in this condition is %f J/(mol*K).\n', S(2));  % J/(mol*K)
SpecifyVolume = Z * 1e3 * R_G * T(l) ./ p(m);
fprintf('Specify Volume of saturated vapor in this condition is %f.\n',SpecifyVolume(1)); % m^3/kmol
fprintf('Specify Volume of saturated vapor in this condition is %f.\n',SpecifyVolume(2)); % m^3/kmol

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

%{
% SubFunction3 ZZroot: Solve the Equation of State.
function ZZ = ZZroot(A,B)
V(1) = 1;
V(2) = -1+B;
V(3) = A-B*(3*B+2);
V(4) = B*(B*B+B-A);
ZZ = roots(V);
% Get rid off the imag root
for 1i = 1:3
    if imag(ZZ(1i)) ~= 0
        ZZ(1i) = 0;
    end
end
% ***************************
ZZ = sort(ZZ);
if abs(ZZ(1)) < 1e-8
    ZZ(1) = ZZ(3);
end
if abs(ZZ(3)) < 1e-8
    ZZ(3) = ZZ(1);
end
end
%}

% SubFunction4 ZZroot2: Solve the Equation of State using Determining Equation.
function [Z_g,Z_l] = ZZroot2(A,B)
c1 = B - 1;
c2 = A - 3 * B^2 - 2 * B;
c3 = B^3 + B^2 - A*B;
q = 2/27 * c1^3 - c1*c2/3 +c3;
r = c2 - c1^2 / 3;
D = q^2 / 4 + r^3 / 27;                                %
OMEGA1 = (-1 + sqrt(3) * 1i) / 2;
OMEGA2 = (-1 - sqrt(3) * 1i) / 2;
Z1 = (-q/2 + sqrt((q/2)^2 + (r/3)^3))^(1/3) + ...
     (-q/2 - sqrt((q/2)^2 + (r/3)^3))^(1/3);
Z2 = OMEGA1 * (-q/2 + sqrt((q/2)^2 + (r/3)^3))^(1/3) + ...
     OMEGA2 * (-q/2 - sqrt((q/2)^2 + (r/3)^3))^(1/3);
Z3 = OMEGA2 * (-q/2 + sqrt((q/2)^2 + (r/3)^3))^(1/3) + ...
     OMEGA1 * (-q/2 - sqrt((q/2)^2 + (r/3)^3))^(1/3);
ZZZ = [Z1 Z2 Z3];                                      % 
ZZZ = ZZZ(imag(ZZZ)==0);                               %
if D > 0
    Z_g = ZZZ(1);
    Z_l = 0;
else
    Z_g = max(ZZZ);
    Z_l = min(ZZZ);
end
end
