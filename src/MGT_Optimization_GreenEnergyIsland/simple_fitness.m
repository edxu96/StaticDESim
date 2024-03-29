% Title: Thermo-economic Optimization of Micro Gas Turbine.
% Based on the theory of nolinear equality and inequality constraints.
% Method: Genetic Algorithm within MATLAB Global Optimization Toolbox.
% Version: 4.5, 2018.4.17, Jie Xu.
% SubTitle: Define Objective Function.
% p_2     -> x(1)     空气压缩机出口压力
% eta_AC  -> x(2)     空气压缩机的等熵效率
% eta_GT  -> x(3)     燃气轮机的等熵效率
% T_3     -> x(4)     回热器空气出口温度
% T_4     -> x(5)     燃气轮机进口温度
function f = simple_fitness(x)
%% 1. 通用辅助数据 --------------------------------------------------------------
global LHV DELTA_p_cc DELTA_p_gRE DELTA_p_aRE DELTA_p_HRSG p_1 ...
       T_1 p_8 p_8p p_9 p_7 T_8 T_8p T_9 c_f A B C ALPHA ...
       K_a K_g CRF phi h_8 h_8p h_9 c_pa c_pg ...
       N W_pMGT m_s ETA_CC;
LHV = 500000 * 1000;                              % 纯甲烷燃料的低热值
DELTA_p_cc = 0.05; DELTA_p_aRE = 0.05; ...
DELTA_p_gRE = 0.03; DELTA_p_HRSG = 0.05;          % 压力损失
p_1 = 101.3*1000; T_1 = 288.15;                   % 环境状态
p_8 = 500*1000; p_8p = p_8; p_9 = p_8; p_7 = p_1; % 余热锅炉水蒸气侧的热力学参数
T_8 = 343.15; T_8p = 410; T_9 = 425;              % 余热锅炉水蒸气侧的热力学参数
c_f = 4 * 10^(-9);                                % 单位能量的燃料价格
A = 837.68; B = -143.22; C = 491.79;              % 回归分析统计的价格方程系数
ALPHA = 0.73;                                     % 费用估价系数
K_a = 1.4;                                        % 空气的等熵系数
K_g = 1.33;                                       % 燃气的等熵系数
CRF = 0.182;                                      % 投资回收系数
phi = 1.06;                                       % 维护因子
h_8p = 640039.2; % (503.97 + (589.30-503.97)/20*16.85) * 1000   回水的焓值
h_8 = 293316;    % (251.56 + (335.29-251.56)/20*10) * 1000      余热回收后的焓值
h_9 = 2748.6 * 1000;                              % 过热水蒸汽焓值
c_pa = 1004;                                      % 空气比热容
c_pg = 1170;                                      % 天然气比热容  
ETA_CC = 0.98877;
%% 2. 负荷参数传递
N = 8000;                              % 装置年运行小时数
W_pMGT = 50*1000;                      % Net Power from MGT
m_s = 0.03;                            % Saturated Steam from HRSG
%% 3. Mathematical Model --------------------------------------------------
T_2 = T_1 .* (1 + 1./x(2) * ...
      ((x(1)./p_1)^((K_a-1)./K_a) - 1));                     % (1)  AC
p_3 = x(1) * (1 - DELTA_p_aRE);                              % (7)  RE 
p_4 = p_3 * (1 - DELTA_p_cc);                                % (5)  CC
p_6 = p_7 ./ (1 - DELTA_p_HRSG);                             % (14) HRSG
p_5 = p_6 ./ (1 - DELTA_p_gRE);                              % (8)  RE
T_5 = x(5) * (1 - x(3) * (1 - (p_4./p_5)^((1-K_g)./K_g)));   % (9)  GT
H = (c_pg * x(5) - ETA_CC * LHV) ./ (c_pa * x(4) - ETA_CC * LHV);
m_g = W_pMGT ./ (c_pg*(x(5)-T_5) - c_pa*(T_2-T_1)*H);
m_a = H * m_g;
m_f = m_g - m_a;                                             % (3)  CC 
T_6 = T_5 - m_a * c_pa * (x(4) - T_2) ./ (m_g * c_pg);       % (6)  RE
T_7p = T_6 - m_s * (h_9 - h_8p) ./ (m_g * c_pg);             % (12) HRSG
T_7 = T_6 - m_s * (h_9 - h_8) ./ (m_g * c_pg);               % (13) HRSG
W_AC = m_a * c_pa * (T_2 - T_1);                             % (2)  AC
W_GT = W_pMGT + W_AC;                                        % (11) GT
%% 4. Economic Model ------------------------------------------------------
ETA_t = W_pMGT ./ (m_f * LHV);                               % 燃气轮机热效率
P = A + B * log(W_pMGT/1000) + C * exp(ETA_t);               % 单位功率价格
% P = 0.9;
z_MT = P * W_pMGT / 1000;                                    % 燃气轮机的购置费
Q_HRSG = (h_9 - h_8) * m_s;                                  % 余热锅炉热负荷
% z_HRSG = z_w * (Q_HRSG ./ Q_w)^ALPHA;                      % 余热锅炉购置费
W = 4.851;
z_HRSG = (Q_HRSG * W)^ALPHA;
%% 5. Define Objective Function -------------------------------------------
f = c_f * m_f * LHV + CRF * phi * (z_MT + z_HRSG) ./ (3600 * N);
end

