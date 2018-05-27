% plot

x(1) = 633186.25; 
x(2) = 0.837; 
x(3) = 0.917; 
x(4) = 769.395751953125; 
x(5) = 1300; 

p_2     = x(1);
eta_AC  = x(2);
eta_GT  = x(3);
T_3     = x(4);
T_4     = x(5);

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

%% Plot the T-s Diagram.
T = [T_1 T_2 T_4 T_6 T_1]; s = [s_1 s_2 s_4 s_6 s_1];
plot(real(s),T);
axis([0 60 200 1400]);
ylabel('温度 T (K)'); xlabel('熵 s (J/kg/K)');
title('布雷顿闭合循环温熵图 T-s Diagram of Brayton Closed Cycle', ...
      'FontSize',20,'FontWeight','bold');
text(s_1, T_1, '1', 'FontSize', 15);
text(s_2, T_2, '2', 'FontSize', 15);
text(s_4, T_4, '4', 'FontSize', 15);
text(s_6, T_6, '6', 'FontSize',15);