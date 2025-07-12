function [] = compareTargetPositionPerformanceCyz()

%% 主测试脚本
clear all; close all; clc;

% 参数设置
fs = 10000;          % 采样频率 (Hz)
N = 1024;            % 采样点数
test_points = 200;   % 测试点数
use_window = true;   % 是否使用窗函数

% 生成测试频率范围 (在两个FFT bin之间扫描)
delta_f = fs / N;    % FFT频率分辨率
f_start = 500;       % 起始频率 (Hz)
f_end = 520;         % 结束频率 (Hz)
f_test = linspace(f_start, f_end, test_points);

% 初始化结果存储
errors_rife = zeros(1, test_points);
errors_jacobsen = zeros(1, test_points);
errors_fft = zeros(1, test_points);
errors_phase_diff = zeros(1, test_points);

% 主测试循环
for i = 1:test_points
    % 生成测试信号 (单频余弦波)
    t = (0:N-1)/fs;
    signal = cos(2*pi*f_test(i)*t);
    
    % 分别使用四种算法估计频率
    f_est_rife = rife_estimate(signal, fs, use_window);
    f_est_jacobsen = jacobsen_estimator(signal, fs, use_window);
    f_est_fft = fft_peak_estimate(signal, fs, use_window);
    f_est_phase_diff = dft_phase_estimation(signal, fs, use_window);
    
    % 计算绝对误差
    errors_rife(i) = abs(f_est_rife - f_test(i));
    errors_jacobsen(i) = abs(f_est_jacobsen - f_test(i));
    errors_fft(i) = abs(f_est_fft - f_test(i));
    errors_phase_diff(i) = abs(f_est_phase_diff - f_test(i));
    
    % 显示进度
    if mod(i, 10) == 0
        fprintf('完成 %.0f%% (真实频率: %.2f Hz)\n', i/test_points*100, f_test(i));
    end
end

% 绘制误差曲线
figure;
plot(f_test, errors_fft, 'b', 'LineWidth', 1.5, 'DisplayName', 'FFT直接估计');
hold on;
plot(f_test, errors_rife, 'r', 'LineWidth', 1.5, 'DisplayName', 'Rife插值');
plot(f_test, errors_jacobsen, 'g', 'LineWidth', 1.5, 'DisplayName', '二次多项式插值');
plot(f_test, errors_phase_diff, 'm', 'LineWidth', 1.5, 'DisplayName', 'DFT Phase');
hold off;

% 图表标注
title(sprintf('频率估计算法性能对比 (N=%d, Fs=%d Hz)', N, fs));
xlabel('真实频率 (Hz)');
ylabel('绝对误差 (Hz)');
legend('show', 'Location', 'northwest');
grid on;
set(gca, 'FontSize', 12);

% 设置y轴范围，确保所有误差可见
max_error = max([errors_fft, errors_rife, errors_jacobsen, errors_phase_diff]);
ylim([0, max_error * 1.1]);

% 显示统计结果
fprintf('\n===== 性能总结 (窗函数: %s) =====\n', string(use_window));
fprintf('算法\t\t最大误差(Hz)\t平均误差(Hz)\t标准差(Hz)\n');
fprintf('FFT Peak\t%.4f\t\t%.4f\t\t%.4f\n', max(errors_fft), mean(errors_fft), std(errors_fft));
fprintf('Rife\t\t%.4f\t\t%.4f\t\t%.4f\n', max(errors_rife), mean(errors_rife), std(errors_rife));
fprintf('Jacobsen\t%.4f\t\t%.4f\t\t%.4f\n', max(errors_jacobsen), mean(errors_jacobsen), std(errors_jacobsen));
fprintf('Phase Diff\t%.4f\t\t%.4f\t\t%.4f\n', max(errors_phase_diff), mean(errors_phase_diff), std(errors_phase_diff));

end

%% FFT直接估计算法
function f_est = fft_peak_estimate(x, fs, use_window)
    N = length(x);
    
    % 窗函数处理
    if use_window
        win = hann(N)';
        x_win = x .* win;
        coherent_gain = sum(win)/N;
    else
        win = ones(1, N);
        x_win = x;
        coherent_gain = 1;
    end
    
    % FFT计算
    X = fft(x_win, N);
    P = abs(X)/N;
    P = P / coherent_gain;
    
    % 单边频谱
    P_single = P(1:floor(N/2)+1);
    P_single(2:end-1) = 2*P_single(2:end-1);
    f_axis = (0:floor(N/2)) * fs / N;
    
    % 找到最大谱线
    [~, k0] = max(P_single);
    f_est = f_axis(k0);
end

%% RIFE插值算法
function f_est = rife_estimate(x, fs, use_window)
    N = length(x);
    
    % 窗函数处理
    if use_window
        win = hann(N)';
        x_win = x .* win;
        coherent_gain = sum(win)/N;
    else
        win = ones(1, N);
        x_win = x;
        coherent_gain = 1;
    end
    
    % FFT计算
    X = fft(x_win, N);
    P = abs(X)/N;
    P = P / coherent_gain;
    
    % 单边频谱
    P_single = P(1:floor(N/2)+1);
    P_single(2:end-1) = 2*P_single(2:end-1);
    f_axis = (0:floor(N/2)) * fs / N;
    
    % 找到最大谱线
    [~, k0] = max(P_single);
    
    % 确定次大谱线方向
    if k0 == 1
        r = 1;
        mag_second = P_single(k0+1);
    elseif k0 == length(P_single)
        r = -1;
        mag_second = P_single(k0-1);
    else
        left_mag = P_single(k0-1);
        right_mag = P_single(k0+1);
        if right_mag > left_mag
            r = 1;
            mag_second = right_mag;
        else
            r = -1;
            mag_second = left_mag;
        end
    end
    
    % 计算幅度比
    mag_max = P_single(k0);
    alpha = mag_second / mag_max;
    
    % Rife算法核心公式 (汉宁窗专用)
    delta = (2*alpha - 1) / (alpha + 1);
    
    % 频率估计
    f_est = f_axis(k0) + r * delta * (fs / N);
end

%% 二次多项式插值算法
function f_est = jacobsen_estimator(x, fs, use_window)
    N = length(x);
    
    % 窗函数处理
    if use_window
        win = hann(N)';
        x_win = x .* win;
        coherent_gain = sum(win)/N;
    else
        win = ones(1, N);
        x_win = x;
        coherent_gain = 1;
    end
    
    % FFT计算
    X = fft(x_win, N);
    P = abs(X)/N;
    P = P / coherent_gain;
    
    % 单边频谱
    P_single = P(1:floor(N/2)+1);
    P_single(2:end-1) = 2*P_single(2:end-1);
    f_axis = (0:floor(N/2)) * fs / N;
    
    % 找到最大谱线
    [~, k0] = max(P_single);
    
    % 提取峰值点及相邻两点幅值
    alpha = P_single(k0);
    
    if k0 > 1
        beta = P_single(k0-1);
    else
        beta = 0;
    end
    
    if k0 < length(P_single)
        gamma = P_single(k0+1);
    else
        gamma = 0;
    end
    
    % Jacobsen插值公式
    if k0 > 1 && k0 < length(P_single)
        delta = (gamma - beta) / (2*(2*alpha - beta - gamma));
    else
        delta = 0; % 边界情况使用基本FFT
    end
    
    % 频率估计
    f_est = f_axis(k0) + delta * (fs / N);
end

%% DFT相位法
function f_est = dft_phase_estimation(signal, fs, use_window)
% 输入:
%   signal - 输入实正弦波信号
%   fs     - 采样频率 (Hz)
% 输出:
%   f0_hat     - 估计频率 (Hz)

N = length(signal);
if mod(N, 2) ~= 0
    error('信号长度N必须是偶数');
end
% 1. 将信号分成两个等长子序列
half = N/2;
s1 = signal(1:half);     % 前一半序列
s2 = signal(half+1:end); % 后一半序列
% 窗函数处理
    if use_window
        win = hann(N/2)';
        x1 = signal(1:N/2) .* win;
        x2 = signal(N/2+1:end) .* win;
    else
        x1 = signal(1:N/2);
        x2 = signal(N/2+1:end);
    end

% 2. 分别计算两个子序列的DFT
S1 = fft(x1);  % N/2点DFT
S2 = fft(x2);  % N/2点DFT

% 3. 只考虑正频率部分（避免负频率干扰）
S1 = S1(1:floor(half/2)+1);
S2 = S2(1:floor(half/2)+1);

% 4. 找到DFT幅度最大谱线的位置（在正频率范围内）
[~, idx1] = max(abs(S1));  % 第一段最大幅度位置
[~, idx2] = max(abs(S2));  % 第二段最大幅度位置

% 5. 检查两段的最大峰值位置是否一致（防止噪声干扰）
if idx1 ~= idx2
    % 如果不一致，选择两段中幅度和更大的位置
    if (abs(S1(idx1)) + abs(S2(idx1))) > (abs(S1(idx2)) + abs(S2(idx2)))
        k0 = idx1;
    else
        k0 = idx2;
    end
else
    k0 = idx1;
end
k0 = k0 - 1;  % MATLAB索引从1开始，转换为从0开始的索引

% 6. 提取最大谱线处的相位
phi1 = angle(S1(k0+1));  % 第一段DFT相位
phi2 = angle(S2(k0+1));  % 第二段DFT相位

% 7. 计算相位差并解模糊（保持在[-π, π]范围内）
delta_phi = phi2 - phi1;
delta_phi = mod(delta_phi + pi, 2*pi) - pi;  % 相位解缠绕

% 8. 计算子序列频率分辨率
delta_f_sub = fs / half;  % 子序列DFT频率分辨率

% 9. 估计相对频偏δ
delta = delta_phi / (2*pi);

% 10. 估计频率
f_est = (k0 + delta) * delta_f_sub;

end