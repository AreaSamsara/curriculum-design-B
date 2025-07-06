function [] = compareTargetPositionPerformance()

%% 主测试脚本
clear all; close all; clc;

% 参数设置
fs = 10240;          % 采样频率 (Hz)
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
    f_est_phase_diff = phase_diff_estimate(signal, fs, use_window);
    
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
plot(f_test, errors_fft, 'b', 'LineWidth', 1.5, 'DisplayName', 'FFT Peak');
hold on;
plot(f_test, errors_rife, 'r', 'LineWidth', 1.5, 'DisplayName', 'Rife');
plot(f_test, errors_jacobsen, 'g', 'LineWidth', 1.5, 'DisplayName', 'Jacobsen');
plot(f_test, errors_phase_diff, 'm', 'LineWidth', 1.5, 'DisplayName', 'Phase Diff');
hold off;

% 图表标注
title(sprintf('频率估计算法性能对比 (N=%d, Fs=%d Hz)', N, fs));
xlabel('真实频率 (Hz)');
ylabel('绝对误差 (Hz)');
legend('show', 'Location', 'northwest');
grid on;
set(gca, 'FontSize', 12);

% 添加FFT bin位置标记
y_lim = ylim;  % 获取当前y轴范围
for k = floor(f_start/delta_f):ceil(f_end/delta_f)
    bin_center = k * delta_f;
    line([bin_center, bin_center], y_lim, 'Color', [0.7 0.7 0.7], 'LineStyle', '--');
    text(bin_center, y_lim(2)*0.95, sprintf('Bin %d', k), ...
        'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', [0.4 0.4 0.4]);
end

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

%% FFT峰值检索算法
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

%% RIFE算法
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

%% Jacobsen算法
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

%% 相位差分法
function f_est = phase_diff_estimate(x, fs, use_window)
    N = length(x);
    
    % 窗函数处理
    if use_window
        win = hann(N)';
        x_win = x .* win;
    else
        win = ones(1, N);
        x_win = x;
    end
    
    % 获取解析信号
    analytic_signal = hilbert(x_win);
    
    % 计算瞬时相位
    phase = unwrap(angle(analytic_signal));
    
    % 计算相位差并解缠绕
    delta_phase = diff(phase);
    delta_phase_unwrap = atan2(sin(delta_phase), cos(delta_phase));
    
    % 计算瞬时频率
    instantaneous_freq = delta_phase_unwrap * fs / (2*pi);
    
    % 计算权重 (信号包络)
    magnitude = abs(analytic_signal);
    weights = (magnitude(1:end-1) + magnitude(2:end)) / 2;
    
    % 最终频率估计 (加权平均)
    f_est = sum(weights .* instantaneous_freq) / sum(weights);
end

end