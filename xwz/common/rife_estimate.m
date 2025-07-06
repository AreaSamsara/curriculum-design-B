% ========================= Rife算法函数封装 =========================
function f_est = rife_estimate(x, fs)
% RIFE_ESTIMATE 使用Rife算法估计信号频率
%   f_est = rife_estimate(x, fs, use_window)
%   输入:
%       x - 输入信号
%       fs - 采样频率 (Hz)
%   输出:
%       f_est - 估计频率 (Hz)

    N = length(x);
    
    % 窗函数处理
    win = hann(N)';
    x_win = x .* win;
    coherent_gain = sum(win)/N;

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