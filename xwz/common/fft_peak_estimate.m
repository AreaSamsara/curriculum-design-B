function f_est = fft_peak_estimate(x,fs)
% FFT峰值检索频率估计
% 输入: x - 输入信号 
% 输出: f_est - 估计频率(Hz) 

    N = length(x);
   
    % 窗函数处理
    win = hann(N)';
    x_win = x .* win;

    % 计算FFT
    X = fft(x_win);
    
    % 找到幅度谱最大值位置
    [~, idx] = max(abs(X));
    
    % 频率估计
    f_est = (idx - 1) * fs / N;  % MATLAB索引从1开始
end