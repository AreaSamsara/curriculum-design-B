function f_est = quinn_estimate(x, fs)
    % Quinn A&M算法
    N = length(x);
    X = fft(x);
    Pxx = abs(X);
    [~, k] = max(Pxx(1:N/2));
    
    if k == 1 || k == N/2
        f_est = (k-1) * fs / N;  % 无法插值，直接返回
        return;
    end
    
    r = Pxx(k+1) / Pxx(k);
    q = Pxx(k-1) / Pxx(k);
    
    if r > q
        % 右侧峰值
        b = (r - 1) / (2 * (2 - r));
    else
        % 左侧峰值
        b = (q - 1) / (2 * (2 - q));
    end
    
    f_est = (k-1 + b) * fs / N;
end