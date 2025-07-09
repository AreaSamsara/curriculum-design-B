function f_est= dft_phase_estimation(signal, fs)
% 基于DFT相位的实正弦波频率和初相高精度估计方法
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

% 2. 分别计算两个子序列的DFT
S1 = fft(s1);  % N/2点DFT
S2 = fft(s2);  % N/2点DFT

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