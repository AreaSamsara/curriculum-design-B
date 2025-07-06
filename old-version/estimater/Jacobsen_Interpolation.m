% 利用Jacobsen改进的直接插值公式
% 输入：
% Fs：采样频率
% x：采样信号序列
% windows：窗函数类型(Hamming/Hanning/Blackman/Black-Harris);不填写时，默认不加窗
% 输出：
% fre_estimate：频率的估计结果
function [fre_estimate] = Jacobsen_Interpolation(Fs,x,windows)
% 计算信号点数
[~,N] = size(x);
% 根据所选窗，选择对应参数
if(strcmp(windows,'none'))
        Q = 1;
    elseif(strcmp(windows,'hamming'))
        Q = 0.6;
    elseif(strcmp(windows,'hanning'))
        Q = 0.55;
    elseif(strcmp(windows,'blackman'))
        Q = 0.55;
    elseif(strcmp(windows,'blackmanharris'))
        Q = 0.56;
    else
        error('Jacobsen改进算法中窗函数输入有误');
    end

% 进行FFT运算，获得频谱
X = fftshift(fft(x));
X = 2*X(N/2+1:end);     %获得单边谱
% 找出最大谱线的序列
k = find(abs(X) == max(abs(X)),1);    
% 根据公式计算出细估计偏差
delta = -Q*real((X(k+1)-X(k-1))/(2*X(k)-X(k-1)-X(k+1)));
% 计算出估计频率
fre_estimate = (k+delta-1)*Fs/N;    % 因为序列是从0开始的，所以要在数组索引的基础上-1
end