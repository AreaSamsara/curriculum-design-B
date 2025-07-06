% 加噪声
% 输入：
% signal：信号采样序列
% N0：信号功率
% type：噪声分布类型
%       Normal：高斯分布
%       Poisson：泊松分布
%       Chisquare：卡方分布
% 输出：
% signal_with_noise：信号叠加噪声
function [signal_with_noise] = noise_add(signal,N0,type)
% 由考虑不同信号类型，设计参数，满足指定的功率
if(strcmp(type,'Normal'))
    A = 0;      % 分布参数1
    B = 1;      % 分布参数2
%     k = sqrt(N0/(B^2+A^2));     % 为使功率达到指定值，乘上的系数
    k = sqrt(N0/B^2); 
    noise = k*random(type,A,B,size(signal));
%     SNR = 0.5/mean((abs(noise)).^2)     %for debug
elseif(strcmp(type,'Poisson'))
    A = 1;      % 分布参数1
%     k = sqrt(N0/(A+A^2));     % 为使功率达到指定值，乘上的系数
    k = sqrt(N0/A);
    noise = k*(random(type,A,size(signal)) - A);
%     SNR = 0.5/mean((abs(noise)).^2)     %for debug
elseif(strcmp(type,'Chisquare'))
    A = 1;      %分布参数1
%     k = sqrt(N0/(2*A+A^2));     % 为使功率达到指定值，乘上的系数
    k = sqrt(N0/(2*A));
    noise = k*(random(type,A,size(signal)) - A);
%     SNR = 0.5/mean((abs(noise)).^2)     %for debug
end
signal_with_noise = signal + noise;
