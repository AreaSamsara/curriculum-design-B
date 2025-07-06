% 利用分段FFT算法
% 输入：
% Fs：采样频率
% x：采样序列
% U：折叠门限（无噪声时，设置为pi/10【自己试出来的】；有噪声时，设置为pi/3【论文给的】）
%   当相位差接近pi/-pi时，容易出现方向判断错误。为解决该问题设置了门限U：
%   当相位差落在[pi-U,pi+U]和[-(pi+U),pi+U]的范围内时，由次大谱的位置来决定估计方向。
% 输出：
% fre_estimate：频率的估计结果
function [fre_estimate] = Segment_FFT(Fs,x,U)
% 计算序列长度
[~,N] = size(x);
% 将原序列分为两部分，分别进行FFT,并转换为单边谱
x1 = x(1:N/2);x2 = x(N/2+1:end);
X1 = fftshift(fft(x1));X1 = 2*X1(N/4+1:end);
X2 = fftshift(fft(x2));X2 = 2*X2(N/4+1:end);
% 找到最大谱线与次大谱线
k1 = find(abs(X1) == max(abs(X1)),1);
% 次大谱线一定在左右相邻谱线中产生
if(abs(X1(k1-1))>abs(X1(k1+1)))
    k_sec1 = k1-1;
else
    k_sec1 = k1+1;
end
% 对于另一段也用相同方法,得到最大谱和左右相邻谱线
k2 = find(abs(X2) == max(abs(X2)),1);
if(abs(X2(k2-1))>abs(X2(k2+1)))
    k_sec2 = k2-1;
else
    k_sec2 = k2+1;
end
theta1 = angle(X1(k1));         %第一段最大谱的相位
theta2 = angle(X2(k2));         %第二段最大谱的相位
delta_theta = diff_theta(theta1,theta2,[k1,k2],[k_sec1,k_sec2],U);  % 根据规则计算出折叠相位差
f_est = delta_theta/(2*pi)*Fs/(N/2);       % 估计的偏差频率
% —————————————————————————————————————————— %
% 加权带来的误差太离谱了，这部分暂时先去掉了。
% % delta = f_est/(Fs/(N/2));                   % 相对频率偏移
% % 计算出最大谱和次大谱的加权平均
% theta1_wight = angle(X1(k1))*abs(X1(k1))/(abs(X1(k1))+abs(X1(k_sec1))) + angle(X1(k_sec1))*abs(X1(k_sec1))/(abs(X1(k1))+abs(X1(k_sec1)));
% theta2_wight = angle(X2(k2))*abs(X2(k2))/(abs(X2(k2))+abs(X2(k_sec2))) + angle(X2(k_sec2))*abs(X2(k_sec2))/(abs(X2(k2))+abs(X2(k_sec2)));
% delta_theta_wight = diff_theta(theta1_wight,theta2_wight,[k1,k2],[k_sec1,k_sec2],U);    % 根据规则计算出折叠相位差
% f_est_wight = delta_theta_wight/(2*pi)*Fs/(N/2);       % 估计的偏差频率
% % 相对频率偏移
% delta = (f_est+f_est_wight)/(Fs/(N/2));
% % 根据维纳滤波，求出最终估计值
% fre_bias = f_est*(1-abs(delta))^2/((1-abs(delta))^2+delta^2) + f_est_wight*delta^2/((1-abs(delta))^2+delta^2);
% —————————————————————————————————————————— %
fre_estimate = (k1-1)*Fs/(N/2) + f_est;
