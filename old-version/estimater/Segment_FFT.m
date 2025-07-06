% ���÷ֶ�FFT�㷨
% ���룺
% Fs������Ƶ��
% x����������
% U���۵����ޣ�������ʱ������Ϊpi/10���Լ��Գ����ġ���������ʱ������Ϊpi/3�����ĸ��ġ���
%   ����λ��ӽ�pi/-piʱ�����׳��ַ����жϴ���Ϊ�������������������U��
%   ����λ������[pi-U,pi+U]��[-(pi+U),pi+U]�ķ�Χ��ʱ���ɴδ��׵�λ�����������Ʒ���
% �����
% fre_estimate��Ƶ�ʵĹ��ƽ��
function [fre_estimate] = Segment_FFT(Fs,x,U)
% �������г���
[~,N] = size(x);
% ��ԭ���з�Ϊ�����֣��ֱ����FFT,��ת��Ϊ������
x1 = x(1:N/2);x2 = x(N/2+1:end);
X1 = fftshift(fft(x1));X1 = 2*X1(N/4+1:end);
X2 = fftshift(fft(x2));X2 = 2*X2(N/4+1:end);
% �ҵ����������δ�����
k1 = find(abs(X1) == max(abs(X1)),1);
% �δ�����һ�����������������в���
if(abs(X1(k1-1))>abs(X1(k1+1)))
    k_sec1 = k1-1;
else
    k_sec1 = k1+1;
end
% ������һ��Ҳ����ͬ����,�õ�����׺�������������
k2 = find(abs(X2) == max(abs(X2)),1);
if(abs(X2(k2-1))>abs(X2(k2+1)))
    k_sec2 = k2-1;
else
    k_sec2 = k2+1;
end
theta1 = angle(X1(k1));         %��һ������׵���λ
theta2 = angle(X2(k2));         %�ڶ�������׵���λ
delta_theta = diff_theta(theta1,theta2,[k1,k2],[k_sec1,k_sec2],U);  % ���ݹ��������۵���λ��
f_est = delta_theta/(2*pi)*Fs/(N/2);       % ���Ƶ�ƫ��Ƶ��
% ������������������������������������������������������������������������������������ %
% ��Ȩ���������̫�����ˣ��ⲿ����ʱ��ȥ���ˡ�
% % delta = f_est/(Fs/(N/2));                   % ���Ƶ��ƫ��
% % ���������׺ʹδ��׵ļ�Ȩƽ��
% theta1_wight = angle(X1(k1))*abs(X1(k1))/(abs(X1(k1))+abs(X1(k_sec1))) + angle(X1(k_sec1))*abs(X1(k_sec1))/(abs(X1(k1))+abs(X1(k_sec1)));
% theta2_wight = angle(X2(k2))*abs(X2(k2))/(abs(X2(k2))+abs(X2(k_sec2))) + angle(X2(k_sec2))*abs(X2(k_sec2))/(abs(X2(k2))+abs(X2(k_sec2)));
% delta_theta_wight = diff_theta(theta1_wight,theta2_wight,[k1,k2],[k_sec1,k_sec2],U);    % ���ݹ��������۵���λ��
% f_est_wight = delta_theta_wight/(2*pi)*Fs/(N/2);       % ���Ƶ�ƫ��Ƶ��
% % ���Ƶ��ƫ��
% delta = (f_est+f_est_wight)/(Fs/(N/2));
% % ����ά���˲���������չ���ֵ
% fre_bias = f_est*(1-abs(delta))^2/((1-abs(delta))^2+delta^2) + f_est_wight*delta^2/((1-abs(delta))^2+delta^2);
% ������������������������������������������������������������������������������������ %
fre_estimate = (k1-1)*Fs/(N/2) + f_est;
