% ������
% ���룺
% signal���źŲ�������
% N0���źŹ���
% type�������ֲ�����
%       Normal����˹�ֲ�
%       Poisson�����ɷֲ�
%       Chisquare�������ֲ�
% �����
% signal_with_noise���źŵ�������
function [signal_with_noise] = noise_add(signal,N0,type)
% �ɿ��ǲ�ͬ�ź����ͣ���Ʋ���������ָ���Ĺ���
if(strcmp(type,'Normal'))
    A = 0;      % �ֲ�����1
    B = 1;      % �ֲ�����2
%     k = sqrt(N0/(B^2+A^2));     % Ϊʹ���ʴﵽָ��ֵ�����ϵ�ϵ��
    k = sqrt(N0/B^2); 
    noise = k*random(type,A,B,size(signal));
%     SNR = 0.5/mean((abs(noise)).^2)     %for debug
elseif(strcmp(type,'Poisson'))
    A = 1;      % �ֲ�����1
%     k = sqrt(N0/(A+A^2));     % Ϊʹ���ʴﵽָ��ֵ�����ϵ�ϵ��
    k = sqrt(N0/A);
    noise = k*(random(type,A,size(signal)) - A);
%     SNR = 0.5/mean((abs(noise)).^2)     %for debug
elseif(strcmp(type,'Chisquare'))
    A = 1;      %�ֲ�����1
%     k = sqrt(N0/(2*A+A^2));     % Ϊʹ���ʴﵽָ��ֵ�����ϵ�ϵ��
    k = sqrt(N0/(2*A));
    noise = k*(random(type,A,size(signal)) - A);
%     SNR = 0.5/mean((abs(noise)).^2)     %for debug
end
signal_with_noise = signal + noise;
