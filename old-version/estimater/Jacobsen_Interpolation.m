% ����Jacobsen�Ľ���ֱ�Ӳ�ֵ��ʽ
% ���룺
% Fs������Ƶ��
% x�������ź�����
% windows������������(Hamming/Hanning/Blackman/Black-Harris);����дʱ��Ĭ�ϲ��Ӵ�
% �����
% fre_estimate��Ƶ�ʵĹ��ƽ��
function [fre_estimate] = Jacobsen_Interpolation(Fs,x,windows)
% �����źŵ���
[~,N] = size(x);
% ������ѡ����ѡ���Ӧ����
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
        error('Jacobsen�Ľ��㷨�д�������������');
    end

% ����FFT���㣬���Ƶ��
X = fftshift(fft(x));
X = 2*X(N/2+1:end);     %��õ�����
% �ҳ�������ߵ�����
k = find(abs(X) == max(abs(X)),1);    
% ���ݹ�ʽ�����ϸ����ƫ��
delta = -Q*real((X(k+1)-X(k-1))/(2*X(k)-X(k-1)-X(k+1)));
% ���������Ƶ��
fre_estimate = (k+delta-1)*Fs/N;    % ��Ϊ�����Ǵ�0��ʼ�ģ�����Ҫ�����������Ļ�����-1
end