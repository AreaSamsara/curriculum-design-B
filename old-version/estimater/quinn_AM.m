% ��������������s=�ź����� N=����Ҷ�任���� fs=����Ƶ��
% �������������f_quinn=ʹ�û���quinn�㷨��A&M�Ľ��㷨��Ƶ�ʹ�ֵ
function f_quinn = quinn_AM(s,N,fs)
    %% quinn�㷨
    s_fft = fft(s,N);       % ����fft����
    [~,index] = max(abs(s_fft));      % Ѱ������׷�
    a1 = real(s_fft(index-1)/s_fft(index));    % ��������׷�����ߵ����߼����м����
    a2 = real(s_fft(index+1)/s_fft(index));    % ��������׷����ұߵ����߼����м����
    deta1 = a1/(1-a1);     % ��������׷���������ߵó���Ƶ�����
    deta2 = a2/(a2-1);     % ��������׷����ұ����ߵó���Ƶ�����
    % �ж�ʹ���ĸ�Ƶ�����
    if deta1>0&&deta2>0    
        deta_1 = deta2;
    else
        deta_1 = deta1;
    end
    
    %% A&M�㷨
    index2 = index + deta_1;  
    p = [0.5,-0.5];
    ss = [0,0];
    for i = 1:2
        for j = 0:N-1
            ss(i) = ss(i) + s(j+1)*exp(-1j*2*pi*(index2+p(i))*j/N);  % ����ʽ
        end
    end
    deta3 = real((ss(1)+ss(2))/(ss(1)-ss(2)))/2;    % �ڶ��ε�����Ƶ�����

    deta_2 = deta_1 + deta3;      % ��Ƶ�����
    f_quinn = (index + deta_2)*fs/N;   % ��Ƶ�ʽ��������õ�Ƶ�ʹ�ֵ
end
