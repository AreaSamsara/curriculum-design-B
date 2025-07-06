% ��������demo������˲���ȥ�˰�ctrl+C��
% �ź�Ƶ�ʡ���situation.f_seq
% ����Ƶ�ʡ���situation.Interfere
% �������͡���situation.noise_type
% ������������situation.Iteration
% �����/�Ÿɱȡ���situation.SNR_dB
% ѡ�񴰺�����ֱ���޸�windows_type
% ѡ����Ʒ�����ֱ���޸�estimate_method

addpath(genpath(pwd));
clear;close all;
%% ��������
stop = false;

while(stop == false)        
    % interactiveUI();
    % break;

    % sel = menu('��ʾ����˵�',...
    %     '������������ƫ��-�ź�Ƶ��',...
    %     '����������������£�RMS���-�ź�Ƶ��',...
    %     '����������������£�RMS���-�ź�Ƶ��',...
    %     '��������RMS���-�����',...
    %     '����������ͬ�������ͣ�RMS���-�����',...
    %     '����������ͬ��������RMS���-�����',...
    %     '�и��ţ����Ÿɱ��£�����ƫ��-�������λ��',...
    %     '�и��ţ����Ÿɱ��£�����ƫ��-�������λ��',...
    %     '�и��ţ�����ƫ��-�Ÿɱ�',...
    %     '�Զ���ģʽ',...
    %     'Exit',...
    %     '����');

    % �滻ԭʼ��menu����Ϊ�����۵�listdlg
    options = {
        '������������ƫ��-�ź�Ƶ��',...
        '����������������£�RMS���-�ź�Ƶ��',...
        '����������������£�RMS���-�ź�Ƶ��',...
        '��������RMS���-�����',...
        '����������ͬ�������ͣ�RMS���-�����',...
        '����������ͬ��������RMS���-�����',...
        '�и��ţ����Ÿɱ��£�����ƫ��-�������λ��',...
        '�и��ţ����Ÿɱ��£�����ƫ��-�������λ��',...
        '�и��ţ�����ƫ��-�Ÿɱ�',...
        '�Զ���ģʽ',...
        'Exit',...
        '����'};
    
    [sel, ok] = listdlg(...
        'PromptString', '��ѡ����ʾģʽ (�ɵ�ѡ)',...
        'ListString', options,...
        'SelectionMode', 'single',...
        'ListSize', [300 200],...  % ����������С
        'Name', '��ʾ����˵�');
    
    if ~ok  % ����û�ȡ��ѡ��
        flag = 0;
        break;
    end

    stop = runEstimate(sel);
    if stop
        break;
    end

    % ��ӵȴ�����
    h = uicontrol('Style', 'pushbutton', 'String', '���ز˵�',...
                 'Position', [10 10 100 30],...
                 'Callback', 'uiresume(gcbf)');
    uiwait(gcf);  % ��ֱͣ�������ť
    delete(h);    % ɾ����ʱ��ť

end

function stop = runEstimate(sel)
    %% ��ʼ������
    fs = 2560;  % ����Ƶ��
    N = 256;    % ��������
    
    stop = false;

    switch(sel)
        case(1)
            %% ������������ƫ��-�ź�Ƶ��
            situation.name = 'ideal';
            situation.f_seq = 210:0.5:220;  %���÷���Ƶ�ʣ���Ҫʱ�����޸�
            situation.SNR_dB = [];
            situation.Iteration = [];
            situation.noise_type = [];
            situation.Interfere = [];
            % windows_type = {'none','hanning','hamming','blackman','blackmanharris'};
            windows_type = {'none'};    % ѡ�񴰺��������޸�
            estimate_method = {'direct','quinn&AM','jacobsen','segment FFT'};   %ѡ������㷨�����޸�
            % estimate_method = {'jacobsen','segment FFT'};
            % estimate_method = {'jacobsen'};
            simlus(fs,N,situation,windows_type,estimate_method);
            case(2)
       %% ����������������£�RMS���-�ź�Ƶ��
            situation.name = 'noise';
            situation.f_seq = 210:0.5:220;  % ���÷���Ƶ�ʣ����޸�
            situation.SNR_dB = 10;          % ��������ȣ��ɸģ���ֻ����һ��ֵ
            situation.Iteration = 1e3;
            situation.Interfere = [];
            situation.noise_type = {'Normal'};      % �������ͣ��ɸģ���ֻ��Ϊһ��
            windows_type = {'none'};                % ���������ͣ��ɸģ���ֻ��Ϊһ��
            estimate_method = {'direct','quinn&AM','jacobsen','segment FFT'};   % �����㷨���ɸģ�������Ϊ���
            simlus(fs,N,situation,windows_type,estimate_method);
        case(3)
        %% ����������������£�RMS���-�ź�Ƶ��
            situation.name = 'noise';
            situation.f_seq = 210:0.5:220;  % ���÷���Ƶ�ʣ����޸�
            situation.SNR_dB = 0;           % ��������ȣ��ɸģ���ֻ����һ��ֵ
            situation.Iteration = 1e3;
            situation.Interfere = [];
            situation.noise_type = {'Normal'};       % �������ͣ��ɸģ���ֻ��Ϊһ��
            windows_type = {'none'};                 % ���������ͣ��ɸģ���ֻ��Ϊһ��
            estimate_method = {'direct','quinn&AM','jacobsen','segment FFT'};   % �����㷨���ɸģ�������Ϊ���
            simlus(fs,N,situation,windows_type,estimate_method);
        case(4)
       %% ��������RMS���-�����
            situation.name = 'noise';   
            situation.f_seq = 215;          % ���÷���Ƶ�ʣ��ɸ�,��ֻ����һ��ֵ
            situation.SNR_dB = 0:0.5:10;    % ��������ȣ��ɸ�
            situation.Iteration = 1e3;
            situation.Interfere = [];
            situation.noise_type = {'Normal'};  % �������ͣ��ɸģ���ֻ��Ϊһ��
            windows_type = {'none'};            % ���������ͣ��ɸģ���ֻ��Ϊһ��
            estimate_method = {'quinn&AM','jacobsen','segment FFT'};    % �����㷨���ɸģ�������Ϊ���
            simlus(fs,N,situation,windows_type,estimate_method);
        case(5)
       %% ����������ͬ�������ͣ�RMS���-�����
            situation.name = 'noise';
            situation.f_seq = 215;          % ���÷���Ƶ�ʣ��ɸ�,��ֻ����һ��ֵ
            situation.SNR_dB = 0:0.5:10;    % ��������ȣ��ɸ�
            situation.Iteration = 1e3;
            situation.Interfere = [];
            situation.noise_type = {'Normal','Poisson','Chisquare'};    % �������ͣ��ɸģ�������Ϊ���
            windows_type = {'none'};                    % ���������ͣ��ɸģ���ֻ��Ϊһ��
            estimate_method = {'segment FFT'};          % �����㷨���ɸģ���ֻ��Ϊһ��
            simlus(fs,N,situation,windows_type,estimate_method);
        case(6)
       %% ����������ͬ��������RMS���-�����
            situation.name = 'noise';
            situation.f_seq = 220;          % ���÷���Ƶ�ʣ��ɸ�,��ֻ����һ��ֵ
            situation.SNR_dB = 0:0.5:10;    % ��������ȣ��ɸ�
            situation.Iteration = 1e3;
            situation.Interfere = [];
            situation.noise_type = {'Normal'};      % �������ͣ��ɸģ���ֻ��Ϊһ��
            windows_type = {'none','hanning','hamming','blackman','blackmanharris'};    % ���������ͣ��ɸģ�������Ϊ���
            %estimate_method = {'quinn&AM'};         % �����㷨���ɸģ���ֻ��Ϊһ��
            %estimate_method = {'jacobsen'};
            estimate_method = {'segment FFT'};
            simlus(fs,N,situation,windows_type,estimate_method);
        case(7)
       %% �и��ţ����Ÿɱ��£�����ƫ��-�������λ��
            situation.name = 'Interfere';
            situation.f_seq = 215;      % ���÷���Ƶ�ʣ��ɸ�,��ֻ����һ��ֵ
            situation.SNR_dB = 10;      % �����Ÿɱȣ��ɸģ���ֻ����һ��ֵ
            situation.Iteration = [];   
            situation.noise_type = [];
            situation.Interfere = 210:0.5:220;  %���õ�Ƶ����Ƶ��
    % windows_type = {'none','hanning','hamming','blackman','blackmanharris'};
            windows_type = {'none'};            % ���������ͣ��ɸ�
            estimate_method = {'quinn&AM','jacobsen','segment FFT'};    % �����㷨���ɸ�
            simlus(fs,N,situation,windows_type,estimate_method);
        case(8)
       %% �и��ţ����Ÿɱ��£�����ƫ��-�������λ��
            situation.name = 'Interfere';
            situation.f_seq = 215;      % ���÷���Ƶ�ʣ��ɸ�,��ֻ����һ��ֵ
            situation.SNR_dB = 0;       % �����Ÿɱȣ��ɸģ���ֻ����һ��ֵ
            situation.Iteration = [];
            situation.noise_type = [];
            situation.Interfere = 210:0.5:220;      %���õ�Ƶ����Ƶ��
    % windows_type = {'none','hanning','hamming','blackman','blackmanharris'};
            windows_type = {'none'};            % ���������ͣ��ɸ�
            estimate_method = {'quinn&AM','jacobsen','segment FFT'};    % �����㷨���ɸ�
            simlus(fs,N,situation,windows_type,estimate_method);
        case(9)
       %% �и��ţ�����ƫ��-�Ÿɱ�
            situation.name = 'Interfere';
            situation.f_seq = 215;      % ���÷���Ƶ�ʣ��ɸ�,��ֻ����һ��ֵ
            situation.SNR_dB = 0:0.5:10;    % �����Ÿɱȣ��ɸ�
            situation.Iteration = [];
            situation.noise_type = [];
            situation.Interfere = 220;      %���õ�Ƶ����Ƶ�ʣ��ɸ�,��ֻ����һ��ֵ
    % windows_type = {'none','hanning','hamming','blackman','blackmanharris'};
            windows_type = {'none'};            % ���������ͣ��ɸ�
            estimate_method = {'quinn&AM','jacobsen','segment FFT'};        % �����㷨���ɸ�
            simlus(fs,N,situation,windows_type,estimate_method);
        case(10)
       %% �Զ���ģʽ
            disp("��������'ideal'/����������'noise'/�и�������'Interfere'")
            situation.name = input('����������');
            if(strcmp(situation.name,'ideal'))
           %% �������
                situation.f_seq = 210:0.5:220;  % ���÷���Ƶ�ʣ��ɸ�
                situation.SNR_dB = [];
                situation.Iteration = [];
                situation.noise_type = [];
                situation.Interfere = [];
                disp("���Ӵ�(���δ�)'none'/������'hanning'/������'hamming'/����������'blackman'/��������-����˹��'blackmanharris'")
                windows_type = input('ѡ�񴰺���(��Ԫ����ʾ)��');
                disp("ֱ�ӷ�'direct'/A&M��������Quinn��'quinn&AM'/Jacobsen�Ľ���ֵ��'jacobsen'/�ֶ�FFT��'segment FFT'")
                estimate_method = input('ѡ������㷨(��Ԫ����ʾ)��');
            elseif(strcmp(situation.name,'Interfere'))
           %% �и������
                disp("�Ÿɱ�'SIR'/���λ��'df'")
                mod_sel = input('ѡ������꣺');
                if(strcmp(mod_sel,'SIR'))
              %% ���ǹ���ƫ��-SIR
                    situation.f_seq = 215;          % ���÷���Ƶ�ʣ��ɸ�,��ֻ����һ��ֵ
                    situation.SNR_dB = 0:0.5:10;    % �����Ÿɱȣ��ɸ�
                    situation.Iteration = [];
                    situation.noise_type = [];
                    situation.Interfere = 220;      % ���õ�Ƶ����Ƶ�ʣ��ɸģ���ֻ����һ��ֵ
                    disp("���Ӵ�(���δ�)'none'/������'hanning'/������'hamming'/����������'blackman'/��������-����˹��'blackmanharris'")
                    windows_type = input('ѡ�񴰺���(��Ԫ����ʾ)��');
                    disp("ֱ�ӷ�'direct'/A&M��������Quinn��'quinn&AM'/Jacobsen�Ľ���ֵ��'jacobsen'/�ֶ�FFT��'segment FFT'")
                    estimate_method = input('ѡ������㷨(��Ԫ����ʾ)��');
                else
              %% ���ǹ���ƫ��-���λ��
                    situation.name = 'Interfere';   
                    situation.f_seq = 215;      % ���÷���Ƶ�ʣ��ɸ�,��ֻ����һ��ֵ
                    situation.SNR_dB = input('�Ÿɱ�Ϊ��');
                    situation.Iteration = [];
                    situation.noise_type = [];
                    situation.Interfere = 210:0.5:220;  % ���õ�Ƶ����Ƶ�ʣ��ɸ�
                    disp("���Ӵ�(���δ�)'none'/������'hanning'/������'hamming'/����������'blackman'/��������-����˹��'blackmanharris'")
                    windows_type = input('ѡ�񴰺���(��Ԫ����ʾ)��');
                    disp("ֱ�ӷ�'direct'/A&M��������Quinn��'quinn&AM'/Jacobsen�Ľ���ֵ��'jacobsen'/�ֶ�FFT��'segment FFT'")
                    estimate_method = input('ѡ������㷨(��Ԫ����ʾ)��');
                end
            elseif(strcmp(situation.name,'noise'))
           %% ���������
                disp("�����'SNR'/�ź�Ƶ��'f'")
                mod_sel = input('ѡ������꣺');
                if(strcmp(mod_sel,'SNR'))
              %% ����RMSƫ��-�����
                    situation.f_seq = 215;       % ���÷���Ƶ�ʣ��ɸ�,��ֻ����һ��ֵ
                    situation.SNR_dB = 0:0.5:10;    % ��������ȣ��ɸ�
                    situation.Iteration = 1e3;
                    disp("��˹�ֲ�'Normal'/���ɷֲ�'Poisson'/�����ֲ�'Chisquare'")
                    situation.noise_type = input('ѡ����������(��Ԫ����ʾ)��');
                    disp("���Ӵ�(���δ�)'none'/������'hanning'/������'hamming'/����������'blackman'/��������-����˹��'blackmanharris'")
                    windows_type = input('ѡ�񴰺���(��Ԫ����ʾ)��');
                    disp("ֱ�ӷ�'direct'/A&M��������Quinn��'quinn&AM'/Jacobsen�Ľ���ֵ��'jacobsen'/�ֶ�FFT��'segment FFT'")
                    estimate_method = input('ѡ������㷨(��Ԫ����ʾ)��');
                else
              %% ����RMSƫ��-�ź�Ƶ��
                    situation.f_seq = 210:0.5:220;      % ���÷���Ƶ�ʣ��ɸ�
                    situation.SNR_dB = input('ѡ������ȣ�(dB)');
                    situation.Iteration = 1e3;
                    disp("��˹�ֲ�'Normal'/���ɷֲ�'Poisson'/�����ֲ�'Chisquare'")
                    situation.noise_type = input('ѡ����������(��Ԫ����ʾ)��');
                    disp("���Ӵ�(���δ�)'none'/������'hanning'/������'hamming'/����������'blackman'/��������-����˹��'blackmanharris'")
                    windows_type = input('ѡ�񴰺���(��Ԫ����ʾ)��');
                    disp("ֱ�ӷ�'direct'/A&M��������Quinn��'quinn&AM'/Jacobsen�Ľ���ֵ��'jacobsen'/�ֶ�FFT��'segment FFT'")
                    estimate_method = input('ѡ������㷨(��Ԫ����ʾ)��');
                end
            end
            simlus(fs,N,situation,windows_type,estimate_method);
        case(11)
       %% �˳�����
            stop = true;
        case(12)
            window();
    end
end