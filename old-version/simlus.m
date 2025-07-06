% ʵ�ָ�������µķ���
% ���룺
% fs������Ƶ��
% N����������
% situation�����������������ṹ�壩
%       .name���������޸��ŵ���������ideal/����������noise/�и�������Interfere
%       .f_seq��������Ƶ������
%       .SNR_dB�������������/�Ÿɱȣ������������������Ϊ��
%       .Iteration�����ؿ����������������ʹ�����ؿ��巨ʱ��ֵ���������Ϊ��
%       .noise_type���������ͣ�Ԫ�����飩
%               ��˹�ֲ�'Normal'/���ɷֲ�'Poisson'/�����ֲ�'Chisquare'
%       .Interfere������Ƶ��
% windows_type������Ĵ��������ͣ�Ԫ�����飩
%       ���Ӵ�(���δ�)'none'/������'hanning'/������'hamming'/����������'blackman'/��������-����˹��'blackmanharris'
% estimate_method��ѡ�õĹ����㷨��Ԫ�����飩
%       ֱ�ӷ�'direct'/A&M��������Quinn��'quinn&AM'/Jacobsen�Ľ���ֵ��'jacobsen'/�ֶ�FFT��'segment FFT'
function [] = simlus(fs,N,situation,windows_type,estimate_method)
%% ��������
delta_f = fs/N;                 %Ƶ�ʷֱ���
%% �������
if(strcmp(situation.name,'ideal'))
    f_seq = situation.f_seq;    %�������Ƶ��
    N_curve = length(windows_type)*length(estimate_method); %ͼ�������ߵĸ���
    legend_str = cell(1,N_curve);    %�洢ͼ��
    % ������������µ�ƫ�����ͼ
    figure;
    hold on;    grid on;
    title('�������޸������뵥Ƶ�źţ�����������ź�Ƶ�ʵĹ�ϵ')
    xlabel(['�ź�Ƶ��/Hz(Ƶ�׷ֱ���',num2str(delta_f),')'])
    ylabel('����ƫ��/Hz')
    for n_wintype = 1:length(windows_type)
        for n_estmethod = 1:length(estimate_method)
            [f_estimate,bias_error] = ideal_estimate(f_seq ,delta_f ,N ,windows_type{n_wintype} ,estimate_method{n_estmethod});
            plot(f_seq,abs(bias_error))
            legend_str{(n_wintype-1)*length(estimate_method)+n_estmethod} = [windows_type{n_wintype},';',estimate_method{n_estmethod}];
        end
    end
    legend(legend_str)
%% ���������
elseif(strcmp(situation.name,'noise'))
    f_seq = situation.f_seq;    % �������Ƶ��
    SNR_dB = situation.SNR_dB;  % ������������
    N_Iteration = situation.Iteration;      % ��������
    noise_type = situation.noise_type;      % ��������
    if(length(SNR_dB)==1 && length(windows_type)==1 && length(noise_type)==1)
        % �̶�����ȡ����������������ͣ��ı��ź�Ƶ�ʣ��������ֹ����㷨
        N_curve = length(estimate_method); %ͼ�������ߵĸ���
        legend_str = cell(1,N_curve);    %�洢ͼ��
        figure;
        hold on;    grid on
        title(['�����Ϊ',num2str(SNR_dB),'dB',...
            '������Ϊ',windows_type{1},...
            '����������Ϊ',noise_type{1}]);
        xlabel(['�ź�Ƶ��/Hz(Ƶ�׷ֱ���',num2str(delta_f),')']);
        ylabel('RMS���')
        for n_estmethod = 1:length(estimate_method)
            [f_estimate,RMSerror] = fre_estimate(f_seq ,delta_f ,SNR_dB ,N ,N_Iteration ,...
                noise_type{1} ,windows_type{1} ,estimate_method{n_estmethod});
            plot(f_seq,RMSerror)
            legend_str{n_estmethod} = [estimate_method{n_estmethod}];
        end
        legend(legend_str)
    elseif(length(f_seq)==1 && length(windows_type)==1 && length(noise_type)==1)
        % �̶��ź�Ƶ�ʡ����������������ͣ��ı�����ȣ�����4�ֹ����㷨
        N_curve = length(estimate_method); %ͼ�������ߵĸ���
        legend_str = cell(1,N_curve);    %�洢ͼ��
        figure;
        hold on;    grid on
        title(['�ź�Ƶ��Ϊ',num2str(f_seq),...
            '������Ϊ',windows_type{1},...
            '����������Ϊ',noise_type{1}]);
        xlabel('�����/dB');
        ylabel('RMS���')
        for n_estmethod = 1:length(estimate_method)
            [f_estimate,RMSerror] = noise_estimate(f_seq ,delta_f ,SNR_dB ,N ,N_Iteration ,...
                noise_type{1} ,windows_type{1} ,estimate_method{n_estmethod});
            plot(SNR_dB,RMSerror)
            legend_str{n_estmethod} = [estimate_method{n_estmethod}];
        end
        legend(legend_str);
    elseif(length(f_seq)==1 && length(windows_type)==1 && length(estimate_method)==1)
        % �̶��ź�Ƶ�ʡ��������������㷨���ı�����ȣ�����3����������
        N_curve = length(noise_type); %ͼ�������ߵĸ���
        legend_str = cell(1,N_curve);    %�洢ͼ��
        figure;
        hold on;    grid on
        title(['�ź�Ƶ��Ϊ',num2str(f_seq),...
            '������Ϊ',windows_type{1},...
            '�����㷨Ϊ',estimate_method{1}]);
        xlabel('�����/dB');
        ylabel('RMS���')
        for n_ntype = 1:length(noise_type)
            [f_estimate,RMSerror] = noise_estimate(f_seq ,delta_f ,SNR_dB ,N ,N_Iteration ,...
                noise_type{n_ntype} ,windows_type{1} ,estimate_method{1});
            plot(SNR_dB,RMSerror)
            legend_str{n_ntype} = [noise_type{n_ntype}];
        end
        legend(legend_str);
    elseif(length(f_seq)==1 && length(noise_type)==1 && length(estimate_method)==1)
        % �̶��ź�Ƶ�ʡ��������͡������㷨���ı�����ȣ�����5�ִ���������
        N_curve = length(windows_type); %ͼ�������ߵĸ���
        legend_str = cell(1,N_curve);    %�洢ͼ��
        figure;
        hold on;    grid on
        title(['�ź�Ƶ��Ϊ',num2str(f_seq),...
            '����������Ϊ',noise_type{1},...
            '�����㷨Ϊ',estimate_method{1}]);
        xlabel('�����/dB');
        ylabel('RMS���')
        for n_wintype = 1:length(windows_type)
            [f_estimate,RMSerror] = noise_estimate(f_seq ,delta_f ,SNR_dB ,N ,N_Iteration ,...
                noise_type{1} ,windows_type{n_wintype} ,estimate_method{1});
            disp(windows_type{n_wintype})
            plot(SNR_dB,RMSerror)
            legend_str{n_wintype} = [windows_type{n_wintype}];
        end
        legend(legend_str);
    end
elseif(strcmp(situation.name,'Interfere'))
    f_seq = situation.f_seq;    % �������Ƶ��
    SIR_dB = situation.SNR_dB;  % ��������Ÿɱ�
    f_int = situation.Interfere;    % ��Ƶ����Ƶ��
    if(length(SIR_dB)==1 && length(f_seq)==1 )
    % �̶��ŸɱȺ��ź�Ƶ�ʣ��ı���ŵ����λ��
        N_curve = length(windows_type)*length(estimate_method); %ͼ�������ߵĸ���
        legend_str = cell(1,N_curve);    %�洢ͼ��
        figure;
        hold on;    grid on;
        title(['��Ƶ�����ź��£��Ÿɱ�Ϊ',num2str(SIR_dB),'dB,�ź�Ƶ��Ϊ',num2str(f_seq),'Hz'])
        xlabel(['���λ��/Hz','(Ƶ�׷ֱ���',num2str(delta_f),')'])
        ylabel('����ƫ��/Hz')
        for n_wintype = 1:length(windows_type)
            for n_estmethod = 1:length(estimate_method)
                [f_estimate,bias_error] = int_estimate_f(f_seq ,f_int ,delta_f ,SIR_dB ,N ,windows_type{n_wintype} ,estimate_method{n_estmethod});
                plot(f_int-f_seq,bias_error)
                legend_str{(n_wintype-1)*length(estimate_method)+n_estmethod} = [windows_type{n_wintype},';',estimate_method{n_estmethod}];
            end
        end
        legend(legend_str);
    elseif(length(f_int)==1 && length(f_seq)==1)
    % �̶�����Ƶ�ʺ��ź�Ƶ�ʣ��ı��Ÿɱ�
        N_curve = length(windows_type)*length(estimate_method); %ͼ�������ߵĸ���
        legend_str = cell(1,N_curve);    %�洢ͼ��
        figure;
        hold on;    grid on;
        title(['��Ƶ�����ź��£�����Ƶ��Ϊ',num2str(f_int),'Hz,�ź�Ƶ��Ϊ',num2str(f_seq),'Hz',...
            '(Ƶ�׷ֱ���',num2str(delta_f),')'])
        xlabel('�Ÿɱ�/dB')
        ylabel('����ƫ��/Hz')
        for n_wintype = 1:length(windows_type)
            for n_estmethod = 1:length(estimate_method)
                [f_estimate,bias_error] = int_estimate_SIR(f_seq ,f_int ,delta_f ,SIR_dB ,N ,windows_type{n_wintype} ,estimate_method{n_estmethod});
                plot(SIR_dB,bias_error)
                legend_str{(n_wintype-1)*length(estimate_method)+n_estmethod} = [windows_type{n_wintype},';',estimate_method{n_estmethod}];
            end
        end
        legend(legend_str);
        
    end
end
