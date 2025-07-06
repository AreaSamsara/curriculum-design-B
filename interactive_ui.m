function interactive_ui()
    % 创建主界面
    fig = uifigure('Name', '参数设置', 'Position', [100 100 400 250]);
    
    % 1. 填空输入框
    uilabel(fig, 'Text', '信号频率 (Hz):', 'Position', [20 180 100 22]);
    freqEdit = uieditfield(fig, 'numeric', 'Position', [130 180 100 22], 'Value', 215);
    
    uilabel(fig, 'Text', '信噪比 (dB):', 'Position', [20 140 100 22]);
    snrEdit = uieditfield(fig, 'numeric', 'Position', [130 140 100 22], 'Value', 10);
    
    % 2. 下拉选择器
    uilabel(fig, 'Text', '窗函数类型:', 'Position', [20 100 100 22]);
    windowDrop = uidropdown(fig, 'Items', {'无窗', '汉宁窗', '汉明窗', '布莱克曼窗'}, ...
                           'Position', [130 100 150 22], 'Value', '汉宁窗');
    
    % 3. 操作按钮
    uibutton(fig, 'push', 'Text', '开始计算', 'Position', [100 40 100 30], ...
             'ButtonPushedFcn', @(btn,event) calculate(freqEdit.Value, snrEdit.Value, windowDrop.Value));
    
    uibutton(fig, 'push', 'Text', '取消', 'Position', [220 40 100 30], ...
             'ButtonPushedFcn', @(btn,event) delete(fig));
end

function calculate(freq, snr, windowType)
    % 处理用户输入
    disp(['频率: ', num2str(freq), ' Hz']);
    disp(['信噪比: ', num2str(snr), ' dB']);
    disp(['窗函数: ', windowType]);
    
    % 在这里调用你的计算函数
    % simlus(fs, N, situation, windows_type, estimate_method);
end