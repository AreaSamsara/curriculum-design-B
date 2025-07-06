function selection = main_menu()
    stop = false;

    while(stop == false)  
        % 替换原始的menu函数为更美观的listdlg
        options = {
            '比较不同窗函数',...
            '比较目标落入不同间隔位置时的性能',...
            '比较不同SNR下的性能',...
            '比较不同干扰类型的性能',...
            '比较噪声/干扰服从不同概率分布函数的性能',...
            '比较不同频率估计算法的性能'};
        
        [sel, ok] = listdlg(...
            'PromptString', '请选择演示模式 (可单选)',...
            'ListString', options,...
            'SelectionMode', 'single',...
            'ListSize', [300 200],...  % 调整弹窗大小
            'Name', '演示界面菜单');
        
        if ~ok  % 如果用户取消选择
            break;
        end

        sel = deal_main_menu(sel)

        stop = old_run_estimate(sel);
        if stop
            break;
        end

        % 添加等待环节
        h = uicontrol('Style', 'pushbutton', 'String', '返回菜单',...
                       'Position', [10 10 100 30],...
                       'Callback', 'uiresume(gcbf)');
        uiwait(gcf);  % 暂停直到点击按钮
        delete(h);    % 删除临时按钮
    end
end

function selection = deal_main_menu(origin_selection)
    selection = 0; 
    switch origin_selection
        case(1)
            selection = 12;
        case(2)
            selection = 2;
        case(3)
            selection = 3;
        case(4)
            selection = 4;
        case(5)
            selection = 5;
        case(6)
            selection = 6;
        otherwise
            selection = 0; 
    end
end