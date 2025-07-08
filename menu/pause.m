function [] = pause()
    % 添加等待环节
    h = uicontrol('Style', 'pushbutton', 'String', '返回菜单',...
                    'Position', [10 10 100 30],...
                    'Callback', 'uiresume(gcbf)');
    uiwait(gcf);  % 暂停直到点击按钮
    delete(h);    % 删除临时按钮
end