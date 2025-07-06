% 需要折叠相位差，保证估计方向不出错
% 输入：
% theta1,theta2：两个相位
% k,k_sec：最大谱和次大谱的位置
% U：折叠门限
% 当相位差接近pi/-pi时，容易出现方向判断错误。为解决该问题设置了门限U：
% 当相位差落在[pi-U,pi+U]和[-(pi+U),pi+U]的范围内时，由次大谱的位置来决定估计方向。
function [delta] = diff_theta(theta1,theta2,k,k_sec,U)
delta = theta2 - theta1;
if((delta<-(pi+U)) || ...
        ((k_sec(1)>k(1)) &&(k_sec(2)>k(2)) && (delta>=-(pi+U)) && (delta<=-(pi-U))))
    delta = delta + 2*pi;
elseif((delta>pi+U) || ...
        ((k_sec(1)<k(1)) && (k_sec(2)<k(2)) && (delta>=pi-U) && (delta<=pi+U)))
    delta = delta - 2*pi;
end