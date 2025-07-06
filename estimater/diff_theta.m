% ��Ҫ�۵���λ���֤���Ʒ��򲻳���
% ���룺
% theta1,theta2��������λ
% k,k_sec������׺ʹδ��׵�λ��
% U���۵�����
% ����λ��ӽ�pi/-piʱ�����׳��ַ����жϴ���Ϊ�������������������U��
% ����λ������[pi-U,pi+U]��[-(pi+U),pi+U]�ķ�Χ��ʱ���ɴδ��׵�λ�����������Ʒ���
function [delta] = diff_theta(theta1,theta2,k,k_sec,U)
delta = theta2 - theta1;
if((delta<-(pi+U)) || ...
        ((k_sec(1)>k(1)) &&(k_sec(2)>k(2)) && (delta>=-(pi+U)) && (delta<=-(pi-U))))
    delta = delta + 2*pi;
elseif((delta>pi+U) || ...
        ((k_sec(1)<k(1)) && (k_sec(2)<k(2)) && (delta>=pi-U) && (delta<=pi+U)))
    delta = delta - 2*pi;
end