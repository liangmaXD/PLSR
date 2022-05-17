%% 预测测试集R2 
function [PRESS,R2] = CalculateR2(Y0,Y_pre)
    Y_m = mean(Y0);
    for i = 1: size(Y_pre,2)
        for k = 1: size(Y0,1)
            if k == 1
                TSS = (Y0(k) - Y_m)^2;
                RSS = (Y0(k) - Y_pre(k,i))^2;
            else
                TSS = (Y0(k) - Y_m)^2 + TSS;
                RSS = (Y0(k) - Y_pre(k,i))^2 + RSS;   % 剩余平方和
            end
        end
        R2(i) = 1 - RSS / TSS;
        PRESS(i) = RSS;
    end
end