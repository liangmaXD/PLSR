%%
% 需要选取特征用PLSRMain函数，使用PLSR方法直接拟合模型用PLSR函数
% PLSRMain函数与PLSR的输入类似，但是需要设置交叉验证数setnum，如五折：setnum = 5
% num_strap为bootstrap的次数
% Times为整个流程的进行的次数，每循环一次Times完成依次特征选择的过程
% num正常情况下请不要赋值，是为防止死机或者断电带来的影响（流程终止后续跑）：
% 根据最后一个生成文件StablePointN.mat文件的后缀数字N，再次启动时令num = N
% -------------------------------------------------------------------------------
% 注意：循环总次数 = num_strap * Times，当特征维度很大的时候，建议：num_strap = 200,
% Times = 100，以保证程序运行时间较短且结果仍然比较客观
% -------------------------------------------------------------------------------
% 注意：本函数调用需要在外独立设置外部验证集！！！以得到模型的泛化能力
% 验证集需要根据训练集进行标准化！！！
%%
%function [F,P,BETA,R2_pre,R2,R2_adjust,StablePoint_result,Y_train_mean,Y_train_std,X_train_mean,X_train_std] = PLSRMain(Y,X,setnum,num_strap,Times,num)
% function [F,P,BETA,R2_pre,R2,R2_adjust,StablePoint_result,Y_train_mean,Y_train_std,X_train_mean,X_train_std] = PLSRMain(Y,X,setnum,Times,num)
function [F,P,BETA,R2_pre,R2,R2_adjust,StablePoint_result,X_train_all,Y_train_all] = PLSRMain(Y,X,setnum,Times,num)
if nargin < 5   %如果输入的参数小于5
    Times = 100;
end
%   Path = cd;
    Path = 'G:\Cortical_spinalcord\相关结果\PLSR\test';
    outputPath = fullfile(Path,'PLSR\result');
    %%
    if ~exist('num','var')
        num = 1;
        if ~exist(outputPath,'dir')
            mkdir(outputPath);
        end
        Group = DivGroup(Y,setnum,Times);
        GroupPath = fullfile(outputPath,'Group');
        DivTrain(Y,X,Group,GroupPath)
        X_train_all = [];
        Y_train_all = [];
        for i=1:Times
            load(fullfile(GroupPath,['Group' num2str(i,'%02d') '.mat']));
            X_train_all = [X_train_all;X_train];
            Y_train_all = [Y_train_all;Y_train];
        end
        elseif num > Times
        error('开始组数超过总组数')        
    else
        GroupPath = [outputPath '\Group'];
        load(fullfile(GroupPath,['Group' num2str(1,'%02d') '.mat']))
        X = [X_test;X_train];
        Y = [Y_test;Y_train];
        clear X_test X_train Y_test Y_train
    end
    for Num = 1: Times
%        [Y_train_mean,Y_train_std,X_train_mean,X_train_std] = PLSRMatch(Num,outputPath,num_strap);
%        [Y_train_mean,Y_train_std,X_train_mean,X_train_std] =  PLSRMatch_myself(Num,outputPath)
        PLSRMatch_myself(Num,outputPath)
    end
        StablePoint_result = zeros(1,size(X,2));
    for Num = 1: Times
        load([outputPath '\StablePoint\StablePoint' num2str(Num,'%02d') '.mat'])
        if(R2_pre>0)
            StablePoint_result(StablePoint) = StablePoint_result(StablePoint) + 1;
        end
    end
%     StablePoint_unique = unique(StablePoint_result);
%     CARS = carspls(X(:,StablePoint_unique),Y,2,set) ;
%     times = zeros(length(StablePoint_unique),1);
%     for i = 1: length(StablePoint_unique)
%         times(i) = sum(StablePoint_result == StablePoint_unique(i));
%     end
 %    StablePoint_result = find(StablePoint_result >= floor((1 / setnum)*Times));
 %   StablePoint_result = find(StablePoint_result >= Times/2);%这里进行了改动
     [Sort,position] = sort(StablePoint_result,'descend');
     StablePoint_result = position(find(Sort(1:30)));%找出现次数排前30的特征   
%     StablePoint_unique(times >= floor((1 / setnum)*Times));
    %% PLSR拟合模型
%     StablePoint_result = StablePoint_unique(CARS.vsel);
    X_result = X(:,StablePoint_result);
    [F,P,BETA,R2_pre,R2,R2_adjust] = PLSR(X_result,Y,setnum);
    save(fullfile(outputPath,'PLSR_result.mat'),'F','P','BETA','R2_pre','R2','R2_adjust','StablePoint_result');
end
