%%
% ��Ҫѡȡ������PLSRMain������ʹ��PLSR����ֱ�����ģ����PLSR����
% PLSRMain������PLSR���������ƣ�������Ҫ���ý�����֤��setnum�������ۣ�setnum = 5
% num_strapΪbootstrap�Ĵ���
% TimesΪ�������̵Ľ��еĴ�����ÿѭ��һ��Times�����������ѡ��Ĺ���
% num����������벻Ҫ��ֵ����Ϊ��ֹ�������߶ϵ������Ӱ�죨������ֹ�����ܣ���
% �������һ�������ļ�StablePointN.mat�ļ��ĺ�׺����N���ٴ�����ʱ��num = N
% -------------------------------------------------------------------------------
% ע�⣺ѭ���ܴ��� = num_strap * Times��������ά�Ⱥܴ��ʱ�򣬽��飺num_strap = 200,
% Times = 100���Ա�֤��������ʱ��϶��ҽ����Ȼ�ȽϿ͹�
% -------------------------------------------------------------------------------
% ע�⣺������������Ҫ������������ⲿ��֤���������Եõ�ģ�͵ķ�������
% ��֤����Ҫ����ѵ�������б�׼��������
%%
%function [F,P,BETA,R2_pre,R2,R2_adjust,StablePoint_result,Y_train_mean,Y_train_std,X_train_mean,X_train_std] = PLSRMain(Y,X,setnum,num_strap,Times,num)
% function [F,P,BETA,R2_pre,R2,R2_adjust,StablePoint_result,Y_train_mean,Y_train_std,X_train_mean,X_train_std] = PLSRMain(Y,X,setnum,Times,num)
function [F,P,BETA,R2_pre,R2,R2_adjust,StablePoint_result,X_train_all,Y_train_all] = PLSRMain(Y,X,setnum,Times,num)
if nargin < 5   %�������Ĳ���С��5
    Times = 100;
end
%   Path = cd;
    Path = 'G:\Cortical_spinalcord\��ؽ��\PLSR\test';
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
        error('��ʼ��������������')        
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
 %   StablePoint_result = find(StablePoint_result >= Times/2);%��������˸Ķ�
     [Sort,position] = sort(StablePoint_result,'descend');
     StablePoint_result = position(find(Sort(1:30)));%�ҳ��ִ�����ǰ30������   
%     StablePoint_unique(times >= floor((1 / setnum)*Times));
    %% PLSR���ģ��
%     StablePoint_result = StablePoint_unique(CARS.vsel);
    X_result = X(:,StablePoint_result);
    [F,P,BETA,R2_pre,R2,R2_adjust] = PLSR(X_result,Y,setnum);
    save(fullfile(outputPath,'PLSR_result.mat'),'F','P','BETA','R2_pre','R2','R2_adjust','StablePoint_result');
end
