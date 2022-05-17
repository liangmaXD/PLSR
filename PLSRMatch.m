function [Y_train_mean,Y_train_std,X_train_mean,X_train_std] = PLSRMatch(Num,outputPath,num_strap)
    fprintf('正在处理第%d组数据\n',Num);
    fprintf('-------------------\n');
    GroupPath = [outputPath '\Group'];
    StablePonitPath = [outputPath '\StablePoint'];
    if ~exist(StablePonitPath,'dir')
        mkdir(StablePonitPath);
    end
    if exist(fullfile(StablePonitPath,['StablePoint' num2str(Num,'%02d') '.mat']),'file')
        delete(fullfile(StablePonitPath,['StablePoint' num2str(Num,'%02d') '.mat']));
    end
    load(fullfile(GroupPath,['Group' num2str(Num,'%02d') '.mat']))    
    if size(X_train,2) > 10000
        for i = 1: size(X_train,2)
            [~,p(i)] = corr(X_train(:,i),Y_train);
        end
        StablePoint = find(p < 0.10); clear p
        X_train = X_train(:,StablePoint);
        X_test = X_test(:,StablePoint);
    end
    %% Z-值化
    Y_train_mean = mean(Y_train);
    Y_train_std = std(Y_train);
    X_train_mean = mean(X_train);
    X_train_std = std(X_train);
    Y_train = (Y_train - Y_train_mean) / Y_train_std;
    Y_test = (Y_test - Y_train_mean) / Y_train_std;
    X_train = (X_train - repmat(X_train_mean,[size(X_train,1) 1])) ...
        ./ repmat(X_train_std,[size(X_train,1) 1]);
    X_test = (X_test - repmat(X_train_mean,[size(X_test,1) 1])) ...
        ./ repmat(X_train_std,[size(X_test,1) 1]);
%    clear Y_train_mean Y_train_std X_train_mean X_train_std
    %% 计算偏最小二乘回归的回归系数
    % 这里套用while循环来实行特征选择
    if ~exist('StablePoint','var')
        StablePoint = 1: size(X_train,2);
    end
%     Times = 1;
    while 1
        Step = [0 1 2 3];
        StepSize = 0.5;
        optimalLVs = GetOptimalLVs(Y_train,X_train,Y_test,X_test);
        while 1
            for i = 1: length(Step)
                PointIndex{i} = PLSRbootstraping(Y_train,X_train,optimalLVs,num_strap,Step(i));
                [~,PRESS(i,1),R2(i,1)] = GetOptimalLVs(Y_train,X_train(:,PointIndex{i}),Y_test,X_test(:,PointIndex{i}));
            end
            [Step,StepSize] = IdentifyStep(PRESS,Step,StepSize);
            fprintf('Step = %d\n',Step)
            fprintf('StepSize = %d\n',-log2(StepSize))
            if log2(StepSize) < -10 || log2(StepSize) > 2
                break
            end
        end
        Loc = find(PRESS == max(PRESS));
        Loc = Loc(1);           % 有两个及以上的结果时取第一个
        if exist(fullfile(StablePonitPath,['StablePoint' num2str(Num,'%02d') '.mat']),'file')
            load(fullfile(StablePonitPath,['StablePoint' num2str(Num,'%02d') '.mat']));
            if length(StablePoint) == length(PointIndex{Loc})
                break;
            end
        end
        StablePoint = StablePoint(PointIndex{Loc});
        R2_pre = R2(Loc);
        X_train = X_train(:, PointIndex{Loc});
        X_test = X_test(:, PointIndex{Loc});
        save(fullfile(StablePonitPath,['StablePoint' num2str(Num,'%02d') '.mat']),'StablePoint','R2_pre')
        clear PointIndex PRESS StablePoint R2 R2_pre Loc
    end
    % 现在已经找到了最佳的潜在因子，下一步进行特征删选。
end


