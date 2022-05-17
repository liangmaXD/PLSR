function [Y_train_mean,Y_train_std,X_train_mean,X_train_std] = PLSRMatch_myself(Num,outputPath)
    fprintf('���ڴ����%d������\n',Num);
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
    %% Z-ֵ��
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
    %% ����ƫ��С���˻ع�Ļع�ϵ��
    % ��������whileѭ����ʵ������ѡ��
    if ~exist('StablePoint','var')
        StablePoint = 1: size(X_train,2);
    end
        CARS = carspls(X_train,Y_train,2,5,'autoscaling',50) %����Ȩ��ɸѡ�������ҵ������Сʱ������      
        Loc =  CARS.vsel    % CARS.vsel%����Ȩ��ɸѡ������
        StablePoint = StablePoint(:,Loc);
        [~,~,~,~,b] = plsregress(X_train(:,StablePoint),Y_train);
        Y_test_pre = [ones(size(X_test(:,StablePoint),1),1) X_test(:,StablePoint)] * b;
       [~,R2_pre] = CalculateR2(Y_test,Y_test_pre);
        X_train = X_train(:,StablePoint);
        X_test = X_test(:,StablePoint);
        save(fullfile(StablePonitPath,['StablePoint' num2str(Num,'%02d') '.mat']),'StablePoint','R2_pre')
        clear  StablePoint  R2_pre Loc CARS
    end
    % �����Ѿ��ҵ�����ѵ�Ǳ�����ӣ���һ����������ɾѡ��
