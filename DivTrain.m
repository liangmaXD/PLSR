function DivTrain(Y,X,Group,outputPath)
%     if exist(outputPath,'dir')
%         warning(['文件目录已存在，请确认是否为当前数据文件！'...
%             '若是，请输入任意键继续。']);
%         pause;
%     else
        mkdir(outputPath);
%     end
    for i = 1: length(Group)
        Train_Med = 1:length(Y);
        l = Group{i};
        Test{i} = Group{i};
        %         for j = [1:i-1 i+1:length(Group)]
        %             Train_Med = [Train_Med,Group{j}];
        %         end
        Train_Med(l) = [];
        Train{i} = Train_Med;
    end
    for i = 1: length(Group)
        if ~exist(fullfile(outputPath,['Group' num2str(i,'%02d') '.mat']),'file')
            clear Y_test X_test Y_train X_train
            test = Test{i};
            train = Train{i};
            Y_test = Y(test,:);
            X_test = X(test,:);
            Y_train = Y(train,:);
            X_train = X(train,:);
            save(fullfile(outputPath,['Group' num2str(i,'%02d') '.mat']),'Y_test','X_test','Y_train','X_train')
        end
    end
end



% %% 留一法
% function DivTrain(Y,X,group,outputPath)
%       
%     mkdir(outputPath);
%     for i = 1: length(group{1,1})
%         Train_Med = 1:length(Y);
%         l = group{1}(1,i);
%         Test{i} = group{1}(1,i);
%         Train_Med(l) = [];
%         Train{i} = Train_Med;
%     end
%     for i = 1: length(group{1,1})
%         if ~exist(fullfile(outputPath,['Group' num2str(i,'%02d') '.mat']),'file')
%             clear Y_test X_test Y_train X_train
%             test = Test{i};
%             train = Train{i};
%             Y_test = Y(test,:);
%             X_test = X(test,:);
%             Y_train = Y(train,:);
%             X_train = X(train,:);
%             save(fullfile(outputPath,['Group' num2str(i,'%02d') '.mat']),'Y_test','X_test','Y_train','X_train')
%         end
%     end
% end
