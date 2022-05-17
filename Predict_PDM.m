%% Ԥ�⹲�����
%% ·���ı�
%������ GroupPath = 'G:\Individual\HC_network\result_predict\indivual_net_PLSR';%�ǵø�·��
%������ ����·��outputPath = 'G:\Individual\HC_network\test';
%PLSRMain   Path = 'G:\Individual\HC_network\test';
%% ��������*����
% clc,clear
% datapath = 'E:\PDM_Network_all\timepoint_2\result_sFNC\sFNC.mat';%ʱ���3�ľ�Ϣ̬
% load(datapath);
% datapath = 'E:\PDM_Network_all\timepoint_2\Analyze_dFNC\dFNC_NO_state.mat';%ʱ���3�Ķ�̬δ��״̬
% load(datapath);
% datapath = 'E:\PDM_Network_all\timepoint_2\Analyze_dFNC\dFNC.mat';%ʱ���3�Ķ�̬��״̬
% load(datapath);
% 
% %-----------------------������
% X_raw = sFNC([1:41],:);% sFNC
% X_raw = std_window_all([1:41],:);% dFNC
% state1_all(find(isnan(state1_all)))=0;% �ҵ�NaN�滻Ϊ0
% state2_all(find(isnan(state2_all)))=0;% �ҵ�NaN�滻Ϊ0
% X_raw = state1_all([1:41],:);
% X_raw = state2_all([1:41],:);
% % X_raw = [sFNC([1:41],:),std_window_all([1:41],:)];
% %-----------------------����
% X_raw = sFNC([42:end],:);% sFNC
% X_raw = std_window_all([42:end],:);% dFNC
% X_raw = state1_all([42:end],:);
% X_raw = state2_all([42:end],:);
% % X_raw = [sFNC([42:end],:),std_window_all([42:end],:)];
% %-----------------------����:�ӷ�����������ѡ������
% % path ='E:\PDM_Network_all\classify\sFNC\models.mat';
% % path ='E:\PDM_Network_all\classify\dFNC_NO_state\models.mat';
% % load(path);
% % sFNC = models.Order10.Data;
% % dFNC = models.Order8.Data;
% % X_raw = sFNC([42:end],:);
% % X_raw = dFNC([42:end],:);
% % X_raw = [sFNC([42:end],:),dFNC([42:end],:)];
% %% ���������Ϊ�����
% load(['E:\PDM_Network\pain_empathy.mat']);
% Y_raw = pain_empathy(1:41,2);%ʱ���2�����˹������
% % Y_raw = pain_empathy(42:end,2);%ʱ���2���˹������
%% ����
% X_raw = [pc1_data_S1_L,pc1_data_S1_R,pc1_data_S2_L,pc1_data_S2_R,pc1_data_acc_L,pc1_data_acc_R,pc1_data_tha_spinal,pc2_data_tha_spinal];%,pc1_data_tha_spinal,pc2_data_tha_spinal
     Y_raw = paintolerance_R;
     X_raw = [S1_FAL,S1_FAR,S1_ADL,S1_ADR,S2_FAL,S2_FAR,S2_ADL,S2_ADR,acc_FAL,acc_FAR,acc_ADL,acc_ADR,spinal_FA,spinal_AD];
%     perm_order = PermuteSort(ones(44,1),1);
%     Y_raw = Y(perm_order,:);
    %% ����㽻����֤
    Group = DivGroup(Y_raw,5); %�������5�۽�����֤
    GroupPath = 'G:\Cortical_spinalcord\��ؽ��\PLSR\test';
    cd G:\Cortical_spinalcord\��ؽ��\PLSR\test;
    DivTrain(Y_raw,X_raw,Group,GroupPath);
    PLSR_result = {};
    Y_pre_all = [];  % Y_pre_all_brain=[];Y_pre_all_brain_spinal=[];
    Y_verify_all = [];  % Y_verify_all_brain = [];Y_verify_all_brain_spinal = [];
    for i = 1:5
%         cd  G:\Cortical_spinalcord\��ؽ��\PLSR\test\spinal;
        load(fullfile(GroupPath,['Group' num2str(i,'%02d') '.mat'])); %��5������ݵ���
        [F,P,BETA,R2_pre,R2,R2_adjust,StablePoint_result,X_train_all,Y_train_all] = PLSRMain(Y_train,X_train,5,100);
        PLSR_result{1,i} = {F,P,BETA,R2_pre,R2,R2_adjust,StablePoint_result};
        save(fullfile(GroupPath,'PLSR_result.mat'),'PLSR_result');
        %% Z-ֵ�� ����ѵ�������б�׼��
        Y_train_mean = mean(Y_train_all);
        Y_train_stdl = std(Y_train_all);
        X_train_mean = mean(X_train_all);
        X_train_std = std(X_train_all);
        Y_test = (Y_test  - Y_train_mean )/Y_train_stdl;
        X_test = (X_test - repmat(X_train_mean,[size(X_test,1) 1])) ...
            ./ repmat(X_train_std,[size(X_test,1) 1]);
        %% ��֤ģ��
        Y_pre = [ones(size(X_test,1),1) X_test(:,StablePoint_result) ] * BETA;
        Y_pre_all = [Y_pre_all;Y_pre];
        Y_verify_all = [Y_verify_all;Y_test];
        [corre,p] = corr(Y_verify_all,Y_pre_all);
        [PRESS,R2_result] = CalculateR2(Y_verify_all,Y_pre_all);
        save(fullfile(GroupPath,['Verify_result.mat']),'Y_verify_all','Y_pre_all','PRESS','R2_result','corre');
        save(fullfile(GroupPath,'corr_p_val'),'p');

    end
%     COUNT_diff(k)=corre_brain-corre_brain_spinal;