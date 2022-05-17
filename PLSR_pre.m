% [F,P,BETA,R2_pre,R2,R2_adjust,stable] = PLSR(X,Y,num_strap,boot_score,set)
% �����һ��Ϊ����X���ڶ���Ϊ��Ӧ����(Y); �������X,Y���б�׼��̬������Z-score��
% num_strapΪ����bootstraping�Ĵ���
% boot_scoreΪbootstraping�趨����ֵ����
% setΪ������֤����,������Ĭ��Ϊ��һ��������֤��
% ������ֻ����������ʱ����������bootstrap���飬���������뽫��ֵ��set
% ��������Ǳ�������µ�ģ�����Fֵ��ģ���������Pֵ��
%       stableֻ���ڽ���bootstraping������²�����
%       ����Ϊ����ɸѡ����������λ�á�
% R2Ϊģ�����ЧӦֵR2��R2_adjustΪģ��У�����R2
% R2_preΪ���Ǳ�������µĽ�����֤Ԥ��R2  
% -------------------------------------------------------------------------
% �����ȡ���ú����ķ�ʽΪ��
% ������bootstrapɸѡ������ֻ���н�����֤ѡǱ������
% [F,P,BETA,R2_pre,R2,R2_adjust] = PLSR(X,Y,set)
% setΪ������֤������ֻ����X,YĬ����һ��
%
% [F,P,BETA,R2_pre,R2,R2_adjust,stable] = PLSR(X,Y,num_strap,boot_score,set)
% ����bootstraping���������ȶ��ԣ�setͬ�����Բ����������һ��������֤ѡȡǱ������
% --------------------------------------------------------------------------
% ע�⣺ boot_scoreΪ2ʱ�����Ӧ�����ȶ���PֵΪ 0.05��3Ϊ 0.01
% �����Ὣÿ������ɸѡ�������������ڽ�����ʾ������PLSR��ϵ�������һ���Ĺ���
% ��ͨ����ؾ���R������ֵ�ֽ�õ��������ܳ���ǰ�����޳��������٣�����һ��ֱ�ӽ���
% �����޳����������ʱ˵����Щ���������ô�С����һ�����Żᵼ�������������ʱ����
% �ᱨ����ʾ�����齫boot_score������ֵ���ͣ����߼�¼�ڵڼ��γ�������״����ͨ��
% �Գ�������޸ģ������к������������ý���n��bootstrap�������Щ����λ�á�
function [F,P,BETA,R2_pre,R2,R2_adjust,stable] = PLSR(X,Y,num_strap,boot_score,set)
    if size(X,1) ~= size(Y,1)
        error('X��Y������������ȣ�');
    end
    if nargin == 3          %�˴�narginΪ�����ĸ���
        set = num_strap;   %setΪ������֤����, num_strapΪ����bootstraping�Ĵ���
    elseif nargin < 5
        set = size(Y,1);   % number of subjects
    end
    if set > size(Y,1)
        error('������֤�������ó�������(%d)��',size(Y,1));
    elseif set == 1
        error('��������Ϊ1����ʱ����֤����')
    end
    X0 = X;
    Y0 = Y;
    X = zscore(X);  % ��׼��̬����Zֵ����
    Y = zscore(Y);  % ��׼��̬����Zֵ����
    PRESS = [];
    for k = 1: min(size(X,1) - ceil(size(Y,1) / set) - 1, size(X,2)) % �������Ǳ��������Ŀ;Ǳ��������Ϊѵ��������������������Сֵ
        [~,~,~,~,BETA] = plsregress(X,Y,k);    %����PLSC��betaΪƫ��С���˻ع�ϵ���� * �ľ���
        Y1 = [ones(size(X,1),1) X] * BETA;    % ���ֵ
        n = size(Y,1);
        Y_m = mean(Y);
        for i = 1: size(Y,1)
            if i == 1
                TSS = (Y(i) - Y_m)^2;
                ESS = (Y1(i) - Y_m)^2;
                RSS = (Y(i) - Y1(i))^2;
            else
                TSS = (Y(i) - Y_m)^2 + TSS;     % �����ƽ���� 
                ESS = (Y1(i) - Y_m)^2 + ESS;    % �ع�ƽ����
                RSS = (Y(i) - Y1(i))^2 + RSS;   % ʣ��ƽ����
            end
        end
        clear Y1
        R2(k) =1 - RSS / TSS;       % ģ��R2
        R2_adjust(k) = 1 - (RSS / (n - k -1)) / (TSS / (n - 1));     % ģ��У�����R2
        F(k) = (ESS / k) / (RSS / (n - k -1));      % ģ��Fֵ
        %% ������֤
        if set ~= size(Y,1)
            set_num = floor(size(X,1) / set);      % ÿ����������
            more_mun = mod(size(X,1),set);        %����/����������
            set_num = repmat(set_num,[set 1]);    %set��1�е�set_num
            for i = 1: more_mun
                set_num(i) = set_num(i) + 1;        % ȷ��ÿ���������
            end
        else
            set_num = ones(size(Y,1),1);
        end
        label = 1 : size(Y,1);
        group = CVset(Y,set,set_num);   % ������֤����,���õ�������������
        array = [];
        %% ���н�����֤����Ԥ��R2
        for m = 1: set
            test = group{m};
            i = 1;
            train = [];
            while i <= set
                if i ~= m
                    train  = [train group{i}];
                end
                i = i + 1;
            end
            x1 = zscore(X(train,:));
            y1 = zscore(Y(train));
            x2 = zscore(X(test,:));
            [~,~,~,~,b] = plsregress(X(train,:),Y(train),k);   % ѵ�������ģ��
            b = roundn(b,-6);
            if m == 1
                Y1 = [ones(size(X(test,:),1),1) X(test,:)]*b;           % Ԥ����֤�� 
            else
                Y1 = [Y1;[ones(size(X(test,:),1),1) X(test,:)]*b];      % Ԥ����֤��
            end
            array = [array test];   % ���ڲ���������飬������Ҫ��¼��ţ��Ա��ϱ��Ա�Ŷ���
        end
        %%
        Y_pre(k) = {Y1};
        array_pre(k) = {array};
        Y2 = Y(array);      % Ԥ����ʵ�ʵı�Ŷ�Ӧ
        for i = 1: size(Y,1)
            if i == 1
                ESS = (Y1(i) - Y_m)^2;
                RSS = (Y2(i) - Y1(i))^2;
            else
                ESS = (Y1(i) - Y_m)^2 + ESS;    % �ع�ƽ����
                RSS = (Y2(i) - Y1(i))^2 + RSS;   % ʣ��ƽ����
            end
        end
        R2_pre(k) = 1 - RSS / TSS;       % ģ��Ԥ��R2
        PRESS(k) = RSS;
        if R2_pre(k) < 0
            R2_pre(k) = 0;
        end
    end
    lv = find(PRESS == min(PRESS));     % �ҵ���ѵ�Ǳ��������/���ܻ����ĵط�
    %% �ع�ϵ������(ʹ��bootstraping�����ȶ���)
    if nargin >= 4
    l0 = 1: size(X,2);
    fprintf('��ʱX��ά��Ϊ%d\n',size(X,2))
    stable = bootstraping(Y,X,lv,num_strap,boot_score);  % ����ϵ������õ��ȶ�������λ��
        if length(stable) ~= size(X,2)
            l = l0(stable);
            X0 = X0(:,stable);
            
            if isempty(X0)
                warning('û��ϵ��ͨ�����趨����ֵ����')
                error('����ֹ���򣡣�')
            end
            [F,P,BETA,R2_pre,R2,R2_adjust,stable1] = PLSR(X0,Y0,num_strap,boot_score,set); % �����ݹ飬�Լ������Լ�����ʵ������ѡ��ĵ�������
            stable = l(stable1);
            return
        else
        end
    end
    %%
    fprintf('Ǳ������Ϊ%dʱ��ģ��Ԥ��Ч�����.\n��ʱ��Ԥ��R2Ϊ%d.\n',lv,R2_pre(lv))
    F = F(lv);
    R2_adjust = R2_adjust(lv);
    array_pre = array_pre{lv};
    Y2 = Y_pre{lv};
    nu = size(X,1) - lv - 1;
    P = fpval(F,lv,nu);
    %% ��ͼ
    [~,~,~,~,BETA] = plsregress(X,Y,lv);
    Y1 = [ones(size(X,1),1) X] * BETA;    % ���ֵ
    Y = Y(array_pre);
    Y1 = Y1(array_pre);
    ero = Y - Y1;   % ���ģ�Ͳв�
    figure(1)
    stem(ero,'b')
    ero = Y - Y2;   % ������֤�в�
    hold on
    stem(ero,'r')
    hold off
    title(['Ǳ������Ϊ' num2str(lv) 'ʱ�Ĳв�'])
    legend('Fitted','Cross-validation')
    
    figure(2)
    plot(1:length(R2),R2,'xb')
    hold on
    plot(1:length(R2_pre),R2_pre,'xr')
    plot(lv,0:0.02:R2(lv),'--k')
    plot(0:0.05:lv,R2_pre(lv),'--k')
    axis([1 length(R2_pre)+1 0 1])
    hold off
    title('PLSģ��ѡ��ͼ')
    xlabel('����')
    ylabel('R-sq')
    legend('Fitted','Cross-validation')
    R2_pre = R2_pre(lv);
    R2 = R2(lv);

    figure(3)
    plot(Y,Y1,'xb')
    hold on
    plot(Y,Y2,'xr')
    legend('Fitted value','Cross-validation value')
    Y_max = max(Y);
    Y_min = min(Y);
    y = Y_min:(Y_max-Y_min)/size(Y,1):Y_max;
    x = Y_min:(Y_max-Y_min)/size(Y,1):Y_max;
    plot(x,y,'--g')
    title(['Ǳ������Ϊ' num2str(lv) '����Ӧͼ'])
    ylabel('Calculate response')
    xlabel('Actual response')
    hold off
end

%% ģ������Pֵ
function p = fpval(x,df1,df2)
% Matlab�Դ�˽�к���
if nargin < 3
    error(message('stats:fpval:TooFewInputs')); 
end

xunder = 1./max(0,x);
xunder(isnan(x)) = NaN;
p = fcdf(xunder,df2,df1);
end

%% ������֤����
function group = CVset(Y,set,set_num)
    Y0 = Y;     
    Y = - sort(-Y);    %���Ӵ�С����
    l = [];
    i = 1;
    while i <= length(Y)
        l1 = find(Y(i) == Y0);
        if length(l1) > 1
            l = [l;l1];
            i = i + length(l1);
        else
            l = [l;l1];
            i = i + 1;
        end
    end
    clear l1
    for i = 1: set
        for j = 1: set_num(i)
            l1(j) = l(i + set*(j-1));
        end
        group{i} = l1;
%         clear l l1
    end
end

%% bootstraping �����������ȶ���
function stable = bootstraping(Y,X,lv,num_strap,boot_score)
    num = length(Y);
    X0 = X;
    Y0 = Y;
    X = zscore(X0);
    Y = zscore(Y0);
    [~,~,~,~,b0] = plsregress(X,Y,lv);
    b0 = b0(2:end);
    parfor i =1: num_strap     % ����bootstrap�õ�ÿһ�ε�bootstrapϵ��
%         for j = 1: length(Y0)
%             l(j) = randperm(num,1);
%         end
        l = ceil(num*rand(num,1));
        Y = Y0(l,:);
        X = X0(l,:);
        X = zscore(X);
        Y = zscore(Y);
        [~,~,~,~,b1] = plsregress(X,Y,lv);
        b(i,:) = b1(2:end);
    end
    b_mean = mean(b);
    b_mean2 = mean(b.^2);
    b_std = sqrt(b_mean2 - b_mean.^2);
    b_std = real(b_std);
    clear b b_mean2 b_mean
    score = b0' ./ b_std;
    l = find(~isfinite(score));
    for i = 1: length(l)
        score(l(i)) = b0(l(i));
    end
    stable = find(abs(score) > boot_score);
end
