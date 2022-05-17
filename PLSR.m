% [F,P,BETA,R2_pre,R2,R2_adjust,stable] = PLSR(X,Y,num_strap,boot_score,set)
% 输入第一项为特征X，第二项为响应变量(Y); 函数会对X,Y进行标准正态化处理（Z-score）
% num_strap为进行bootstraping的次数
% boot_score为bootstraping设定的阈值分数，为二值矩阵。eg.[1 3];
% set为交叉验证折数,不输入默认为留一法交叉验证。
% 当函数只有三个输入时：将不进行bootstrap检验，第三个输入将赋值给set
% 输出：最佳潜在因子下的模型拟合F值，模型拟合显著P值；
%       stable只能在进行bootstraping的情况下才能输
%       出，为最终筛选出的特征的位置。
% R2为模型拟合效应值R2，R2_adjust为模型校正后的R2
% R2_pre为最佳潜在因子下的交叉验证预测R2
% --------------------------------------------------------------------------
% 建议采取调用函数的方式为：
% 不进行bootstrap筛选特征，只进行交叉验证选潜在因子
% [F,P,BETA,R2_pre,R2,R2_adjust] = PLSR(X,Y,set)
% set为交叉验证折数，只输入X,Y默认留一法
% 
% [F,P,BETA,R2_pre,R2,R2_adjust,stable] = PLSR(X,Y,num_strap,boot_score,set)
% 采用bootstraping测试特征稳定性，set同样可以不输入采用留一法交叉验证选取潜在因子
% --------------------------------------------------------------------------
% 注意： boot_score为2时，相对应特征稳定性P值为 0.05；3为 0.01
% --------------------------------------------------------------------------
% 程序改进2018-9-27
% 函数会将每次特征筛选后的数量在命令窗口进行提示。相较于之前按固定阈值进行筛选进行
% 了改进，修改为动态阈值。从第一个（低）阈值开始依次增加0.1，直到递增到第二个（高）
% 阈值且在高阈值下所有特征系数稳定（大于高阈值），程序进行结果输出。
% --------------------------------------------------------------------------
% 改进理由：PLSR的特征之间存在影响，使得在噪声较多的情况下影响目标特征的稳定性。导
% 致始终在同一阈值下工作时，若阈值较低，不能有效的去除噪声；若阈值较高，使得刚开始
% 在进行特征筛选时，高阈值将一些真实有效的特征剔除。
function [F,P,BETA,R2_pre,R2,R2_adjust,stable] = PLSR(X,Y,num_strap,boot_score,set)
    if size(X,1) ~= size(Y,1)
        error('X与Y的行数必须相等！');
    end
    if nargin == 3
        set = num_strap;
    elseif nargin < 5
        set = size(Y,1);
    end
    if set > size(Y,1)
        error('交叉验证折数不得超过人数(%d)！',size(Y,1));
    elseif set == 1
        error('折数不得为1，此时无验证集！')
    end
    X0 = X;
    Y0 = Y;
    X = zscore(X);  % 标准正态化（Z值化）
    Y = zscore(Y);  % 标准正态化（Z值化）
    PRESS = [];
    for k = 1: min(min(size(X,1) - ceil(size(Y,1) / set) - 1, size(X,2)),15)% 计算最大潜在因子数目
%         fprintf('潜在因子为%d时\n',k)
        [~,~,~,~,BETA] = plsregress(X,Y,k);
        Y1 = [ones(size(X,1),1) X] * BETA;    % 拟合值
        n = size(Y,1);
        Y_m = mean(Y);
        for i = 1: size(Y,1)
            if i == 1
                TSS = (Y(i) - Y_m)^2;
                ESS = (Y1(i) - Y_m)^2;
                RSS = (Y(i) - Y1(i))^2;
            else
                TSS = (Y(i) - Y_m)^2 + TSS;     % 总离差平方和 
                ESS = (Y1(i) - Y_m)^2 + ESS;    % 回归平方和
                RSS = (Y(i) - Y1(i))^2 + RSS;   % 剩余平方和
            end
        end
        clear Y1
        R2(k) =1 - RSS / TSS;       % 模型R2
        R2_adjust(k) = 1 - (RSS / (n - k -1)) / (TSS / (n - 1));     % 模型校正后的R2
        F(k) = (ESS / k) / (RSS / (n - k -1));      % 模型F值
        %% 交叉验证
        if set ~= size(Y,1)
            set_num = floor(size(X,1) / set);      % 每组最少人数
            more_mun = mod(size(X,1),set);
            set_num = repmat(set_num,[set 1]);
            for i = 1: more_mun
                set_num(i) = set_num(i) + 1;        % 确定每组具体人数
            end
        else
            set_num = ones(size(Y,1),1);
        end
        label = 1 : size(Y,1);
        group = CVset(Y,set,set_num);   % 交叉验证分组,调用第三个函数代码
        array = [];
        %% 进行交叉验证计算预测R2
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
%             x1 = zscore(X(train,:));
%             y1 = zscore(Y(train));
%             x2 = zscore(X(test,:));
            [~,~,~,~,b] = plsregress(X(train,:),Y(train),k);   % 训练集拟合模型
%             b = roundn(b,-6);%这里进行了改动
            if m == 1
                Y1 = [ones(size(X(test,:),1),1) X(test,:)]*b;           % 预测验证集 
            else
                Y1 = [Y1;[ones(size(X(test,:),1),1) X(test,:)]*b];      % 预测验证集
            end
            array = [array test];   % 由于采用随机分组，所以需要记录编号，以保障被试编号对其
        end
        %%
        Y_pre(k) = {Y1};
        array_pre(k) = {array};
        Y2 = Y(array);      % 预测与实际的编号对应
        for i = 1: size(Y,1)
            if i == 1
                ESS = (Y1(i) - Y_m)^2;
                RSS = (Y2(i) - Y1(i))^2;
            else
                ESS = (Y1(i) - Y_m)^2 + ESS;    % 回归平方和
                RSS = (Y2(i) - Y1(i))^2 + RSS;   % 剩余平方和
            end
        end
        R2_pre(k) = 1 - RSS / TSS;       % 模型预测R2
        for i = 1: length(R2_pre)
            R2_pre_adjust(i) = R2_pre(i) / ((log(i^i) / 100) + 1);
        end
        PRESS(k) = RSS;
        if R2_pre(k) < 0
            R2_pre(k) = 0;
        end
    end
    lv = find(R2_pre_adjust == max(R2_pre_adjust));     % 找到最佳的潜在因子数/可能会出错的地方
    lv = lv(1);
    %% 回归系数检验(使用bootstraping检验稳定性)
    if nargin >= 4
    l0 = 1: size(X,2);
    fprintf('此时X的维度为%d\n',size(X,2))
    stable = bootstraping(Y,X,lv,num_strap,boot_score(1));  % 进行系数检验得到稳定的特征位置
        if length(stable) ~= size(X,2) || boot_score(1) ~= boot_score(2)
            l = l0(stable);
            X0 = X0(:,stable);
            if isempty(X0);
                warning('没有系数通过所设定的阈值！！')
                error('已终止程序！！')
            end
            if boot_score(1) < boot_score(2)
                boot_score(1) = boot_score(1) + 0.1;
            else
                boot_score(1) = boot_score(2);
            end
            [F,P,BETA,R2_pre,R2,R2_adjust,stable1] = PLSR(X0,Y0,num_strap,boot_score,set); % 函数递归，自己调用自己，以实现特征选择的迭代过程
            stable = l(stable1);
            return
        else
        end
    end
    %%
    fprintf('潜在因子为%d时，模型预测效果最佳.\n此时的预测R2为%d.\n',lv,R2_pre(lv))
    F = F(lv);
    R2_adjust = R2_adjust(lv);
    array_pre = array_pre{lv};
    Y2 = Y_pre{lv};
    nu = size(X,1) - lv - 1;
    P = fpval(F,lv,nu);
    
    %% 作图
    [~,~,~,~,BETA] = plsregress(X,Y,lv);
    Y1 = [ones(size(X,1),1) X] * BETA;    % 拟合值
    Y = Y(array_pre);
    Y1 = Y1(array_pre);
    ero = Y - Y1;   % 最佳模型残差
    figure(1)
    stem(ero,'b')
    ero = Y - Y2;   % 交叉验证残差
    hold on
    stem(ero,'r')
    hold off
    title(['潜在因子为' num2str(lv) '时的残差'])
    legend('Fitted','Cross-validation')
    
    figure(2)
    plot(1:length(R2),R2,'xb')
    hold on
    plot(1:length(R2_pre),R2_pre,'xr')
    plot(lv,0:0.02:R2(lv),'--k')
    plot(0:0.05:lv,R2_pre(lv),'--k')
    axis([1 length(R2_pre)+1 0 1])
    hold off
    title('PLS模型选择图')
    xlabel('分量')
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
    title(['潜在因子为' num2str(lv) '的响应图'])
    ylabel('Calculate response')
    xlabel('Actual response')
    hold off
end

%% 模型显著P值
function p = fpval(x,df1,df2)
% Matlab自带私有函数
if nargin < 3
    error(message('stats:fpval:TooFewInputs')); 
end

xunder = 1./max(0,x);
xunder(isnan(x)) = NaN;
p = fcdf(xunder,df2,df1);
end

%% 交叉验证分组
function group = CVset(Y,set,set_num)
    Y0 = Y;
    Y = - sort(-Y);
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
        clear l1
    end
end

%% bootstraping 检验特征的稳定性
function stable = bootstraping(Y,X,lv,num_strap,boot_score)
    num = length(Y);
    X0 = X;
    Y0 = Y;
    X = zscore(X0);
    Y = zscore(Y0);
    [~,~,~,~,b0] = plsregress(X,Y,lv);
    b0 = b0(2:end);
    parfor i =1: num_strap     % 进行bootstrap得到每一次的bootstrap系数
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
    if isempty(stable) %|| length(stable) / length(score) > 0.99
        stable = 1 : length(score);
        if isempty(stable)
        fprintf(['没有通过所设定阈值的特征，此时最稳定特征的分数为%.2f\n'...
            '输出前一次筛选出稳定特征得到的结果。\n'],max(abs(score)))
        end
    end
end
