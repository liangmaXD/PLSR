function stable = PLSRbootstraping(Y,X,lv,num_strap,boot_score)
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
    if isempty(stable)
        stable = 1 : length(score);
    end
end