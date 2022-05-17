function [stable,score] = bootstraping(Y,X,lv,num_strap,boot_score)
    num1 = length(Y);
    X2 = X;
    Y2 = Y;
    X = zscore(X2);
    Y = zscore(Y2);
    [~,~,~,~,b0] = plsregress(X,Y,lv);
    b0 = b0(2:end);
    parfor j =1: num_strap     % 进行bootstrap得到每一次的bootstrap系数
%         for j = 1: length(Y0)
%             l(j) = randperm(num,1);
%         end
        l1 = ceil(num1*rand(num1,1));
        Y = Y2(l1,:);
        X = X2(l1,:);
        X = zscore(X);
        Y = zscore(Y);
        [~,~,~,~,b1] = plsregress(X,Y,lv);
        b2(j,:) = b1(2:end);
    end
    b_mean = mean(b2);
    b_mean2 = mean(b2.^2);
    b_std = sqrt(b_mean2 - b_mean.^2);
    b_std = real(b_std);
    clear b b_mean2 b_mean
    score = b0' ./ b_std;
    l1 = find(~isfinite(score));
    for j = 1: length(l1)
        score(l1(j)) = b0(l1(j));
    end
    stable = find(abs(score) > boot_score);
end