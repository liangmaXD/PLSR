function [optimalLVs,PRESS,R2] = GetOptimalLVs(Y_train,X_train,Y_test,X_test)
    for i = 1: min(min(size(Y_train,1) - 1, size(X_train,2)),20)
        [~,~,~,~,b] = plsregress(X_train,Y_train,i);
        Y_test_pre(:,i) = [ones(size(X_test,1),1) X_test] * b;
        Y_train_pre(:,i) = [ones(size(X_train,1),1) X_train] * b;
        Beta(i,:) = b;
    end
    [~,R2_pre] = CalculateR2(Y_test,Y_test_pre);
    [~,R2_match] = CalculateR2(Y_train,Y_train_pre);
    l1 = R2_pre < 0;
    l2 = R2_match < 0;
    l = l1 & l2;
    R2 = R2_pre .* abs(R2_pre) .* R2_match; %
    R2(l) = - R2(l);
    clear l1 l2 l
%     for i = 1: length(R2)
%         R2_adjust(i) = R2(i) / ((log(i^i) / 100) + 1);
%     end
%     if ~isempty(find(R2>0))
%         optimalLVs = find(R2 == max(R2));
%         optimalLVs = optimalLVs(1);
%     else
%         optimalLVs = find(PRESS == min(PRESS));
%         PRESS = min(PRESS);
%         PRESS  = PRESS(1);
%     end
    optimalLVs = find(R2 == max(R2));
    optimalLVs = optimalLVs(1);
%     optimalLVs = find(PRESS == min(PRESS));
    PRESS = max(R2);
    R2 = R2_pre(optimalLVs);
end