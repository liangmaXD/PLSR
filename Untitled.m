disp('... Permutations ...')
NUM_PERMS=500;
for iter_perm = 1: NUM_PERMS
    
    % Display number of permutations (every 50 permuts)
    if mod(iter_perm,50) == 0, disp(num2str(iter_perm)); end
    
    % Leave X unchanged (no need to permute both X and Y matrices)
    Xp = X; % X is already normalized
    
    % Permute Y by shuffling rows (subjects) within groups
    perm_order = PermuteSort(diagnosis_grouping,Time); 
    Yp = Y0(perm_order,:);

    %====================PLSC_head==============================
    % Normalize permuted Y
    Yp = myPLS_norm(Yp,NUM_GROUPS,diagnosis_grouping,CONST_NORM);
    
    % Cross-covariance matrix between X and permuted Y
    Rp = myPLS_cov(Xp,Yp,CONST_NORM,diagnosis_grouping);
    
    % SVD of Rp
    [Up,Sp,Vp] = svd(Rp,'econ');
    
    % Procrustas transform (to correct for axis rotation/reflection)
    rotatemat = rri_bootprocrust(U, Up);
    Up = Up * Sp * rotatemat; 
    Sp = sqrt(sum(Up.^2)); 
    %====================PLSC_end==============================
    
    % Keep singular values for sample distribution of singular values
%     Sp = diag(Sp')';
    permsamp(:,iter_perm) = Sp';
    
    % 计数，有多少次随机情况的潜在因子>自身结果
    if iter_perm == 1
        sp = (Sp' >= diag(S));
    else
        sp = sp + (Sp' >= diag(S));
    end

end