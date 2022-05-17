function Result = PermuteSort(diagnosis_grouping,Time)
    NumGroup = length(unique(diagnosis_grouping));
    Num = length(diagnosis_grouping) / Time;
    Mid = randperm(Num);
    for i = 1: NumGroup
        GroupNum(i) = length(find(diagnosis_grouping == i)) / Time;
        MidResult{i} = Mid(1: GroupNum(i));
        Mid(1: GroupNum(i)) = [];
    end
    for i = 1: length(MidResult)
        Group_Mid = MidResult{i};
        for k = 1: NumGroup-1
            l = find(Group_Mid > sum(GroupNum(1:k)));
            Group_Mid(l) = Group_Mid(l) + (GroupNum(k))*(Time-1);
        end
        MidResult{i} = Group_Mid';
    end
    Result = [];
    for k = 1: NumGroup
        Result = [Result;MidResult{k}];
        for i = 1: Time - 1
            Group_Mid = MidResult{k};
            l = [];
            for j = 1: NumGroup
                l = find( (Group_Mid <= sum(GroupNum(1:j) .* Time)) &...
                    (Group_Mid > sum(GroupNum(1:j-1).* Time)));
                Result = [Result;Group_Mid(l) + GroupNum(j)*i];
            end
        end
    end
end