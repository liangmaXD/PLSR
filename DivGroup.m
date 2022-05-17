function group = DivGroup(Y,setnum,Times)
    if setnum < size(Y,1)
        set_num = floor(size(Y,1) / setnum);      % 每组最少人数
        more_mun = mod(size(Y,1),setnum);
        set_num = repmat(set_num,[setnum 1]);
        for i = 1: more_mun
            set_num(i) = set_num(i) + 1;        % 确定每组具体人数
        end
    else
        set_num = ones(size(Y,1),1);
    end
    Y0 = Y;
    Y = - sort(-Y);
    %     fprintf('%d\n',nargin);
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
    if nargin < 3
        for i = 1: setnum
            for j = 1: set_num(i)
                l1(j) = l(i + setnum*(j-1));
            end
            group{i} = l1;
            clear l1
        end
    else
        for i = 1: ceil(size(Y,1) / setnum)
            if length(l) > setnum
                store{i} = l(1:setnum);
                l(1:setnum) = [];
            else
                store{i} = l;
            end
        end
        if Times > setnum^(length(store) - 1) * length(store{end})
            Times = setnum^(length(store) - 1) * length(store{end});
        end
        for i = 1: Times
            l = [];
            for j = 1: length(store)
                l = randperm(length(store{j}),1);
                group_num = store{j};
                Sub(j) = group_num(l);
            end
            group{i} = Sub;
        end
    end
end

% %% 留一法
% function group = DivGroup(Y,setnum,Times)
%     if setnum < size(Y,1)
%         %set_num = floor(size(Y,1) / setnum);      % 每组最少人数
%         set_num =1;
%         more_mun = mod(size(Y,1),setnum);
%         set_num = repmat(set_num,[setnum 1]);
%         for i = 1: more_mun
%             set_num(i) = set_num(i) + 1;        % 确定每组具体人数
%         end
%     else
%         set_num = ones(size(Y,1),1);
%     end
% %      Y0 = Y;
% %      Y = - sort(-Y);
%     %     fprintf('%d\n',nargin);
%     l = [];
%     i = 1;
%     while i <= length(Y)
% %         l1 = find(Y(i) == Y0);
%          l1 = i;
% %         if length(l1) > 1
% %             l = [l;l1];
% %             i = i + length(l1);
% %         else
%             l = [l;l1];
%             i = i + 1;
%         end
%     end
%     clear l1
%     if nargin < 3
%         for i = 1: setnum
%             for j = 1: set_num(i)
%                 l1(j) = l(i + setnum*(j-1));
%             end
%             group{i} = l1;
%             clear l1
%         end
%     else
%         for i = 1: ceil(size(Y,1) / setnum)
%             if length(l) > setnum
%                 store{i} = l(1:setnum);
%                 l(1:setnum) = [];
%             else
%                 store{i} = l;
%             end
%         end
%         if Times > setnum^(length(store) - 1) * length(store{end})
%             Times = setnum^(length(store) - 1) * length(store{end});
%         end
%         for i = 1: Times
%             l = [];
%             for j = 1: length(store)
%                 l = randperm(length(store{j}),1);
%                 group_num = store{j};
%                 Sub(j) = group_num(l);
%             end
%             group{i} = Sub;
%         end
%     end
% end