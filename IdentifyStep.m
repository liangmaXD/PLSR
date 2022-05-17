function [Step,StepSize] = IdentifyStep(PRESS,Step,StepSize)
    if Step(3) > 4
        StepSize = 2^(-20);
        return
    end
    if sum(PRESS == sort(PRESS)) == length(PRESS)
        if (PRESS(1) == PRESS(end) && Step(1) ~= 0)
            StepSize = 2^(-20);
        elseif PRESS(1) == PRESS(end) && Step(1) == 0
            if StepSize > 2^(-4)
                StepSize = StepSize / 2;
                Step = 0:StepSize:StepSize*4;
            else
                StepSize = 2^(-20);
            end
        else
            StepSize = StepSize * 2;
            Step = Step(end) - StepSize:StepSize:Step(end) + 3*StepSize;
        end
        return
    elseif (sum(PRESS == -sort(-PRESS)) == length(PRESS)) ||...
            (PRESS(1) == PRESS(end) && PRESS(1) == max(PRESS))
        if Step(1) == 0
            if StepSize > 2^(-4)
                Step = Step(1):StepSize:Step(1)+StepSize*4;
                StepSize = StepSize / 2;
            else
                StepSize = 2^(-20);
            end
        else 
            StepSize = StepSize / 2;
            Step = Step(1) - StepSize:StepSize:Step(1)+StepSize*3;
        end
        return
    end
    PRESS_sort = -sort(-PRESS);
    i = 1;
    SortNum = [];
    while i <= length(PRESS_sort)
        SortNum = [SortNum;find(PRESS == PRESS_sort(i))];
        i = i + length(find(PRESS == PRESS_sort(i)));
    end
    %     SortNum = SortNum;
    if abs(SortNum(1) - SortNum(2)) <= StepSize*2 ||...
            abs(SortNum(1) - SortNum(2)) > StepSize*4 ||...
            Step(end) == Step(SortNum(1))
        if Step(SortNum(1)) - 2*StepSize >= 0
            Step = Step(SortNum(1)) - StepSize: StepSize:...
                Step(SortNum(1)) + 3*StepSize;
        else
            Step = Step(SortNum(1)): StepSize : Step(SortNum(1)) + 4 * StepSize;
        end
        StepSize = StepSize / 2;
    else
        Step = min(Step(SortNum(1)),Step(SortNum(2))):StepSize:...
            min(Step(SortNum(1)),Step(SortNum(2)))+4*StepSize;
        StepSize = StepSize / 2;
    end
end