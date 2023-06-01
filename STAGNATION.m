function [i_s] = STAGNATION(np,vt)
    i_s = 1;
    if vt(1) >= 0
        sign0 = 1;
    else
        sign0 = -1;
    end
    for i=2:np
        if vt(i) >= 0
            sign = 1;
        else
            sign = -1;
        end
        if sign == sign0
            continue;
        else
            i_s = i;
            break;
        end
    end
end



