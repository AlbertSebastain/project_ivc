function [code] = CAVLC_encode(source)
    nc = 1;
    %source = zeros(4,4,3);
    %b = [6,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0];
    %source(:,:,1) = b;
    
    load('CAVLC_table.mat');
    x_zig1 = ZigZag4x4(source);
    code_final = strings(1,3);
    for dim = 1:3
    x_zig = x_zig1(:,dim);
    if sum(x_zig == 0) == 16
        code = 1;
        return
    end
    total_coeffcient = sum(x_zig ~= 0);
    ind = find(x_zig ~= 0);
    trail_ones = 0;
    ind_trail = 0;
    sign_trail = 0;
    total_zeros = sum(x_zig(1:ind(end)) == 0);
    for i = ind(end):-1:1
        if x_zig(i) == 0
            continue;
        elseif x_zig(i) ~= 1 && x_zig(i) ~= -1
            break;
        else
            trail_ones = trail_ones+1;
            ind_trail(trail_ones) = i;
            sign_trail(trail_ones) = max(0,x_zig(i));
        end
    end
    trail_ones = min(trail_ones,3);
    if total_coeffcient > 10 && trail_ones < 3
        suffixLength = 1;
    else
        suffixLength = 0;
    end
    ind_trail = ind_trail(1:trail_ones);
    sign_trail = sign_trail(1:trail_ones);
    sign_trail = ~sign_trail;
    first_ind = find(index_total_coe_trail_ones(:,1) == trail_ones);
    ind_table1 = first_ind(find(index_total_coe_trail_ones(first_ind,2) == total_coeffcient));
    code_1 = total_coe_trail_ones(ind_table1,nc);
    code_2 = "";
    if ~isempty(sign_trail)
        code_2 = join(string(double(sign_trail)),'');
    end
    code = join([code_1,code_2],'');
    level_without_trail = setdiff(ind,ind_trail);
    for i = size(level_without_trail):-1:1
        ind_level = level_without_trail(i);
        if x_zig(ind_level) > 0
            level = bitshift(x_zig(ind_level),1)-2;
        else
            level = bitshift(x_zig(ind_level),1)-1;
        end
        level_prefix = level/(bitshift(1,suffixLength));
        level_suffix = mod(level,bitshift(1,suffixLength));
        code_3 = join(string([zeros(1,level_prefix),1]),'');
        code = join([code,code_3],'');
        code_4 = dec2bin(level_suffix,suffixLength);
        if suffixLength == 0
            suffixLength = suffixLength+1;
            code_4 = "";
        elseif level > (bitshift(3,suffixLength)-1) && suffixLength < 6
            suffixLength = suffixLength+1;
        end
        code = join([code,code_4],'');
    end
    code_5 = total_zeros_total_coeffic(total_zeros+1,total_coeffcient);
    code = join([code,code_5],'');
    zero_left = total_zeros;
    for i = size(ind):-1:2
        ind1 = ind(i);
        run_before = ind1-ind(i-1)-1;
        zero_left_ind = min(zero_left,7);
        code_6 = table_run_before(run_before+1,zero_left_ind);
        code = join([code,code_6],'');
        zero_left = zero_left-run_before;
    end
    code_final(dim) = code;
    end
end