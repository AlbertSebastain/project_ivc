function code = exponential_golomb_code(k,x)
    k = floor(log2(mode(x)));
    len = length(x);
    code = strings(1,len);
    for i = 1:len
        temp_str = dec2bin(x(i));
        y = dec2bin(bin2dec(temp_str(1:end-k))+1);
        len_y = strlength(y);
        % uniary code
        u_y = join(string([ones(1,len_y-1),0]),'');
        append = dec2bin(mod(x(i),(2^k)),k);
        temp_code = [u_y,y(2:end),append];
        code(i) = join(temp_code,'');
        
    end
    code = join(code,'');
end
    