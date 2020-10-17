function [code,len_code,k] = exponential_golomb_code(x,k)
    if nargin  == 1
        k = ceil(log2(mode(x)));
        if k == -Inf
            k = 0;
        end
    end
    len = length(x);
    code = strings(1,len);
    len_code = 0;
    code_book = cell(1,2);
    code_book_c = strings(len,1);
    for i = 1:len
        temp = x(i);
        y = bitshift(temp,-k)+1;
        len_y = strlength(dec2bin(y));
        u_y = dec2bin(sum(bitset(0,2:len_y)));
        s_y = dec2bin(y);
        s_y = s_y(2:end);
        append = dec2bin(mod(x(i),(2^k)),k);
        temp_code = [u_y,s_y,append];
        code(i) = join(temp_code,'');
        len_code = len_code+strlength(code(i));
        code_book_c(i) = code(i);
    end
    %code_book{1} = x;
    %code_book{2} = code_book_c;
    %code = join(code,'');
end
        