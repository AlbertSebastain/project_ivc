function x = exponential_golomb_decode(code,k)
%     source = codebook{1};
%     code_book_c = codebook{2};
%     len = numel(source);
%     x = zeros(1,len);
%     for i = 1:len
%         ind = source(code_book_c == code(i));
%         x(i) = ind(1);
%code = code{1,1};
len = numel(code);
x = zeros(1,len);
    for i =  1:len
        x_temp = convertStringsToChars(code(i));
        xr = x_temp(end-k+1:end);
%         if isempty(xr)
%             xr = 0;
%         end
        temp = x_temp(1:end-k);
        if strlength(temp) ~= 1
            ind = strfind(temp,'0');
            n = ind(1)-1;
            m = bin2dec(temp(ind(1)+1:ind(1)+n));
        else
            n = 0;
            m = 0;
        end
        xq = dec2bin(2^n+m-1);
        x_rec = [xq,xr];
        x_rec = bin2dec(x_rec);
        x(i) = x_rec;
    end
end