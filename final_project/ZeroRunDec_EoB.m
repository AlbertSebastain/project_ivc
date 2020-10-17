function dst = ZeroRunDec_EoB(src,EOB)
%EoB = 1000;
%  Function Name : ZeroRunDec1.m zero run level decoder
%  Input         : src (zero run encoded sequence 1xM with EoB signs)
%                  EoB (end of block sign)
%
%  Output        : dst (reconstructed zig-zag scanned sequence 1xN)
k = 0;
len_src = length(src);
m = 1;
dst = [];
while k < len_src
    k = k+1;
    if m == 65
        m = 1;
    end
    if src(k) == EOB
        dst(end+1:end+64-m+1) = 0;
        m = 1;
        continue
    end
    if src(k) == 0
        dst(end+1:end+src(k+1)+1) = 0;
        m = m+src(k+1)+1;
        k = k+1;
       
        continue
    end
    dst(end+1) = src(k);
    m = m+1;
end
size(dst);
size(src);
end