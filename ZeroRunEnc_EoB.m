function zze = ZeroRunEnc_EoB(zz,EOB)
%EOB = 1000;
%  Input         : zz (Zig-zag scanned sequence, 1xN)
%                  EOB (End Of Block symbol, scalar)
%
%  Output        : zze (zero-run-level encoded sequence, 1xM)
len = length(zz);
zze = [];
for rec = 1:len/64
    if rec ~= len/64
        zz_t = zz((rec-1)*64+1:rec*64);
    else
        zz_t = zz((rec-1)*64+1:len);
    end

[lens] = length(zz_t);
k = 1;
zze_t = [];
flag = 0;
m = 0;
while k <= lens
    if zz_t(k) == 0
        count = -1;
        m = k;
        while (zz_t(m) == 0 )
            if m == 1
                flag1 = 1;
            end
            count = count+1;
            m = m+1;
            if m > lens
                flag = 1;
                break
            end
            
        end
        zze_t = [zze_t,0,count];
        k = m-1;
    else
        zze_t = [zze_t,zz_t(k)];
    end
    k = k+1;
end
len2 = length(zze_t);
if flag == 1
    zze_t = zze_t(1:len2-1);
    zze_t(end) = EOB;
end
zze = [zze,zze_t];
end
end