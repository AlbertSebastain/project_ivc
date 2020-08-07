% [len,col] = size(total_zeors_total_coeffic);
% for i = 1:len
%     for j = 1:col
%         if total_zeors_total_coeffic(i,j) ~= "";
%             total_zeors_total_coeffic(i,j) = strrep(total_zeors_total_coeffic(i,j),' ','');
%         end
%     end
% end
n = 4;
for i = 12:62
    index_total_coe_trail_ones(i,2) = n;
    if index_total_coe_trail_ones(i,1) == 0
        n = n+1;
    end
end