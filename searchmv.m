function mv = searchmv(search,block)
best_ssd = inf;
    for i = 1:9
        for j = 1:9
            search_block = search(i:i+7,j:j+7);
            ssd_current = sum((search_block-block).^2,'all');
            if ssd_current<best_ssd
                best_ssd = ssd_current;
                mv = [i,j]-[5,5];
            end
        end
    end
end