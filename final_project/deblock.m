function new_im = deblock(img,level)
    alpha = 10;
    beta =11;
    [len,col,dim] = size(img);
    dl = len/8;
    dc = col/8;
    new_im = img;
    wi = [0.1,0.4,0.4,0.1];
    for d = 1:dim
        for i = 1:dl-1
            for j = 1:dc
                block_up = new_im((i-1)*8+1:i*8,(j-1)*8+1:j*8,d);
                block_down = new_im(i*8+1:(i+1)*8,(j-1)*8+1:j*8,d);
                for ind = 1:8
                    p0 = block_up(8,ind);
                    q0 = block_down(1,ind);
                    p1 = block_up(7,ind);
                    q1 = block_up(2,ind);
                    p2 = block_up(6,ind);
                    p3 = block_up(5,ind);
                    q2 = block_down(3,ind);
                    q3 = block_down(3,ind);
                    if (abs(p0-q0) <= (alpha/4)+2 && abs(p1-q1) <= beta && level == 2)
                        new_im(i*8,(j-1)*8+ind,d) = (p2+2*p1+2*p0+2*q0+q1+4)/8;
                        new_im(i*8-1,(j-1)*8+ind,d) = (p2+p1+p0+q0)/4;
                        new_im(i*8-2,(j-1)*8+ind,d) = (2*p3+3*p2+p1+p0+q0+4)/8;
                        new_im(i*8+1,(j-1)*8+ind,d) = (p1+2*p0+2*q0+2*q1+q2+4)/8;
                        new_im(i*8+2,(j-1)*8+ind,d) = (p0+q0+q1+p2+2)/4;
                        new_im(i*8+3,(j-1)*8+ind,d) = (2*q3+3*q2+q1+q0+p0+4)/8;
                    elseif (abs(p0-q0) <= alpha && abs(p1-q1) <= beta)
                        
                        vector = [p1;p0;q0;q1];
                        filters = wi*vector;
                        %delta = ((q-p)*4+(p2-q2)+4)/8;
                        %delta_p = max(-2,delta);
                        %delta_p = min(2,delta_p);
                        new_im(i*8,(j-1)*8+ind,d) = filters;
                        new_im(i*8+1,(j-1)*8+ind,d) = filters;
                    end
                end
            end
        end
    end
    for d = 1:dim
        for i = 1:dl
            for j = 1:dc-1
                block_left = img((i-1)*8+1:i*8,(j-1)*8+1:j*8,d);
                block_right = img((i-1)*8+1:i*8,j*8+1:(j+1)*8,d);
                for ind = 1:8
                    p0 = block_left(ind,8);
                    q0 = block_right(ind,1);
                    p1 = block_left(ind,7);
                    q1 = block_right(ind,2);
                    p2 = block_left(ind,6);
                    p3 = block_left(ind,5);
                    q2 = block_right(ind,3);
                    q3 = block_right(ind,4);
                     if (abs(p0-q0) <= (alpha/4)+2 && abs(p1-q1) <= beta && level == 2)
                         new_im((i-1)*8+ind,j*8,d) = (p2+2*p1+2*p0+2*q0+q1+4)/8;
                         new_im((i-1)*8+ind,j*8-1,d) = (p2+p1+p0+q0+2)/4;
                         new_im((i-1)*8+ind,j*8-2,d) = (2*p3+3*p2+p1+p0+q0+4)/8;
                         new_im((i-1)*8+ind,j*8+1,d) = (p1+2*p0+2*q0+2*q1+q2+4)/8;
                         new_im((i-1)*8+ind,j*8+2,d) = (p0+q0+q1+q2+2)/4;
                         new_im((i-1)*8+ind,j*8+3,d) = (2*p3+3*p2+p1+p0+q0+4)/8;
                     elseif (abs(p0-q0) <= alpha && abs(p1-q1) <= beta)
                        vector = [p2;p0;q0;q2];
                        filters = wi*vector;
                        %delta = ((q-p)*4+(p2-q2)+4)/8;
                        %delta_p = max(-2,delta);
                        %delta_p = min(2,delta_p);
                        new_im((i-1)*8+ind,j*8,d) = filters;
                        new_im((i-1)*8+ind,j*8+1,d) = filters;
                    end
                end
            end
        end
    end
end