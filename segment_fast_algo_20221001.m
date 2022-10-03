img = imread("lena.png");
img = double(rgb2gray(img));

% fast algorithm
th = 25;
R = zeros(size(img));
mean_num_list = zeros(1,2); %A: 1, B: 2

sz = size(img);
r_cnt = 0;
for m = 1:sz(1)
    for n = 1:sz(2)
        row_check = false; col_check = false;
        % row_check = 1 -> check m -1
        % col_check = 1 -> check n -1
        if n ~= 1
            col_check = true;
        end
        if m ~= 1
            row_check = true;
        end
        
        if (row_check == false) && (col_check == false)
            % first element
            r_cnt = r_cnt + 1;
            R(m,n) = r_cnt;
            mean_num_list(r_cnt,:) = [img(m,n), 1];
        else
            % others, compare by rc and cc
            row_combine = false; col_combine = false;
            if row_check
                if abs(img(m,n) - mean_num_list(R(m-1,n),1)) <= th
                    row_combine = true;
                end
            end
            if col_check
                if abs(img(m,n) - mean_num_list(R(m,n-1),1)) <= th
                    col_combine = true;
                end
            end
            % priority: row && col -> row -> col
            if row_combine && col_combine   % top and left are the same region
                % update mean, num and region combine
                % if regions are the same
                if R(m-1,n) == R(m,n-1)
                    
                    new_num = mean_num_list(R(m-1,n),2) + 1;
                    new_mean = (prod(mean_num_list(R(m-1,n),:))+ img(m,n))/ new_num;
                    R(m,n) = R(m-1,n);
                    mean_num_list(R(m-1,n),:) = [new_mean, new_num];
                else
                    new_num = mean_num_list(R(m-1,n),2) + mean_num_list(R(m,n-1),2) + 1;
                    new_mean = (prod(mean_num_list(R(m-1,n),:)) + prod(mean_num_list(R(m,n-1),:)) + img(m,n))/new_num;
                    if R(m-1,n) > R(m,n-1)  % R(m-1,n) combine to R(m,n-1)
                        R(m,n) = R(m,n-1);
                        mean_num_list(R(m,n-1),:) = [new_mean, new_num];
                        R(R == R(m-1,n)) = R(m,n-1);
                    else    % R(m,n-1) combine to R(m-1,n)
                        R(m,n) = R(m-1,n);
                        mean_num_list(R(m-1,n),:) = [new_mean, new_num];
                        R(R == R(m,n-1)) = R(m-1,n);
                    end
                end    
            elseif col_combine
                R(m,n) = R(m,n-1);
                new_num = mean_num_list(R(m,n-1),2) + 1;
                new_mean = (prod(mean_num_list(R(m,n-1),:)) + img(m,n))/ new_num;
                mean_num_list(R(m,n-1),:) = [new_mean, new_num];
            elseif row_combine
                R(m,n) = R(m-1,n);
                new_num = mean_num_list(R(m-1,n),2) + 1;
                new_mean = (prod(mean_num_list(R(m-1,n),:)) + img(m,n))/ new_num;
                mean_num_list(R(m-1,n),:) = [new_mean, new_num];
            else
                % add new region
                r_cnt = r_cnt + 1;
                R(m,n) = r_cnt;
                mean_num_list = cat(1,mean_num_list,[img(m,n), 1]);
            end
        end
    end
end

for i = 1:r_cnt
    R(R == i) = mean_num_list(i,1);
end
R = uint8(R);
imshow(R);