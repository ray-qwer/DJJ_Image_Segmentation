img = imread("lena.png");
img = double(rgb2gray(img));
tic
% fast algorithm
th = 25;
min_pixel = ceil(512*512/1000);
R = zeros(size(img));
mean_num_list = zeros(1,2); %A: 1, B: 2
list_init = true;

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
        flag = 0; ref_point = 0; replace_point = 0;
        % flag =0 -> new region; flag =1 -> add a point to a region;
        % flag =2 -> combine one region to another
        if row_combine && col_combine   % top and left are the same region
            % update mean, num and region combine
            % if regions are the same
            if R(m-1,n) == R(m,n-1)
                flag = 1; ref_point = R(m-1,n);
            else
                if R(m-1,n) > R(m,n-1)  % R(m-1,n) combine to R(m,n-1)
                    flag = 2; ref_point = R(m,n-1); replace_point = R(m-1,n);
                else    % R(m,n-1) combine to R(m-1,n)
                    flag = 2; ref_point = R(m-1,n); replace_point = R(m,n-1);
                end
            end    
        elseif col_combine
            flag = 1; ref_point = R(m,n-1);
        elseif row_combine
            flag = 1; ref_point = R(m-1,n);
        else
        end
        switch flag
            case 0  % add new region
                r_cnt = r_cnt + 1;
                R(m,n) = r_cnt;
                if list_init
                    mean_num_list(r_cnt,:) = [img(m,n),1];
                    list_init = false;
                else
                    mean_num_list = cat(1,mean_num_list,[img(m,n),1]);
                end
            case 1  % add one new point
                R(m,n) = ref_point;
                new_num = mean_num_list(ref_point,2) + 1;
                new_mean = (prod(mean_num_list(ref_point,:)) + img(m,n))/ new_num;
                mean_num_list(ref_point,:) = [new_mean, new_num];
            case 2  % combine one region to another
                R(m,n) = ref_point;
                new_num = mean_num_list(ref_point,2) + mean_num_list(replace_point,2) + 1;
                new_mean = (prod(mean_num_list(ref_point,:)) + prod(mean_num_list(replace_point,:)) + img(m,n)) / new_num;
                mean_num_list(ref_point,:) = [new_mean, new_num];
                mean_num_list(replace_point,:) = [0, 0];
                R(R == replace_point) = ref_point;
            otherwise
                disp("error");
        end
        
    end
end
toc
t = mean_num_list(:,2);
index = find(t<min_pixel & t >0);
neighbor_matrix = [0,1,0;1,0,1;0,1,0];
tic
for i = 1:length(index)
    M = zeros(sz);
    di_M = zeros(sz);
    [r,c] = find(R == index(i));
    M((c-1)*sz(1)+r) = 1;
    di_M(conv2(M,neighbor_matrix,'same') > 0) = 1;
    [n_r,n_c] =  find((di_M - M) == 1); % get neighbor
    u = unique(R((n_c-1)*sz(1)+n_r));
    if ~isempty(u(mean_num_list(u,2) > min_pixel))
        u = u(mean_num_list(u,2) > min_pixel);
    end
    [min_dist,u_index]=min((abs(mean_num_list(index(i),1) - mean_num_list(u,1))./sqrt(mean_num_list(u,2))));
%     [min_dist,u_index]=min((abs(mean_num_list(index(i),1) - mean_num_list(u,1))));
    % update mean, num, region
    region_index = u(u_index(1));
    new_num = mean_num_list(region_index,2) + mean_num_list(index(i),2);
    mean_num_list(region_index,2) = new_num;
    R(R == index(i)) = region_index;
    mean_num_list(index(i),:) = [0, 0];
end
toc
for i = 1:r_cnt
    R(R == i) = mean_num_list(i,1);
end
R = uint8(R);
imshow(R);