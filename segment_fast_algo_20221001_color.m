img = double(imread("lena_color.png"));

tic
% fast algorithm
th = 25;
sz = size(img);
min_pixel = ceil(512*512/1000);
R = zeros(sz(1:2));
m_n_list = zeros(0,4); %A: r,g,b, B: 4
list_init = true;


r_cnt = 0;
for m = 1:sz(1)
    for n = 1:sz(2)
        row_combine = false; col_combine = false;
        % compare by rc and cc
        pixel = squeeze(img(m,n,:))';
        if n ~= 1
            if mean(abs(pixel - m_n_list(R(m,n-1),1:3))) <= th
                col_combine = true;
            end
        end
        if m ~= 1
            if mean(abs(pixel - m_n_list(R(m-1,n),1:3))) <= th
                row_combine = true;
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
                
                m_n_list = cat(1,m_n_list,[pixel,1]);
                
            case 1  % add one new point
                R(m,n) = ref_point;
                new_num = m_n_list(ref_point,4) + 1;
                new_mean = (m_n_list(ref_point,1:3).* m_n_list(ref_point,4) + pixel)./ new_num;
                m_n_list(ref_point,:) = [new_mean, new_num];
            case 2  % combine one region to another
                R(m,n) = ref_point;
                new_num = m_n_list(ref_point,4) + m_n_list(replace_point,4) + 1;
                new_mean = (m_n_list(ref_point,1:3).*m_n_list(ref_point,4) + m_n_list(replace_point,1:3).*m_n_list(replace_point,4) + pixel) / new_num;
                m_n_list(ref_point,:) = [new_mean, new_num];
                m_n_list(replace_point,:) = zeros(1,4);
                R(R == replace_point) = ref_point;
            otherwise
                disp("error");
        end
        
    end
end
toc
t = m_n_list(:,4);
index = find(t<min_pixel & t >0);
neighbor_matrix = [0,1,0;1,0,1;0,1,0];
tic
for i = 1:length(index)
    M = zeros(sz(1:2));
    di_M = zeros(sz);
    [r,c] = find(R == index(i));
    M((c-1)*sz(1)+r) = 1;
    di_M(conv2(M,neighbor_matrix,'same') > 0) = 1;
    [n_r,n_c] =  find((di_M - M) == 1); % get neighbor
    u = unique(R((n_c-1)*sz(1)+n_r));
    if ~isempty(u(m_n_list(u,4) > min_pixel))
        u = u(m_n_list(u,4) > min_pixel);
    end
%     [min_dist,u_index]=min((abs(mean_num_list(index(i),1) - mean_num_list(u,1))./sqrt(mean_num_list(u,2))));
    [min_dist,u_index]=min(mean(abs(m_n_list(index(i),1:3) - m_n_list(u,1:3)),2));
    % update mean, num, region
    region_index = u(u_index(1));
    new_num = m_n_list(region_index,4) + m_n_list(index(i),4);
    m_n_list(region_index,4) = new_num;
    R(R == index(i)) = region_index;
    m_n_list(index(i),:) = zeros(1,4);
end
toc
red = zeros(size(R));
blue = zeros(size(R));
green = zeros(size(R));
for i = 1:r_cnt
    red(R == i) = m_n_list(i,1);
    blue(R == i) = m_n_list(i,2);
    green(R == i) = m_n_list(i,3);
end
imgout = cat(3,red,blue);
imgout = cat(3,imgout, green);
imgout = uint8(imgout);
imshow(imgout);