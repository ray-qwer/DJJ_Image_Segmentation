img = imread("lena.png");
img = double(rgb2gray(img));
tic
% fast algorithm
th = 30;
min_pixel = ceil(512*512/1000);
R = zeros(size(img));
m_n_list = zeros(0,2); %A: 1, B: 2
init = true;
% edge finding with sign gaussian filter
sigma = 10; ch=1; cv=100;
t = [-10:10];
sgf = sign(t).*exp(-sigma*abs(t));
gh = conv2(img, sgf,'same').*ch;
gv = conv2(img, sgf','same').*cv;
lambda = 0.1;

sz = size(img);
r_cnt = 0;
for m = 1:sz(1)
    for n = 1:sz(2)
        row_combine = false; col_combine = false;
        % compare by rc and cc
        pixel = img(m,n); 
        updiff = 0; leftdiff = 0; gradUpDiff = 0; gradLeftDiff = 0;
        if n > 1
            leftdiff = pixel - m_n_list(R(m,n-1),1);  gradLeftDiff = abs(leftdiff)+abs(gh(m,n))*lambda;
            if gradLeftDiff <= th
                col_combine = true;
            end
        end
        if m > 1
            updiff = pixel - m_n_list(R(m-1,n),1); gradUpDiff = abs(updiff)+abs(gv(m,n))*lambda;
            if gradUpDiff <= th
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
                if abs(updiff - leftdiff) <= th
                    flag = 2; ref_point = min(R(m-1,n),R(m,n-1)); replace_point = max(R(m-1,n),R(m,n-1));
                elseif gradUpDiff > gradLeftDiff
                    flag = 1; ref_point = R(m,n-1);
                else
                    flag = 1; ref_point = R(m-1,n);
                end
            end    
        elseif col_combine
            flag = 1; ref_point = R(m,n-1);
        elseif row_combine
            flag = 1; ref_point = R(m-1,n);
        else
            r_cnt = r_cnt + 1; ref_point = r_cnt; flag = 0;
        end

        R(m,n) = ref_point; 
        switch flag
            case 0  % add new region
                m_n_list = cat(1,m_n_list,[pixel,1]);
            case 1  % add one new point
                new_num = m_n_list(ref_point,2) + 1;
                new_mean = (prod(m_n_list(ref_point,:)) + pixel)/ new_num;
                m_n_list(ref_point,:) = [new_mean, new_num];
            case 2  % combine one region to another
                new_num = m_n_list(ref_point,2) + m_n_list(replace_point,2) + 1;
                new_mean = (prod(m_n_list(ref_point,:)) + prod(m_n_list(replace_point,:)) + pixel) / new_num;
                m_n_list(ref_point,:) = [new_mean, new_num];
                m_n_list(replace_point,:) = [0, 0];
                R(R == replace_point) = ref_point;
            otherwise
                disp("error");
        end
        
    end
end
toc
t = m_n_list(:,2);
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
    if ~isempty(u(m_n_list(u,2) > min_pixel))
        u = u(m_n_list(u,2) > min_pixel);
    end
%     [min_dist,u_index]=min((abs(mean_num_list(index(i),1) - mean_num_list(u,1))./sqrt(mean_num_list(u,2))));
    [min_dist,u_index]=min((abs(m_n_list(index(i),1) - m_n_list(u,1))));
    % update mean, num, region
    region_index = u(u_index(1));
    new_num = m_n_list(region_index,2) + m_n_list(index(i),2);
    m_n_list(region_index,2) = new_num;
    R(R == index(i)) = region_index;
    m_n_list(index(i),:) = [0, 0];
end
toc
for i = 1:r_cnt
    R(R == i) = m_n_list(i,1);
end
R = uint8(R);
imshow(R);