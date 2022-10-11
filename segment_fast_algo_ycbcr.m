img = imread("lena_color.png");
img = double(rgb2ycbcr(img));
% threshold need to be modified
% diff equation is not the same
th = 25; sigma = 1; lambda = [0.1, 0.1, 0.1]; diff_factor = [1, 0.5, 1.5];
image_name = "lena_cr1.5_1";
genSegImage(img,lambda, sigma, diff_factor, th, image_name);
% for th = [17:33]
%      image_name = "lena_th"+th;
%      genSegImage(img,lambda, sigma, diff_factor, th, image_name );
% end
% th = 25;
% for diff = [0.5:0.1:1.5]
%     diff_factor = [diff, 1+(1-diff)/2, 1+(1-diff)/2 ];
%     image_name = "lena_diff"+diff;
%     genSegImage(img, lambda, sigma, diff_factor, th, image_name);
% end
% diff = [0.5, 1.25, 1.25];
% for l = [0.1:0.2:2]
%     lambda = [l, l, l];
%     image_name = "lena_l"+l;
%     genSegImage(img, lambda, sigma, diff_factor, th, image_name);
% end

function genSegImage(img, lambda, sigma, diff_factor, th, image_name)
% fast algorithm
% th = 25;
sz = size(img);
min_pixel = ceil(512*512/500);
R = zeros(sz(1:2));
m_n_list = zeros(0,4); %A: 1, B: 2
% edge finding with sign gaussian filter
% sigma = 1; 
ch=1; cv=1;
t = [-10:10];
sgf = sign(t).*exp(-sigma.*abs(t));
ghrgb = convn(img, sgf,'same').*ch;
gvrgb = convn(img, sgf','same').*cv;
% lambda = [0.1, 0.1, 0.1]*8;

% diff_factor = [0.5,1,1];

r_cnt = 0;
for m = 1:sz(1)
    for n = 1:sz(2)
        row_combine = false; col_combine = false;
        % compare by rc and cc
        pixel = squeeze(img(m,n,:))'; up_region = 0; left_region = 0;
        updiff = 0; leftdiff = 0; gradUpDiff = 0; gradLeftDiff = 0;
        if n > 1
            left_region = R(m,n-1);
            leftdiff = (pixel - m_n_list(left_region,1:3)).*diff_factor;  gradLeftDiff = sum(abs(leftdiff))*1+lambda*abs(squeeze(ghrgb(m,n,:)));
            if gradLeftDiff <= th
                col_combine = true;
            end
        end
        if m > 1
            up_region = R(m-1,n);
            updiff = (pixel - m_n_list(up_region,1:3)).*diff_factor; gradUpDiff = sum(abs(updiff))*1+lambda*abs(squeeze(gvrgb(m,n,:)));
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
            if up_region == left_region
                flag = 1; ref_point = up_region;
            else
                if sum(abs(updiff - leftdiff)) <= th
                    flag = 2; ref_point = min(up_region,left_region); replace_point = max(up_region,left_region);
                elseif gradUpDiff > gradLeftDiff
                    flag = 1; ref_point = left_region;
                else
                    flag = 1; ref_point = up_region;
                end
            end    
        elseif col_combine
            flag = 1; ref_point = left_region;
        elseif row_combine
            flag = 1; ref_point = up_region;
        else
            r_cnt = r_cnt + 1; ref_point = r_cnt; flag = 0;
        end

        R(m,n) = ref_point; 
        switch flag
            case 0  % add new region
                m_n_list = cat(1,m_n_list,[pixel,1]);
            case 1  % add one new point
                new_num = m_n_list(ref_point,end) + 1;
                new_mean = (m_n_list(ref_point,1:3).* m_n_list(ref_point,end) + pixel)/ new_num;
                m_n_list(ref_point,:) = [new_mean, new_num];
            case 2  % combine one region to another
                new_num = m_n_list(ref_point,end) + m_n_list(replace_point,end) + 1;
                new_mean = (m_n_list(ref_point,1:3).*m_n_list(ref_point,end) + m_n_list(replace_point,1:3).*m_n_list(replace_point,end) + pixel) / new_num;
                m_n_list(ref_point,:) = [new_mean, new_num];
                m_n_list(replace_point,:) = zeros(1,4);
                R(R == replace_point) = ref_point;
            otherwise
                disp("error");
        end
    end
end

t = m_n_list(:,end);
index = find(t<min_pixel & t >0);
neighbor_matrix = [0,1,0;1,0,1;0,1,0];
tic
for i = 1:length(index)
    M = zeros(sz(1:2));
    di_M = zeros(sz(1:2));
    [r,c] = find(R == index(i));
    M((c-1)*sz(1)+r) = 1;
    di_M(conv2(M,neighbor_matrix,'same') > 0) = 1;
    [n_r,n_c] =  find((di_M - M) == 1); % get neighbor
    u = unique(R((n_c-1)*sz(1)+n_r));
    if ~isempty(u(m_n_list(u,2) > min_pixel))
        u = u(m_n_list(u,2) > min_pixel);
    end
%     [min_dist,u_index]=min((abs(mean_num_list(index(i),1) - mean_num_list(u,1))./sqrt(mean_num_list(u,2))));
    [~,u_index]=min(sum(abs(((m_n_list(index(i),1:3) - m_n_list(u,1:3)).*diff_factor(1:3))),2));
    % update mean, num, region
    region_index = u(u_index(1));
    new_num = m_n_list(region_index,end) + m_n_list(index(i),end);
    m_n_list(region_index,end) = new_num;
    R(R == index(i)) = region_index;
    m_n_list(index(i),:) = zeros(1,4);
end
toc
rlayer = zeros(size(R));
blayer = zeros(size(R));
glayer = zeros(size(R));
for i = 1:r_cnt
    rlayer(R == i) = m_n_list(i,1);
    blayer(R == i) = m_n_list(i,2);
    glayer(R == i) = m_n_list(i,3);
end
imgout = cat(3, rlayer,blayer);
imgout = cat(3, imgout,glayer);
imgout = uint8(imgout);
imgout = ycbcr2rgb(imgout);
imwrite(imgout,"./fast_algo_image/"+image_name+".png");
end
% imshow(imgout);