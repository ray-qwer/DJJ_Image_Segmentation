clear;
img = imread("lena_color.png");
img = double(rgb2ycbcr(img));

% initial k point 
sz = size(img);
k = 500;
Knums = randperm(prod(sz(1:2)),k);
Knums_x = mod(Knums-1,sz(1))+1; Knums_y = ceil(Knums/sz(1));

% edge detection
sigma = 1;
t = [-10:10];
sgf = sign(t).*exp(-sigma.*abs(t));
gx = convn(img, sgf','same');
gy = convn(img, sgf,'same');
g = zeros(sz(1:2));
% norm2
g = (abs(gx(:,:,1)).^2+abs(gx(:,:,2)).^2+abs(gx(:,:,3).^2) + abs(gy(:,:,1)).^2+abs(gy(:,:,2)).^2+abs(gy(:,:,3)).^2).^0.5;

% find new knums with lowest gradient
L = round(sz(1)/sqrt(k));
for i = 1:k
    c_x = Knums_x(i) - (-L:L); c_x = c_x(c_x > 0 & c_x <= sz(1));
    c_y = Knums_y(i) - (-L:L); c_y = c_y(c_y > 0 & c_y <= sz(2));
    [~,index] = min(g(c_x,c_y));
    Knums_x(i) = c_x(mod(index(1)-1,length(c_x))+1); Knums_y(i) = c_y(ceil(index(1)/length(c_x)));
end

lambda1 = 0.2; lambda2 = 0.6;
[M,N] = meshgrid([1:sz(1)],[1:sz(2)]);
kmeans = zeros(k,5);    % 5: m, n, y, cb, cr
for c = 1:k
    kmeans(c,:) = [Knums_x(c), Knums_y(c), squeeze(img(Knums_x(c),Knums_y(c),:))'];
end
R = zeros(sz(1:2));
Y = img(:,:,1); Cb = img(:,:,2); Cr = img(:,:,3);

% size for checking area
dist_x = round(sz(1)/(k^0.25));
dist_y = round(sz(2)/(k^0.25));
% dist_x = 100; dist_y = 100;
d = ones([sz(1:2),k])*inf;
for repeat = 1:15
    % dist
    for i = 1:k
        c_x = round(kmeans(i,1)) - (-dist_x:dist_x); c_x = c_x(c_x > 0 & c_x <= sz(1));
        c_y = round(kmeans(i,2)) - (-dist_y:dist_y); c_y = c_y(c_y > 0 & c_y <= sz(2));
        d(c_x,c_y,i) = sqrt(lambda1.*((N(c_x,c_y)-kmeans(i,1)).^2+(M(c_x,c_y)-kmeans(i,2)).^2)+ lambda2.*( (Y(c_x,c_y) - kmeans(i,3)).^2+ ...
            (Cb(c_x,c_y)- kmeans(i,4)).^2 + (Cr(c_x,c_y)-kmeans(i,5)).^2));
    end
    
    % distribute region 
    for i = 1:sz(1)
        for j = 1:sz(1)
            [~,r] = min(d(i,j,:));
            R(i,j) = r;
        end
    end
    
    % update mean
    for i = 1:k
        loc = find(R == i);
        kmeans(i,:) = [mean(N(loc)), mean(M(loc)), mean(Y(loc)), mean(Cb(loc)), mean(Cr(loc))];
    end
end

% paintiing

Y_out = zeros(sz(1:2)); Cb_out = zeros(sz(1:2)); Cr_out = zeros(sz(1:2));
% bwlabel
tic
addition_region = 0;
for i = 1:k
    
    bin_fig = zeros(size(R));
    bin_fig(R == i) = 1;
    [b, n] = bwlabel(bin_fig,4);
    if n == 1
        loc = find(R == i);
        Y_out(loc) = kmeans(i,3); Cb_out(loc) = kmeans(i,4); Cr_out(loc) = kmeans(i,5);
    else
        for j = 1:n
            loc = find(b == j);
            m_Y = mean(Y(loc)); m_Cb = mean(Cb(loc)); m_Cr = mean(Cr(loc));
            Y_out(loc) = m_Y; Cb_out(loc) = m_Cb; Cr_out(loc) = m_Cr;
            if j == 1
                kmeans(i,:) = [mean(N(loc)), mean(M(loc)), m_Y, m_Cb, m_Cr];
            else
                kmeans = cat(1,kmeans,[mean(N(loc)), mean(M(loc)), m_Y, m_Cb, m_Cr]);
                addition_region = addition_region + 1;
            end
        end
    end
end
toc

% fast algo
tic
f_R = zeros(sz(1:2));
r_cnt = 0;
for m = 1:sz(1)
    for n = 1:sz(2)
        row_combine = false; col_combine = false;
        % compare by rc and cc
        pixel = R(m,n);
        up_region = 0; left_region = 0;
        if n > 1
            left_region = f_R(m,n-1);
            if pixel == R(m,n-1)
                col_combine = true;
            end
        end
        if m > 1
            up_region = f_R(m-1,n);
            if pixel == R(m-1,n)
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
                flag = 2; ref_point = min(up_region, left_region); replace_point = max(up_region, left_region);
            end
        elseif col_combine
            flag = 1; ref_point = left_region;
        elseif row_combine
            flag = 1; ref_point = up_region;
        else
            r_cnt = r_cnt + 1; ref_point = r_cnt; flag = 0;
        end

        f_R(m,n) = ref_point; 
        switch flag
            case 0  % add new region
                
            case 1  % add one new point
                
            case 2  % combine one region to another
                f_R(f_R == replace_point) = ref_point;
            otherwise
                disp("error");
        end
    end
end
for i = 1:r_cnt
    loc = find(f_R == i);
    m_Y = mean(Y(loc)); m_Cb = mean(Cb(loc)); m_Cr = mean(Cr(loc));
    Y_out(loc) = m_Y; Cb_out(loc) = m_Cb; Cr_out(loc) = m_Cr;
end
toc
imgout = ycbcr2rgb(uint8(cat(3,Y_out,Cb_out,Cr_out)));
imshow(imgout);
% hold on;
% plot(round(kmeans(:,1)),round(kmeans(:,2)),"*r");
% plot(Knums_x,Knums_y,"*b");
% hold off;