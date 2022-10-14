img = imread("lena_color.png");
img = double(rgb2ycbcr(img));

% initial k point 
sz = size(img);
k = 30;
Knums = randperm(prod(sz(1:2)),k);
Knums_x = mod(Knums,sz(1)); Knums_y = ceil(Knums/sz(1));
lambda1 = 0.3; lambda2 = 0.6;
[M,N] = meshgrid([1:sz(1)],[1:sz(2)]);
kmeans = zeros(k,5);    % 5: m, n, y, cb, cr
for i = 1:k
    kmeans(i,:) = [Knums_x(i), Knums_y(i), squeeze(img(Knums_x(i),Knums_y(i),:))'];
end
R = zeros(sz(1:2));
Y = img(:,:,1); Cb = img(:,:,2); Cr = img(:,:,3);

for repeat = 1:15   
    % dist
    d = zeros([sz(1:2),k]);
    for i = 1:k
        d(:,:,i) = sqrt(lambda1.*((M-kmeans(i,1)).^2+(N-kmeans(i,2)).^2)+ lambda2.*( (Y - kmeans(i,3)).^2+ ...
            (Cb- kmeans(i,4)).^2 + (Cr-kmeans(i,5)).^2));
    end
    
    % distribute region 
    for i = 1:sz(1)
        for j = 1:sz(1)
            [~,r] = min(d(i,j,:));
            R(i,j) = r(1);
        end
    end
    
    % update mean
    for i = 1:k
        loc = find(R == i);
        kmeans(i,:) = [mean(M(loc)), mean(N(loc)), mean(Y(loc)), mean(Cb(loc)), mean(Cr(loc))];
    end
end
% paintiing
Y_out = zeros(sz(1:2)); Cb_out = Y_out; Cr_out = Y_out;
for i = 1:k
    loc = find(R == i);
    Y_out(loc) = kmeans(i,3); Cb_out(loc) = kmeans(i,4); Cr_out(loc) = kmeans(i,5);
end
imgout = ycbcr2rgb(uint8(cat(3,Y_out,Cb_out,Cr_out)));
imshow(imgout);
imwrite(imgout,"./kmeans/lena_iter15_k30.png");

