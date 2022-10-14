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
g = ((abs(gx(:,:,1))+abs(gx(:,:,2))+abs(gx(:,:,3))).^2 + (abs(gy(:,:,1))+abs(gy(:,:,2))+abs(gy(:,:,3))).^2).^0.5;

% find new knums with lowest gradient
L = 20;
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

d = ones([sz(1:2),k])*inf;
for repeat = 1:15
    % dist
    for i = 1:k
        c_x = round(kmeans(i,1)) - (-dist_x:dist_x); c_x = c_x(c_x > 0 & c_x <= sz(1));
        c_y = round(kmeans(i,2)) - (-dist_y:dist_y); c_y = c_y(c_y > 0 & c_y <= sz(2));
        d(c_x,c_y,i) = sqrt(lambda1.*((M(c_x,c_y)-kmeans(i,1)).^2+(N(c_x,c_y)-kmeans(i,2)).^2)+ lambda2.*( (Y(c_x,c_y) - kmeans(i,3)).^2+ ...
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