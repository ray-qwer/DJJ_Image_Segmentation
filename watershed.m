img = imread("lena_color.png");
img = double(rgb2ycbcr (img));

% some var
sz = size(img);
Q = 3;

% gradient 
%   filter initial
sigma = 1;
t = [-10:10];
sgf = sign(t).*exp(-sigma.*abs(t));

gx = convn(img, sgf','same').^2;
gy = convn(img, sgf,'same').^2;
g = zeros(sz(1:2));

g = sqrt(sum(gx,3) + sum(gy,3));
L = round(g/Q);
[R,region] = bwlabel(L == 0);

% level = 1;
for level = 1:max(L(:))
    level_find = find(L == level);
    
    % fast algo
    addi_region = 0;
    R_sub = zeros(sz(1:2));
    for i = level_find'
        
        NB = findNeighbor(i, sz(1:2));
        if all(R(NB) == 0)
            if all(R_sub(NB) == 0)
                addi_region = addi_region + 1;
                R_sub(i) = region + addi_region;
            else
                k = R_sub(NB);
                NB1 = NB(k~=0 & k<= region);
                if ~isempty(NB1)
                    [~,rNB] = max(g(i) - g(NB1));
                    R_sub(i) = R_sub(NB1(rNB)); 
                    sub = R_sub(i);
                else
                    % combine all with connection
                    k = k(k~=0);
                    R_sub(i) = min(k);
                    sub = R_sub(i);
                end
                for j = k
                    if j ~= sub && j > region
                        R_sub(R_sub == j) = sub;
                    end
                end
            end
        else
            [~,rNB] = max(g(i)-g(NB));
            R_sub(i) = R(NB(rNB));
            k = R_sub(NB);
            for j = k
                if j > region
                    R_sub(R_sub == j) = R_sub(i); 
                end
            end
        end   
    end
region = addi_region + region;
R = R_sub + R;
end
imgout1 = label2rgb(R);
figure(1);
imshow(imgout1);
r_layer = zeros(sz(1:2)); g_layer = r_layer; b_layer = r_layer;
r_img = img(:,:,1); g_img = img(:,:,2); b_img = img(:,:,3);
for i = 1:max(R(:))
    k = find(R == i);
    r_layer(k) = mean(r_img(k)); g_layer(k) = mean(g_img(k)); b_layer(k) = mean(b_img(k));
end
figure(2);
imgout2 = uint8(cat(3,r_layer,g_layer,b_layer));
imshow(ycbcr2rgb (imgout2));
function Neighbor =findNeighbor(i,sz)
    Neighbor = [];
    x = mod(i-1,sz(1))+1;
    y = ceil(i/sz(1));
    switch x
        case 1
            Neighbor = [Neighbor, i+1];
        case sz(1)
            Neighbor = [Neighbor, i-1];
        otherwise
            Neighbor = [Neighbor, i-1, i+1];
    end
    switch y
        case 1
            Neighbor = [Neighbor, i+sz(1)];
        case sz(2)
            Neighbor = [Neighbor, i-sz(1)];
        otherwise
            Neighbor = [Neighbor, i-sz(1), i+sz(1)];
    end
end
