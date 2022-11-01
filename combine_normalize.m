img = imread("./peppers.jpg");
gray_img = double(rgb2gray(img));
sz = size(img);
[h,w] = size(gray_img);
nC = ceil(w*h/200);
seg = mex_ers(gray_img,nC);
seg1 = seg;
% edge=(seg~=seg(:,[1,1:w-1])) | (seg~=seg([1,1:h-1],:));

% parameter to modify

% weight order: y, cb, cr, h, s, v, g1, g2, g3
% 0 <= g <= 255**(1/2)
% 0 <= Y <= 255;    0 <= Cb <= 255; 0 <= Cr <= 255
% 0 <= H < 360;     0 <= S <= 1;    0 <= L <= 255
weight = [1, 0.8, 0.8,  1/20, 1, 1, 1.2, 0.8, 0.8]; 

YCbCrImg = double(rgb2ycbcr(img));
% gradient for edge detection
sigma = 1;
t = [-10:10];
sgf = sign(t).*exp(-sigma.*abs(t));
gx = convn(YCbCrImg, sgf','same');
gy = convn(YCbCrImg, sgf,'same');
g = (sqrt(gx.^2 + gy.^2));
% seperate g to 3 layer
g1 = g(:,:,1); g2 = g(:,:,2); g3 = g(:,:,3);
g1 = normalize(g1); g2 = normalize(g2); g3 = normalize(g3);
g1 = g1 - min(g1(:)); g2 = g2 - min(g2(:)); g3 = g3 - min(g3(:));
% Y, cb, cr
YCbCr = double(rgb2ycbcr(img));
Y = YCbCr(:,:,1); Cb = YCbCr(:,:,2); Cr = YCbCr(:,:,3);
Y = normalize(Y); Cb = normalize(Cb); Cr = normalize(Cr);
% HSL
[H,S,L] = rgb2hsl(img);
S = normalize(S); L = normalize(L);
score = 3+mean(g1(:))+mean(g2(:))+mean(g3(:));

maxLabel = max(seg(:)); % the label is start from zero
for i = 0:maxLabel
    % combination
    [region_edge, region_adj] = findEdgeRegion(seg, i, 4);
    combine_region = [];
    for j = region_adj
        edge_adj = (region_edge & seg == j);
        % factor to check if combine
        e_g1 = mean(g1(edge_adj)); e_g2 = mean(g2(edge_adj)); e_g3 = mean(g3(edge_adj));    % -255 <= g <= 255
        e_y  = abs(mean(Y(seg == i)) -  mean(Y(seg == j)));     % 0 <= Y <= 255
        e_Cb = abs(mean(Cb(seg == i)) - mean(Cb(seg == j)));    % 0 <= Cb <= 255
        e_Cr = abs(mean(Cr(seg == i)) - mean(Cr(seg == j)));    % 0 <= Cr <= 255
        e_H  = min(abs(mean(H(seg == i)) -  mean(H(seg == j))),abs(mean(H(seg == i)) + 360 -  mean(H(seg == j)))); 
        % 0 <= H < 360, but 1 is similar to 359 => min(H1-H2, H1+360-H2)
        e_S  = abs(mean(S(seg == i)) -  mean(S(seg == j)));     % 0 <= S <= 1
        e_L  = abs(mean(L(seg == i)) -  mean(L(seg == j)));     % 0 <= L <= 255
        if score >= sum(weight.*[e_y,e_Cb, e_Cr, e_H, e_S, e_L, e_g1, e_g2, e_g3])
            combine_region = [combine_region, j];
        end
    end
    seg(ismember(seg,combine_region)) = i;
end

edge=(seg~=seg(:,[1,1:w-1])) | (seg~=seg([1,1:h-1],:));
% figure(1);
imgout = uint8(255*edge + double(img)*0.7);
% imshow(imgout);

% save file
fileID = fopen("./combine/recordNum_n.txt","r");
if fileID ~= -1
    picNum = fscanf(fileID, "%d")+1;
    fclose(fileID);
else
    picNum = 1;
end


imwrite(imgout, "./combine/n"+picNum+".png");
fileID = fopen("./combine/recordNum_n.txt","w");
fprintf(fileID,"%d\n",picNum);
fclose(fileID);
writematrix([weight,score],"./combine/record_n.csv",'WriteMode','append');


function [edge, adj] = findEdgeRegion(seg, region, channel)
    switch channel
        case 4
            kernel =[0,1,0; 1,1,1; 0,1,0];
        case 8
            kernel = [1,1,1; 1,1,1; 1,1,1];
        otherwise
            error("kernel should be 4 or 8");
    end
    if region > max(seg(:))
        edge = []; adj = []; return;
    end
    d = zeros(size(seg));
    d(seg == region) = 1;
    edge = (conv2(d,kernel,'same')~=0) & (d == 0);
    adj = unique(seg(edge))';
end