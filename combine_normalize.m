img = imread("./baboon.png");
gray_img = double(rgb2gray(img));
sz = size(img);
[h,w] = size(gray_img);
nC = ceil(w*h/200);
seg = mex_ers(gray_img,nC);
seg1 = seg;
% edge=(seg~=seg(:,[1,1:w-1])) | (seg~=seg([1,1:h-1],:));

% parameter to modify

% weight order: y, cb, cr, h, s, v, g1, g2, g3, 
% 0 <= g <= 255**(1/2)
% 0 <= Y <= 255;    0 <= Cb <= 255; 0 <= Cr <= 255
% 0 <= H < 360;     0 <= S <= 1;    0 <= L <= 255
texture_factor = 0.9;
lap_factor = 1.2;
alpha = 0.2;
score = 6.5;
weight = [1, 0.8, 0.8,  1/15, 1, 1, 1.2, 1.2, 1.2, ones(1,3)*lap_factor, ones(1,6)*texture_factor]; 

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
gx = abs(gx); gy = abs(gy);
gx1 = normalize(gx(:,:,1)); gx2 = normalize(gx(:,:,2)); gx3 = normalize(gx(:,:,3));
gy1 = normalize(gy(:,:,1)); gy2 = normalize(gy(:,:,2)); gy3 = normalize(gy(:,:,3));

% g1 = g1 - min(g1(:)); g2 = g2 - min(g2(:)); g3 = g3 - min(g3(:));
% Y, cb, cr
YCbCr = double(rgb2ycbcr(img));
Y = YCbCr(:,:,1); Cb = YCbCr(:,:,2); Cr = YCbCr(:,:,3);
Y = normalize(Y); Cb = normalize(Cb); Cr = normalize(Cr);
% HSL
[H,S,L] = rgb2hsl(img);
S = normalize(S); L = normalize(L);
% lap
lap = abs(laplacian(img, 2.5));
Lap1 = normalize(lap(:,:,1)); Lap2 = normalize(lap(:,:,2)); Lap3 = normalize(lap(:,:,3));


% score = 3+mean(g1(:))+mean(g2(:))+mean(g3(:));


maxLabel = max(seg(:)); % the label is start from zero
for i = 0:maxLabel
    % combination
    [region_edge, region_adj] = findEdgeRegion(seg, i, 4);
    combine_region = [];
    for j = region_adj
        edge_adj = (region_edge & seg == j);
        % factor to check if combine
        % e: edge, rd: region diff
        r_i = (seg == i); r_j = (seg == j);
        e_g1 = mean(g1(edge_adj)); e_g2 = mean(g2(edge_adj)); e_g3 = mean(g3(edge_adj));    % -255 <= g <= 255
        e_lap1 = mean(Lap1(edge_adj)); e_lap2 = mean(Lap2(edge_adj)); e_lap3 = mean(Lap3(edge_adj));
        rd_y  = abs(mean(Y(r_i)) -  mean(Y(r_j)));     % 0 <= Y <= 255
        rd_Cb = abs(mean(Cb(r_i)) - mean(Cb(r_j)));    % 0 <= Cb <= 255
        rd_Cr = abs(mean(Cr(r_i)) - mean(Cr(r_j)));    % 0 <= Cr <= 255
        rd_H  = min(abs(mean(H(r_i)) -  mean(H(r_j))),abs(mean(H(r_i)) + 360 -  mean(H(r_j)))); 
        % 0 <= H < 360, but 1 is similar to 359 => min(H1-H2, H1+360-H2)
        rd_S  = abs(mean(S(r_i)) -  mean(S(r_j)));     % 0 <= S <= 1
        rd_L  = abs(mean(L(r_i)) -  mean(L(r_j)));     % 0 <= L <= 255
        rd_gx1 = abs(mean(gx1(r_i))^alpha - mean(gx1(r_j))^alpha);      % texture feature
        rd_gx2 = abs(mean(gx2(r_i))^alpha - mean(gx2(r_j))^alpha);
        rd_gx3 = abs(mean(gx3(r_i))^alpha - mean(gx3(r_j))^alpha);
        rd_gy1 = abs(mean(gy1(r_i))^alpha - mean(gy1(r_j))^alpha);
        rd_gy2 = abs(mean(gy2(r_i))^alpha - mean(gy2(r_j))^alpha);
        rd_gy3 = abs(mean(gy3(r_i))^alpha - mean(gy3(r_j))^alpha);
        if score > sum(weight.*[rd_y,rd_Cb, rd_Cr, rd_H, rd_S, rd_L, e_g1, e_g2, e_g3, e_lap1, e_lap2, e_lap3, ...
                rd_gx1, rd_gx2, rd_gx3, rd_gy1, rd_gy2, rd_gy3])
            combine_region = [combine_region, j];
        end
    end
    seg(ismember(seg,combine_region)) = i;
end

edge=(seg~=seg(:,[1,1:w-1])) | (seg~=seg([1,1:h-1],:));
figure(1);
imgout = uint8(255*edge + double(img)*0.7);
imshow(imgout);

% save file
% fileID = fopen("./combine/recordNum_n.txt","r");
% if fileID ~= -1
%     picNum = fscanf(fileID, "%d")+1;
%     fclose(fileID);
% else
%     picNum = 1;
% end
% 
% 
% imwrite(imgout, "./combine/n"+picNum+".png");
% fileID = fopen("./combine/recordNum_n.txt","w");
% fprintf(fileID,"%d\n",picNum);
% fclose(fileID);
% writematrix([weight,score],"./combine/record_n.csv",'WriteMode','append');


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