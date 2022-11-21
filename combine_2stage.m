clear;
load("lena_color.mat");

sz = size(img);
h = sz(1); w = sz(2);
nC = ceil(w*h/200);

texture_factor = 1;
lap_factor = 1.3;
edge_factor = 0.6;

score = 7;
% score_rough = 6;

w_general = [1, 0.8, 0.8,  1/30, 0.8, 1,];
weight = [w_general,ones(1,3)*edge_factor, ones(1,3)*lap_factor, ones(1,6)*texture_factor];


YCbCrImg = double(rgb2ycbcr(img));
sf_sigma = 1;   % small filter sigma
lf_sigma = 0.5  % long filter sigma, choose lower than 0.5, i think
sf_t = [-10:10];
lf_t = [-20:20];

[gx, gy] = gradient_img(YCbCrImg, sf_sigma, sf_t);

g = (sqrt(gx.^2 + gy.^2));
% seperate g to 3 layer
g1 = g(:,:,1); g2 = g(:,:,2); g3 = g(:,:,3);
g1 = normalize(g1); g2 = normalize(g2); g3 = normalize(g3);
gx = abs(gx); gy = abs(gy);
% gx1 = normalize(gx(:,:,1)); gx2 = normalize(gx(:,:,2)); gx3 = normalize(gx(:,:,3));
% gy1 = normalize(gy(:,:,1)); gy2 = normalize(gy(:,:,2)); gy3 = normalize(gy(:,:,3));
gx = gx .* (10/255); gy = gy.* (10/255);
gx1 = gx(:,:,1); gx2 = gx(:,:,2); gx3 = gx(:,:,3);
gy1 = gy(:,:,1); gy2 = gy(:,:,2); gy3 = gy(:,:,3);

% YCbCr
Y = YCbCr(:,:,1); Cb = YCbCr(:,:,2); Cr = YCbCr(:,:,3);
Y = normalize(Y); Cb = normalize(Cb); Cr = normalize(Cr);

% HSL
[H, S, L] = rgb2hsl(img);
S = normalize(S); L = normalize(L);

% small filter lap
lap = abs(laplacian(img, 2.5, sf_t));
lap = lap.* (10/255);
Lap1 = lap(:,:,1); Lap2 = lap(:,:,2); Lap3 = lap(:,:,3);

maxLabel = max(seg(:));

img_size = prod(sz(1:2))

% first stage
for i = 0: maxLabel
    [region_edge, region_adj] = findEdgeRegion(seg, i, 4);
    combine_region = [];
    r_i = (seg == i);
    area_i = nnz(r_i);
    if (area_i == 0)
        continue;
    end
    boundary_i = getBoundariesLength2D(r_i);
    for j = region_adj
        edge_adj = (region_edge & seg == j);
        r_j = (seg == j);
        
        % calculate score
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

        % adjust threshold
        boundary_j = getBoundariesLength2D(r_j);
        area_j = nnz(r_j);
        boundary_ij = nnz(edge_adj);
        t1 = min(area_i, area_j);
        t2 = boundary_ij / min(boundary_i, boundary_j);
        gA = mean(gx1(r_i))^2 + mean(gy1(r_i))^2;
        gB = mean(gx1(r_j))^2 + mean(gy1(r_j))^2;
        t3 = min(gA^ alpha, gB^ alpha);
        score_r = 1 + (exp(-t1/img_size))/5  + (t2)/10 +((t3-b3)+ abs(t3-b3))*3 -(exp(-rd_H/180))/2;
        score_adj = score* score_r;
        if score_adj > sum(weight.*[rd_y,rd_Cb, rd_Cr, rd_H, rd_S, rd_L, e_g1, e_g2, e_g3, e_lap1, e_lap2, e_lap3, ...
            rd_gx1, rd_gx2, rd_gx3, rd_gy1, rd_gy2, rd_gy3])
            combine_region = [combine_region, j];
        end
    end
    seg(ismember(seg, combine_region)) = i;
end

% long filter applied before 2nd stage
% updating variables
score = 6;
w_general = [1, 0.8, 0.8,  1/30, 0.8, 1,];
weight = [w_general,ones(1,3)*edge_factor, ones(1,3)*lap_factor, ones(1,6)*texture_factor];

[gx, gy] = gradient_img(YCbCrImg, lf_sigma, lf_t);
g = (sqrt(gx.^2 + gy.^2));
% seperate g to 3 layer
g1 = g(:,:,1); g2 = g(:,:,2); g3 = g(:,:,3);
g1 = normalize(g1); g2 = normalize(g2); g3 = normalize(g3);
gx = abs(gx); gy = abs(gy);
% gx1 = normalize(gx(:,:,1)); gx2 = normalize(gx(:,:,2)); gx3 = normalize(gx(:,:,3));
% gy1 = normalize(gy(:,:,1)); gy2 = normalize(gy(:,:,2)); gy3 = normalize(gy(:,:,3));
gx = gx .* (10/255); gy = gy.* (10/255);
gx1 = gx(:,:,1); gx2 = gx(:,:,2); gx3 = gx(:,:,3);
gy1 = gy(:,:,1); gy2 = gy(:,:,2); gy3 = gy(:,:,3);

lap = abs(laplacian(img, 8, lf_t));
lap = lap.* (10/255);
Lap1 = lap(:,:,1); Lap2 = lap(:,:,2); Lap3 = lap(:,:,3);

% second stage
for i = 0: maxLabel
    [region_edge, region_adj] = findEdgeRegion(seg, i, 4);
    combine_region = [];
    r_i = (seg == i);
    area_i = nnz(r_i);
    if (area_i == 0)
        continue;
    end
    boundary_i = getBoundariesLength2D(r_i);
    for j = region_adj
        edge_adj = (region_edge & seg == j);
        r_j = (seg == j);
        
        % calculate score
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

        % adjust threshold
        boundary_j = getBoundariesLength2D(r_j);
        area_j = nnz(r_j);
        boundary_ij = nnz(edge_adj);
        t1 = min(area_i, area_j);
        t2 = boundary_ij / min(boundary_i, boundary_j);
        gA = mean(gx1(r_i))^2 + mean(gy1(r_i))^2;
        gB = mean(gx1(r_j))^2 + mean(gy1(r_j))^2;
        t3 = min(gA^ alpha, gB^ alpha);
        score_r = 1 + (exp(-t1/img_size))/5  + (t2)/10 +((t3-b3)+ abs(t3-b3))*3 -(exp(-rd_H/180))/2;
        score_adj = score* score_r;
        if score_adj > sum(weight.*[rd_y,rd_Cb, rd_Cr, rd_H, rd_S, rd_L, e_g1, e_g2, e_g3, e_lap1, e_lap2, e_lap3, ...
            rd_gx1, rd_gx2, rd_gx3, rd_gy1, rd_gy2, rd_gy3])
            combine_region = [combine_region, j];
        end
    end
    seg(ismember(seg, combine_region)) = i;
end
