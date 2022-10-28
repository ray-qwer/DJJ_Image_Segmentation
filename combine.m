img = imread("./lena_color.png");
gray_img = double(rgb2gray(img));
sz = size(img);
[h,w] = size(gray_img);
nC = ceil(w*h/200);
seg = mex_ers(gray_img,nC);
% edge=(seg~=seg(:,[1,1:w-1])) | (seg~=seg([1,1:h-1],:));

% gradient for edge detection
sigma = 1;
t = [-10:10];
sgf = sign(t).*exp(-sigma.*abs(t));
gx = convn(img, sgf','same');
gy = convn(img, sgf,'same');
g = (sqrt(gx.^2 + gy.^2));

% adjacent region and edge with 4 connection. 
kernel = [0,1,0; 1,1,1; 0,1,0];
d = zeros(size(seg));
d(seg == 20) = 1;
dilate = (conv2(d, kernel,'same')~=0)&(d == 0);
edge_7 = (dilate) & (seg == 7);
adj = unique(seg(dilate));

% Y, cb, cr
YCbCr = rgb2ycbcr(img);
Y = YCbCr(:,:,1); Cb = YCbCr(:,:,2); Cr = YCbCr(:,:,3);
% HSL
[H,S,L] = rgb2hsl(img);




