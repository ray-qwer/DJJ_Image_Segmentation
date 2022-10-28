name = 'lena_color';
path = '';
imgPath = strcat(name, '.png');
oriImg = imread(imgPath);
grey_img = double(rgb2gray(oriImg));
[h,w]=size(grey_img);
nC = floor(w*h/200);
% nC = 40;
t = cputime;
segments = mex_ers(grey_img,nC);
edge=(segments~=segments(:,[1,1:w-1])) | (segments~=segments([1,1:h-1],:));
image(edge*255+grey_img*0.7)
colormap(gray(256))
