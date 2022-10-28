function [H,S,L] =rgb2hsl(img)
    if length(size(img))~=3 && size(img,3) ~= 3
        error("please input a color image.");
    end
    img = double(img);
    [h,w,~] = size(img);
    H = zeros(h,w); S = H; L = H;
    rgbMax = max(img,[],3);
    rgbMin = min(img,[],3);
    L = (rgbMax + rgbMin)./2;
    MmM = rgbMax - rgbMin;
    for i = 1:h
        for j = 1:w
            r = img(i,j,1); g = img(i,j,2); b = img(i,j,3);
            Max = rgbMax(i,j); Min = rgbMin(i,j);
            mmm = MmM(i,j); l = L(i,j);
            if Max == Min
                H(i,j) = 0;
            elseif Max == r 
                if g >= b
                    H(i,j) = 60*(g-b)/ mmm;
                else
                    H(i,j) = 60*(g-b)/ mmm + 360;
                end
            elseif Max == g
                H(i,j) = 60*(b-r)/ mmm + 120;
            else
                H(i,j) = 60*(r-g)/ mmm + 240;
            end
            
            if Max == Min || l == 0
                S(i,j) = 0;
            elseif l <= 1/2*255
                S(i,j) = mmm/ (2*l);
            else
                S(i,j) = mmm/ (2*255-2*l);
            end
        end
    end
end