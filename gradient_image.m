function [gx, gy] =gradient_img(img, sigma, t)
    sgf = sign(t) .*exp(-sigma.*abs(t));
    gx = convn(img, sgf', 'same');
    gy = convn(img, sgf, 'same');
end