function [gx, gy] = gradient_image(img, sigma, t)
    sgf = sign(t) .*exp(-sigma.*abs(t));
    gx = convn(img, sgf', 'same');
    gy = convn(img, sgf, 'same');
end