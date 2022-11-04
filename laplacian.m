function lap = laplacian(img,L)
    % lap_filter
    sigma = log(10)/(pi*L^2);
    t = -10:10;
    lap_filter = -(2*pi*sigma - (2*pi*sigma)^2 .* (t.^2)) .*exp(-pi*sigma.*(t.^2));
    
    img_x = convn(img, lap_filter,'same');
    img_y = convn(img, lap_filter','same');
    lap = img_x + img_y;
end