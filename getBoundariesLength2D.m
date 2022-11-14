function l = getBoundariesLength2D(img)
    img = padarray(img,[1,1], 0, 'both');
    kernel = [1,-1];
    l = nnz(conv2(img, kernel)) + nnz(conv2(img, kernel')) -1;
end