img = imread("lena.png");
img = double(rgb2gray(img));

% fast algorithm
th = 25;
region = zeros(size(img));
mean_num_list = zeros(1,2); %A: 1, B: 2

sz = size(img);
r_cnt = 0;
for m = 1:sz(1)
    for n = 1:sz(2)
        row_check = 0; col_check = 0;
        % row_check = 1 -> check m -1
        % col_check = 1 -> check n -1
        if n ~= 1
            col_check = 1;
        end
        if m ~= 1
            row_check = 1;
        end
        
        if (row_check == col_check) == 0
            % first element
            r_cnt = r_cnt + 1;
            region(m,n) = r_cnt;
            mean_num_list(r_cnt) = [img(m,n), 1];
        else
            % others, compare by rc and cc
            
        end
    end
end