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