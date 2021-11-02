function m = mom(img, p,q)

[row col] = size(img);
    for i = 1:row
        for j = 1:col
            m(i,j) = i.^p * j.^q * img(i,j);
        end
    end
    m = sum(sum(m));
end 