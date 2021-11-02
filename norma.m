function n = norma(f)

xmax = max(f);
xmin = min(f);
n = (f-xmin)/(xmax-xmin);

end 