clear

% Cards
% x1=880;
% y1=350;
% rect=[x1 y1 500 500];
x1=1500;
y1=1550;
%300 - 1000
rect=[x1 y1 500 500];

% %Cactus
% x1=60;
% y1=370;
% rect=[x1 y1 90 90];

% % cactus loop
% for n=100:2:150
%     img{n} = imread(sprintf('I%d.jpg',n));

% cards loop    
for n=584:1:610
   img{n} = imread(sprintf('IMG_4%d.jpg',n));
   img{n}= rgb2gray(img{n});
    
    img{n}=imcrop(img{n},rect);
    [M N] = size(img{n});
    img{n} = im2double(img{n});
    
    %Variance 
    i_hat{n} = (1/(M*N))*sum(sum(img{n}));
    f_var{n} = (1/(M*N))*sum(sum((img{n}-i_hat{n}).^2));
    
    % EIG
    fxy{n} = padarray(img{n},[1 1],'symmetric','both');
    fxy_row{n}= circshift(fxy{n}, [0 1]); % right 
    fxy_col{n}= circshift(fxy{n}, [1 0]); %down
    fx{n}= fxy_row{n} - fxy{n}; % fx= f(x+1,y)-f(x,y)
    fy{n}= fxy_col{n} - fxy{n}; % fy= f(x,y+1)-f(x,y)
    
    f_eig{n}= sum(sum(fx{n}.^2+fy{n}.^2));
    
    %EIL
    fxy_d{n}= circshift(fxy{n}, [1 0]); %f(x,y+1)
    fxy_r{n}= circshift(fxy{n}, [0 1]); %f(x+1,y)
    fxy_u{n}= circshift(fxy{n}, [-1 0]); %f(x,y-1)
    fxy_l{n}= circshift(fxy{n}, [0 -1]); %f(x-1,y)
    fxy_ul{n}= circshift(fxy{n}, [-1 -1]); %f(x-1,y-1)
    fxy_ur{n}= circshift(fxy{n}, [-1 1]); %f(x+1,y-1)
    fxy_dr{n}= circshift(fxy{n}, [1 1]); % f(x+1,y+1)
    fxy_dl{n}= circshift(fxy{n}, [1 -1]); %f(x-1,y+1)
    
      % fxx+fyy= -f(x-1,y-1)-4f(x-1,y)-f(x-1,y+1)-4f(x,y-1)
      %          +20f(x,y)-4f(x,y+1)-f(x+1,y-1)-4f(x+1,y)-f(x+1,y+1)
    fxx{n}=-fxy_ul{n}-4*fxy_l{n}-fxy_dl{n}-4*fxy_u{n};
    fyy{n}=20*fxy{n}-4*fxy_d{n}-fxy_ur{n}-4*fxy_r{n}-fxy_dr{n};
  
    f_eil{n} = sum(sum((fxx{n}+fyy{n}).^2));
    
    % M22
    m00= mom(img{n},0,0);
    m10= mom(img{n},1,0);
    m01= mom(img{n},0,1);
    
    xhat = m10/m00;
    yhat = m01/m00;
    
    x = 1:M;
    y = 1:N;
    f_m22{n}= sum(sum(((x-xhat)'.^2*(y-yhat).^2).* img{n} ));
    
    % FT1 & FT2
    A=fftshift(fft2(img{n}));
    logA=log(1+abs(A)); 
    out=logA/max(logA(:));
    logA(floor(N/2)+1, floor(M/2)+1) = 0;
    A = abs(logA);
    %imshow(out)
        
    ci = [M/2, M/2, 5];     % center and radius of circle ([c_row, c_col, r])
    [xx,yy] = ndgrid((1:M)-ci(1),(1:N)-ci(2));
    mask = logical((xx.^2 + yy.^2)<ci(3)^2);
        
    circleim = double(zeros(M,N));
    circleim = A.*mask;
    invcircleim = A.*(1-mask);
        
    uv12 = sum(sum(circleim));
    uv23 = sum(sum(invcircleim));

    f_ft1{n} = uv23/uv12;
    f_ft2{n} = uv23;
    
end

var = cell2mat(f_var);
var = norma(var);
eig = cell2mat(f_eig);
eig = norma(eig);
eil = cell2mat(f_eil);
eil = norma(eil);
m22 = cell2mat(f_m22);
m22 = norma(m22);
ft1 = cell2mat(f_ft1);
ft1 = norma(ft1);
ft2 = cell2mat(f_ft2);
ft2 = norma(ft2);

plot(var,'-o','MarkerSize',3, 'LineWidth',1);
hold on;
plot(eig, '-o','MarkerSize',3, 'LineWidth',1);
plot(eil, '-o','MarkerSize',3, 'LineWidth',1);
plot(m22, '-o','MarkerSize',3, 'LineWidth',1);
plot(ft1, '-o','MarkerSize',3, 'LineWidth',1);
plot(ft2, '-o','MarkerSize',3, 'LineWidth',1);

hold off;
lgd = legend(["VAR","EIG","EIL","M22","FT1","FT2"]);
title(lgd,'Focus measures');
ax = gca;
ax.YGrid = 'off';
ax.XGrid = 'on';
xticks(0:1:30)
