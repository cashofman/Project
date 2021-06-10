%% start
clear all
close all
clc
%% load image
I = imread('./triangle1.png');
%I = double(I);
%I = rgb2gray(I);
%I = I(100:end-50,1:end-80);
s=4;
I = I(1:s:end,1:s:end,:);
I = double(I);
imshow(I,gray);

[col1, row1] = size(I);
%% selecting area
grayImage = I;
imshow(I, []);
axis on;
title('Original Grayscale Image');
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.

% Ask user to draw freehand mask.
message = sprintf('Left click and hold to begin drawing.\nSimply press enter to finish');
uiwait(msgbox(message));
hFH = imfreehand(); % Actual line of code to do the drawing.
% Create a binary image ("mask") from the ROI object.
binaryImage = hFH.createMask();
xy = hFH.getPosition;


%% deleting area
P1 = [320/s, 204/s];
P2 = [320/s, 1003/s];
P3 = [1012/s, 602/s];
s = det([P1-P2;P3-P1]);
J2 = I;
C = ones(col1, row1); % Confidence term
F = ones(col1, row1); % Filled term


for j=1:row1
    for i=1:col1
        P = [i,j];
        if s*det([P3-P;P2-P3])>=0 && s*det([P1-P;P3-P1])>=0 && s*det([P2-P;P1-P2])>=0
            J2(i,j) = NaN;
            C(i,j) = 0;
            F(i,j) = 0;
        end
    end
end 

figure
R = J2;
G = J2;
R(isnan(R)) = 255;
G(isnan(R)) = 0;
B = G;
Im(:,:,1) = R;
Im(:,:,2) = G;
Im(:,:,3) = B;
Im = uint8(Im);
x = imshow(Im);

saveas(x,'./Results/triangle_start.tif')
im = J2;



%% first iteration
n_range = [7, 9, 11];
for k = 1:numel(n_range)
    n = n_range(k);
    iter = 0;
    n1 = floor(n/2);
    alpha = 255;
    
    %read image
    I = imread('./triangle1.png');
    s=4;
    I = I(1:s:end,1:s:end,:);
    I = double(I);
    imshow(I,gray);
    [col1, row1] = size(I);
    
    % deleting area
    P1 = [320/s, 204/s];
    P2 = [320/s, 1003/s];
    P3 = [1012/s, 602/s];
    s = det([P1-P2;P3-P1]);
    J2 = I;
    C = ones(col1, row1); % Confidence term
    F = ones(col1, row1); % Filled term


    for j=1:row1
        for i=1:col1
            P = [i,j];
            if s*det([P3-P;P2-P3])>=0 && s*det([P1-P;P3-P1])>=0 && s*det([P2-P;P1-P2])>=0
                J2(i,j) = NaN;
                C(i,j) = 0;
                F(i,j) = 0;
            end
        end
    end 
    im = J2;
%    for jn = 1:100
    while sum(sum(1-F))>1
        extra =1;
        iter = iter+1;
        iter
        E = edge2(F);
        E2 = edge(F)-E;
        dx = [[1, 0, -1];[1, 0, -1];[1, 0, -1]];
        dy = [[1, 1, 1];[0, 0, 0];[-1, -1, -1]];
        Fx1 = conv2(dx, E2);
        Fy1 = conv2(dy, E2);
        Fx1 = Fx1(2:col1+1,2:row1+1);
        Fy1 = Fy1(2:col1+1,2:row1+1);

        [Ix, Iy] = gradient(im);
        Ix = Ix / alpha; Iy = Iy / alpha;
        Gx = Fy1.*Ix; Gy = Fx1.*Iy;
        Gx(isnan(Gx))=0; Gy(isnan(Gy))=0;

        Cp = confidence(C, n);
        Dp = sqrt(Gx.^2 + Gy.^2);
        Dp(isnan(Dp))=0;
        priority = Cp.*Dp;

        if sum(sum(priority)) == 0
            priority = Cp.*E;
            [y1, x1] = find(priority==max(max(priority)));
            res = randperm(length(y1),1);
            y = y1(res);
            x = x1(res);
            extra = 0;
        else
            [y1, x1] = find(priority==max(max(priority)));
            res = randperm(length(y1),1);
            y = y1(res);
            x = x1(res);

        end

        phi_p = phi(im, y, x, n);
        fill = phi(F, y, x, n); 
        antifill = 1 - fill;
        indices = find(isnan(phi_p) == 0);
        p_arr = phi_p(indices);


        err_start = 10^10;
        y2 = [];
        x2 = [];
        dist = [];
        if extra == 1
            for j = n1+1:row1-n1
                for i = n1+1:col1-n1
                    if all(all(phi(F, i, j, n)==1))
                        phi_q = phi(im, i, j, n); 
                        q_arr = phi_q(indices);
                        err = sum((q_arr-p_arr).^2,'all');
                        if err < err_start %&& (i-13)^2+(j-118)^2 > (n+2)^2 && (i-161)^2+(j-33)^2 > (n+2)^2 && (i-161)^2+(j-204)^2 > (n+2)^2 && (i-y)^2+(j-x)^2 > n^2
                            % up 25,235, left 322, 64 right 322,405
                            err_start = err;
                            y2 = i;
                            x2 = j;
                            dist = sqrt((y-i)^2 + (x-j)^2);
                        elseif err==err_start
                            y2 = [y2, i];
                            x2 = [x2, j];
                            dist = [dist, sqrt((y-i)^2 + (x-j)^2)]; 
                        end
                    end
                end
            end
            idy = find(dist==min(min(dist)));
            idy = idy(1);
            idy = randperm(length(y2),1);

            x1 = x2(idy);
            y1 = y2(idy);
            phi_q = phi(im, y1, x1, n); 

            phi_p(isnan(phi_p)) = 0;
            to_fill = phi_q.*antifill + phi_p.*fill;
            im(y-n1:y+n1,x-n1:x+n1) = to_fill;
            F(y-n1:y+n1,x-n1:x+n1) = 1;
            C(y-n1:y+n1,x-n1:x+n1) = C(y-n1:y+n1,x-n1:x+n1).*fill + Cp(y,x)*antifill;
        elseif extra == 0
             phi_p = phi(im, y, x, n);
             phi_p(isnan(phi_p)) = 0;

             fill = phi(F, y, x, n);
             antifill = 1 - fill;

             value = im(y,x);
             to_fill = phi_p.*fill + value * antifill;

             im(y-n1:y+n1,x-n1:x+n1) = to_fill;
             F(y-n1:y+n1,x-n1:x+n1) = 1;
             C(y-n1:y+n1,x-n1:x+n1) = C(y-n1:y+n1,x-n1:x+n1).*fill + Cp(y,x)*antifill;
        end


        if mod(iter,75) == 0
            close all
            R = im;
            G = im;
            R(isnan(R)) = 255;
            G(isnan(R)) = 0;
            B = G;
            Im(:,:,1) = R;
            Im(:,:,2) = G;
            Im(:,:,3) = B;
            Im = uint8(Im);
            fig_name = strcat('./Results/triangle_',num2str(n),'/triangle_no_dist',num2str(iter),'.tif');
            imwrite(Im,fig_name)
        end
        sum(sum(1-F))
    end
    close all
    R = im;
    G = im;
    R(isnan(R)) = 255;
    G(isnan(R)) = 0;
    B = G;
    Im(:,:,1) = R;
    Im(:,:,2) = G;
    Im(:,:,3) = B;
    Im = uint8(Im);
    figure
    imshow(Cp*255,gray)
    %imshow(Im)
    fig_name = strcat('./Results/triangle_',num2str(n),'/triangle_no_dist',num2str(iter),'.tif');
    imwrite(Im,fig_name)
end
%%



%%
R = im;
G = im;
R(isnan(R)) = 255;
G(isnan(R)) = 0;
B = G;
Im(:,:,1) = R;
Im(:,:,2) = G;
Im(:,:,3) = B;
Im = uint8(Im);
figure
imshow(Im)
fig_name = strcat('./Results/triangle_',num2str(11),'/triangle_',num2str(iter),'.tif');
saveas(gcf,fig_name)
%% Functions
function kernel = phi(Im, y, x, n)
    n1 = floor(n/2);
    kernel = Im(y-n1:y+n1,x-n1:x+n1);
end


function SSD = err_value(phi_p, phi_q, fill)
    phi_p = double(phi_p); phi_q = double(phi_q);
    SSD = sum(sum((phi_p.*fill - phi_q.*fill).^2));
end
function Cp = confidence(C, n)
    [l, k] = size(C);
    kernel = ones(n,n);
    n1 = floor(n/2);
    a = conv2(kernel,C);
    Cp = a(n1+1:n1+l, n1+1:n1+k)/sum(sum(kernel));
    Cp(1:4,1:end) = 0;
    Cp(end-4:end,1:end) = 0;
    Cp(1:end,1:4) = 0;
    Cp(1:end,end-4:end) = 0;
end
function [ I, theta ] = isophote(L)
% source: https://nl.mathworks.com/matlabcentral/fileexchange/48811-image-isophotes?s_tid=mwa_osa_a
    L = double(L);
    [Lx,Ly] = gradient(L);
    I = sqrt(Lx.^2+Ly.^2);
    theta = atan(Ly./Lx);
    %theta(isnan(theta))=0;
    I(isnan(I))=0;
end


function [x, y] = is_edge(f)
    [row, col] = size(f);
    x = [];
    y = [];
    for i = 1:row
        for j = 1:col
            if f(i, j) == 0 && abs(f(i-1,j)) + abs(f(i+1,j))+ abs(f(i,j-1)) + abs(f(i,j+1)) > 0
                x = [x, j];
                y = [y, i];
            end
        end    
    end
end

function result = edge2(f)
    [col, row] = size(f);
    edge = conv2(ones(3,3), abs(f));
    edge = edge(2:col+1,2:row+1);
    edge = edge ~=9;
    result = edge.*f;
    result = result./result;
    result(1,:) = NaN; result(:,1) = NaN; result(end,:) = NaN; result(:,end) = NaN;
    result(isnan(result))=0;
end

function result = edge(f)
    [col, row] = size(f);
    edge = conv2(ones(3,3), abs(f));
    edge = edge(2:col+1,2:row+1);
    edge1 = edge ~=9;
    filled2 = f - edge1;
    edge2 = conv2(ones(3,3), abs(filled2));
    edge2 = edge2(2:col+1,2:row+1);
    edge2 = edge2 ~=9;
    result = edge2.*f;
    result = result./result;
    result(1:3,:) = NaN; result(:,1:3) = NaN; result(end-2:end,:) = NaN; result(:,end-2:end) = NaN;
    result(isnan(result))=0;
end