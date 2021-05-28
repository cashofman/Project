%% start
clear all
close all
clc
%% load image
I = imread('./mountaineer.jpg');
%I = double(I);
I = rgb2gray(I);
%I = I(100:end-50,1:end-80);
s=6;
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
P1 = [112, 54];
P2 = [112, 418];
P3 = [432, 234];
s = det([P1-P2;P3-P1]);
J2 = I;
C = ones(col1, row1); % Confidence term
F = ones(col1, row1); % Filled term


for j=1:row1
    for i=1:col1
         P = [i,j];
        if binaryImage(i,j) == 1
            %s*det([P3-P;P2-P3])>=0 && s*det([P1-P;P3-P1])>=0 && s*det([P2-P;P1-P2])>=0
            %binaryImage(i,j) == 1
            J2(i,j) = NaN;
            C(i,j) = 0;
            F(i,j) = 0;
        end
    end
end 


imshow(J2)
im = J2;
iter = 0;
%%
E0 = edge2(F);
%% first iteration
n = 9;
n1 = floor(n/2);
alpha = 255;
E1 = edge2(F);

while sum(sum(1-F))>1
    iter = iter+1;
    E = edge2(F);
    [F1, theta1] = isophote(E); % normal vector of the edge
    [F2, theta2] = isophote(im); % isophote of image
    theta1 = theta1 * -1; 
    F2 = F2/alpha;
    [Fx, Fy] = gradient(E);
    [Ix, Iy] = gradient(im);
    Ix = Ix / alpha; Iy = Iy / alpha;
    Gx = Fx.*Iy; Gy = Fy.*Ix;
    Gx(isnan(Gx))=0; Gy(isnan(Gy))=0;
    
    Cp = confidence(C, n);
    Dp = sqrt(Gx.^2 + Gy.^2);%F2.*abs(cos(theta2).*cos(theta1)+sin(theta2).*sin(theta1));
    Dp(isnan(Dp))=0;
    priority = Cp.*Dp;
    
    if sum(sum(priority)) == 0
        [y1, x1] = find(E==max(max(E)));
        res = randperm(length(y1),1);
        y = y1(res);
        x = x1(res);
    else
        [y1, x1] = find(priority==max(max(priority)));
        res = randperm(length(y1),1);
        y = y1(res);
        x = x1(res);
    end

    phi_p = phi(im, y, x, n); %im(y-n1:y+n1,x-n1:x+n1);
    fill = phi(F, y, x, n); 
    antifill = 1 - fill;
    indices = find(isnan(phi_p) == 0);
    p_arr = phi_p(indices);
    
    
    err_start = 10^10;
    y2 = [];
    x2 = [];
    for j = n1+1:row1-n1
        for i = n1+1:col1-n1
            if all(all(phi(F, i, j, n)==1))
                phi_q = phi(im, i, j, n); 
                q_arr = phi_q(indices);
                err = sum((q_arr-p_arr).^2);
                if err < err_start
                    err_start = err;
                    y2 = i;
                    x2 = j;
                elseif err==err_start
                    y2 = [y2, i];
                    x2 = [x2, j];
                end
            end
        end
    end
    idy = randperm(length(y2),1);
    x1 = x2(idy);
    y1 = y2(idy);
    phi_q = phi(im, y1, x1, n); %im(y1-n1:y1+n1,x1-n1:x1+n1);
    
    phi_p(isnan(phi_p)) = 0;
    to_fill = phi_q.*antifill + phi_p.*fill;
    im(y-n1:y+n1,x-n1:x+n1) = to_fill;
    F(y-n1:y+n1,x-n1:x+n1) = 1;
    C(y-n1:y+n1,x-n1:x+n1) = Cp(y,x);
    sum(sum(1-F))
end
sum(sum(1-F))
close all   
E = E0 + E1 + edge2(F);
figure()
imshow(im,gray)
figure()
imshow(E*255,gray)
figure()
imshow(priority*255,gray)
%%
E = edge2(F);
imshow(F2);
sum(sum(priority))

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
function result = edge(f)
    [row, col] = size(f);
    edge = conv2(ones(3,3), abs(f));
    edge = edge(2:row+1,2:col+1);
    f1 = f==0;
    result = edge.*f1;
    result = result./result;
end

function result = edge2(f)
    [col, row] = size(f);
    edge = conv2(ones(3,3), abs(f));
    edge = edge(2:col+1,2:row+1);
    edge = edge ~= 9;
    result = edge.*f;
    result = result./result;
    result(1,:) = NaN; result(:,1) = NaN; result(end,:) = NaN; result(:,end) = NaN;
    result(isnan(result))=0;
end

