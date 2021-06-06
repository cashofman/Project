%% start
clear all
close all
clc
%% load image
I = imread('./waterval.jpeg');
s=1;
I = I(1:s:end,1:s:end,1:3);
DI = double(I);
imshow(I);
%I = double(I);
[col1, row1, dim1] = size(I);
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
J2 = DI;
C = ones(col1, row1); % Confidence term
F = ones(col1, row1); % Filled term
for j=1:row1
    for i=1:col1
        if binaryImage(i,j) == 1
            J2(i,j,:) = NaN;
            C(i,j) = 0;
            F(i,j) = 0;
        end
    end
end 
figure
imshow(uint8(J2))
im = J2;
%% first iteration
%im = double(im);
E0 = edge2(F);
n = 9;
n1 = floor(n/2);
alpha = 255;
E1 = edge2(F);
priority = 1;
% for i = 1:50
%     i
while sum(sum(1-F))>1
%     iter = iter+1;
%     [F1, theta1] = isophote(E); % normal vector of the edge
%     [F2, theta2] = isophote(im); % isophote of image
%     theta1 = theta1 * -1; 
%     F2 = F2/alpha;

    E = edge2(F);
    E2 = edge(F)-E;
    dx = [[1, 0, -1];[1, 0, -1];[1, 0, -1]];
    dy = [[1, 1, 1];[0, 0, 0];[-1, -1, -1]];
    Fx1 = conv2(dx, E2);
    Fy1 = conv2(dy, E2);
    Fx1 = Fx1(2:col1+1,2:row1+1);
    Fy1 = Fy1(2:col1+1,2:row1+1);
    
    [Fx, Fy] = gradient(E);
    [Ix, Iy] = gradient(double(im),100);
    Ix = Ix / alpha; Iy = Iy / alpha;
    Gx = Fy1.*Ix; Gy = Fx1.*Iy;
    Gx(isnan(Gx))=0; Gy(isnan(Gy))=0;

    I1 = sqrt(Ix.^2 + Iy.^2);

    Cp = confidence(C, n);
    Dp = sqrt(Gx.^2 + Gy.^2);
    Dp(isnan(Dp))=0;??
    priority = Cp.*Dp;
    priority = sum(priority,3);
    
    if sum(sum(priority)) == 0
        priority = E2.*Cp;
        [y1, x1] = find(priority==max(max(priority)));
        res = randperm(length(y1),1);
        y = y1(res);
        x = x1(res);
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
    for j = n1+1:row1-n1
        for i = n1+1:col1-n1
            if all(all(phi(F, i, j, n)==1))
                phi_q = phi(im, i, j, n); 
                q_arr = phi_q(indices);
                err = sum((q_arr-p_arr).^2,'all');
                if err < err_start 
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
    phi_q = phi(im, y1, x1, n); %im(y1-n1:y1+n1,x1-n1:x1+n1);
    
    phi_p(isnan(phi_p)) = 0;
    to_fill = phi_q.*antifill + phi_p.*fill;
    im(y-n1:y+n1,x-n1:x+n1, :) = to_fill;
    F(y-n1:y+n1,x-n1:x+n1) = 1;
    C(y-n1:y+n1,x-n1:x+n1) = Cp(y,x);
    sum(sum(1-F))
end
sum(sum(1-F))
close all   
E = E0 + E1 + edge2(F);
figure()
imshow(uint8(im))
figure()
imshow(E*255,gray)
figure()
imshow(priority*255)
%%
figure
imshow(uint8(im))

%% Functions
function kernel = phi(Im, y, x, n)
    n1 = floor(n/2);
    kernel = Im(y-n1:y+n1,x-n1:x+n1,:);
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
    result(1,:) = NaN; result(:,1) = NaN; result(end,:) = NaN; result(:,end) = NaN;
    result(isnan(result))=0;
end