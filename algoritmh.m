%% initialise dip
addpath('/Users/cashofman/Desktop/Master/ADIP/dip/common/dipimage')
dip_initialise
addpath '/Users/cashofman/Desktop/Master/ADIP/dip/common/dipimage'
%% start
clear all
close all
clc
%% load image
I = imread('/Users/cashofman/Desktop/Master/ADIP/Project/triangle_4.png');
I = I(:,:,1);
s=5;
I = I(1:s:end,1:s:end);
J = mat2im(I);
dipshow(I);
[col1, row1] = imsize(J);
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
% P1 = [180*3/4, 138*3/4];
% P2 = [180*3/4, 360*3/4];
% P3 = [415*3/4, 250*3/4];
% s = det([P1-P2;P3-P1]);
%s*det([P3-P;P2-P3])>=0 && s*det([P1-P;P3-P1])>=0 && s*det([P2-P;P1-P2])>=0
J2 = I;
C = ones(row1, col1); % Confidence term
F = ones(row1, col1); % Filled term


for j=1:col1
    for i=1:row1
%         P = [i,j];
        if binaryImage(i,j) == 1
            J2(i,j) = NaN;
            C(i,j) = 0;
            F(i,j) = 0;
        end
    end
end 

I3= mat2im(J2);
dipshow(I3)
im = J2;
%%
n = 7;
alpha = 255;
%[X, Y] = is_edge(F); % all edge points
[F1, theta1] = isophote(F); % normal vector of the edge
[F2, theta2] = isophote(im); % isophote of image
theta1 = theta1*-1;
theta2 = theta2;
Cp = confidence(C, n);
E2 = edge2(F);


%%
priority = E2.*Cp.*abs(F2/alpha.*((cos(theta2).*cos(theta1)+sin(theta2).*sin(theta1))));

dipshow(priority,'base')
%dipshow(E2.*theta2,'angle')
%% first iteration
n = 13;
n1 = floor(n/2);
alpha = 255;


%for m = 1:10
while sum(sum(1-F))>1

    
    %[X, Y] = is_edge(F); % all edge points
    [F1, theta1] = isophote(F); % normal vector of the edge
    [F2, theta2] = isophote(im); % isophote of image
    E = edge2(F);
    Cp = confidence(C, n);
    priority = E.*Cp.*abs(F2/alpha.*((cos(theta2).*cos(theta1)+sin(theta2).*sin(theta1))));

    %P = diag(priority(Y,X));
    %dipshow(mat2im(priority),'base')
    % res1 = find(P==max(P));
    % res = randperm(length(res1),1);
    [y1, x1] = find(priority==max(max(priority)));
    res = randperm(length(y1),1);
    y = y1(res);
    x = x1(res);
    % x = X(res);
    % y = Y(res);

    phi_p = im(y-n1:y+n1,x-n1:x+n1);
    fill = F(y-n1:y+n1,x-n1:x+n1);
    antifill = 1 - fill;

    err_start = 10^10;

    y2 = [];
    x2 = [];
    for j = n1+1:col1-n1
        for i = n1+1:row1-n1
            if all(all(F(i-n1:i+n1,j-n1:j+n1))) == 1
                phi_q = im(i-n1:i+n1,j-n1:j+n1);
                err = err_value(phi_p, phi_q, fill);
                if err < err_start
                    err_start = err;
                    y2 = i;
                    x2 = j;
                elseif err==0 
                    y2 = [y2, i];
                    x2 = [x2, j];

                end
            end
        end
    end
    idy = randperm(length(y2),1);
    x1 = x2(idy);
    y1 = y2(idy);
    phi_q = im(y1-n1:y1+n1,x1-n1:x1+n1);


    to_fill = double(phi_q).*antifill + double(im(y-n1:y+n1,x-n1:x+n1)).*fill;
    im(y-n1:y+n1,x-n1:x+n1) = uint8(to_fill);
    F(y-n1:y+n1,x-n1:x+n1) = 1;
    C(y-n1:y+n1,x-n1:x+n1) = Cp(y-n1:y+n1,x-n1:x+n1);
    sum(sum(1-F))
end
sum(sum(1-F))
close all   
% E = edge(F);
% dipshow(E,'base')
% dipshow(priority,'base')
dipshow(mat2im(im))

%%
subplot(2, 2, 1);
imshow(I);
subplot(2, 2, 2);
imshow(J2);
subplot(2, 2, 3);
imshow(im);
subplot(2, 2, 4);
imshow(C);


%% Functions
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
    theta(isnan(theta))=0;
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
    [row, col] = size(f);
    edge = conv2(ones(3,3), abs(f));
    edge = edge(2:row+1,2:col+1);
    edge = edge ~= 9;
    result = edge.*f;
    result = result./result;
    result(1,:) = NaN; result(:,1) = NaN; result(end,:) = NaN; result(:,end) = NaN;
end
