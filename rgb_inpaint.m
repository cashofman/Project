%% start
clear all
close all
clc
%% load image
i = './crimsi_images/jelle.jpeg';
I = imread(i);
s=3;
I = I(1:s:end,1:s:end,1:3);
DI = double(I);
imshow(I);
[col1, row1, dim1] = size(I);
%% selecting area
% source: https://nl.mathworks.com/matlabcentral/answers/155561-how-to-replace-manually-selected-part-of-the-image-with-background-color
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
%% deleting cropped area
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
%% Save cropped image
R = im(:,:,1);
G = im(:,:,2);
B = im(:,:,3);

R(isnan(R)) = 255;
G(isnan(G)) = 0;
B(isnan(B)) = 0;

Im(:,:,1) = R;
Im(:,:,2) = G;
Im(:,:,3) = B;

%imwrite(uint8(Im),'baby_cropped.tif')
imshow(uint8(Im))
%% first iteration
iter = 0;
n = 9;
n1 = floor(n/2);
alpha = 255;
priority = 1;
im2 = zeros(2*n+col1, 2*n+row1, dim1);
F2 = ones(2*n+col1, 2*n+row1);
C2 = ones(2*n+col1, 2*n+row1);
im2(n+1:end-n,n+1:end-n,:) = im;
F2(n+1:end-n,n+1:end-n) = F;
C2(n+1:end-n,n+1:end-n) = C;
C = C2;
F = F2;
im = im2;

imshow(uint8(im2))
%%

while sum(sum(1-F))>1
    iter = iter+1
    
    % calculate the edge and the normal vector of the edge
    E = edge2(F);
    E2 = edge(F)-E;
    dx = [[1, 0, -1];[1, 0, -1];[1, 0, -1]];
    dy = [[1, 1, 1];[0, 0, 0];[-1, -1, -1]];
    Fx = conv2(dx, E2);
    Fy = conv2(dy, E2);
    Fx1 = Fx(2:end-1,2:end-1);
    Fy1 = Fy(2:end-1,2:end-1);
    
    % calculate the image isophote
    [Ix, Iy] = gradient(double(im),100);
    Ix = Ix / alpha; Iy = Iy / alpha;
    Gx = Fy1.*Ix; Gy = Fx1.*Iy;
    Gx(isnan(Gx))=0; Gy(isnan(Gy))=0;
    
    % calculate the priority 
    Cp = confidence(C, n);
    Dp = sqrt(Gx.^2 + Gy.^2);
    Dp(isnan(Dp))=0;
    priority = Cp.*Dp;
    priority = sum(priority,3);
    
    % find the highest priority pixel
    [y1, x1] = find(priority==max(max(priority)));
    res = randperm(length(y1),1);
    y = y1(res);
    x = x1(res);
  
    % extract the patch which needs to be filled from the image
    phi_p = phi(im, y, x, n); 
    fill = phi(F, y, x, n); 
    antifill = 1 - fill;
    indices = find(isnan(phi_p) == 0);
    p_arr = phi_p(indices);
    
    % find the patch(es) with the lowest error
    err_start = 10^10;
    y2 = [];
    x2 = [];
    for j = n+n1+1:row1-n-n1
        for i = n+n1+1:col1-n-n1
            if all(all(phi(F, i, j, n)==1))
                phi_q = phi(im, i, j, n); 
                q_arr = phi_q(indices);
                err = sum((q_arr-p_arr).^2,'all');
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
    
    % chose the exemplar patch
    idy = randperm(length(y2),1);
    x1 = x2(idy);
    y1 = y2(idy);
    phi_q = phi(im, y1, x1, n); 
    
    % fill the image, update the fill and confidence matrix
    phi_p(isnan(phi_p)) = 0;
    to_fill = phi_q.*antifill + phi_p.*fill;
    im(y-n1:y+n1,x-n1:x+n1, :) = to_fill;
    F(y-n1:y+n1,x-n1:x+n1) = 1;
    C(y-n1:y+n1,x-n1:x+n1) = C(y-n1:y+n1,x-n1:x+n1)*fill + Cp(y,x)*antifill;
    sum(sum(1-F))
    
    % save intermediat result
    if mod(iter,75) == 0
        fig_name = strcat('./Results/jelle',num2str(n),'iter_',num2str(iter),'.tif');
        R = im(:,:,1);
        G = im(:,:,2);
        B = im(:,:,3);
        
        R(isnan(R)) = 255;
        G(isnan(G)) = 0;
        B(isnan(B)) = 0;
        
        R = R(n+1:end-n,n+1:end-n);
        G = G(n+1:end-n,n+1:end-n);
        B = B(n+1:end-n,n+1:end-n);
        
        Im(:,:,1) = R;
        Im(:,:,2) = G;
        Im(:,:,3) = B;
        imwrite(uint8(Im),fig_name)
        imshow(uint8(Im))
    end
end
im = im(n+1:end-n,n+1:end-n, :);
figure
imshow(uint8(im))
%% Functions
function kernel = phi(Im, y, x, n)
% extracts the area from an image around a given x,y coordinate given a
% patch size n
    n1 = floor(n/2);
    kernel = Im(y-n1:y+n1,x-n1:x+n1,:);
end

function Cp = confidence(C, n)
% calculates the confidence value of a pixel given the confidence
% matrix of the image
    [l, k] = size(C);
    n1 = floor(n/2);
    kernel = ones(n,n);
    kernel(n1+1,n1+1) = 0;
    a = conv2(kernel,C);
    Cp = a(n1+1:n1+l, n1+1:n1+k)/sum(sum(kernel));
    Cp(1:4,1:end) = 0;
    Cp(end-4:end,1:end) = 0;
    Cp(1:end,1:4) = 0;
    Cp(1:end,end-4:end) = 0;
end

function result = edge2(f)
% calculates the edge directly next to the unfilled pixels
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
% calculates the edge directly next to the unfilled pixels, and the pixels
% next to the edge
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