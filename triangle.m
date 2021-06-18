%% start
clear all
close all
clc
%% load image
I = imread('./triangle1.png');
s=6;
I = I(1:s:end,1:s:end,:);
I = double(I);
imshow(I,gray);

[col1, row1] = size(I);
%% inpainting
n_range = [9, 13, 15];
alpha = 255;
for k = 1:numel(n_range)
    n = n_range(k);
    iter = 0;
    n1 = floor(n/2);
    
    %read image
    I = imread('./triangle1.png');
    s=3;
    I = I(1:s:end,1:s:end,:);
    I = double(I);
    imshow(I,gray);
    [col1, row1] = size(I);
    
    % Corner points of triangle
    P1 = [320/s, 204/s];
    P2 = [320/s, 1003/s];
    P3 = [1018/s, 602/s];
    s = det([P1-P2;P3-P1]);
    
    % create image for reconstruction, confidence matrix and filled matrix
    J2 = I;
    C = ones(col1, row1); % Confidence term
    F = ones(col1, row1); % Filled term
    
    % delete image inside triangle
    for j=1:row1
        for i=1:col1
            P = [i,j];
            if s*det([P3-P;P2-P3])>=0 && s*det([P1-P;P3-P1])>=0 && s*det([P2-P;P1-P2])>=0 % inside triangle
                J2(i,j) = NaN; % set pixel value to NaN
                C(i,j) = 0; % set confidence value to 0
                F(i,j) = 0; % set filled value to 0
            end
        end
    end
    
    im = J2;
    while sum(sum(1-F))>1
        extra = 1;
        iter = iter+1;
        iter
        
        % calculate the edge and the normal vector of the edge
        E = edge2(F);
        E2 = edge(F)-E;
        dx = [[1, 0, -1];[1, 0, -1];[1, 0, -1]];
        dy = [[1, 1, 1];[0, 0, 0];[-1, -1, -1]];
        Fx1 = conv2(dx, E2);
        Fy1 = conv2(dy, E2);
        Fx1 = Fx1(2:col1+1,2:row1+1);
        Fy1 = Fy1(2:col1+1,2:row1+1);

        % calculate the image isophote
        [Ix, Iy] = gradient(im);
        Ix = Ix / alpha; Iy = Iy / alpha;
        Gx = Fy1.*Ix; Gy = Fx1.*Iy;
        Gx(isnan(Gx))=0; Gy(isnan(Gy))=0;

        % calculate the priority 
        Cp = confidence(C, n);
        Dp = sqrt(Gx.^2 + Gy.^2);
        Dp(isnan(Dp))=0;
        priority = Cp.*Dp;
        
        % find the highest priority, if priority is 0 only look at the
        % confidence term
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

        % extract the patch which needs to be filled from the image
        phi_p = phi(im, y, x, n);
        fill = phi(F, y, x, n); 
        antifill = 1 - fill;
        indices = find(isnan(phi_p) == 0);
        p_arr = phi_p(indices);

        % find the exemplar patch(es) with the lowest error
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
            % chose the closest patch
            idy = find(dist==min(min(dist)));
            idy = idy(1);
            
            % construct the exemplar patch
            x1 = x2(idy);
            y1 = y2(idy);
            phi_q = phi(im, y1, x1, n); 
            
            % fill the image, update the fill and confidence matrix
            phi_p(isnan(phi_p)) = 0;
            to_fill = phi_q.*antifill + phi_p.*fill;
            im(y-n1:y+n1,x-n1:x+n1) = to_fill;
            F(y-n1:y+n1,x-n1:x+n1) = 1;
            C(y-n1:y+n1,x-n1:x+n1) = C(y-n1:y+n1,x-n1:x+n1).*fill + Cp(y,x)*antifill;
            
        elseif extra == 0
            % if priority is 0, fill the target pixel with its surrounding
            % pixels
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

        % save intermediate results
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
            fig_name = strcat('./Results/triangle_',num2str(n),'/test',num2str(iter),'.tif');
            imwrite(Im,fig_name)
            imshow(Im)
        end
        sum(sum(1-F))
    end
    
    close all
    Im(:,:,1) = R;
    Im(:,:,2) = G;
    Im(:,:,3) = B;
    Im = uint8(Im);
    figure
    fig_name = strcat('./Results/triangle_',num2str(n),'/test',num2str(iter),'.tif');
    imwrite(Im,fig_name)
end

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