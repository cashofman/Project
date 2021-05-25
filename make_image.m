clear all
col = 2^10;
row = col;
Im = 255 * ones(row,col);

%% outer triangle
x1 = floor(col/5);
x2 = 4*floor(col/5);
x3 = floor(col/2);
y1 = 4*floor(row/5)-80;
y2 = 4*floor(row/5)-80;
y3 = floor(y1 - sqrt((x2-x1)^2-(x3-x1)^2));

P1 = [x1, y1];
P2 = [x2, y2];
P3 = [x3, y3];
s = det([P1-P2;P3-P1]);
for i=1:col
    for j=1:row
        P = [i, j];
        if s*det([P3-P;P2-P3])>=0 && s*det([P1-P;P3-P1])>=0 && s*det([P2-P;P1-P2])>=0
            Im(j,i) = 0;
        end
    end
end
%% inner triangle
x = 50;
y = x*cos(pi/3);
z = x*sin(pi/3);
dy = sqrt(z^2+(x+y)^2);
x1 = x1 + (x + y);
x2 = x2 - (x + y);
x3 = x3;
y1 = y1 - z;
y2 = y2 - z;
y3 = y3 + dy;

P1 = [x1, y1];
P2 = [x2, y2];
P3 = [x3, y3];
s = det([P1-P2;P3-P1]);
for i=1:col
    for j=1:row
        P = [i, j];
        if s*det([P3-P;P2-P3])>=0 && s*det([P1-P;P3-P1])>=0 && s*det([P2-P;P1-P2])>=0
            Im(j,i) = 255;
        end
    end
end

%% circle

x = x1;
y = y3+166/2;
r = 80;

for i=1:col
    for j=1:row
        if sqrt((x-i)^2+(y-j)^2)<r
            Im(j,i) = 0;
        end
    end
end
dipshow(Im)
%%
x = x2;
y = y3+166/2;


for i=1:col
    for j=1:row
        if sqrt((x-i)^2+(y-j)^2)<r
            Im(j,i) = 0;
        end
    end
end

%%
x = x3;
y = y1+156;

for i=1:col
    for j=1:row
        if sqrt((x-i)^2+(y-j)^2)<r
            Im(j,i) = 0;
        end
    end
end

dipshow(Im)