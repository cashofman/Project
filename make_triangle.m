a = 1;
b = a * cos(60*pi/180);
c = a * sin(60*pi/180);
lim = 1.5 * a;

x = [0, c, -c, 0, c];
y = [a, -b, -b, a, -b];
x2 = [-c, c, 0];
y2 = [b, b, -a];

plot(x,y,'k','LineWidth',3)
hold on
plot(x2, y2, '.k', 'MarkerSize',200)
set(gca,'Color','white','xticklabel',[],'yticklabel',[])%setting background black and removing labels
xlim([-lim lim])
ylim([-lim lim])
daspect([1 1 1]);%equal data unit length along x and y axis
