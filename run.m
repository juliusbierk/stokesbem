clc
phi = 0;
tic()
if system('g++ bem.cpp -larmadillo -fopenmp')==0
    disp('Calculating mesh and boundary conditions.')
    create_system(phi);
    system('./a.out');

    filename = 'velocity.txt';
    delimiter = ' ';
    formatSpec = '%f%f%f%f%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, ...
                'MultipleDelimsAsOne', true,  'ReturnOnError', false);
    fclose(fileID);
    X = dataArray{:, 1};
    Y = dataArray{:, 2};
    Z = dataArray{:, 3};
    UX = dataArray{:, 4};
    UY = dataArray{:, 5};
    UZ = dataArray{:, 6};
    clearvars filename delimiter formatSpec fileID dataArray ans;
end

%% Plot
load system
omega = 50; % Hz
addpath '../'
Xo = X;
Yo = Y;
Zo = Z;
U = UX;
V = UY;
W = UZ;

% Meaner
x = [];
y = [];
u = [];
v = [];
while numel(X)>0
    logger = (X==X(1)) & (Y==Y(1)) & (Z>-1) & (Z<1);
    logger2 = (X==X(1)) & (Y==Y(1));
    x(end+1) = X(1);
    y(end+1) = Y(1);
    u(end+1) = mean(U(logger));
    v(end+1) = mean(V(logger));
    X(logger2) = [];
    Y(logger2) = [];
    Z(logger2) = [];
    U(logger2) = [];
    V(logger2) = [];
    W(logger2) = [];
end

%%
clf()
colormap(flipud(parula))
xx = min(x):0.2:max(x);
yy = min(y):0.2:max(y);
[X,Y] = meshgrid(xx,yy);
U = griddata(x,y,u,X,Y);
V = griddata(x,y,v,X,Y);
radii = logspace(-0.8,0.8,15);
contourf(X,Y,2./sqrt(omega*sqrt(U.^2+V.^2)),radii) % 50 hertz
colorbar()
hold on

xx = min(x):0.5:max(x);
yy = min(y):0.5:max(y);
[X,Y] = meshgrid(xx,yy);
U = griddata(x,y,u,X,Y);
V = griddata(x,y,v,X,Y);
patch('faces',mt,'vertices',mp,'facecol',cont_colour(0.5),'edgecol','k','facealpha',1)
plot3(sx(:,1),sx(:,2),sx(:,3),'.','Color',[0,0,0])
quiver3(X,Y,X*0,U,V,U*0,3)

%quiver3(Xo,Yo,Zo,UX,UY,UZ,0)
%quiver3(sx(:,1),sx(:,2),sx(:,3),sv(:,1),sv(:,2),sv(:,3),0)
axis([-5,15,-10,10])
axis equal
%print('-dpng',sprintf('imgs/%05i.png',index))
title('Radius at which Pe = 1')
xlabel('x [\mum]')
ylabel('y [\mum]')

X = Xo;
Y = Yo;
Z = Zo;
