close all
clear all

global rn Kn alphain alphani  ri Ki


% model parameters
rn = 1.5;
Kn = 100;
alphain = 0.01;
alphani = 0.01;
ri = 1.5;
Ki = 100;

%Non-squared
f = @eq_ni;
t1 = 160;
t2 = 160;
na = 15;
y1 = linspace(0,t1,na); % from zero to t1, 30 arrows
y2 = linspace(0,t2,na);
[x,y] = meshgrid(y1, y2);
size(x)
size(y)

u = zeros(size(x));
v = zeros(size(y));

t = 0;
for i = 1:numel(x)
    Yprime = f(t,[x(i);y(i)]);
    u(i) = Yprime(1);
    v(i) = Yprime(2);
end

syms un vn

un = u./sqrt(u.^2 + v.^2);
vn = v./sqrt(v.^2 + v.^2);

%plot trajectories

figure(1);
subplot(2,2,1)
quiver(x, y, un, vn, 0.7,'Color',  [0.8 0.8 0.8],'linewidth',1.5 )
%axis equal tight;
syms n i0;
hold on
fimplicit(@(n,i0) ri*(1-i0/Ki)-alphani*n, [0 t1 0 t2],'Color','#0072BD', 'linewidth',2)
hold on
fimplicit(@(n,i0) rn*(1-n/Kn)-alphain*i0, [0 t1 0 t2],'Color','#EDB120','linewidth',2)
hold on
axis([ 0 t1 0 t2])
hold on
y10 = 10;
x10 = 20;
[ts,ys] = ode45(@eq_ni,0:0.1:1000,[x10;y10]);
%    plot(ys(:,1),ys(:,2),'r-','linewidth',1.5)
%plot(x10,y10,'r.','MarkerSize',20)
plot(ys(end,1),ys(end,2),'s','MarkerSize',15,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6])
set(gca,'fontsize',14)
title('(a) Coexistence')
set(gca,'xtick',[])
set(gca,'ytick',[])
xlabel('N')
ylabel('I_0')


%% FIGURE 2, COMPETITIVE EXCLUSION
% model parameters
rn = 1.5;
Kn = 100;
alphain = 0.025;
alphani = 0.025;
ri = 1.5;
Ki = 100;

%Non-squared
f = @eq_ni;
t1 = 110;
t2 = 110;
y1 = linspace(0,t1,na); % from zero to t1, 30 arrows
y2 = linspace(0,t2,na);
[x,y] = meshgrid(y1, y2);
size(x)
size(y)

u = zeros(size(x));
v = zeros(size(y));

t = 0;
for i = 1:numel(x)
    Yprime = f(t,[x(i);y(i)]);
    u(i) = Yprime(1);
    v(i) = Yprime(2);
end

syms un vn

un = u./sqrt(u.^2 + v.^2);
vn = v./sqrt(v.^2 + v.^2);

%plot trajectories

figure(1);
subplot(2,2,2)
quiver(x, y, un, vn, 0.7,'Color',  [0.8 0.8 0.8],'linewidth',1.5 )
%axis equal tight;
syms n i0;
hold on
fimplicit(@(n,i0) ri*(1-i0/Ki)-alphani*n, [0 t1 0 t2],'Color','#0072BD', 'linewidth',2)
hold on
fimplicit(@(n,i0) rn*(1-n/Kn)-alphain*i0, [0 t1 0 t2],'Color','#EDB120','linewidth',2)
hold on
axis([ 0 t1 0 t2])
hold on
y10 = 10;
x10 = 20;
[ts,ys] = ode45(@eq_ni,0:0.1:1000,[x10;y10]);
%    plot(ys(:,1),ys(:,2),'r-','linewidth',1.5)
%plot(x10,y10,'r.','MarkerSize',20)
plot(ys(end,1),ys(end,2),'s','MarkerSize',15,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6])
y10 = 20;
x10 = 10;
[ts,ys] = ode45(@eq_ni,0:0.1:1000,[x10;y10]);
%    plot(ys(:,1),ys(:,2),'r-','linewidth',1.5)
%plot(x10,y10,'r.','MarkerSize',20)
plot(ys(end,1),ys(end,2),'s','MarkerSize',15,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6])
set(gca,'fontsize',14)
xlabel('N')
ylabel('I_0')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('(b) Competitive exclusion (bistability)')


%% FIGURE 3, NATIVE ONLY
% model parameters
rn = 1.5;
Kn = 90;
alphain = 0.02;
alphani = 0.02;
ri = 1.3;
Ki =50;

%Non-squared
f = @eq_ni;
t1 = 100;
t2 = 100;
y1 = linspace(0,t1,na); % from zero to t1, 30 arrows
y2 = linspace(0,t2,na);
[x,y] = meshgrid(y1, y2);
size(x)
size(y)

u = zeros(size(x));
v = zeros(size(y));

t = 0;
for i = 1:numel(x)
    Yprime = f(t,[x(i);y(i)]);
    u(i) = Yprime(1);
    v(i) = Yprime(2);
end

syms un vn

un = u./sqrt(u.^2 + v.^2);
vn = v./sqrt(v.^2 + v.^2);

%plot trajectories

figure(1);
subplot(2,2,3)
quiver(x, y, un, vn, 0.7,'Color',  [0.8 0.8 0.8],'linewidth',1.5 )
%axis equal tight;
syms n i0;
hold on
fimplicit(@(n,i0) ri*(1-i0/Ki)-alphani*n, [0 t1 0 t2],'Color','#0072BD', 'linewidth',2)
hold on
fimplicit(@(n,i0) rn*(1-n/Kn)-alphain*i0, [0 t1 0 t2],'Color','#EDB120','linewidth',2)
hold on
axis([ 0 t1 0 t2])
hold on
y10 = 10;
x10 = 20;
[ts,ys] = ode45(@eq_ni,0:0.1:1000,[x10;y10]);
%    plot(ys(:,1),ys(:,2),'r-','linewidth',1.5)
%plot(x10,y10,'r.','MarkerSize',20)
plot(ys(end,1),ys(end,2),'s','MarkerSize',15,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6])
y10 = 20;
x10 = 10;
[ts,ys] = ode45(@eq_ni,0:0.1:1000,[x10;y10]);
%    plot(ys(:,1),ys(:,2),'r-','linewidth',1.5)
%plot(x10,y10,'r.','MarkerSize',20)
plot(ys(end,1),ys(end,2),'s','MarkerSize',15,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6])
set(gca,'fontsize',14)
xlabel('N')
ylabel('I_0')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('(c) Natives only')




%% FIGURE 4, INVASIVE ONLY
% model parameters
rn = 1.2;
Kn = 80;
alphain = 0.015;
alphani = 0.03;
ri = 3.4;
Ki =110;

%Non-squared
f = @eq_ni;
t1 = 120;
t2 = 120;
y1 = linspace(0,t1,na); % from zero to t1, 30 arrows
y2 = linspace(0,t2,na);
[x,y] = meshgrid(y1, y2);
size(x)
size(y)

u = zeros(size(x));
v = zeros(size(y));

t = 0;
for i = 1:numel(x)
    Yprime = f(t,[x(i);y(i)]);
    u(i) = Yprime(1);
    v(i) = Yprime(2);
end

syms un vn

un = u./sqrt(u.^2 + v.^2);
vn = v./sqrt(v.^2 + v.^2);

%plot trajectories

figure(1);
subplot(2,2,4)
quiver(x, y, un, vn, 0.7,'Color',  [0.8 0.8 0.8],'linewidth',1.5 )
%axis equal tight;
syms n i0;
hold on
fimplicit(@(n,i0) ri*(1-i0/Ki)-alphani*n, [0 t1 0 t2],'Color','#0072BD', 'linewidth',2)
hold on
fimplicit(@(n,i0) rn*(1-n/Kn)-alphain*i0, [0 t1 0 t2],'Color','#EDB120','linewidth',2)
hold on
axis([ 0 t1 0 t2])
hold on
y10 = 10;
x10 = 20;
[ts,ys] = ode45(@eq_ni,0:0.1:1000,[x10;y10]);
%    plot(ys(:,1),ys(:,2),'r-','linewidth',1.5)
%plot(x10,y10,'r.','MarkerSize',20)
plot(ys(end,1),ys(end,2),'s','MarkerSize',15,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6])
y10 = 20;
x10 = 10;
[ts,ys] = ode45(@eq_ni,0:0.1:1000,[x10;y10]);
%    plot(ys(:,1),ys(:,2),'r-','linewidth',1.5)
%plot(x10,y10,'r.','MarkerSize',20)
plot(ys(end,1),ys(end,2),'s','MarkerSize',15,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6])
set(gca,'fontsize',14)
xlabel('N')
ylabel('I_0')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('(d) Invaders only')
