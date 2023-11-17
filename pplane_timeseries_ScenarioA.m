close all
clear all

global rn Kn alphain alphani  ri Ki


% model parameters
rn = 1.5;
Kn = 100;
alphain = 0.02;
alphani = 0.01;
ri = 1.5;
Ki = 20;
Ki_i = Ki;
Ki_1 = 70;
Ki_2 = 90;

%Non-squared
f = @eq_ni;
t1 = 160;
t2 = 120;
na = 15;
y1 = linspace(0,t1,na); % from zero to t1, 30 arrows
y2 = linspace(0,t2,na);
[x,y] = meshgrid(y1, y2);
%size(x)
%size(y)

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


%% pre-invasion
figure(1);
subplot(3,3,1)
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
title('(a) ')
set(gca,'xtick',[])
set(gca,'ytick',[])
xlabel('N')
ylabel('I_0')



%% FIGURE 2, INCREASE A BIT THE CARRYING CAPACITY


Ki = Ki_1;

%Non-squared
f = @eq_ni;
t1 = 160;
t2 = 120;
y1 = linspace(0,t1,na); % from zero to t1, 30 arrows
y2 = linspace(0,t2,na);
[x,y] = meshgrid(y1, y2);

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
subplot(3,3,2)
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
ylabel('I_m')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('(b) ')




%% Increase more the carrying capacity: Invaders displace natives

% model parameters
Ki = Ki_2;

%Non-squared
f = @eq_ni;
t1 = 160;
t2 = 120;
y1 = linspace(0,t1,na); % from zero to t1, 30 arrows
y2 = linspace(0,t2,na);
[x,y] = meshgrid(y1, y2);

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
subplot(3,3,3)
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
ylabel('I_m')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('(c) ')


%% Look at timeseries

global rm  Km  alphanm alphamn eps epsi 


% if choice = 1, Coexistence, but native in low density
% if choice = 2, Competitive exclusion

choice = 1 ;

Ki = Ki_i;

rm = ri;
alphanm = alphani;
alphamn = alphain;

if choice == 1
Km = Ki_1;
else
Km = Ki_2;
end


eps1 = 1e-6;
epsi1 = 0;

Tfin = 60;

n0 = 20;
i00 = 2;
im0 = 0;

% Choose ss for initial conditions

eps = 0;
epsi = 0;

options = odeset('RelTol',1e-4,'AbsTol',1e-6);
[T,Yss] = ode45(@eq_nondim_K, 0:.1:Tfin, [n0; i00;0], options);

nss = Yss(end,1);
i0ss = Yss(end,2);

eps = eps1;
epsi = epsi1;


n0 = nss;
i00 = i0ss;
im0 = 0;

options = odeset('RelTol',1e-4,'AbsTol',1e-6);
[T,Y1] = ode45(@eq_nondim_K, 0:.1:Tfin, [nss; i0ss;im0], options);

figure(1)
subplot(3,3,[4,5,6])
plot(T,Y1(:,1),'Color','#EDB120', 'linewidth',2)
hold on
plot(T,Y1(:,2),'--','Color','#0072BD','linewidth',1)
hold on
plot(T,Y1(:,3),':','Color','#0072BD','linewidth',1)
hold on
Ym = Y1(:,3)+Y1(:,2);
plot(T,Ym,'Color','#0072BD','linewidth',2)
hold on
legend({'N','I_0','I_m'}, 'Location','west', 'FontSize',12)
xlabel('time')
ylabel('Population')
%title('1 plant, X fungi')
axis([0 Tfin 0 110])
set(gca,'fontsize',14)




%% Complete displacement

choice = 2 ;

if choice == 1
Km = Ki_1;
else
Km = Ki_2;
end


eps1 = 1e-6;
epsi1 = 0;



n0 = 20;
i00 = 2;

% Choose ss for initial conditions

eps = 0;
epsi = 0;

options = odeset('RelTol',1e-4,'AbsTol',1e-6);
[T,Yss] = ode45(@eq_nondim_K, 0:.1:Tfin, [n0; i00;0], options);

nss = Yss(end,1);
i0ss = Yss(end,2);

eps = eps1;
epsi = epsi1;


n0 = nss;
i00 = i0ss;
im0 = 0;

options = odeset('RelTol',1e-4,'AbsTol',1e-6);
[T,Y1] = ode45(@eq_nondim_K, 0:.1:Tfin, [nss; i0ss;im0], options);

figure(1)
subplot(3,3,[7,8,9])
plot(T,Y1(:,1),'Color','#EDB120', 'linewidth',2)
hold on
plot(T,Y1(:,2),'--','Color','#0072BD','linewidth',1)
hold on
plot(T,Y1(:,3),':','Color','#0072BD','linewidth',1)
hold on
Ym = Y1(:,3)+Y1(:,2);
plot(T,Ym,'Color','#0072BD','linewidth',2)
hold on
legend({'N','I_0','I_m'}, 'Location','west', 'FontSize',12)
xlabel('time')
ylabel('Population')
%title('1 plant, X fungi')
%axis([0 Tfin 0 max(max(Ym))+10])
axis([0 Tfin 0 110])
set(gca,'fontsize',14)



