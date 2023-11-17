%close all
clear all

global rn Kn alphain alphani  ri Ki


% model parameters
rn = 1.5;
Kn = 100;
alphain = 0.02;
alphani = 0.02;
ri = 1.5;
Ki = 60;
Ki_i = Ki;
Ki_1 = 60; 
Ki_2 = 90;

r_i = ri;
r_i_1 = 3.5;

Tfin = 30;

%original, for simulations
eps_o = 1e-6;
epsi_o = 0;


%Non-squared
f = @eq_ni;
t1 = 180;
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
figure(2);
subplot(2,3,1)
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
%set(gca,'xtick',[])
%set(gca,'ytick',[])
xlabel('N')
ylabel('I_0')



%% FIGURE 2, Increase in carrying capacity which increases the fertility as well


Ki = Ki_1; 
ri = r_i_1; 

%Non-squared
f = @eq_ni;
t1 = 180;
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

figure(2);
subplot(2,3,2)
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
%set(gca,'xtick',[])
%set(gca,'ytick',[])
title('(b) ')




%% More increase in carrying capacity

% model parameters
Ki = Ki_2;
ri = r_i_1;

%Non-squared
f = @eq_ni;
t1 = 180;
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

figure(2);
subplot(2,3,3)
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
%set(gca,'xtick',[])
%set(gca,'ytick',[])
title('(c) ')


%% Look at timeseries

global rm  Km  alphanm alphamn eps epsi 


% if choice = 1, Coexistence, but native in low density
% if choice = 2, Competitive exclusion

eps = eps_o;
epsi = epsi_o;


choice = 1 ;

ri = r_i;
Ki = Ki_i;

rm = ri;
alphanm = alphani;
alphamn = alphain;

if choice == 1
Km = Ki_1;
rm = r_i_1;
else
Km = Ki_2;
rm = r_i_1;
end



% To find initial conditions
eps = 0;
i0 = 8;
im0 = 0;
options = odeset('RelTol',1e-4,'AbsTol',1e-6);
[T,Yss] = ode45(@eq_nondim_K, 0:.1:5, [Kn; i0;0], options);
n0 = min(Yss(:,1));


eps = eps_o;
options = odeset('RelTol',1e-4,'AbsTol',1e-6);
[T,Y1] = ode45(@eq_nondim_K, 0:.1:Tfin, [n0; i0;im0], options);

figure(2)
subplot(2,2,3)
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
title('(d)')
axis([0 Tfin 0 110])
set(gca,'fontsize',14)




%% Complete displacement

choice = 2 ;

if choice == 1
Km = Ki_1;
else
Km = Ki_2;
end


% To find initial conditions
options = odeset('RelTol',1e-4,'AbsTol',1e-6);
[T,Yss] = ode45(@eq_nondim_K, 0:.1:5, [Kn; i0;0], options);
n0 = min(Yss(:,1));
eps = eps_o;

options = odeset('RelTol',1e-4,'AbsTol',1e-6);
[T,Y1] = ode45(@eq_nondim_K, 0:.1:Tfin, [n0; i0;im0], options);

figure(2)
subplot(2,2,4)
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
title('(e)')
axis([0 Tfin 0 110])
set(gca,'fontsize',14)



