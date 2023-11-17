close all
clear all

global rn Kn alphain alphani  ri Ki alphamn alphanm rm eps epsi Km


% model parameters
rn = 1.5;
Kn = 100;
alphain = 0.2;
alphani = 0.01;
alphanm = alphani;
alphamn = alphain;
ri = rn;
rm = ri;
Ki = 20;
Km = 80;




%Non-squared
f = @eq_ni;
t1 = 150;
t2 = 100;
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
subplot(1,2,1)
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





%% More increase in carrying capacity

% model parameters
Ki = Km;

%Non-squared
f = @eq_ni;
t1 = 150;
t2 = 100;
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
subplot(1,2,2)
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


%% Look at timeseries: microbiome exchange occurs versus not.


% no microbiome exchange

Ki = 20;
Tfin = 20;
n0 = Kn;
i00 = 1;
im0 = 0;

eps = 0;
epsi = 0;

options = odeset('RelTol',1e-4,'AbsTol',1e-6);
[T,Yext] = ode45(@eq_nondim_K, 0:.1:Tfin, [n0; i00;im0], options);
%microbiome exchange
eps = 1e-6;
[T,Ymic] = ode45(@eq_nondim_K, 0:.1:Tfin, [n0; i00;im0], options);


figure(2)
%subplot(1,2,1)
plot(T,Ymic(:,1),'-','Color','#EDB120','linewidth',2)
hold on
%plot(T,Yext(:,1),'-','Color','#EDB120','linewidth',2)
%hold on
plot(T,Ymic(:,3)+Ymic(:,2),'-','Color','#0072BD','linewidth',2)
hold on
plot(T,Yext(:,2)+Yext(:,3),'--','Color','#0072BD','linewidth',2)
hold on
%legend({'I_0','I_m'}, 'Location','west', 'FontSize',12)
xlabel('time [arbitrary units]')
ylabel('Population')
title('(a)')
%axis([0 Tfin 0 110])
legend({'N','no microbiome exchange','microbiome exchange'},'Location','west','FontSize',12)
set(gca,'fontsize',14)
axis([0 Tfin 0 110])


