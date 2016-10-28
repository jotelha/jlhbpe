%% set everything up, plot 1
n = 1;
f = -2.82;
c = 1; % mol
DeltaG0m = -26.28e3; % J / mol

theta = 0:0.001:1;

%plot(bockrisSwinkelTerm(theta,n));
%plot(bockrisSwinkelRHS(theta,c,f,DeltaG0m));
figure(1);
plot(theta,log(bockrisSwinkelTerm(theta,n)),theta,log(bockrisSwinkelRHS(theta,c,f,DeltaG0m)));

%% plot 2
x0 = [0,1-eps];

c0 = 0:1:25;
theta0 = zeros(1,numel(c0));
for i = 1:numel(c0)
    theta0(i) = fzero( @(th) bockrisSwinkelTerm(th,n) - bockrisSwinkelRHS(th,c0(i),f,DeltaG0m), x0 );
end

figure(2);
plot(c0,theta0);

%% import one
files = dir('*1d*');
file = files(end).name;

[dat,delOut,headerOut] = importdata(file,' ',9);

X1 = dat.data(:,2);
Y1 = dat.data(:,1);
cDSm1 = dat.data(:,6);

%% plot 3
minX = min(X1);
maxX = max(X1);
spanX = maxX - minX;
stepX = spanX / 1000;

X2 = minX:stepX:maxX;
Y2 = 0;

selection = (Y1 == 0);
X2 = X1(selection);
cDSm2 = cDSm1(selection);

%cDSm2 = interp2(X1,Y1,cDSm1,X2,Y2);

figure(3)
plot(X2,cDSm2);

%% plot 4
figure(4);
n = 10;
x0 = [0,1-eps];
theta0 = zeros(1,numel(cDSm2));
for i = 1:numel(cDSm2)
    theta0(i) = fzero( @(th) bockrisSwinkelTerm(th,n) - bockrisSwinkelRHS(th,cDSm2(i),f,DeltaG0m), x0 );
end

figure(4);
plot(X2,theta0);

%% import all
x0 = [0,1-eps];
n = 1;
f = -2.82;
DeltaG0m = -26.28e3; % J / mol

files = dir('*1d*');
N = numel(files);
dat = cell(N,1);
%X1 = 
for i = 1:N
    file = files(i).name;

    [d,delOut,headerOut] = importdata(file,' ',9);
    dat{i} = d.data;
    

    X1 = dat{i}(:,2);
    Y1 = dat{i}(:,1);
    cDSm1 = dat{i}(:,6);
    
    selection = (Y1 == 0);
    X2 = X1(selection);
    cDSm2 = cDSm1(selection);
    
    theta0 = zeros(1,numel(cDSm2));
    for j = 1:numel(cDSm2)
        theta0(j) = fzero( @(th) bockrisSwinkelTerm(th,n) - bockrisSwinkelRHS(th,cDSm2(j),f,DeltaG0m), x0 );
    end
    figure(6);
    plot(X2,theta0);
    hold on;
end
hold off;