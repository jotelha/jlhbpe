%% prepare figure
hFig1 = figure(1);
pm = 2;
pn = 4;

ax = cell(1,pm*pn);
for i = 1:(pn*pm)
    ax{i} = subplot(pm,pn,i);
end

%% set everything up, plot 1
n = 1;
f = -2.82;
c = 1; % mol
DeltaG0m = -26.28e3; % J / mol

theta = 0:0.001:1;

%plot(bockrisSwinkelTerm(theta,n));
%plot(bockrisSwinkelRHS(theta,c,f,DeltaG0m));

plot(ax{1},theta,log(bockrisSwinkelTerm(theta,n)),theta,log(bockrisSwinkelRHS(theta,c,f,DeltaG0m)));

%% plot 2
x0 = [0,1-eps];

c0 = 0:1:25;
theta0BockrisSwinkel = zeros(1,numel(c0));
for i = 1:numel(c0)
    theta0BockrisSwinkel(i) = fzero( @(th) bockrisSwinkelTerm(th,n) - bockrisSwinkelRHS(th,c0(i),f,DeltaG0m), x0 );
end

%figure(2);
plot(ax{2},c0,theta0BockrisSwinkel);

%% import one
files = dir('*1d*');
file = files(end).name;

[dat,delOut,headerOut] = importdata(file,' ',9);

X1 = dat.data(:,2);
Y1 = dat.data(:,1);
cDSm1 = dat.data(:,6);

%% plot 3, surface concentrations
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

%figure(3)
plot(ax{3},X2,cDSm2);

%% plot 4, Bockris Swinkel surface coverage
n = 10;
x0 = [0,1-eps];
theta0BockrisSwinkel = zeros(1,numel(cDSm2));
for i = 1:numel(cDSm2)
    theta0BockrisSwinkel(i) = fzero( @(th) bockrisSwinkelTerm(th,n) - bockrisSwinkelRHS(th,cDSm2(i),f,DeltaG0m), x0 );
end

plot(ax{4},X2,theta0BockrisSwinkel);

%% plot 5, Langmuir surface coverage
beta = 1;
% first import one
theta0Langmuir = beta.*cDSm2./(1+beta.*cDSm2);

plot(ax{5},X2,theta0BockrisSwinkel);

%% plot 6, adsorption reference
files = dir('zhang2015control_adsorption_mass*');
N = numel(files);
datRef = cell(N,1);

xRef = cell(N,1);
mRef = cell(N,1);

for i = 1:N
    file = files(i).name;

    [d,delOut,headerOut] = importdata(file,'\t',2);
    datRef{i} = d.data;
    
    xRef{i} = datRef{i}(:,1);
    mRef{i} = datRef{i}(:,2);
    
    plot(ax{6},xRef{i},mRef{i});
    hold(ax{6},'on');
end
hold(ax{6},'off');


%% concentrations, import all
x0 = [0,1-eps];
n = 1;
f = -2.82;
DeltaG0m = -26.28e3; % J / mol

beta = 1; % beta = e ^ (-DeltaG0i/RT)

files = dir('*1d*');
N = numel(files);
dat = cell(N,1);

xNum = cell(N,1);
cDSmNum = cell(N,1);

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
    
    theta0BockrisSwinkel = zeros(1,numel(cDSm2));
    for j = 1:numel(cDSm2)
        theta0BockrisSwinkel(j) = fzero( @(th) bockrisSwinkelTerm(th,n) - bockrisSwinkelRHS(th,cDSm2(j),f,DeltaG0m), x0 );
    end
    theta0Langmuir = beta.*cDSm2./(1+beta.*cDSm2);
    
    xNum{i} = X2;
    cDSmNum{i} = cDSm2;

    plot(ax{7},X2,theta0BockrisSwinkel);
    plot(ax{8},X2,theta0Langmuir);

    hold(ax{7},'on');
    hold(ax{8},'on');

end
hold(ax{7},'off');
hold(ax{8},'off');

%% set titles
t = {   'Bockris-Swinkel terms, LHS and RHS',...
        'Bockris-Swinkel surface coverage, linear concentration',...
        'Imported surface concentrations, one set',...
        'Bockris-Swinkel surface coverage, one set',...
        'Langmuir surface coverage, one set',...
        'Adsorption mass reference, zhang2015control',...
        'Bockris-Swinkel surface coverage, all',...        
        'Langmuir surface coverage, all' };


for i = 1:(pn*pm)
    title(ax{i},t{i});
end

%% second figure
hFig2 = figure(2);
pm = 1;
pn = 2;

ax = cell(1,pm*pn);
for i = 1:(pn*pm)
    ax{i} = subplot(pm,pn,i);
end
%% fit
N = numel(xNum);
selection = cell(N,1);
cDSmSel = cell(N,1);
xSel = cell(N,1);
mSel = cell(N,1);

for i=1:N
    selection{i} = (xNum{i} >= min(xRef{i})) & (xNum{i} <= max(xRef{i}));
    xSel{i} = xNum{i}(selection{i});
    cDSmSel{i} = cDSmNum{i}(selection{i});
    mSel{i} = interp1(xRef{i},mRef{i},xSel{i});
end

% ft = fittype('m0*beta*c/(1+beta*c)','independent','c');
% fo = fitoptions(ft);
fo = fitoptions('Method','NonlinearLeastSquares');
fo.TolFun=1e-12;
fo.TolX=1e-12;
fo.DiffMinChange=1e-12;
fo.StartPoint = [0.5,5e-7];
%fo.StartPoint = [1,1];
fo.Lower = [0 0];
fo.Robust = 'LAR';
fo.Display='iter';
ft = fittype('m0*beta*c/(1+beta*c)','independent','c','options',fo);

%fo.Robust = 'LAR';
langmuirFit = cell(N,1);
for i=1:N
    langmuirFit{i} = fit(cDSmSel{i},mSel{i},ft,fo);
    plot(ax{1},xSel{i},langmuirFit{i}(cDSmSel{i}));
    hold(ax{1},'on');
end
hold(ax{1},'off');

% xTot = xSel{
% for i=1:N
%     
%     langmuirFit{i} = fit(cDSmSel{i},mSel{i},ft,fo);
%     plot(ax{1},xSel{i},langmuirFit{i}(cDSmSel{i}));
%     hold(ax{1},'on');
% end
%% plot 2, adsorption reference, again
for i = 1:N
    plot(ax{2},xRef{i},mRef{i});
    hold(ax{2},'on');
end
hold(ax{2},'off');

