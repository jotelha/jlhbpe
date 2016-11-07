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
pm = 3;
pn = 6;

ax = cell(1,pm*pn);
for i = 1:(pn*pm)
    ax{i} = subplot(pm,pn,i);
end
%% fit Langmuir, single
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
%fo.Display='iter';
fo.Display='final';
ft = fittype('m0*beta*c/(1+beta*c)','independent','c','options',fo);

%fo.Robust = 'LAR';
langmuirFit = cell(N,1);
for i=1:N
    langmuirFit{i} = fit(cDSmSel{i},mSel{i},ft,fo);
%     plot(ax{1},xSel{i},langmuirFit{i}(cDSmSel{i}));
%     hold(ax{1},'on');
end
% hold(ax{1},'off');

% xTot = xSel{
% for i=1:N
%     
%     langmuirFit{i} = fit(cDSmSel{i},mSel{i},ft,fo);
%     plot(ax{1},xSel{i},langmuirFit{i}(cDSmSel{i}));
%     hold(ax{1},'on');
% end

%% fit Langmuir, all
tmp = cellfun(@(c) c' ,cDSmSel,'UniformOutput',false);
cDSmTot = [tmp{:}]';

tmp = cellfun(@(c) c' ,mSel,'UniformOutput',false);
mTot = [tmp{:}]';

langmuirTotFit = fit(cDSmTot,mTot,ft,fo);
% for i=1:N
%     plot(ax{2},xSel{i},langmuirTotFit(cDSmSel{i}));
%     hold(ax{2},'on');
% end

%% Langmuir plots
hFig1 = figure(1);
pm = 1;
pn = 2;

ax = cell(1,pm*pn);
for i = 1:(pn*pm)
    ax{i} = subplot(pm,pn,i);
end


for i=1:N
    plot(ax{1},xSel{i}(1:20:end),mSel{i}(1:20:end),'o');
    hold(ax{1},'on');
    plot(ax{1},xSel{i},langmuirFit{i}(cDSmSel{i}));
end
hold(ax{1},'off');


for i=1:N
    plot(ax{2},xSel{i}(1:20:end),mSel{i}(1:20:end),'o');
    hold(ax{2},'on');
    plot(ax{2},xSel{i},langmuirTotFit(cDSmSel{i}));
end
hold(ax{2},'off');

%% plot 3, adsorption reference, again
for i = 1:N
    plot(ax{3},xRef{i},mRef{i});
    hold(ax{3},'on');
end
hold(ax{3},'off');

%% fit Tamkin, single

% ft = fittype('m0*beta*c/(1+beta*c)','independent','c');
% fo = fitoptions(ft);
fo = fitoptions('Method','NonlinearLeastSquares');
fo.TolFun=1e-12;
fo.TolX=1e-12;
fo.DiffMinChange=1e-12;
fo.StartPoint = [max(mTot)/2,max(mTot)/(2*max(cDSmTot))];
%fo.StartPoint = [1,1];
fo.Lower = [-Inf eps];
fo.Robust = 'LAR';
% fo.Display='iter';
fo.Display='final';
ft = fittype('-A+log(c)/F','independent','c','options',fo);

%fo.Robust = 'LAR';
tamkinFit = cell(N,1);
for i=1:N
    tamkinFit{i} = fit(cDSmSel{i},mSel{i},ft,fo);
    plot(ax{4},xSel{i},tamkinFit{i}(cDSmSel{i}));
    hold(ax{4},'on');
end
hold(ax{4},'off');

%% fit Tamkin, all

tamkinTotFit = fit(cDSmTot,mTot,ft,fo);
for i=1:N
    plot(ax{5},xSel{i},tamkinTotFit(cDSmSel{i}));
    hold(ax{5},'on');
end

%% fit Temkin full, single

fo = fitoptions('Method','NonlinearLeastSquares');
% fo.TolFun=1e-12;
% fo.TolX=1e-12;
% fo.DiffMinChange=1e-12;
fo.StartPoint = [1/max(mTot),1,1];
%fo.StartPoint = [1,1];
fo.Lower = [-Inf -Inf -Inf];
fo.Upper = [Inf Inf Inf];
fo.Robust = 'LAR';
fo.Display='iter';
%fo.Display='final';
ft = fittype('-1/f*log((1+beta1*c)/(1+beta0*c))','independent','c','options',fo);

%fo.Robust = 'LAR';
fullTemkinFit = cell(N,1);
for i=1:N
    fullTemkinFit{i} = fit(cDSmSel{i},mSel{i},ft,fo);
%     plot(ax{7},xSel{i},unnamedFit{i}(cDSmSel{i}));
%     hold(ax{7},'on');
end
% hold(ax{7},'off');

fullTemkinTotFit = fit(cDSmTot,mTot,ft,fo);


%% full Temkin plots

hFig2 = figure(2);
pm = 1;
pn = 2;

ax = cell(1,pm*pn);
for i = 1:(pn*pm)
    ax{i} = subplot(pm,pn,i);
end


for i=1:N
    plot(ax{1},xSel{i}(1:20:end),mSel{i}(1:20:end),'o');
    hold(ax{1},'on');
    plot(ax{1},xSel{i},fullTemkinFit{i}(cDSmSel{i}));
end
hold(ax{1},'off');


for i=1:N
    plot(ax{2},xSel{i}(1:20:end),mSel{i}(1:20:end),'o');
    hold(ax{2},'on');
    plot(ax{2},xSel{i},fullTemkinTotFit(cDSmSel{i}));
end
hold(ax{2},'off');

%% fit Bockris-Swinkel, single
beta0 = exp(26280/(jlh.Constants.R*jlh.Constants.T))/55.5;
betaInv0 = 1/beta0;
f0 = -2.82;
n0 = 10;

% ft = fittype('betaInv*exp(-f*theta/theta0)*(theta/theta0/(1-theta/theta0)^n)*((theta/theta0+n*(1-theta/theta0))^(n-1)/n^n)','independent','theta','problem','theta0','options',fo);
bockrisSwinkelModel = @(betaInv,f,n,theta) betaInv*exp(-f*theta).*(theta./(1-theta).^n).*((theta+n*(1-theta)).^(n-1)/n^n);

%fo.Robust = 'LAR';
bockrisSwinkelFit = cell(N,1);
for i=1:N
    fo = fitoptions('Method','NonlinearLeastSquares');
    fo.StartPoint = [betaInv0 2*max(mSel{i})];
    fo.Lower = [eps max(mSel{i})];
    fo.Upper = [max(cDSmSel{i}) Inf];
    fo.Robust = 'LAR';
    fo.Display='iter';
    %fo.Display='final';
    ft = fittype(@(betaInv,m0,m) bockrisSwinkelModel(betaInv,f0,n0,m/m0),'independent','m');
    bockrisSwinkelFit{i} = fit(mSel{i},cDSmSel{i},ft,fo);
end    
bockrisSwinkelFit2 = cell(N,1);
for i=1:N
    fo = fitoptions('Method','NonlinearLeastSquares');
    % fo = fitoptions;
    % fo.TolFun=1e-12;
    % fo.TolX=1e-12;
    % fo.DiffMinChange=1e-12;
    % fo.StartPoint = [max(mTot)/2,1/(2*max(cDSmTot)),1/(2*max(cDSmTot))];
    fo.StartPoint = [-2.82,1];
    fo.Lower = [-10 0];
    fo.Upper = [5 20];
    fo.Robust = 'LAR';
    fo.Display='iter';
    %fo.Display='final';
    ft = fittype(@(f,n,m) bockrisSwinkelModel(bockrisSwinkelFit{i}.betaInv,f,n,m/bockrisSwinkelFit{i}.m0),'independent','m');
    bockrisSwinkelFit2{i} = fit(mSel{i},cDSmSel{i},ft,fo);
end
for i=1:N
    fo = fitoptions('Method','NonlinearLeastSquares');
    fo.StartPoint = [bockrisSwinkelFit{i}.betaInv bockrisSwinkelFit{i}.m0 bockrisSwinkelFit2{i}.f bockrisSwinkelFit2{i}.n];
    fo.Lower = [eps, max(mSel{i}), -10, 0];
    fo.Upper = [Inf Inf 5 20];
    fo.Robust = 'LAR';
    fo.Display='iter';
    %fo.Display='final';
    ft = fittype(@(betaInv,m0,f,n,m) bockrisSwinkelModel(betaInv,f,n,m/m0),'independent','m');
    bockrisSwinkelFit3{i} = fit(mSel{i},cDSmSel{i},ft,fo);
end 


m0BockrisSwinkel = cell(3,1);
for i=1:N
    x0 = [0+eps,bockrisSwinkelFit3{i}.m0-eps];
    for j = 1:numel(cDSmSel{i})
        m0BockrisSwinkel{i}(j) = fzero( @(m0) bockrisSwinkelFit3{i}(m0)-cDSmSel{i}(j), x0 );
    end 
end 
%% Bockris figure
hFig2 = figure(3);
pm = 1;
pn = 3;

ax = cell(1,pm*pn);
for i = 1:(pn*pm)
    ax{i} = subplot(pm,pn,i);
end

for i=1:N
    plot(ax{1},xSel{i},bockrisSwinkelFit{i}(mSel{i}));
    hold(ax{1},'on');
end
hold(ax{1},'off');


for i=1:N
    plot(ax{2},xSel{i}(1:20:end),mSel{i}(1:20:end),'o');
        hold(ax{2},'on');
    plot(ax{2},xSel{i},m0BockrisSwinkel{i});
end
hold(ax{2},'off');


%% fit unknown, all

bockrisSwinkelTotFit = fit(0.9*mTot/max(mTot),cDSmTot,ft,fo);
for i=1:N
    plot(ax{11},xSel{i},bockrisSwinkelTotFit(cDSmSel{i}));
    hold(ax{11},'on');
end
hold(ax{11},'off');
%% Flory-Huggins

%fo.Robust = 'LAR';
floryFit = cell(N,1);
for i=1:N

    fo = fitoptions('Method','NonlinearLeastSquares');
    %fo.TolFun=min(cDSmTot)/100;
    %fo.TolX=min(mTot);
    %fo.DiffMinChange=min(cDSmTot)/100;
    %fo.StartPoint = [0.5,5e-7];
    fo.StartPoint = [max(cDSmSel{i}), max(mSel{i}), 2];
    fo.Lower = [eps max(mSel{i})+eps 0];
    fo.Upper = [Inf Inf 20];

    %fo.Robust = 'LAR';
    %fo.MaxIter = 100;
    %fo.MaxFunEvals = 500;
    fo.Display='iter';
    %fo.Display='final';
    ft = fittype('betaInv*m/(m0*(1-m/m0)^(n))*exp(1-n)','independent','m','options',fo);

    floryFit{i} = fit(mSel{i},cDSmSel{i},ft,fo);
    %plot(ax{4},xSel{i},floryFit{i}(mSel{i}));
    %hold(ax{4},'on');
end

m0Flory = cell(3,1);
for i=1:N
    x0 = [0+eps,floryFit{i}.m0-eps];
    for j = 1:numel(cDSmSel{i})
        m0Flory{i}(j) = fzero( @(m0) floryFit{i}(m0)-cDSmSel{i}(j), x0 );
    end 
end 
%hold(ax{4},'off');

%% Fory-Huggins, all

%fo.Robust = 'LAR';
floryTotFit = fit(mTot,cDSmTot,ft,fo);
m0FloryTot = cell(3,1);
x0 = [0+eps,floryTotFit.m0-eps];
for i=1:N
    for j = 1:numel(cDSmSel{i})
        m0FloryTot{i}(j) = fzero( @(m0) floryTotFit(m0)-cDSmSel{i}(j), x0 );
    end 
end 

%% Flory figure
hFig4 = figure(4);
pm = 1;
pn = 3;

ax = cell(1,pm*pn);
for i = 1:(pn*pm)
    ax{i} = subplot(pm,pn,i);
end

for i=1:N
    plot(ax{1},xSel{i},floryFit{i}(mSel{i}));
    hold(ax{1},'on');
end
hold(ax{1},'off');


for i=1:N
    plot(ax{2},xSel{i}(1:20:end),mSel{i}(1:20:end),'o');
    hold(ax{2},'on');
    plot(ax{2},xSel{i},m0Flory{i});
end
hold(ax{2},'off');

for i=1:N
    plot(ax{3},xSel{i}(1:20:end),mSel{i}(1:20:end),'o');
    hold(ax{3},'on');
    plot(ax{3},xSel{i},m0FloryTot{i});
end
hold(ax{2},'off');


%% plot 6, adsorption reference, again
for i = 1:N
    plot(ax{6},mSel{i},cDSmSel{i});
    hold(ax{6},'on');
end
hold(ax{3},'off');

%% plot 6, adsorption reference, again
for i = 1:N
    plot(ax{9},mSel{i},log(cDSmSel{i}));
    hold(ax{9},'on');
end
hold(ax{3},'off');