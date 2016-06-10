function obj = mesh1D(obj)
    % use 1st species to estimate analytical solution
%     z1 = obj.z(1,1); 
    z1 = max(obj.z);
    % RT = obj.RT;
    % UT = obj.UT; % thermal voltage
    accuracy = 1028;
%     ionicStrength = 1/2* sum(obj.z.^2.*obj.c_bulk); 
%     lambdaD = sqrt(constants.VacuumPermittivity*epsilon_r*RT/(2*F^2*ionicStrength));

    fprintf('  Ionic strength I = %g \n', obj.ionicStrength);
    fprintf('  Debye length lambdaD = %g m \n', obj.lambdaD);
    fprintf('  Domain L = %g m \n', obj.L);



    % those equations are NOT corrected by width of Stern layer
    % phiPB = @(x) 4*RT/(abs(z1)*F)*atanh( tanh( abs(z1)*F * E_m/(4*RT) )*exp(-x/lambdaD));
    % phiPBx = @(x) ( - 4*RT/( lambdaD*abs(z1)*F)* ...
    %         sinh( abs(z1) * F* E_m/(4*RT) ) * ...
    %         cosh( abs(z1) * F* E_m/(4*RT) ) ) ./ ...
    %     ( exp(x/lambdaD)* cosh( abs(z1) * F * E_m/(4*RT) )^2 ...
    %         - exp(-x/lambdaD) * sinh( abs(z1) * F* E_m/(4*RT) )^2);

    % dimensionless, with correction for Stern layer:
    % (x*L+lambdaS)/lambdaD = x/epsilon+delta
%     phi0 = obj.delta_phi;
    phi0 = obj.phi_bpe;
%     gamma = tanh( abs(z1) * phi0/4);
    oldDigits = digits(accuracy);
    gamma = tanh(vpa(phi0/4));
    offset = obj.delta;
%     offset = 0;

%     expTerm = @(x) exp(-((x+offset)/obj.epsilon)); % correction for position
    
%     phiPB = @(x) 4/(abs(z1))*atanh( gamma * expTerm);
    phiPB = @(x) 4*atanh( gamma * exp(-(x/obj.epsilon+offset)) );

%     phiPBx = @(x) ( - 4/obj.epsilon*abs(z1))* ...
%             sinh( abs(z1) * phi0/4 ) * ...
%             cosh( abs(z1) * phi0/4 )  ./ ...
%         ( exp((x/obj.epsilon+obj.delta))* cosh( abs(z1) * phi0/4 )^2 ...
%             - exp(-((x+obj.delta)/obj.epsilon)) * sinh( abs(z1) * phi0/4 )^2);
    phiPBx = @(x) ( - 4/obj.epsilon)* ...
            sinh( vpa(phi0/4) ) * ...
            cosh( vpa(phi0/4) )  ./ ...
        ( exp(x/obj.epsilon+offset)* cosh( vpa(phi0/4) )^2 ...
            - exp(-(x/obj.epsilon+offset)) * sinh( vpa(phi0/4) )^2);

    % estimate convection velocity vector    
    % consider convection-diffusion pde
    % ut + div*(-c*grad u) + beta* grad u = f
    % compare with NP equation
    % identify ||beta|| = |z| F D / (RT) * |grad phi|
    % and convection coefficient c = D
    % betaEstimate = @(x)abs(z1)*F*D1/RT * abs(phiPBx(x));

    % dimensionless approach:
    % div( - grad c - z*phix*c) = 0
    % ||beta|| = |z|*|phix| = |z|*|grad phi|
    % convection coefficient c = 1
    betaEstimate = @(x) abs(z1)*abs(phiPBx(x));

    %estimate Peclet number Pe = ||beta||*h/2c < 1 for stability 
    PeEstimate = @(x,h) betaEstimate(x).*h/2;

    %% create mesh according to Pe estimates
    % settings
    vertexSpacing = obj.epsilon/1000; % only on plots
    plotInterval = obj.epsilon;
    plotIntervalPe = obj.epsilon;
    hMaxFactor = obj.hMaxFactor; %1e-2; % how much smaller than the theoretical element size allowed by Peclet condition should the elements actually be?

    x0 = (0:vertexSpacing:1)';
    xnum = size(x0,1);
    hnum = 5;
    initExponent = log(obj.epsilon/100);
    h0 = exp(initExponent:-1:(initExponent-hnum+1));
    xarr = repmat(x0,1,hnum);
    harr = repmat(h0,xnum,1);
    % PeArr = PeEstimate(xarr,harr);

    xCur = 0;
    hCur = 0;
    v = [0];
    int = [];
    nv = 1;
    constantLength = false;
    failN = 1;
    % adjusting mesh intervals such that Peclet number does not exceed 1 anywhere
    opt = optimset('TolX',10^(-accuracy));
    while(xCur<1)
        PeStabilityCriterion = @(h) double( PeEstimate(xCur,h) ) - 1;
%         PeStabilityCriterionSquare = @(h) PeStabilityCriterion(h)^2;
        if ~constantLength
            try
                hMax = fzero(PeStabilityCriterion,[hCur,1],opt);
%                 hMax = fminbnd(PeStabilityCriterionSquare,hCur,1,opt);

                hCur = hMax*hMaxFactor; % choose tenth of actual maximum h to be safe
                % hCur = zero;
            catch
                fprintf('  Stopped growth at h = %g \n', hCur);
                constantLength = true;
                failN = nv;
            end 
        end
        xCur = xCur + hCur;
        nv = nv+1;
        v(nv) = xCur;
        int(nv-1) = hCur;
    end
%     v(nv) = 1; % correct last entry
    
    digits(oldDigits);

    fprintf('    Element size at surface h0 = %g m \n', int(1) );
    fprintf('    Maximum element size hL = %g m \n', hCur );
    fprintf('    Maximum element size reached at vertex %d \n', failN );
    fprintf('    Total number of vertices n = %d \n', nv );




    % select first Debye length
    obj.vFirstDebyeLength = v(v <= obj.epsilon);
    obj.nvFirstDebyeLength = numel(obj.vFirstDebyeLength);
    obj.vFirstDebyeLength = [obj.vFirstDebyeLength, v(obj.nvFirstDebyeLength+1)];
    adjustmentFactor = obj.epsilon / obj.vFirstDebyeLength(end);
    obj.intFirstDebyeLength = int(1:(obj.nvFirstDebyeLength-1))*adjustmentFactor;
    obj.vFirstDebyeLength = obj.vFirstDebyeLength*adjustmentFactor;
    fprintf('    Number of vertices within first Debye length above surface nD = %d \n', obj.nvFirstDebyeLength );

    %distributionFirstDebyeLength = int(2:nvFirstDebyeLength)./int(1:(nvFirstDebyeLength-1));
    %obj.distributionFirstDebyeLength = obj.intFirstDebyeLength(2:end)./obj.intFirstDebyeLength(1:(end-1));
    obj.distributionFirstDebyeLength = obj.vFirstDebyeLength(2:end)/obj.epsilon;
%     precision = ceil(abs(log10(min (obj.distributionFirstDebyeLength(2:end)-obj.distributionFirstDebyeLength(1:(end-1))))));
    precision = ceil(abs(log10(min(obj.distributionFirstDebyeLength)/100)));
    format = sprintf('%%%d.%df ',precision,precision);
%     obj.distributionFirstDebyeLengthStr = ['0 ',sprintf(format,obj.distributionFirstDebyeLength)];
    obj.distributionFirstDebyeLengthStr = ['0 ',sprintf(format,obj.distributionFirstDebyeLength)];


    %distributionRemaining = int((nvRemaining+1):end)./int(nbRemaining:(end-1));
    % distributionRemaining = intRemaining(2:end)./intRemaining(1:(end-1));
    % precision = ceil(abs(log10(min (distributionRemaining(2:end)-distributionRemaining(1:(end-1))))));
    % format = sprintf('%%%d.%df; ',precision,precision);
    % distributionRemainingStr = ['0; ',sprintf(format,distributionRemaining)];
    
    % selection for extended ddl
    obj.vExtendedDdl = v( (obj.epsilon < v) & (v <= (1+obj.extendedDdlFactor)*obj.epsilon) );
    obj.nvExtendedDdl = numel(obj.vExtendedDdl);
    obj.vExtendedDdl = [ obj.vExtendedDdl, v(obj.nvFirstDebyeLength+obj.nvExtendedDdl+1) ];
    adjustmentFactor = obj.extendedDdlFactor*obj.epsilon / (obj.vExtendedDdl(end)-obj.vExtendedDdl(1));
    obj.intExtendedDdl = int(obj.nvFirstDebyeLength:(obj.nvExtendedDdl+obj.nvFirstDebyeLength-1))*adjustmentFactor;
    obj.vExtendedDdl = (obj.vExtendedDdl - obj.vExtendedDdl(1))*adjustmentFactor+obj.epsilon;
    fprintf('    Number of vertices within extended DDL = %d \n', obj.nvExtendedDdl );
    obj.distributionExtendedDdl = (obj.vExtendedDdl(2:end)-obj.epsilon)/(obj.extendedDdlFactor*obj.epsilon);
    precision = ceil(abs(log10(min(obj.distributionExtendedDdl)/100)));
    format = sprintf('%%%d.%df ',precision,precision);
    obj.distributionExtendedDdlStr = ['0 ',sprintf(format,obj.distributionExtendedDdl)];
    
    obj.vRemaining = v(v > (obj.extendedDdlFactor+1)*obj.epsilon);
    obj.nvRemaining = numel(obj.vRemaining);
    adjustmentFactor = (1-(obj.extendedDdlFactor+1)*obj.epsilon ) / (obj.vRemaining(end)-obj.vRemaining(1));
    obj.intRemaining = int(obj.nvFirstDebyeLength:end)*adjustmentFactor;
    obj.vRemaining = (obj.vRemaining - obj.vRemaining(1))*adjustmentFactor+(1+obj.extendedDdlFactor)*obj.epsilon;

    fprintf('    Number of vertices outside first Debye length above surface nD = %d \n', obj.nvRemaining );  
%     obj.distributionRemaining = obj.intRemaining(2:end)./obj.intRemaining(1:(end-1));
    obj.distributionRemaining = (obj.vRemaining(2:end)-(obj.extendedDdlFactor+1)*obj.epsilon )/(1-(obj.extendedDdlFactor+1)*obj.epsilon );
%     precision = ceil(abs(log10(min (obj.distributionRemaining(2:end)-obj.distributionRemaining(1:(end-1))))));
    precision = ceil(abs(log10(min(obj.distributionRemaining)/100)));
    format = sprintf('%%%d.%df ',precision,precision);
    obj.distributionRemainingStr = ['0 ',sprintf(format,obj.distributionRemaining)];

     %% store in object
    % find entry closest to lambdaD and replace
%     v(obj.nvFirstDebyeLength) = obj.epsilon; % maybe better to insert another element?

%     obj.vMesh1D = v;
%     obj.distributionFirstDebyeLength = distributionFirstDebyeLength;
    
    %% add mesh to comsol model
    % 
%     edg = [0:(nv-2);1:(nv-1)]; % matlab starts to count at 1, COMSOL at 0
%     vtx = [0,nvFirstDebyeLength-1,nv-1];

%     % manualGeom = m.geom.create('manualGeom',1); % 1 for one dimension
%     % manualMesh = m.mesh().create('manualMesh','manualGeom');
%     %manualMesh.clearMesh();
%     manualMesh.data().clearData();
%     manualMesh.feature().clear(); % clear everything just in case
%     % manualMesh = m.mesh().create('manualMesh','geom');
%     manualMesh.data().setElem('edg',edg);
%     manualMesh.data().setElem('vtx',vtx); % geometry vertices, indices
%     manualMesh.data().setVertex(v); % mesh vertices, coordinates
%     manualMesh.data().createMesh();
% 
%     stationaryStudyStep1.set('mesh', {'geom' 'manualMesh'});

    %% show
    f = figure('visible',obj.plotsVisible);
    %f = figure();
    nRows = 3;
    nCols = 2;
    set(f,'Position',400*[0 0 nCols nRows/obj.widthToHeight]);

    for J = 1:nRows
        for K = 1:nCols
            i = (J-1)*nCols + K;
            subplot(nRows,nCols,i);
            switch J
                case 1
                    switch K
                        case 1
                            plot(x0(x0 < plotInterval),phiPB(x0(x0 < plotInterval))); % potential
                            t = 'Poisson-Boltzmann potential ditribution';
                            ylabel('phi / V');
                            xlabel('x / m');
                        case 2
                            plot(x0(x0 < plotInterval),phiPBx(x0(x0 < plotInterval))); %concentrations
                            t = 'Potential gradient';
                            ylabel('phi_x / V m^{-1}');
                            xlabel('x / m');
                    end
                case 2
                    switch K
                        case 1
                            plot(x0(x0 < plotInterval),log(abs(phiPB(x0(x0 < plotInterval))))); % potential
                            t = 'Poisson-Boltzmann log potential';
                            ylabel('log( abs( phi ) )');
                            xlabel('x / m');
                        case 2
                            plot(x0(x0 < plotInterval),log(abs(phiPBx(x0(x0 < plotInterval))))); %concentrations
                            t = 'Log potential gradient';
                            ylabel('log( abs( phi_x ) )');
                            xlabel('x / m');
                    end
                case 3
                    switch K
                        case 1
                            %hnum = 5;
                            %initExponent = log(lambdaD/1000);
                            %h0 = exp(initExponent:-1:(initExponent-hnum+1));
                            %xarr = repmat(x0,1,hnum);
                            %harr = repmat(h0,xnum,1);
                            %Pe = cell(1,hnum);
                            %[Pe{:}] = arrayfun(PeArray,h0);
                            plot(xarr(x0 < plotIntervalPe,:),PeEstimate(xarr(x0 < plotIntervalPe,:),harr(x0 < plotIntervalPe,:))); % potential
                            t = 'local Peclet numbers for constant element size h';
                            ylabel('Pe');
                            xlabel('x / m');
                            lgd = [repmat('h = ',hnum,1), num2str(h0',3), repmat(' m',hnum,1)];
                            legend( lgd );
                        case 2
                            plot(v(1:failN),int(1:failN),'LineStyle','none','Marker','*'); % mesh
                            t = 'safe mesh meeting Peclet criterion';
                            ylabel('h / m');
                            xlabel('x / m');
                    end
            end
            title(t);
        end
    end

%     set(f,'PaperPosition',5*[0 0 nCols nRows/widthToHeight]);
%     print(f,'-r300','-dpng','img/pecletStabilityMesh.png');
%     %saveas(f, 'img/convergencePlot.png'
    obj.savePlot(f,'pecletStabilityMesh.png',nRows,nCols);

end