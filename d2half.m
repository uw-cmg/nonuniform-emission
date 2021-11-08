function [TVJRCTL,result] = d2half(I)
    tic;
    S = inputstandardize(I);
    Ar = S.Ar;
    Vcathode = S.Vcathode;
    zlist = S.zlist;
    Schottky = S.Schottky;
    VFermi = S.VFermi;

    % constants
    % h = 6.62607015e-34;
    kB = 1.380649e-23;
    m = 9.10938356e-31;
    epsilon0 = 8.854187817e-12;
    elec = 1.602176634e-19;
    A = 1.20173e6;

    % Mesh
    xnum = S.meshsize(1);
    ynum = S.meshsize(2);
    znum = length(zlist);
    Z=repmat(permute(zlist,[3 2 1]),[xnum ynum 1]);

    % Initial U as U0
    Vanode = S.TVanode(1,2);
    if isfield(S,'U0')
        U = S.U0;
    else
        U = repmat(Vcathode,[1 1 znum])+(Vanode-Vcathode)./zlist(end).*Z;
    end
    U(:,:,end)= Vanode;
    
    if Schottky
        USchottkyz = permute(elec/(16*pi*epsilon0)./zlist,[3 2 1]);
        U = U + USchottkyz;
    end

    z_multip = [0;(zlist(2:end-1)-zlist(1:end-2)).*(zlist(3:end)-zlist(2:end-1));0];
    Thomas_upper  = [0;2*(zlist(2:end-1)-zlist(1:end-2))./(zlist(3:end)-zlist(1:end-2));0];
    Thomas_lower  = [0;2*(zlist(3:end)  -zlist(2:end-1))./(zlist(3:end)-zlist(1:end-2));0];
    halfxnum = floor(xnum/2);
    z_multip = repmat(permute(z_multip,[3 2 1]),[halfxnum+1 ynum 1]);

    xterm = 2/(S.step(1)^2)*(cos(2*pi/xnum*(0:halfxnum)')-1); % x half range
    yterm = 2/(S.step(2)^2)*(cos(2*pi/ynum*(0:ynum-1)')-1);
    Thomas_diag = -2+z_multip.*(repmat(xterm,[1 ynum znum])+repmat(yterm',[halfxnum+1 1 znum]));

    fftU = nan(size(U));

    TVJRCTL = nan(size(S.TVanode,1),7);
    % 1, T: temperature
    % 2, V: Vanode
    % TVJRCTL(:,1:2) = S.TVanode;
    % 3, J: 2.5-D calculated current density.
    % 4, R: Richardson-Dushman equation.
    % 5, C: Space charge: 1-D Child-Langmuir law.
    % TVJRCTL(:,5) = 4*epsilon0/9*sqrt(2*elec/m)/(zlist(end)^2)*S.TVanode(:,2).^(3/2);
    % 6, T: 1-D Child-Langmuir law with finite temperature correction.
    % 7, L: Number of iteration of loopi

    newconv_cri1 = min(S.conv_cri(1), Vanode*S.conv_cri(2)); % in V.

    for rowi = 1:size(S.TVanode,1)
        T = S.TVanode(rowi,1);
        Vanode = S.TVanode(rowi,2);
        JCL0 = 4*epsilon0/9*sqrt(2*elec/m)/(zlist(end)^2)*mean((Vanode-Vcathode).^(3/2),[1,2]);
        JCLT = 4*epsilon0/9*sqrt(2*elec/m)/(zlist(end)^2)*mean((Vanode-Vcathode).^(3/2).*TcalCL(T,Vanode-Vcathode,S.CLTintegral),[1,2]);
        
        kBT = kB*T;
        eoverkBT = elec/kBT;
        c = Ar.*sqrt(pi/2)/epsilon0/sqrt(kBT/m)*A*T^2;
        
        uminmin_est = min(min(Vcathode,[],'all'), log(JCLT/(A*T^2))/eoverkBT);
        newconv_cri2 = min(S.conv_cri(2), newconv_cri1/(-uminmin_est));
        
        U(U>0) = Vanode/mean(U(:,:,end),'all') * U(U>0);
        
        loopconv = false;
        for loopi = 1:S.looplimit
            if Schottky
                Ecathode = (U(:,:,2)-USchottkyz(1,1,2)-Vcathode)/(zlist(2)-zlist(1));
                U1Schottky  = Vcathode + sqrt(max(0,Ecathode*elec/(4*pi*epsilon0)));
                zSchottky   = zlist(1) + sqrt(max(0,elec/(16*pi*epsilon0)./Ecathode));
                U2 = U(:,:,2);
                U1Schottky(zSchottky>zlist(2)) = U2(zSchottky>zlist(2));
                U1Schottky(Ecathode<=0) = Vcathode(Ecathode<=0)+USchottkyz(1,1,2);
                U(:,:,1) = U1Schottky;
            end
            
            Uminl = U;       
            for zi = 2:znum
                Uminl(:,:,zi) = min(U(:,:,zi),Uminl(:,:,zi-1));
            end
            Uminr = U;
            for zi = znum-1:-1:1
                Uminr(:,:,zi) = min(U(:,:,zi),Uminr(:,:,zi+1));
            end
            
            % Solve \nabla^2u=f
            % where
            % erfcl  = erfc (sqrt(eoverkBT*(U-Uminl)));
            % erfcxl = erfcx(sqrt(eoverkBT*(U-Uminl)));
            % f = c.*exp(eoverkBT*(Uminl-VFermi)).*(erfcxl+max(0,erfcxl-erfcx(sqrt(eoverkBT*(U-Uminr))).*exp(eoverkBT*(Uminr-Uminl))));
            %   = c.*exp(eoverkBT*(U    -VFermi)).*(erfcl +max(0,erfcl -erfc (sqrt(eoverkBT*(U-Uminr)))));
            % Remember to avoid NaN, Inf.
            
            erfcxl = erfcx(sqrt(eoverkBT*(U-Uminl)));
            f = c.*exp(eoverkBT*(Uminl-VFermi)).*(erfcxl+max(0,erfcxl-erfcx(sqrt(eoverkBT*(U-Uminr))).*exp(eoverkBT*(Uminr-Uminl))));
             
            change_erfc = ~isfinite(f); % To avoid NaN, Inf
            if sum(change_erfc,'all')
                erfcl  = erfc (sqrt(eoverkBT*(U-Uminl)));
                ferfc = c.*exp(eoverkBT*(U    -VFermi)).*(erfcl +max(0,erfcl -erfc (sqrt(eoverkBT*(U-Uminr)))));
                
                f(change_erfc)=ferfc(change_erfc);
                if isfinite(f)
                else
                    break
                end
            end

            % FFT method here
            fftf = fft2(f); 

            % boundary condition scaled
            fftf = z_multip.*fftf(1:halfxnum+1,:,:);
            
            if Schottky
                fftf(1,1,end)=xnum*ynum*(Vanode-USchottkyz(end));
            else
                fftf(1,1,end)=xnum*ynum*Vanode; 
            end
            
            ff2U1 = fft2(Vcathode);      
            fftf(:,:,1)=ff2U1(1:halfxnum+1,:);

            v = zeros(halfxnum+1,ynum,znum);
            for zi=2:znum-1
                w = Thomas_diag(:,:,zi) - Thomas_lower(zi).*v(:,:,zi-1);
                fftf(:,:,zi) = ( fftf(:,:,zi) - Thomas_lower(zi).*fftf(:,:,zi-1) )./w;
                v(:,:,zi) = Thomas_upper(zi)./w;
            end
            for zi=znum-1:-1:2
                fftf(:,:,zi) = fftf(:,:,zi) - v(:,:,zi).*fftf(:,:,zi+1);
            end

            fftU(1:halfxnum+1,:,:) = fftf;
            fftU(halfxnum+2:xnum,1,:) = conj(flip(fftf(2:halfxnum,1,:)));
            fftU(halfxnum+2:xnum,2:end,:) = conj(flip(flip(fftf(2:halfxnum,2:end,:)),2));

            Unew = real(ifft2(fftU));

            if sum(isnan(Unew),'all')
                break
            end

            if Schottky
                Unew = Unew + USchottkyz;
            end
            
            deltaU = Unew - U;     
          
            % step size choices
            if isfield(S,'iteroption') && isnumeric(S.iteroption)
                iter_step = S.iteroption;
            else
                iter_step = 1/sqrt(loopi);
            end

            if Schottky
                if abs(deltaU(:,:,2:end)) < max(newconv_cri1, newconv_cri2 * abs(U(:,:,2:end)))
                    U = U + iter_step*deltaU;
                    loopconv= true;
                    break
                end
            else
                if abs(deltaU) < max(newconv_cri1, newconv_cri2 * abs(U))
                    U = U + iter_step*deltaU;
                    loopconv= true;
                    break
                end
            end
            U = U + iter_step*deltaU;
        end
        
        if ~loopconv
            loopi = -loopi;
        end

        JRD = mean(Ar.*A.*T^2.*exp(eoverkBT*(Vcathode-VFermi)),[1 2]);
        
        if Schottky
            Ecathode = (U(:,:,2)-USchottkyz(1,1,2)-Vcathode)/(zlist(2)-zlist(1));
            U1Schottky  = Vcathode + (Ecathode>0).*sqrt(Ecathode*elec/(4*pi*epsilon0));
            zSchottky = (Ecathode>0).*sqrt(elec/(16*pi*epsilon0)./Ecathode);
            U2 = U(:,:,2);
            U1Schottky(zSchottky>zlist(2)) = U2(zSchottky>zlist(2));
            U1Schottky(Ecathode<=0) = Vcathode(Ecathode<=0)+USchottkyz(1,1,2);
            umin = min(cat(3,U1Schottky,U(:,:,2:end)),[],3);           
        else
            umin = min(U,[],3);
        end
        
        Jmap = Ar.*A.*T^2.*exp(eoverkBT*(umin-VFermi));
        J2halfD = mean(Jmap,[1 2]); 
  
        TVJRCTL(rowi,:) = [T,Vanode,J2halfD,JRD,JCL0,JCLT,loopi];

        if isfield(S,'Tstopcriterion') && rowi > S.Tstopcriterion(2) && J2halfD/TVJRCTL(rowi-S.Tstopcriterion(2),3) < 1+S.Tstopcriterion(1)
            TVJRCTL = TVJRCTL(1:rowi,:);
            break
        end

    end
    
    % variables in struct 'result'
    if isfield(S,'outputoption') && ~isempty(S.outputoption)
        if ismember("TVJRCTL",S.outputoption) || ismember("all",S.outputoption)
            result.TVJRCTL = TVJRCTL;
        end
        if ismember("I",S.outputoption) || ismember("all",S.outputoption)
            result.I = I;
        end
        if ismember("S",S.outputoption) || ismember("all",S.outputoption)
            result.S = S;
        end
        if ismember("umin",S.outputoption) || ismember("all",S.outputoption)
            result.umin = umin;
        end
        if ismember("Jmap",S.outputoption) || ismember("all",S.outputoption)
            result.Jmap = Jmap;
        end
        if ismember("U",S.outputoption) || ismember("all",S.outputoption)
            result.U = U;                
        end
        if ismember("zlist",S.outputoption) || ismember("all",S.outputoption)
            result.zlist = zlist;                
        end
        if ( ismember("U1Schottky",S.outputoption) || ismember("all",S.outputoption) ) && exist('U1Schottky','var')
            result.U1Schottky = U1Schottky;
        end
        if ismember("runtime",S.outputoption) || ismember("all",S.outputoption)
            result.runtime = toc;                
        end
    end
    
    % output to .mat file
    if isfield(S,'outputmatname')
        if ~isfield(S,'outputoption') || ismember("J",S.outputoption) || ~exist('result','var')
            save(strcat(S.outputmatname,'.mat'),'TVJRCTL');
        elseif isempty(S.outputoption) || ismember("N",S.outputoption) % No output
        else
            save(strcat(S.outputmatname,'.mat'),'TVJRCTL','result');
        end
    end
end

%% Standardize input

function S = inputstandardize(I)
    % Ar,
    % Vcathode, = -2 for work function 2eV.
    % VFermi, = 2 for Fermi level energy = -2eV.
    % Get cathode mesh size
    S.meshsize = [];
    if isfield(I,'Ar') && ~isscalar(I.Ar)
        if isempty(S.meshsize)
            S.meshsize = size(I.Ar);
        elseif ~isequal(S.meshsize,size(I.Ar))
            error('Error: Sizes of Ar, Vcathode, VFermi do not match.');
        end
    end
    if isfield(I,'Vcathode') && ~isscalar(I.Vcathode)
        if isempty(S.meshsize)
            S.meshsize = size(I.Vcathode);
        elseif ~isequal(S.meshsize,size(I.Vcathode))
            error('Error: Sizes of Ar, Vcathode, VFermi do not match.');
        end
    end
    if isfield(I,'VFermi') && ~isscalar(I.VFermi)
        if isempty(S.meshsize)
            S.meshsize = size(I.VFermi);
        elseif ~isequal(S.meshsize,size(I.VFermi))
            error('Error: Sizes of Ar, Vcathode, VFermi do not match.');
        end
    end
    if isempty(S.meshsize)
        S.meshsize = [1,1];
    end
    
    % Resize scalar inputs
    if ~isfield(I,'Ar')
        S.Ar = ones(S.meshsize);
    elseif isscalar(I.Ar)
        S.Ar = repmat(I.Ar,S.meshsize);
    else
        S.Ar = I.Ar;
    end

    if ~isfield(I,'Vcathode')
        S.Vcathode = zeros(S.meshsize);
    elseif isscalar(I.Vcathode)
        S.Vcathode = repmat(I.Vcathode,S.meshsize);
    else
        S.Vcathode = I.Vcathode;
    end
    
    if ~isfield(I,'VFermi')
        S.VFermi = zeros(S.meshsize);
    elseif isscalar(I.VFermi)
        S.VFermi = repmat(I.VFermi,S.meshsize);
    else
        S.VFermi = I.VFermi;
    end
    
    % Fill not finite input for Vcathode
    nf = ~isfinite(S.Vcathode);
    if sum(nf,'all')
        S.Ar(nf) = 0;
        S.VFermi(nf) = Inf;
        if isfield(I,'wfnanLaplacecrit') && isnumeric(I.wfnanLaplacecrit) && isscalar(I.wfnanLaplacecrit)
            S.Vcathode = laplace5point(S.Vcathode,wfnanLaplacecrit);
        else
            S.Vcathode = laplace5point(S.Vcathode,1e-12);
        end
    end
    
    % step = [xstep, ystep]
    if isscalar(I.step)
        S.step = I.step*[1,1];
    else
        S.step = vec2hor(I.step);
    end
    
    % TVanode = [T, Vanode];
    if isfield(I,'TVanode')
        S.TVanode = I.TVanode;
    elseif isfield(I,'T') && isfield(I,'Vanode')
        if isscalar(I.T)
            S.TVanode = [repmat(I.T,[length(I.Vanode),1]),hor2vec(I.Vanode)];
        elseif isscalar(I.Vanode)
            S.TVanode = [hor2vec(I.T),repmat(I.Vanode,[length(I.T),1])];
        else
            S.TVanode = [hor2vec(I.T),hor2vec(I.Vanode)];
        end
    else
        error('Error: Invalid input of T or Vanode.');
    end
    
    % If consider Schottky effect
    if isfield(I,'Schottky')
        S.Schottky = logical(I.Schottky);
    else
        S.Schottky = true;
    end
    
    % zlist
    if isfield(I,'dAK')
        if isfield(I,'zlist')
            S.zlist = I.dAK/I.zlist(end) * hor2vec(I.zlist);
        elseif isfield(I,'fixedzlist') && logical(I.fixedzlist)
            S.zlist = I.dAK * [0:2^-12:2^-7,(2^-7+2^-9):2^-9:2^-3,(2^-3+2^-7):2^-7:1]';
        else                  
            if isfield(I,'zlistset')
                zlistset = I.zlistset;
            else
                zlistset = [2^-11,2^-9,2^-7,64,64];
            end
            
            S.zlist = I.dAK * defaultunitzlist(max(S.Vcathode,[],'all'),min(S.Vcathode,[],'all'),min(S.TVanode(:,2)),max(S.TVanode(:,1)),I.dAK,zlistset);
        end
    else
        S.zlist = hor2vec(I.zlist);
    end
    
    % 
    if isfield(I,'looplimit')
        S.looplimit = I.looplimit;
    else
        S.looplimit = 1000;
    end
    
    % Tstopcriterion
    if isfield(I,'Tstopcriterion') && isnumeric(I.Tstopcriterion)
        if isscalar(I.Tstopcriterion)
            S.Tstopcriterion = [I.Tstopcriterion,1];
        elseif length(I.Tstopcriterion) == 2
            S.Tstopcriterion = vec2hor(I.Tstopcriterion);
        end
    end
    
    % iteroption
    if isfield(I,'iteroption') && isnumeric(I.iteroption) && isscalar(I.iteroption) && I.iteroption > 0 && I.iteroption <= 1
        S.iteroption = I.iteroption;
    end
    
    % outputmatname
    if isfield(I,'outputmatname')
        S.outputmatname = char(I.outputmatname);
    end
    % outputoption
    if isfield(I,'outputoption')
        S.outputoption = I.outputoption;
    end
    
    % conv_cri = [conv_cri1, conv_cri2]
    % Defaults:
    % conv_cri1 = 1e-4; % in V.
    % conv_cri2 = 1e-4; % in fraction.
    if ~isfield(I,'conv_cri')
        S.conv_cri = [1e-4,1e-4];
    elseif isscalar(I.conv_cri)
        S.conv_cri = I.conv_cri*[1,1];
    else
        S.conv_cri = vec2hor(I.conv_cri);
    end   
    
    % U0
    if isfield(I,'U0') && isequal(size(I.U0),[size(S.wf),length(S.zlist)])
        S.U0 = I.U0;
    end
    
    if isfield(I,'CLTintegral') && isscalar(I.CLTintegral) && ~logical(I.CLTintegral)
        S.CLTintegral = false;
    else
        S.CLTintegral = true; 
    end
    
end

%% Default zlist
function zlist = defaultunitzlist(Vkmax,Vkmin,Vanode,Tmax,dAK,zlistset)
    % constants
    kB = 1.380649e-23;
    m = 9.10938356e-31;
    epsilon0 = 8.854187817e-12;
    elec = 1.602176634e-19;
    A = 1.20173e6;

    JFSCL = 4*epsilon0/9*sqrt(2*elec/m)*Vanode^(3/2)/(dAK^2);    
    umin = min(Vkmin, log(JFSCL/(A*Tmax^2))*kB*Tmax/elec);

    % where U(z) = Ucathode
    zUk = 2*(Vkmax-umin)^(3/4)/((Vkmin-umin)^(3/4)+(Vanode-umin)^(3/4));
    if logical(zUk)
        zstep = min(zlistset(1),zUk/zlistset(4));
    else
        zstep = zlistset(1);
    end
    zlist = [0:zstep/8:0.95*zstep,zstep:zstep/4:8.1*zstep,9*zstep:zstep:zUk+zstep];
    
    % where U(z) = 0
    zU0 = ((-umin)^(3/4)+(Vkmax-umin)^(3/4))/((Vkmin-umin)^(3/4)+(Vanode-umin)^(3/4));
    zstep = min(zlistset(2),zU0/zlistset(5));    
    zlist = [zlist,zlist(end)+zstep:zstep:zU0+zstep];
    zlist = [zlist,zlist(end)+zlistset(3):zlistset(3):1];
    zlist(end) = 1;
    zlist = zlist';
end

%% 5-point Laplace's equation
function [B,iter]=laplace5point(A,meanerrcri)
    B = A;
    nfA = ~isfinite(A);
    Amean = mean(A(~nfA));
    B(nfA) = Amean;
    for iter = 1:1e4
        wfconv = cconv2(B);
        meanerr = max(abs(B(nfA)-wfconv(nfA)));
        if meanerr < Amean * meanerrcri
            break
        end
        B(nfA) = wfconv(nfA);
    end
end

function B=cconv2(A)
    C = conv2(A,[0,1,0;1,0,1;0,1,0]/4);
    B = C(2:end-1,2:end-1);
    B(:,1) = B(:,1) + C(2:end-1,end);
    B(:,end) = B(:,end) + C(2:end-1,1);
    B(1,:) = B(1,:) + C(end,2:end-1);
    B(end,:) = B(end,:) + C(1,2:end-1);
end

%% Conversion between horizontal and vertical vectors
function v = hor2vec(h)
    if size(h,1) == 1
        v = h';
    else
        v = h;
    end
end

function h = vec2hor(v)
    if size(v,2) == 1
        h = v';
    else
        h = v;
    end
end

%% Finite temperature correction on Child-Langmuir law
function result = TcalCL(T,VAK,ifintegration)
    elec = 1.602176634e-19;
    kB = 1.380649e-23;
    eVoverkT = elec.*VAK./kB./T;
    
    if nargin < 2.5 || ifintegration
        fun = @(u) u./sqrt(erfcx(u)-1+2/sqrt(pi).*u);
        result = nan(length(eVoverkT),1);
        for Vi = 1:length(eVoverkT)
            eVoverkTi = eVoverkT(Vi);
            result(Vi) = 9/(2*sqrt(pi)).*eVoverkTi.^(-3/2)*(integral(fun,0,sqrt(eVoverkTi),'RelTol',1e-11,'AbsTol',1e-11)).^2;
        end
    else
        result = 1+3/2.*sqrt(pi./eVoverkT);
    end    
end