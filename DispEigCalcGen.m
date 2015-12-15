function [eigVecs,eigVals] = ...
    DispEigCalcGen(DiffMobObj,GridObj,ParamObj,...
    kxHolder,kyHolder,Interactions)


% keyboard

% Vectors are column vectors
Nm  = ParamObj.Nm;
kx0holder = ParamObj.Nx/2 +1;
ky0holder = ParamObj.Ny/2 +1;

Nc = 10;
% Calculate coeff
[Coeff_best, ~] = CoeffCalcExpCos2D(Nc,GridObj.phi,ParamObj.bcE); 
% Build equil distribution
f = DistBuilderExpCos2Dsing(Nc,GridObj.phi,Coeff_best);        
f_FT = fftshift( fft( f ) );
rhoEqlFT = zeros(ParamObj.Nx,ParamObj.Ny,Nm);
c = pi * ParamObj.bcP / (ParamObj.L_rod ^ 2);

% rhoEqlFT(ParamObj.Nx/2+1,ParamObj.Ny/2+1,:) = c .* f_FT;


rhoEqlFT(ParamObj.Nx/2+1,ParamObj.Ny/2+1,:) = ...
    c * ParamObj.Nx * ParamObj.Ny  .* f_FT;    
    
% rhoEqlTemp = c * ParamObj.Nx * ParamObj.Ny  .* f_FT(Nm/2+1);
% rhoEqlFT = zeros(Nm,Nm,Nm);
% rhoEqlFT(ParamObj.Nx/2+1,ParamObj.Ny/2+1,Nm/2+1) = rhoEqlTemp;


% keyboard

if 1
    rhoEqlFT =  real(rhoEqlFT);   
end
% keyboard

% Prefactors
MdiagVec = -( DiffMobObj.D_par + DiffMobObj.D_perp ) / 2 .* ...
    ( GridObj.kx(kxHolder) .^ 2 + GridObj.ky(kyHolder) .^ 2 )...
    - DiffMobObj.D_rot * GridObj.km .^2 ;
MdiagVec = MdiagVec';


% rho_0 = c / (2*pi)

%     Fm = MayerFncDiffBtwPntsCalc(...
%         ParamObj.Nx, ParamObj.Ny, ParamObj.Nm, ...
%         ParamObj.Lx, ParamObj.Ly, GridObj.dx,...
%         GridObj.dy, GridObj.dphi, ParamObj.L_rod);
  Fm = MayerFncDiffBtwPntsCalc2(...
        ParamObj.Nx, ParamObj.Ny, ParamObj.Nm, ...
        ParamObj.Lx, ParamObj.Ly, ParamObj.L_rod);
    
Fm_FT = fftshift(fftn( Fm ) );
    %     toc
%         keyboard
    % Divide by 1/N^3 from FT normalization
    IntPreFac = 2 * pi * ParamObj.Lx * ParamObj.Ly / ...
        ( ParamObj.Nx * ParamObj.Ny * ParamObj.Nm ) ^ 2;
    
   
%% Build the Matrix
posksqr   = ( DiffMobObj.D_par + DiffMobObj.D_perp ) / 2 .* ...
    ( GridObj.kx(kxHolder) .^ 2 + GridObj.ky(kyHolder).^ 2 ) ;
MintAlpha = zeros(Nm,Nm);
MintAlphaP2 = zeros(Nm,Nm);
MintAlphaM2 = zeros(Nm,Nm);
MintBeta  = zeros(Nm,Nm);

alphatempP2 = IntPreFac .* (DiffMobObj.D_par - DiffMobObj.D_perp) / 4 * ...
    ( GridObj.kx(kxHolder) + sqrt(-1) * GridObj.ky(kyHolder) ) ^ 2;

alphatempM2 = IntPreFac .* (DiffMobObj.D_par - DiffMobObj.D_perp) / 4 * ...
    ( GridObj.kx(kxHolder) - sqrt(-1) * GridObj.ky(kyHolder) ) ^ 2;

% keyboard
km   = GridObj.km;
for i = 1:Nm;

    for j = 1:Nm;
        % alpha = k^2 + m^2 + m*n
        Alphaholder = Nm / 2 + 1 - j + i;
        AlphaholderP2 = Nm / 2 + 1 - j + i + 2; % exp(m + n - 2)
        AlphaholderM2 = Nm / 2 + 1 - j + i - 2; % exp(m + n + 2)
        
        if Alphaholder < 1
            Alphaholder = Alphaholder + Nm;
        end
        if Alphaholder > Nm
            Alphaholder = Alphaholder - Nm;
        end
        
        if AlphaholderP2 < 1
            AlphaholderP2 = AlphaholderP2 + Nm;
        end
        if AlphaholderP2 > Nm
            AlphaholderP2 = AlphaholderP2 - Nm;
        end
        
        if AlphaholderM2 < 1
            AlphaholderM2 = AlphaholderM2 + Nm;
        end
        if AlphaholderM2 > Nm
            AlphaholderM2 = AlphaholderM2 - Nm;
        end
        
        alphatemp = IntPreFac .* ...
            (posksqr + DiffMobObj.D_rot .* ...
            (km(j).^2 + km(Alphaholder) .* km(j) ) );
      
        if alphatempP2 ~= 0
        MintAlphaP2(i,j) = alphatempP2 .* ...
            rhoEqlFT(ParamObj.Nx/2+1,ParamObj.Ny/2+1,AlphaholderP2) .* ...
            Fm_FT(kxHolder,kyHolder,j);
        end
        if  alphatempM2 ~= 0
        MintAlphaM2(i,j) = alphatempM2 .* ...
            rhoEqlFT(ParamObj.Nx/2+1,ParamObj.Ny/2+1,AlphaholderM2) .* ...
            Fm_FT(kxHolder,kyHolder,j);
        end
        
        if Alphaholder == Nm/2 +1
%             keyboard
%             IntPreFacTmp = (pi * ParamObj.Lx * ParamObj.Ly * ParamObj.bcP) / ...
%         ( ParamObj.L_rod^2 * ParamObj.Nx * ParamObj.Ny * ParamObj.Nm );
%     
%        InsIntCntr = IntPreFacTmp .*  DiffMobObj.D_rot * GridObj.km(i) .^2 .* ...
%              Fm_FT(kxHolder,kyHolder,j);
         
%          IntPreFacTmp .*  DiffMobObj.D_rot * GridObj.km(i) .^2  - ...
%              alphatemp .* rhoEqlFT(ParamObj.Nx/2+1,ParamObj.Ny/2+1,Alphaholder);
        end
        
        if i == 207 && j == 207
%             keyboard
        end
        
        MintAlpha(i,j) = alphatemp .* ...
            rhoEqlFT(ParamObj.Nx/2+1,ParamObj.Ny/2+1,Alphaholder) .* ...
            Fm_FT(kxHolder,kyHolder,j);

        Betaholder = Nm / 2 + 1 - j + i;
        
        if Betaholder < 1
            Betaholder = Betaholder + Nm;
        end
        if Betaholder > Nm
            Betaholder = Betaholder - Nm;
        end
   
        betatemp = IntPreFac .* DiffMobObj.D_rot .* ...
            (km(Betaholder).^2 + km(Betaholder) .* km(j) );
        
        MintBeta(i,j) = betatemp .* ...
            rhoEqlFT(kx0holder,ky0holder,Betaholder) .* ...
            Fm_FT(kx0holder,ky0holder,Betaholder );
    end
    
%     keyboard
end

% keyboard
%% Driving

Mp1PreFac = -ParamObj.vD ./ 2 .* ( ...
    sqrt(-1) .* GridObj.kx(kxHolder) - GridObj.ky(kyHolder) );
Mm1PreFac = -ParamObj.vD ./ 2 .* ( ...
    sqrt(-1) .* GridObj.kx(kxHolder) + GridObj.ky(kyHolder) );

Mp1vec = Mp1PreFac .* ones(Nm,1);
Mm1vec = Mm1PreFac .* ones(Nm,1);


%% Aniso Just diff
Mp2PreFac =  - ( DiffMobObj.D_par - DiffMobObj.D_perp ) / 4 .* ...
    ( GridObj.kx(kxHolder) + sqrt(-1) .* GridObj.ky(kyHolder)  ) .^2;

Mm2PreFac = - ( DiffMobObj.D_par - DiffMobObj.D_perp ) / 4 .* ...
    ( GridObj.kx(kxHolder) - sqrt(-1) .* GridObj.ky(kyHolder)  ) .^2;

Mp2vec = Mp2PreFac .* ones(Nm,1);
Mm2vec = Mm2PreFac .* ones(Nm,1);
%%

% keyboard

if Interactions
    Mint = MintAlpha + MintAlphaP2 + MintAlphaM2 + MintBeta;
else
    Mint = 0;
end
% keyboard 
    M =  diag( MdiagVec') +...
        diag(  Mp1vec(1:Nm-1) ,1) +...
        diag(  Mm1vec(1:Nm-1) ,-1) +...
        diag(  Mp2vec(1:Nm-2) ,2) +...
        diag(  Mm2vec(1:Nm-2) ,-2) ;
    % PBC driving
    M(1,ParamObj.Nm) = Mm1PreFac;
    M(ParamObj.Nm,1) = Mp1PreFac;
    
    % PBC aniso. Just diffusion now
    M(1,ParamObj.Nm-1) = Mm2PreFac;
    M(ParamObj.Nm-1,1) = Mp2PreFac;
    M(2,ParamObj.Nm)   = Mm2PreFac; 
    M(ParamObj.Nm,2)   = Mp2PreFac;

    M = M + Mint;
   % PBC driving
    
%   keyboard 
        
    [eigVecs,eigVals] = eig(M);

% keyboard
end
% toc