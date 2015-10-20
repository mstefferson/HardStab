%% Just isotropic diffusion


function [M,MintAlpha,eigVals] = ...
    DispEigCalcGenIsDiff(DiffMobObj,GridObj,ParamObj,Diffusion,...
    kxHolder,kyHolder)


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
rhoEqlFT = zeros(Nm,Nm,Nm);
c = pi * ParamObj.bcP / (ParamObj.L_rod ^ 2);
% keyboard
% rhoEqlFT(ParamObj.Nx/2+1,ParamObj.Ny/2+1,:) = c .* f_FT;

% I need to work through the factor of Nm, but if should be there
rhoEqlFT(ParamObj.Nx/2+1,ParamObj.Ny/2+1,:) = ...
    c * ParamObj.Nx * ParamObj.Ny  .* f_FT;

% keyboard

% Prefactors
DiagPreFac = -( DiffMobObj.D_par + DiffMobObj.D_perp ) / 2 .* ...
    ( GridObj.kx(kxHolder) .^ 2 + GridObj.ky(kyHolder) .^ 2 )...
    - DiffMobObj.D_rot * GridObj.km .^2 ;
DiagPreFac = DiagPreFac';

% Build diagonal
if Diffusion
    DiffContrib = DiagPreFac;
else
    DiffContrib = zeros(Nm,1);
end

% rho_0 = c / (2*pi)

    Fm = MayerFncDiffBtwPntsCalc(...
        ParamObj.Nx, ParamObj.Ny, ParamObj.Nm, ...
        ParamObj.Lx, ParamObj.Ly, GridObj.dx,...
        GridObj.dy, GridObj.dphi, ParamObj.L_rod);
    Fm_FT = fftshift(fftn( Fm ) );
    %     toc
    %     keyboard
    % Divide by 1/N^3 from FT normalization
    IntPreFac = 2 * pi * ParamObj.Lx * ParamObj.Ly / ...
        ( ParamObj.Nx * ParamObj.Ny * ParamObj.Nm ) ^ 2;
    
   
%% Build the Matrix
posksqr   = ( DiffMobObj.D_par + DiffMobObj.D_perp ) / 2 .* ...
    GridObj.kx(kxHolder) .^ 2 + GridObj.ky(kyHolder) .^ 2;
MintAlpha = zeros(Nm,Nm);
MintBeta  = zeros(Nm,Nm);
km   = GridObj.km;
for i = 1:Nm;

    for j = 1:Nm;
        % alpha = k^2 + m^2 + m*n
        Alphaholder = Nm / 2 + 1 - j + i;
        if Alphaholder < 1
            Alphaholder = Alphaholder + Nm;
        end
        if Alphaholder > Nm
            Alphaholder = Alphaholder - Nm;
        end
        alphatemp = IntPreFac .* ...
            (posksqr + DiffMobObj.D_par .* ...
            (km(j).^2 + km(Alphaholder) .* km(j) ) );
      
        MintAlpha(i,j) = alphatemp .* ...
            rhoEqlFT(Nm/2+1,Nm/2+1,Alphaholder) .* ...
            Fm_FT(kxHolder,kyHolder,j);
        
        Betaholder = Nm / 2 + 1 - j + i;
        
        if Betaholder < 1
            Betaholder = Betaholder + Nm;
        end
        if Betaholder > Nm
            Betaholder = Betaholder - Nm;
        end
                
%         keyboard
        betatemp = IntPreFac .* ...
            DiffMobObj.D_par .* ...
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

%%

% keyboard
% keyboard 
    Mint = MintAlpha + MintBeta;
    M =  diag( DiffContrib') + Mint +...
        diag( Mm1vec(1:Nm-1) ,-1) + diag(  Mp1vec(1:Nm-1) ,1) ;
           
   % PBC driving
    M(1,ParamObj.Nm) = Mm1PreFac;
    M(ParamObj.Nm,1) = Mp1PreFac;
%   keyboard 
        
    [eigVals] = eig(M);

% keyboard
end
% toc