%% Isotropic!
% close all
function DrivenDispWrapBcNGen

CurrentDir = pwd;
addpath( genpath( CurrentDir) );

Interactions = 1;
Diffusion    = 1;
AnisoDiff    = 0;
PerturbGen   = 1;
SparseMat    = 0;
SaveMe       = 1;

xMode    = 0;
yMode    = 0;

NVaryStr = 'VaryAll';% VaryNx VaryNy VaryNxNy VaryNm VaryAll
NVec  = [32 64 128 256];
Nx    = 2^8;
Ny    = Nx;
Nm    = 2^8;

kxHolder = Nx/2+1 + xMode;
kyHolder = Ny/2+1 + yMode;

vD           = 0;
% bc           = 1.48;
dbc   = 0.05;
bcVec = [ 1.20:dbc:1.60 ];
bcE   = 1.4;

maxRealEigVal = zeros( length(NVec),length(bcVec) ); %include extra ind
maxImagEigVal = zeros( length(NVec),length(bcVec)  );
MaxEig = 10; % Guess at Max eigenvalue

L_rod   = 1;
Lx      = 10;
Ly      = Lx;
Mob_0   = 1;

Mob_par  = 2 * Mob_0;
Mob_perp = Mob_0;
Mob_rot = 6 *Mob_par/(L_rod^2);
Mob_rot  = Mob_par;
D_par = Mob_par;

if AnisoDiff;
    D_perp = Mob_perp;
else
    D_perp = D_par;
end
D_rot = Mob_rot;
D_rot =  Mob_par;

DiffMobObj = struct('Mob_par', Mob_par,'D_par',D_par,'D_perp',D_perp,...
    'Mob_rot', Mob_rot,'D_rot',D_rot);


ParamObj = struct('Nx',Nx,'Ny',Ny,'Nm',Nm,'Lx',Lx,'Ly',Ly,'L_rod',L_rod,...
    'bcP',bcVec(1),'bcE',bcE,'vD',vD);
[GridObj] = DispGridMaker(...
    ParamObj.Nx,ParamObj.Ny,ParamObj.Nm,ParamObj.Lx,ParamObj.Ly);
% keyboard

for i = 1:length(NVec)
    % Nx
    if ( strcmp('VaryNx',NVaryStr) ) || ( strcmp('VaryNxNy',NVaryStr) ) ...
            || ( strcmp('VaryAll',NVaryStr) )
        Nx = NVec(i);
        ParamObj.Nx = Nx;
        kxHolder = Nx/2+1 + xMode;
    end
    % Ny
    if ( strcmp('VaryNy',NVaryStr) ) || ( strcmp('VaryNxNy',NVaryStr) ) ...
            || ( strcmp('VaryAll',NVaryStr) )
        Ny = NVec(i);
        ParamObj.Ny = Ny;
        kyHolder = Ny/2+1 + yMode;
    end
    
    % Nm
    if ( strcmp('VaryNm',NVaryStr) ) || ( strcmp('VaryAll',NVaryStr) )
        Nm = NVec(i);
        ParamObj.Nm = Nm;
        
    end

    %     keyboard
    
    [GridObj] = DispGridMaker(...
        ParamObj.Nx,ParamObj.Ny,ParamObj.Nm,ParamObj.Lx,ParamObj.Ly);
    
    for j = 1:length(bcVec)
        ParamObj.bcP = bcVec(j);
        if PerturbGen
             [eigVecs,eigVals] = ...
                DispEigCalcGen(DiffMobObj,GridObj,ParamObj,...
                kxHolder,kyHolder,Interactions);       
        else
             [eigVecs,eigVals] = DispEigCalcIsoSS(DiffMobObj,GridObj,...
                ParamObj,Interactions,Diffusion,...
                SparseMat ,kxHolder,kyHolder);
        end
        maxRealEigVal(i,j) = max( real( diag(eigVals) ) );
        maxImagEigVal(i,j) = max( imag( diag(eigVals) ) );
        %         keyboard
    end
%     keyboard
end

% keyboard
ParamStr1 = ...
sprintf('N = %d\nLx = %.1f\nvD = %.1f\nbcE  = %.2f\n',...
      Nm, Lx,vD,bcE);
ParamStr2 = ...
    sprintf('kx = %d\nky = %d\nAnisoDiff = %d\nPerturbGen = %d ',...
      xMode,yMode,AnisoDiff,PerturbGen);
 
EigPlotBcN(bcVec,NVec,bcE,maxRealEigVal,maxImagEigVal,...
    ParamStr1,ParamStr2,SaveMe,xMode,yMode,AnisoDiff,PerturbGen,...
    NVaryStr,MaxEig);

function GridObj = DispGridMaker(Nx,Ny,Nm,Lx,Ly)

dx   = Lx/Nx;
dy   = Ly/Ny;
dphi = 2*pi/Nm;
% Make vectors and grids
x                = ( -Lx/2 : dx: Lx/2 - dx);
y                = ( -Ly/2 : dy: Ly/2 - dy);
phi              = ( 0: dphi: (2*pi - dphi) );
% Make k-space spacings
dkx          = 2*pi/Lx;
dky          = 2*pi/Ly;
dkm          = 1;
% Make k vectors and grids
kx_max           = pi/dx;            %Maximum spatial k-vector allowed by grid
ky_max           = pi/dy;            %Maximum spatial k-vector allowed by grid
km_max           = pi / dphi;        %Maximum angular k-vector
kx               = ( -kx_max: dkx: (kx_max - dkx) );
ky               = ( -ky_max: dky: (ky_max - dky) );
km               = ( -km_max: dkm: (km_max - dkm) );

%Put it all in an object
GridObj = struct('x',x,'y',y,'phi',phi,'dx',dx,'dy',dy,'dphi',dphi,...
    'kx',kx,'ky',ky,'km',km);

