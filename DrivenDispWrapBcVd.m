%% Isotropic!
% close all
function DrivenDispWrapBcVd

CurrentDir = pwd;
addpath( genpath( CurrentDir) );

Interactions = 1;
Diffusion    = 1;
AnisoDiff    = 1;
SparseMat    = 0;
SaveMe       = 1;

xMode    = 1;
yMode    = 0;

Nx    = 2^8;
Ny    = Nx;
Nm    = 2^8;

dbc   = 0.01;
bcVec = [ 1.40:dbc:1.60 ];
dvD   = 5;
vDVec = [0:dvD:50];

maxRealEigVal = zeros( length(vDVec)+1,length(bcVec)+1 ); %include extra ind
maxImagEigVal = zeros( length(vDVec)+1,length(bcVec)+1  );

kxHolder = Nx/2+1 + xMode;
kyHolder = Ny/2+1 + yMode;
    
L_rod   = 1;
Lx      = 10;
Ly      = Lx;
Mob_0   = 1;

Mob_par  = 2 * Mob_0;
Mob_perp = Mob_0;
Mob_rot = 6 *Mob_par/(L_rod^2);

D_par = Mob_par;

if AnisoDiff;
    D_perp = Mob_perp;
else
    D_perp = D_par;
end
D_rot = Mob_rot;

DiffMobObj = struct('Mob_par', Mob_par,'D_par',D_par,'D_perp',D_perp,...
    'Mob_rot', Mob_rot,'D_rot',D_rot);

ParamObj = struct('Nx',Nx,'Ny',Ny,'Nm',Nm,'Lx',Lx,'Ly',Ly,'L_rod',L_rod,...
    'bc',bcVec(1),'vD',vDVec(1));
[GridObj] = DispGridMaker(...
        ParamObj.Nx,ParamObj.Ny,ParamObj.Nm,ParamObj.Lx,ParamObj.Ly);
% keyboard

for i = 1:length(vDVec)
    ParamObj.vD = vDVec(i);
    
    for j = 1:length(bcVec)
        ParamObj.bc = bcVec(j);
        [eigVals] = DispEigCalcIsoSS(DiffMobObj,GridObj,...
            ParamObj,Interactions,Diffusion,...
            SparseMat ,kxHolder,kyHolder);
        maxRealEigVal(i,j) = max( real( eigVals ) );
        maxImagEigVal(i,j) = max( imag( eigVals ) );
%         keyboard
    end       
end

% keyboard

ParamStr1 = ...
    sprintf('N = %d\nLrod = %d\nLx = %.1f\nbcE  = %.2e\n',...
      Nm, L_rod,Lx,bcE);
ParamStr2 = ...
    sprintf('kx = %d\nky = %d\nAnisoDiff = %d\nPerturbGen = %d ',...
      xMode,yMode,AnisoDiff,PerturbGen);

EigPlotBcVdIso(bcVec,dbc,vDVec,dvD,maxRealEigVal,maxImagEigVal,...
    ParamStr1,ParamStr2,SaveMe,xMode,yMode,AnisoDiff)
% 
% figure

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

