%% Isotropic!
% close all
function DrivenDispCheckWrapper

CurrentDir = pwd;
addpath( genpath( CurrentDir) );

Interactions = 1;
Diffusion    = 1;
AnisoDiff    = 1;
SparseMat    = 0;

vD           = 0;
bc           = 1.48;
bcVec = [ 1.4:0.01:1.6 ] ;
NVec  = [256]  ;
maxRealEigVal = zeros( length(NVec),length(bcVec) );
maxImagEigVal = zeros( length(NVec),length(bcVec) );

Nx    = 2^7;
Ny    = Nx;
Nm    = 2^7;

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
    'bc',bc,'vD',vD);

% keyboard

for i = 1:length(NVec)
    %
    %     Nx = NVec(i);
    %     ParamObj.Nx = NVec(i);
    %     Ny = NVec(i);
    %     ParamObj.Ny = NVec(i);
    Nm = NVec(i);
    ParamObj.Nm = Nm;
    
    kxHolder = Nx/2+2;
    kyHolder = Nx/2+1;
    
    [GridObj] = DispGridMaker(...
        ParamObj.Nx,ParamObj.Ny,ParamObj.Nm,ParamObj.Lx,ParamObj.Ly);
    
    for j = 1:length(bcVec)
        ParamObj.bc = bcVec(j);
        [eigVals] = DispEigCalc(DiffMobObj,GridObj,ParamObj,Interactions,Diffusion,...
            SparseMat ,kxHolder,kyHolder);
        maxRealEigVal(i,j) = max( real( eigVals ) );
        maxImagEigVal(i,j) = max( imag( eigVals ) );
    end
    
    
end

% keyboard
figure
plot(bcVec,maxRealEigVal)
legendCell = cellstr(num2str(NVec', 'N=%-d'));
xlabel('bc');ylabel('Maximum eigenvalue');
title('Nx = Ny = 128')
savefig(gcf,'MaxEigVsbcNxNy128')
legend(legendCell)

figure
plot(bcVec,maxImagEigVal)
legendCell = cellstr(num2str(NVec', 'N=%-d'));
xlabel('bc');ylabel('Maximum eigenvalue');
title('Nx = Ny = 128')
savefig(gcf,'MaxEigVsbcNxNy128')
legend(legendCell)



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

