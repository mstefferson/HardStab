%% Isotropic!
% close all
function DrivenDispCheck

CurrentDir = pwd;
addpath( genpath( CurrentDir) );

Interactions = 1;
Diffusion    = 1;
SparseMat    = 0;
vD           = 0;
bc           = 1.48;

Nx    = 2^8;
Ny    = Nx;
Nm    = 2^8;

kxHolder = Nx/2+1;
kyHolder = Nx/2+1;

L_rod = 1;
Lx    = 10;
Ly    = Lx;
Mob_pos = 1;
Mob_rot = 6 *Mob_pos/(L_rod^2);

D_pos = Mob_pos;
D_rot = Mob_rot;

% calculate bc
b     = L_rod^2/pi;               % Average excluded volume per particle
c     = bc / b;                   % Concentration

DiffMobObj = struct('Mob_pos', Mob_pos,'D_pos',D_pos,'Mob_rot', Mob_rot,'D_rot',D_rot);

GridObj = DispGridMaker(Nx,Ny,Nm,Lx,Ly);

ParamObj = struct('Nx',Nx,'Ny',Ny,'Nm',Nm,'Lx',Lx,'Ly',Ly,'L_rod',L_rod,...
    'bc',bc,'vD',vD);

[eigVals] = DispEigCalc(DiffMobObj,GridObj,ParamObj,Interactions,Diffusion,...
    SparseMat ,kxHolder,kyHolder);

maxEigVal = max( real( eigVals ) );
disp(maxEigVal)

% disp(maxEigVal)
% %%%%%%%%%%%%%%Plot Eigenvalues%%%%%%%%%%%%%%%%%%%%
% figure
% subplot(2,1,1)
% plot(1:Nm, sort( real(eigVals) ),1:Nm, sort( real( MdiagVec ) ) )
% legend('eigenvalues','eigenvalues no driving')
% title('real part eigenvals')
% subplot(2,1,2)
% plot(1:Nm, sort( imag(eigVals) ) )
% title('imaginary part eigenvals')
% 
% fprintf('Max eigenval = %f\n', max( real( eigVals) ) )
% % figure
% % subplot(2,1,1)
% % plot( GridObj.km,GridObj.km,real( ResortEigs ), GridObj.km, DiffContrib,...
% %     GridObj.km, real(IntContrib))
% % legend('Eigs','Diff','Inter','Sum')
% % title('real')
% % subplot(2,1,2)
% % plot( imag( diag(eigVals)  ) )
% % title('imag')
function GridObj = DispGridMaker(Nx,Ny,Nm,Lx,Ly)

dx   = Lx/Nx;              
dy   = Ly/Ny;              
dphi = 2*pi/Nm;
% Make vectors and grids
% x                = ( 0: dx : Lx - dx );
x                = ( -Lx/2 : dx: Lx/2 - dx); 
% y                = ( 0: dy : Ly - dy );
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

