% Input creater for FT2DrotHRDiffMainExe
%

% Add path
MatHomeDir = sprintf('/home/mws/Documents/MATLAB/');
addpath( genpath( strcat(MatHomeDir,'Programs') ) );
addpath( genpath( strcat(MatHomeDir,'Research/BG/DDFT') ) );
addpath(strcat(MatHomeDir,'Research/BG/INtransEq'))

% Now can change number of grid points in the x, y, phi direction
Run  = 1; % Run main from here
Move = 0; % Move files to a nice location

%%%%%%%% Trial %%%%%%%%%%%%
trial    = 53;

%%%%%% Turn on/off interactions%%%%%%%%%
Interactions = 1;
SaveMe       = 1;   
Movies       = 1; % Movies won't run if save is zero

%%%%%%%%%%%%% Box and Rod Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx      = 32;
Ny      = 32;
Nm      = 32;
L_rod   = 1;                  % Length of the rods
Lx      = 2*pi;               % Box length
Ly      = 2*pi;               % Box length
AspctRt = 8;
v0      = 0;                  %Driving velocity

%%%%%%%%%%%%%%%Time recording %%%%%%%%%%%%%%%%%%%%%%%%%%
delta_t     = 1e-3;  %time step
t_record    = 1e-2;  %time interval for recording dynamics
t_tot       = 1e-1;  %total time
ss_epsilon  = 1e-6;   %steady state condition

% The number of k-modes above and below k = 0 added as a perturbation 
% to a uniform density.
NumModesX   = 10;
NumModesY   = 10;
NumModesM   = 10;


%%%%%%%%%%%%%%%%%%%%% Physical Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tmp      = 1;            % Temperature
% mobility
Mob_par  = 1; 
Mob_perp = 1; 
Mob_rot  = 1;

%%%%%%%%% Initial density parameters%%%%%%%%%%%%%%%%%%
% Weight of the spatial sinusoidal perturbation. 
% Perturbations added to rho(i,j,k) = 1. Must be small
WeightPos      = 1e-2;    
WeightAng      = 1e-2;  
Random         = 1;       % Random perturbation coeffs
bc             = 1.0;    % Dimensionless  scaled concentration

% Type of initial Condition
IntGauss   = 0;
IntPw      = 0;
IntSepPw   = 0;
IntEqPw    = 0;    % Distribution from perturbing equil. dist.
IntEqSepPw = 0;
IntLoad    = 0;
IntNemPw   = 1;    % Distribution from perturbing equil. dist.

% Save a string saying what you want
[IntDenType, IntDenIndicator] = IntDenIndicatorMaker(...
    IntGauss,IntPw,IntSepPw,IntEqPw,IntEqSepPw,IntLoad,IntNemPw);
     
if IntLoad
    IntDenName = sprintf('DenBlow2');
else
    IntDenName = sprintf('Irrelevant');
end

% Concentration and rod stuff
b       = L_rod^2/pi;               % Average excluded volume per particle
c       = bc / b;                   % Concentration
Norm    = c * Lx * Ly;              % number of particles
Diam    = L_rod / AspctRt;              % "Diameter" aka width of the rods

% Turn movies off is Save is off
if SaveMe == 0; Movies = 0;end
% Turn off Drive parameter  is v0 = 0;
if v0  == 0; Drive = 0; else Drive = 1;end

% Parameter vectors
Paramtmp = [trial Interactions Drive Movies SaveMe Nx Ny Nm Lx Ly...
    L_rod Tmp Norm WeightPos WeightAng Random ...
    NumModesX NumModesY NumModesM bc Mob_par Mob_perp Mob_rot v0];
Timetmp  = [delta_t t_record t_tot ss_epsilon];

% Make the output directory string and input file
FileDir    = sprintf('HRdiff_N%i%i%i_bc%.1f_Int%i_vd%.1f_%s%i%i%i_t%i',...
    Nx,Ny,Nm,bc,Interactions,v0,IntDenIndicator,NumModesX,NumModesY,NumModesM,trial);
FileInpt   = sprintf('Inpt_N%i%i%i_bc%.1f_Int%i_vd%.1f_%s%i%i%i_HEstep_t%i.txt', ...
    Nx,Ny,Nm,bc,Interactions,v0,IntDenIndicator,NumModesX,NumModesY,NumModesM,trial);

Where2SavePath    = sprintf('/home/mws/Documents/Research/BG/DDFT/Outputs/%s',FileDir);

fid = fopen(FileInpt,'wt');  % Note the 'wt' for writing in text mode
fprintf(fid,'%s\n%s\n%s\n',FileDir,Where2SavePath,IntDenType);
fprintf(fid,'Param\t');
fprintf(fid,'%f\t',Paramtmp);
fprintf(fid,'\n');
fprintf(fid,'Time\t');
fprintf(fid,'%f\t',Timetmp);

fclose(fid);

% keyboard
if Move == 1
    movefile(FileInpt,'C:/Users/MWS/Documents/Research/BG/Files2Run/HRddft')
end

% keyboard
if Run == 1
    % Add path
    % Make a temporary dir to run from
    RndInt     = round( 1000 * rand());
    TempDirStr = ...
        sprintf('/home/mws/Documents/MATLAB/Temps/HRddft/temp_%i_t%i',...
        RndInt,trial);
    WhereAmI = pwd;
    mkdir(TempDirStr)
    % Move input into temporary directory and move there
    FilePathOld = sprintf('%s/%s',WhereAmI,FileInpt);
    FilePathNew = sprintf('%s/%s',TempDirStr,FileInpt);
    movefile(FilePathOld,FilePathNew);
    cd(TempDirStr)
    [DenFinal, DenFTFinal, GridObj, DidIBreak,SteadyState] = ...
        HR2DrotDrMain(FileInpt);
end




