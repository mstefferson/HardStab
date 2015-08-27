% DdftPrgmRunnerIDCube
%
% Runs a bunch text files or an individual program 
% created by into FT2Drot_HRDiffMain_exe_v9

RunAll = 1;     % Run all inputs
RunOne = 0;     % Run just one file
Clear  = 0;     % Clear temp directories
% Add path

MatHomeDir = sprintf('/home/mws/Documents/MATLAB/');
addpath( genpath( strcat(MatHomeDir,'Programs') ) );
addpath( genpath( strcat(MatHomeDir,'Research/BG/DiffFT/HRddft') ) );
addpath(strcat(MatHomeDir,'Research/BG/INtransEq'))

% Clear Temp Folders
if Clear
    fclose('all');
    cd C:/Users/MWS/Documents/MATLAB/Temps/HRddft
    DelDir = ls;
    [n, m] = size(DelDir);
    for ii = 3:n
        DelDirStrTemp = sprintf('%s', DelDir(ii,:) );
        keyboard
        rmdir(DelDirStrTemp)
    end
end
% Where to find the files
SrcStr = sprintf('C:/Users/MWS/Documents/Research/BG/Files2Run/HRddft');

%%%%%%%%%%%%% Running 1 file File %%%%%%%%%%%%%%%%%%%%%%%%%%
if RunOne
    % Enter a file name
    FileName = sprintf('blah.txt');
    Temp = 1234;         % Just a random integer
    % Run from a temp folder
    TempDirStr = sprintf('C:/Users/MWS/Documents/MATLAB/Temps/HRtemp_v%i',Temp);
    mkdir(TempDirStr)
    cd(TempDirStr)
    % Move input into temporary directory
    FilePathOld = sprintf('%s/%s',SrcStr,FileName);
    FilePathNew = sprintf('%s/%s',TempDirStr,FileName);
    movefile(FilePathOld,FilePathNew);
    FT2Drot_HRDiffMain_exe(FilePathNew)  
end

%%%%%%%%%%%% Running many files %%%%%%%%%%%%%%%%%%%%%%%
if RunAll
    RndInt = round( 100* rand());
    % Move to directory to grab files
    cd C:/Users/MWS/Documents/Research/BG/Files2Run/HRddft
    files = ls('*.txt');
    [NumFiles, void] = size(files);
    
    for ii = 1:NumFiles
        % Run from a temp folder
        TempDirStr = sprintf('C:/Users/MWS/Documents/MATLAB/Temps/HRddft/temp_v%i',RndInt + ii);
        mkdir(TempDirStr)
        cd(TempDirStr)
        % Move input into temporary directory
        FilePathOld = sprintf('%s/%s',SrcStr,files(ii,:));
        FilePathNew = sprintf('%s/%s',TempDirStr,files(ii,:));
        movefile(FilePathOld,FilePathNew);
%         keyboard
        HR2DrotExeMainIDCube(FilePathNew)
    end  
end

cd C:/Users/MWS/Documents/Research/BG/Results/HRddft/IsoDiffCube/