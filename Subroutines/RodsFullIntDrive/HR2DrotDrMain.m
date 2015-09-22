% HR2DrotDrMain.m
%
% Program is the main() for running the diffusion of 2D hard rods with
% orientation. Program handles interactions using DDFT

function [DenFinal, DenFTFinal, GridObj, DidIBreak,SteadyState] = ...
    HR2DrotDrMain(InputFile)
% Add paths (this should already be added, but just to be careful)
% Save error messages in file
try
   
    %     keyboard
    tMainID  = tic;
    %Grab the parameters
    % keyboard
    DataTemp    = importdata(InputFile);
    ParamVec    = DataTemp.data(1,:);
    TimeVec     = DataTemp.data(2,~isnan(DataTemp.data(2,:)));    %Pull out the time stuff
    FileNameMat = DataTemp.textdata(1);
    Path2Save   = DataTemp.textdata(2);
    IntDenType  = DataTemp.textdata(3);
    
    % Make some  objects
    %     ParamNmVec = {'trial' 'Interactions' 'Nx' 'Ny' 'Nm' 'Lx' 'Ly' 'L_rod' 'Diam' 'Eta_visc'...
    %         'Tmp' 'Norm' 'WeightPos' 'WeightAng' 'NumModesX' 'NumModesY' 'NumModesM' 'bc'};
    ParamNmVec = {'trial' 'Interactions' 'Drive' 'Movies' 'SaveMe'...
        'Nx' 'Ny' 'Nm' 'Lx' 'Ly' 'L_rod' 'Diam' 'Eta_visc'...
        'kB' 'Tmp' 'Norm' 'WeightPos' 'WeightAng' 'NumModesX' 'NumModesY' ...
        'NumModesM' 'bc' 'Mob_pos','Mob_rot'  };
    TimeNmVec  = {'delta_t' 't_record' 't_tot' 'ss_epsilon'};
    
    ParamObj   = struct('NmVec',{ParamNmVec},'ValVec',ParamVec,'trial',ParamVec(1),...
        'Interactions',ParamVec(2), 'Drive',ParamVec(3),...
        'Movies',ParamVec(4),'SaveMe',ParamVec(5),...
        'Nx', ParamVec(6),'Ny', ParamVec(7),'Nm', ParamVec(8),...
        'Lx', ParamVec(9),'Ly', ParamVec(10),'L_rod', ParamVec(11), ...
        'Tmp',ParamVec(12), 'Norm',ParamVec(13), 'WeightPos',ParamVec(14), ...
        'WeightAng',ParamVec(15), 'Random',ParamVec(16),...
        'NumModesX',ParamVec(17), 'NumModesY',ParamVec(18), ...
        'NumModesM', ParamVec(19),'bc',ParamVec(20), ...
        'Mob_par',ParamVec(21),'Mob_perp',ParamVec(22), ...
        'Mob_rot',ParamVec(23), 'v0', ParamVec(24) );
    
    
    %     keyboard
    % Create a file that holds warning print statements
    WarningStmtString = sprintf('WarningStmts_%i.txt',ParamObj.trial);
    wfid  = fopen(WarningStmtString,'a+');    % a+ allows to append data
    
    LocString = sprintf('Location_%i.txt',ParamObj.trial);
    lfid      = fopen(LocString,'a+');    % a+ allows to append data
    fprintf(lfid,'Starting main, current code\n');
    
    %Time Recording
    N_time   = ceil(TimeVec(3)/TimeVec(1)); %number of time steps
    N_record = ceil(TimeVec(3)/TimeVec(2)); %number of time points to record. Does not include initial density
    N_count  = ceil(TimeVec(2)/TimeVec(1)); %spacing between times to record
    
    
    TimeObj = struct('TimeNmVec',{TimeNmVec},'TimeVec',{TimeVec},...
        'delta_t',TimeVec(1), 't_record',TimeVec(2), 't_tot',TimeVec(3), 'ss_epsilon',TimeVec(4),...
        'N_time', N_time, 'N_record',N_record,'N_count',N_count);
    
    
    %%%Make all the grid stuff%%%%%%%%%%%%%%
    [GridObj] = GridMakerPBCxk(...
        ParamObj.Nx,ParamObj.Ny,ParamObj.Nm,ParamObj.Lx,ParamObj.Ly);
    fprintf(lfid,'Made grid\n');
    
    %Make diffusion coeff (send smallest dx dy for stability
    [DiffMobObj] =  DiffDrCoupCoeffCalcGivenMob( wfid,ParamObj.Tmp,...
        ParamObj.Mob_par,ParamObj.Mob_perp,ParamObj.Mob_rot,...
        TimeObj.delta_t, min(GridObj.dx,GridObj.dy),...
        GridObj.dphi,GridObj.kx2D, GridObj.ky2D,ParamObj.v0);
    
    fprintf(lfid,'Made diffusion object\n');
    
    % Initialize density
    [rho] = MakeConcFromInd(GridObj,ParamObj,IntDenType);
    
    % Run the main code
    tBodyID      = tic;
    [DenRecObj]  = HR2DrotDrDenEvolverFTBody(wfid,lfid,rho,ParamObj, ...
        TimeObj,GridObj,DiffMobObj);
     BodyRunTime  = toc(tBodyID);
    fprintf(lfid,'Made density object\n');
    fprintf(lfid,'Body Run Time = %f\n\n', BodyRunTime);
    % keyboard
    % Store final density and transform
    DenFinal   = DenRecObj.Density_rec(:,:,:,end);
    DenFTFinal = DenRecObj.DensityFT_rec(:,:,:,end);
    DidIBreak  = DenRecObj.DidIBreak;
    SteadyState = DenRecObj.SteadyState;
    %     keyboard
    
    if ParamObj.SaveMe
        MemObj = 0;
        % Save all parameters
        PTMGDObj = struct('ParamObj',ParamObj,'TimeObj',TimeObj,'MemObj',MemObj,...
            'GridObj',GridObj,'D_par',DiffMobObj.D_par,...
            'D_perp',DiffMobObj.D_perp,'D_rot',DiffMobObj.D_rot);
        % Save everything. Save seperately for big files
        DenStr = sprintf('DenRec_%i',ParamObj.trial);
        PTMGDStr = sprintf('PTMGD_%i',ParamObj.trial);
        save(DenStr,'DenRecObj','-v7.3')
        save(PTMGDStr,'PTMGDObj','-v7.3')
    end
    
    % Run movies if you want
    if ParamObj.Movies == 1
        % Build OP records
        tOpID           = tic ;
%                 keyboard
        if  DenRecObj.DidIBreak == 0
            [OrderParamObj] = ...
                CPNrecMaker(ParamObj.Nx,ParamObj.Ny,...
                DenRecObj.TimeRecVec,GridObj,DenRecObj.Density_rec);
        else %Don't incldue the blowed up denesity for movies. 
            TimeRecVecTemp = DenRecObj.TimeRecVec(1:end-1);
            [OrderParamObj] = ...
                CPNrecMaker(ParamObj.Nx,ParamObj.Ny,TimeRecVecTemp,...
                GridObj,....
                DenRecObj.Density_rec(:,:,:,1:length(TimeRecVecTemp)));
        end
        OpRunTime       = toc(tOpID);
        fprintf(lfid,'Made interaction order paramater object\n');
        
        % Make matlab movies
        tMovID       = tic;
        
        [MovieObj]   = OPMatMovieMakerTgthr(ParamObj.Nx,ParamObj.Ny,ParamObj.Nm,...
            GridObj.x,GridObj.y,GridObj.phi,...
            OrderParamObj,DenRecObj.Density_rec,ParamObj.trial);
        %         keyboard
        MovRunTime   = toc(tMovID);
        
        fprintf(lfid,'Made movies\n');
        % Record how long it took
        fprintf(lfid,'OrderParam Run time = %f\n', OpRunTime);
        fprintf(lfid,'Make Mov Run Time = %f\n',  MovRunTime);
    end
    
    
    % Save how long everything took
    fprintf(lfid,'Everything saved\n');
    TotRunTime = toc(tMainID);
    fprintf(lfid,'Total Run time = %f\n', TotRunTime);
    
    fclose('all');
    % keyboard
    if ParamObj.SaveMe
        % Move everything
        MoveStrTxt = sprintf('*%i.txt', ParamObj.trial);
        MoveStrMat = sprintf('*%i.mat', ParamObj.trial);
        movefile(MoveStrTxt,Path2Save{1});
        movefile(MoveStrMat,Path2Save{1});
%         cd /home/mws/Documents/Research/BG/DDFT/Outputs
    end
    
catch err %Catch errors and move them
    %         keyboard
    
    ErrFileNmStr = sprintf('logFile%i.txt',ParamObj.trial);
    efid         = fopen(ErrFileNmStr,'a+');
    % write the error to file and to screen
    % first line: message
    fprintf(efid,'%s', err.getReport('extended', 'hyperlinks','off')) ;
    fprintf('%s', err.getReport('extended')) ;
    fclose(efid);
    fclose('all');
    
    keyboard
    if ParamObj.SaveMe
        DenStr = sprintf('DenRec_%i',ParamObj.trial);
        PTMGDStr = sprintf('PTMGD_%i',ParamObj.trial);
        PTMGDObj = struct('ParamObj',ParamObj,'TimeObj',TimeObj,...
            'GridObj',GridObj,'D_par',DiffMobObj.D_par,...
            'D_perp',DiffMobObj.D_perp,...
            'D_rot',DiffMobObj.D_rot);
        save(DenStr,'DenRecObj','-v7.3')
        save(PTMGDStr,'PTMGDObj','-v7.3')
        
        MoveStrTxt = sprintf('*%i.txt', ParamObj.trial);
        MoveStrMat = sprintf('*%i.mat', ParamObj.trial);
        
        movefile(MoveStrTxt,Path2Save{1});
        movefile(MoveStrMat,Path2Save{1});
    end
    
end %End try and catch

% clc
close all

end % end main