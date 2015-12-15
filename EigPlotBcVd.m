function EigPlotBcVd(bcVec,dbc,vDVec,dvD,bcE,maxRealEigVal,maxImagEigVal,...
    ParamStr1,ParamStr2,SaveMe,xMode,yMode,AnisoDiff,PerturbGen,Cmax,N)
figure
subplot(2,1,1)
pcolor( [bcVec bcVec(end)+dbc] ,[vDVec vDVec(end) + dvD],...
    maxRealEigVal)
Ax = gca;
Ax.CLim = [-6 3];
% Ax.YTick = [1:length(vDVec)+1];
% Ax.YTickLabel = num2cell( [vDVec vDVec(end) + dvD]);
colorbar
xlabel('bc');ylabel('vD');
if PerturbGen
    titstr =  sprintf(' Max Real Eigenvalue about Gen bc = %.2f ',bcE);
else
    titstr =  ' Max Real Eigenvalue about Iso';
end
title(titstr)
textbp(ParamStr1,'color','white');

subplot(2,1,2)
% pcolor( [bcVec bcVec(end)+dbc] ,1:length(vDVec)+1,maxImagEigVal)
pcolor( [bcVec bcVec(end)+dbc] ,[vDVec vDVec(end) + dvD],...
    maxImagEigVal)
Ax = gca;
% Ax.YTick = [1:length(vDVec)+1];
% Ax.YTickLabel = num2cell( [vDVec vDVec(end) * dvD]);
colorbar
xlabel('bc');ylabel('vD');
titstr =  sprintf(' Max Imag Eigenvalue');
title(titstr)
textbp(ParamStr2,'color','white');

if SaveMe
    if PerturbGen
        if bcE >= 1.5
        savestr = sprintf('EigsVsbcvDkx%dky%dAiD%dGenN%d',xMode,yMode,AnisoDiff,N);
        else
        savestr = sprintf('EigsVsbcvDkx%dky%dAiD%dGenI%d',xMode,yMode,AnisoDiff,N);    
        end
    else
        savestr = sprintf('EigsVsbcvDkx%dky%dAiD%dIso%d',xMode,yMode,AnisoDiff,N);
    end
    savefig(gcf,savestr)
end


NemSol = zeros( length(vDVec)+1,length(bcVec)+1 );
Unstab = maxRealEigVal > 0;
NemSol(Unstab) = maxRealEigVal(Unstab);
figure
% pcolor( [bcVec bcVec(end)+dbc] ,1:length(vDVec)+1,NemSol)
pcolor( [bcVec bcVec(end)+dbc] ,[vDVec vDVec(end) + dvD],NemSol)
Ax = gca;
Ax.CLim = [0 Cmax];
% Ax.YTick = [1:length(vDVec)+1];
% Ax.YTickLabel = num2cell( [vDVec vDVec(end) * dvD]);
colorbar
xlabel('bc');ylabel('vD');
if PerturbGen
    titstr =  sprintf(' Max Real Eigenvalue about Gen bc = %.2f ',bcE);
else
    titstr =  ' Max Real Eigenvalue about Iso';
end
title(titstr)
textbp(ParamStr1,'color','white');
textbp(ParamStr2,'color','white');


if SaveMe
     if PerturbGen
         if bcE >= 1.5    
             savestr = sprintf('MaxEigsVsbcvDkx%dky%dAiD%dGenN%d',xMode,yMode,AnisoDiff,N);
         else
             savestr = sprintf('MaxEigsVsbcvDkx%dky%dAiD%dGenI%d',xMode,yMode,AnisoDiff,N);
         end
        
    else
        savestr = sprintf('MaxEigsVsbcvDkx%dky%dAiD%dIso%d',xMode,yMode,AnisoDiff,N);
     end
    savefig(gcf,savestr)
end

