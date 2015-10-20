function EigPlotBcVd(bcVec,dbc,vDVec,dvD,bcE,maxRealEigVal,maxImagEigVal,...
    ParamStr1,ParamStr2,SaveMe,xMode,yMode,AnisoDiff,PerturbGen)
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
    titstr =  sprintf(' Max Real Eigenvalue about Gen bc = %.2e ',bcE);
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
        savestr = sprintf('EigsVsbcvDkx%dky%dAiD%dGen',xMode,yMode,AnisoDiff);
    else
        savestr = sprintf('EigsVsbcvDkx%dky%dAiD%dIso',xMode,yMode,AnisoDiff);
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
Ax.CLim = [0 3];
% Ax.YTick = [1:length(vDVec)+1];
% Ax.YTickLabel = num2cell( [vDVec vDVec(end) * dvD]);
colorbar
xlabel('bc');ylabel('vD');
if PerturbGen
    titstr =  sprintf(' Max Real Eigenvalue about Gen bc = %.2e ',bcE);
else
    titstr =  ' Max Real Eigenvalue about Iso';
end
title(titstr)
textbp(ParamStr1,'color','white');
textbp(ParamStr2,'color','white');


if SaveMe
     if PerturbGen
        savestr = sprintf('MaxEigsVsbcvDkx%dky%dAiD%dGen',xMode,yMode,AnisoDiff);
    else
        savestr = sprintf('MaxEigsVsbcvDkx%dky%dAiD%dIso',xMode,yMode,AnisoDiff);
     end
    savefig(gcf,savestr)
end

