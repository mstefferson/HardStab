function EigPlotBcN(bcVec,NVec,bcE,maxRealEigVal,maxImagEigVal,...
    ParamStr1,ParamStr2,SaveMe,xMode,yMode,AnisoDiff,PerturbGen,NVaryStr,MaxEig)
figure
subplot(2,1,1)
plot(bcVec,maxRealEigVal)
Ax = gca;
Ax.YLim = [0 MaxEig];
legendCell = cellstr(num2str(NVec', 'N=%-d'));
legend(legendCell)

xlabel('bc');ylabel('Max Real Eig Val');
if PerturbGen
    str2 = sprintf(': Gen bcE = %.2e', bcE);
    titstr = [NVaryStr str2];
else
    titstr = [NVaryStr ': about Iso'];
end

title(titstr)
textbp(ParamStr1,'color','black');

subplot(2,1,2)
plot(bcVec,maxImagEigVal)
Ax = gca;
legendCell = cellstr(num2str(NVec', 'N=%-d'));
legend(legendCell)
xlabel('bc');ylabel('Imag Eig Val');
title(titstr)
textbp(ParamStr2,'color','black');

if SaveMe
    if PerturbGen
        savestr1 = sprintf('EigsbcNkx%dky%dAiD%dGen',xMode,yMode,AnisoDiff);
    else
        savestr1 = sprintf('EigsbcNkx%dky%dAiD%dIso',xMode,yMode,AnisoDiff);
    end
    savestr  = [savestr1 NVaryStr];
    savefig(gcf,savestr)
end