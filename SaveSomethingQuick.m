MaxTemp = maxRealEigVal(1:end-1,1:end-1)';
ImagTemp = maxImagEigVal(1:end-1,1:end-1)';

ParamStr1 = ...
    sprintf('N = %d\nLrod = %d\nLx = %.1f\ndx = %.2e\n',...
      Nm, L_rod,Lx,GridObj.dx);
ParamStr2 = ...
    sprintf('kx = %d\nky = %d',...
      xMode,yMode);

figure
subplot(2,1,1)
plot(bcVec,MaxTemp)
Ax = gca;
legendcell = cellstr( num2str( NVec','N=%d' ) );
legend(legendcell,'location','best')
xlabel('bc');ylabel('Max Real Eig Val');
textbp(ParamStr1)
titstr =  sprintf(' Real: Loop Nx'); 
title(titstr)

% Ax.YTick = [1:length(vDVec)+1];
% Ax.YTickLabel = num2cell( [vDVec vDVec(end) + dvD]);
subplot(2,1,2)
plot(bcVec,ImagTemp)
Ax = gca;
legendcell = cellstr( num2str( NVec','N=%d' ) );
legend(legendcell,'location','best')
xlabel('bc');ylabel('Max Real Eig Val');
textbp(ParamStr2)
titstr =  sprintf(' Imag'); 
title(titstr)
