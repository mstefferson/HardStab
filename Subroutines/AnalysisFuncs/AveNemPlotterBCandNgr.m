% AveNemPlotterBCandNgr
%
% Plots the average nematic parameter vs scaled concentration bc
% for various grid points

function AveNemPlotterBCandNgr(AveNemMatFit,AveNemMatPde,StdNemMatPde,...
    bcVec,NgrVec,Lx,L_rod,v0, SaveMe,trial)

%PLot it up

h(1) = figure;

xlabel('bc'); ylabel('<N>');
set(gca,'YLim',[0 1])
titstr    = sprintf('Average Nematic order');
Paramstr = sprintf('Lx/L_{rod}  = %.2f\nL_{rod} = %.2f\nv0 = %.2f',...
    Lx/L_rod, L_rod, v0);

title(titstr);

hold all

for i=1:length(NgrVec)
    errorbar(gca,bcVec,AveNemMatPde(i,:),StdNemMatPde(i,:),'-')
end

plot(gca,bcVec,AveNemMatFit,'--')
%         keyboard

legendCell1 = cellstr(num2str(NgrVec','N pde %d'));
legendCell2 = cellstr(num2str(NgrVec','N fit %d'));
legendCell  = [legendCell1; legendCell2];
legend(legendCell,'location','best')
% keyboard
textbp(Paramstr)

if SaveMe
    filename = sprintf('PlotAveNem%.2f_Lrod%.2f_t%d.fig',v0,L_rod,trial);
    savefig(gcf,filename);
end

figure()
titstr    = sprintf('Average Nematic order PDE');
subplot(2,1,1)
surf(bcVec,NgrVec,AveNemMatPde)
xlabel('bc'); ylabel('Ngr');zlabel('<N>');
set(gca,'ZLim',[0 1],'CLim',[0 1])
colorbar
title(titstr);


titstr    = sprintf('Average Nematic order PDE');
subplot(2,1,2)
surf(bcVec,NgrVec,AveNemMatFit)
xlabel('bc'); ylabel('Ngr');zlabel('<N>');
set(gca,'ZLim',[0 1],'CLim',[0 1])
colorbar
title(titstr);

if SaveMe
    filename = sprintf('SurfAveNem%.2f_Lrod%.2f_t%d.fig',v0,L_rod,trial);
    savefig(gcf,filename);
end
end % End PlotMe
