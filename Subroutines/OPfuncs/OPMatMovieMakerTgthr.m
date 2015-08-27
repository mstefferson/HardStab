function [MovieObj] = ...
    OPMatMovieMakerTgthr(Nx,Ny,Nm,x,y,phi,OrderParamObj,Density_rec,SaveMe,trial)

%%%% Concentration %%%%%
nFrames = OrderParamObj.nFrames;

%Initialize the movie structure array
M_All(nFrames)  = struct('cdata',zeros(Nx,Ny,3,'int8'), 'colormap',[]);

% keyboard
xPos  = Nx/2;
yPos  = Ny/2;
eps = 0.000001;

% Set up figure
figure()
% keyboard
set(gcf,'renderer','zbuffer')

axh1 = subplot(2,2,1); % Save the handle of the subplot
colorbar('peer',axh1)
axpos1 = get(axh1,'position'); % Save the position as ax
set(axh1,'NextPlot','replaceChildren',...
    'CLim',[min(min(min(OrderParamObj.C_rec))) max(max(max(OrderParamObj.C_rec)))],'YDir','normal');
set(axh1,'position',axpos1,'NextPlot','replaceChildren'); % Manually setting this holds the position with colorbar

axh2 = subplot(2,2,2); % Save the handle of the subplot
axpos2 = get(axh2,'position'); % Save the position as ax
set(axh2,'position',axpos2); % Manually setting this holds the position with colorbar

axh3 = subplot(2,2,3); % Save the handle of the subplot
axpos3 = get(axh3,'position'); % Save the position as ax
set(axh3,'position',axpos3); % Manually setting this holds the position with colorbar
set(axh3,'NextPlot','replaceChildren',...
        'YLim',[ min(min(min(min( Density_rec( xPos, yPos, :, : ) ) ) ) ) - eps ...
        max(max(max(max( Density_rec( xPos, yPos, :, : ) ) ) ) )  + eps] )
    
axh4 = subplot(2,2,4); % Save the handle of the subplot
axpos4 = get(axh4,'position'); % Save the position as ax
set(axh4,'position',axpos4); % Manually setting this holds the position with colorbar
    

for ii = 1:nFrames
    
    %     set(axh1,'position',axpos1,'NextPlot','replaceChildren'); % Manually setting this holds the position with colorbar
    pcolor(axh1,x,y,OrderParamObj.C_rec(:,:,ii)')
    shading(axh1,'interp');
    TitlStr = sprintf('Concentration t = %.2f', OrderParamObj.TimeRec(ii));
    title(axh1,TitlStr)
     
    pcolor(axh2,x,y,OrderParamObj.POP_rec(:,:,ii)');
    shading(axh2,'interp'); colorbar('peer',axh2);
    set(axh2,'NextPlot','replaceChildren',...
    'CLim',[0 max(max(max(OrderParamObj.POP_rec)))],'YDir','normal');
    %     colorbar;
    %     set(gcf,'renderer','zbuffer')
    hold(axh2,'on')
    quiver(axh2,x,y,OrderParamObj.nx_POP_rec(:,:,ii)',OrderParamObj.ny_POP_rec(:,:,ii)','color',[1,1,1]);
    hold(axh2,'off')
    TitlStr = sprintf('Polar Order t = %.2f', OrderParamObj.TimeRec(ii));
    title(axh2,TitlStr)

    plot( axh3,phi, reshape( Density_rec(xPos,yPos,:,ii), 1, Nm ) );
    TtlStr = sprintf('Distribution (%d,%d) t = %.2f',...
                      xPos,yPos,OrderParamObj.TimeRec(ii));
    title(axh3,TtlStr);
   
    
    set(axh4,'NextPlot','replaceChildren',...
        'CLim',[0 max(max(max(OrderParamObj.NOP_rec)))],'YDir','normal');
    pcolor(axh4,x,y,OrderParamObj.NOP_rec(:,:,ii)');
    shading(axh4,'interp'); colorbar('peer',axh4)
    hold(axh4,'on')
    quiver(axh4,x,y,OrderParamObj.NADx_rec(:,:,ii)',OrderParamObj.NADy_rec(:,:,ii)',...
        'color',[1,1,1],'ShowArrowHead','off');
    hold(axh4,'off')
    TitlStr = sprintf('Nematic Order t = %.2f', OrderParamObj.TimeRec(ii));
    %     TitlStr = sprintf('Nematic Order');
    title(axh4,TitlStr)
%     keyboard
    
    M_All(ii) = getframe(gcf); %Store the frame
end


MovieObj = struct('M_All',M_All);
if SaveMe
MovStr = sprintf('MovieAll_%i',trial);
save(MovStr,'M_All','-v7.3')
end
% keyboard
close all
