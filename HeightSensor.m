function dataout=HeightSensor(name)
    NSMU = NSMatlabUtilities();

    %open a NanoScope image file 
    
    %name='Test2.spm';
    NSMU.Open(which(name));
    
    %numChan = NSMU.GetNumberOfChannels();
    numChan = 1;
    
%     f = figure();
%     movegui(f,'north');
%     pos = get(f,'Position');
%     set(f,'Position', pos + [-150 0 300 0]);
    
    %-------------------------------------------
    %Display Image for each channel
    for ichan = 1:numChan
        
        %Get data and some descriptive info
        [data, scale_units, type_desc] = NSMU.GetImageData(ichan, NSMU.METRIC);
        lineDir = NSMU.GetLineDirection(ichan);
        scanLine = NSMU.GetScanLine(ichan);
        AspectRatio = NSMU.GetImageAspectRatio(ichan);
        
        %compute planefit coefficients and display on command line
        [a, b, c, fitTypeStr] = NSMU.GetPlanefitSettings(ichan);
        disp(['Chan # ',num2str(ichan),': Aspect Ratio = ',num2str(AspectRatio),...
            ', PlaneFit Coeffs = (',num2str(a),',',num2str(b),',',num2str(c),')']);
       
%         sp = subplot(1,double(numChan),double(ichan));
% 
%         %this code spaces things nicely, but I got here iteratively, not
%         %logically :=}
%         pos = get(sp, 'Position'); % gives the position of current sub-plot
%         if ichan == 1
%             new_pos = pos + [-.04 0 .02 .02]; 
%         elseif ichan == 2
%             new_pos = pos + [-.015 0 .02 .02];
%         elseif ichan == 3
%             new_pos = pos + [0 0 .02 .02];
%         end        
%         set(sp, 'Position',new_pos );
        
%         % now plot and annotate
         dataout=flipud(data);
%         image(flipud(data),'CDataMapping','scaled');
%         set(gca,'YDir','normal');
%         axis('tight', 'square'); 
%         colorbar;
%         title(type_desc)  
%         xlabel(strcat(lineDir,'; ',scanLine));
%         axis('tight', 'square');
%         drawnow;
    end

    NSMU.Close();
%end
