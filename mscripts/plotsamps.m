function hplot=plotsamps(handles);
% function to plot samples from sampgui
% slb 4apr04

% say what's going on
%hwait=commentgui_devel('Havenotea','Position',[20 20 100 21]);


% set the current figure to be the one we started plotting in
figure(handles.figno);
handles.res=str2double(get(handles.ResEdit,'String'));

% things probs won't want to change:
dotted1d=handles.clevels;
conf2dcol='qbwr-'; % default plot type (change color after plotting)
margmeth='marg';
uselike = 0; % I never really set this up properly
npars=2;

if (handles.isheld==0)
    hold off;
    clflogo
else
    hold on;
end

% sort out priors and cuts
sz=size(handles.samples);
isamps=1:sz(2); % by default use them all
for i=1:sz(1)
    if (handles.cuts(1,i)==1)
        isamps=isamps(find(handles.samples(i,isamps)>=handles.cuts(2,i)));
        isamps=isamps(find(handles.samples(i,isamps)<=handles.cuts(3,i)));
    else
        fprintf(1,'ERROR, havent planned for this type of cut yet!\n')
    end
end

count=handles.samples(1,isamps); % just the weights (cut out for historical reasons)


if ((handles.ycol==1)&(handles.xcol~=1)) % have probability up y axis
    % make a 1d plot
    % it shouldn't be possible to get to this point if zcol>1 but
    % if it happens anyway, then zcol is just ignored
    
    z=handles.samples(handles.xcol,isamps);
    % .. if there is anything to plot.. was this parameter varied?
    if (min(z)==max(z)) 
        if (~isempty(which('cell2mat'))) % for older matlab versions
            fprintf(1,['Error: ',protect(cell2mat(handles.names(handles.xcol))), ...
                    ' was fixed']);
        else
            fprintf(1,'Error: parameter %i was fixed\n',handles.xcol);
        end
        return; 
    end
    
    [p x]=samp2prob1d(z,count,handles.res,'marg');
    if (length(x)>1)
        xstep=x(2)-x(1);
        pnorm=p/max(p);
        htmp=plot(x,pnorm,'r-');
        set(htmp,'color',handles.colors);
        hold on
        [pmean xmean]=samp2prob1d(z,handles.samples(2,isamps),handles.res,'mean');
        htmp=plot(xmean,pmean/max(pmean),'r:');
        set(htmp,'color',handles.colors);        
        for i=1:length(handles.clevels)
            [handles.lower(i),handles.upper(i)]= ...
                confreg(x,p,handles.clevels(i),100,'n');
            bfp=(handles.upper(i)+handles.lower(i))/2;
            err=(handles.upper(i)-handles.lower(i))/2;
            if (~isempty(which('cell2mat'))) % for older matlab versions
                fprintf(1,'%12.6f < %s < %12.6f at %3.1f percent confidence\n', ...
                    handles.lower(i),cell2mat(handles.names(handles.xcol)),handles.upper(i), ...
                    handles.clevels(i)*100)
            else
                fprintf(1,'%12.6f < Parameter %i < %12.6f at %3.1f percent confidence\n', ...
                    handles.lower(i),handles.xcol,handles.upper(i), ...
                    handles.clevels(i)*100)
            end
            v=axis;
            htmp=plot([1 1]*handles.lower(i),[v(3) v(4)],'r--');
            set(htmp,'color',handles.colors);
            htmp=plot([1 1]*handles.upper(i),[v(3) v(4)],'r--');
            set(htmp,'color',handles.colors);
           % again.. need to sort the style of this line..
        end
        %axis([axlims(ipar,1) axlims(ipar,2) 0 1.05])
        if (handles.isheld==0) hold off; end
        xlabel(handles.names(handles.xcol))
        ylabel('Probability')
        set(gca,'FontSize',handles.fontsize)
        set(get(gca,'XLabel'),'FontSize',handles.fontsize)
        set(get(gca,'YLabel'),'FontSize',handles.fontsize)
    end
    
elseif (handles.xcol==1) % have sample number along x axis
    % make a 1d plot
    % it shouldn't be possible to get to this point if zcol>1 but
    % if it happens anyway, then zcol is just ignored
    
    if (handles.ycol==1) % just plot prob against line number
        htmp=plot(handles.samples(2,isamps));
        set(htmp,'color',handles.colors);
        ylabel('-log(Likelihood)')
        tmp=size(handles.samples);
        axis tight
        if (max(handles.samples(2,isamps))>10*mean(handles.samples(2,isamps)))
            h=title('If line drops sharply maybe need to increase No. to burn')
        end      
    else 
        z=handles.samples(handles.ycol,isamps);
        % .. if there is anything to plot.. was this parameter varied?
        if (min(z)==max(z)) 
            
            if (~isempty(which('cell2mat'))) % for older matlab versions
                fprintf(1,['Error: ',protect(cell2mat(handles.names(handles.xcol))), ...
                    ' was fixed\n']);
        else
            fprintf(1,'Error: Parameter %i was fixed\n',handles.xcol);
        end
        return; 
    end
    
    htmp=plot(z);
    set(htmp,'color',handles.colors);
    ylabel(handles.names(handles.ycol))
    tmp=size(handles.samples);
%    axis([0 tmp(2) min(z) max(z)])
    axis tight
    
  end
  xlabel('Line number')
  set(gca,'FontSize',handles.fontsize)
  set(get(gca,'XLabel'),'FontSize',handles.fontsize)
  set(get(gca,'YLabel'),'FontSize',handles.fontsize)
  
else
    
    if (handles.zcol<=2) % deal with zcol=--- or zcol=prob here  
        
        % plot a 2d plot (or 3d surface)
        x=handles.samples(handles.xcol,isamps);
        y=handles.samples(handles.ycol,isamps);
        
        % do some checks:  (shouldn't be possible to get inside here from cosmogui..)
        if (handles.xcol==handles.ycol)
            fprintf(1,'Error: x and y axes are identical');
            return;
        end
        if (min(x)==max(x))               
            if (~isempty(which('cell2mat'))) % for older matlab versions
                fprintf(1,['Error: ',protect(cell2mat(handles.names(handles.xcol))), ...
                        ' was fixed\n']);
            else
                fprintf(1,'Error: Parameter %i was fixed\n',handles.xcol);
            end
            return; 
        end
        if (min(y)==max(y)) 
            if (~isempty(which('cell2mat'))) % for older matlab versions
                fprintf(1,['Error: ',protect(cell2mat(handles.names(handles.ycol))), ...
                        ' was fixed\n']);
            else
                fprintf(1,'Error: Parameter %i was fixed\n',handles.ycol);
            end
            return; 
        end
        
        switch handles.plottype
            
            case 1 % contours plus shading
                [xim,x1vals,x2vals]=samp2imw([x' y']',count,handles.res);
                if (sum(xim)==0)
                    fprintf(1,'ERROR: Probability is zero everywhere')
                else
                    % contour(x1vals,x2vals,xim')
                    if (handles.zcol==1)
                        handles.contours=plot2dconfg_exp(x1vals,x2vals,xim,...
                            handles.clevels,10,100,'qcr-');
 %                      set(handles.contours,'color',handles.colors)
                    elseif (handles.zcol==2)
                        % scatter plot
                        z=exp(-(handles.samples(2,isamps)-...
                            min(handles.samples(2,isamps))));
                        fprintf(1,'Making scatter plot...\n')
                        % handles.scatter=plot3(x,y,z,'b.','markersize',handles.markersize);
                        handles.scatter=scatter(x,y,5,z,'filled');
                        view(2); box on
                    else
                        fprintf(1,'ERROR: zcol>2 in wrong bit of code\n');
                    end              
                end
                
            case 2 
                [xim,x1vals,x2vals]=samp2imw([x' y']',count,handles.res);
                if (sum(xim)==0)
                    fprintf(1,'ERROR: Probability is zero everywhere\n')
                else
                    if (handles.zcol==1)
                        % contours only if zcol=1
                        % contour(x1vals,x2vals,xim')
                        handles.contours=plot2dconfg_exp(x1vals,x2vals,xim,...
                            handles.clevels,4,100,'qbwr-');
                      set(handles.contours,'color',handles.colors)
                    elseif (handles.zcol == 2)
                        fprintf(1,'Making surface plot ...\n')
                        handles.surface=surf2dconfg_exp(x1vals,x2vals,xim,...
                            handles.clevels,4,100,'qcr-');
                        view(3); axis tight; axis vis3d;
                    else
                        fprintf(1,'ERROR: zcol>2 in wrong bit of code\n');
                    end              
                end              
                
            case 3
                if (handles.zcol==1)
                    % shading only 
                    clflogo; % clear fig first, as overwrites anyway, and otherwise messes up cmap
                    [xim,x1vals,x2vals]=samp2imw([x' y']',count,handles.res);
                    % contour(x1vals,x2vals,xim')
                    handles.contours=plot2dconfg_exp(x1vals,x2vals,xim,...
                        handles.clevels,4,100,'nc');
%                      set(handles.contours,'color',handles.colors)
                elseif (handles.zcol==2)
                    % histogram
                    [xim,x1vals,x2vals]=samp2imw([x' y']',count,handles.res);
                    handles.image=surf(x1vals,x2vals,xim');
                    shading flat; box on; 
                    view(3); grid off; axis tight; axis vis3d
                    set(gca,'ydir','normal')
                else
                    fprintf(1,'ERROR: zcol>2 in wrong bit of code\n');
                end
                
            case 4
                if (handles.zcol==1)
                    % scatter plot
                    handles.scatter=plot(x,y,'b.','markersize',handles.markersize);
                    %          handles.scatter=scatter(x,y,20,handles.samples(2,isamps),'x');
                      set(handles.scatter,'color',handles.colors)
                else
                    fprintf(1,'ERROR: zcol>2 in wrong bit of code\n');
                end
                
            case 5 % histogram
                if (handles.zcol==1)
                    clflogo
                    [ximl,x1vals,x2vals]=samp2imw([x' y']', ...
                        exp(-handles.samples(2,isamps)+min(handles.samples(2,isamps))).*count,handles.res);
                    [xim,x1vals,x2vals]=samp2imw([x' y']',count,handles.res);
                    ximl(find(xim>0))=ximl(find(xim>0))./xim(find(xim>0));
                    handles.image=imagesc(x1vals,x2vals,ximl');
                    %          handles.image=surf(x1vals,x2vals,ximl');
                    %          shading flat; box on; view(2); grid off
                    set(gca,'ydir','normal')
                else
                    fprintf(1,'ERROR: zcol>2 in wrong bit of code\n');
                end
                
        end
        set(gca,'layer','top')
        box on
        axis tight
        % set(handles.contours,'color',handles.colour) % sim linestyle...
        % axis([axlims(ipar2,1) axlims(ipar2,2) axlims(ipar1,1) axlims(ipar1,2)])
        xlabel(handles.names(handles.xcol))
        ylabel(handles.names(handles.ycol))
        zlabel('Unnormalised Probability')
        set(gca,'FontSize',handles.fontsize)
        set(get(gca,'XLabel'),'FontSize',handles.fontsize)
        set(get(gca,'YLabel'),'FontSize',handles.fontsize)
        set(get(gca,'ZLabel'),'FontSize',handles.fontsize)
        
    else
        % do a 3d plot
        x=handles.samples(handles.xcol,isamps);
        y=handles.samples(handles.ycol,isamps);
        z=handles.samples(handles.zcol,isamps);
        % should do some checks on this really.
        
        switch handles.plottype
            
            case 1 % scatter
                fprintf(1,'Making scatter plot...\n')
                % zcol>2 so do a 3d plot
                handles.scatter=scatter(x,y,5,z,'filled');
                box on
                axis tight
                xlabel(handles.names(handles.xcol))
                ylabel(handles.names(handles.ycol))
                handles.colorbar=colorbar;
                set(get(handles.colorbar,'Ylabel'),'String',...
                    handles.names(handles.zcol));
                set(gca,'FontSize',handles.fontsize)
                set(get(gca,'XLabel'),'FontSize',handles.fontsize)
                set(get(gca,'YLabel'),'FontSize',handles.fontsize)
                set(handles.colorbar,'FontSize',handles.fontsize)
                
            case 2 % 68 % isosurface
                % zcol>2 so do a 3d plot
                
                % make a 3d histogram
                [xim,xn]=samp3imw([x' y' z']',handles.samples(1,isamps),handles.res);
                % plot 
                handles.iso=plot3dconf(xn(1,:),xn(2,:),xn(3,:),xim,0.68,5,1);
%               set(handles.iso,'color',handles.colors)
                % sort out the axis labels etc
                xlabel(handles.names(handles.xcol))
                ylabel(handles.names(handles.ycol))
                zlabel(handles.names(handles.zcol))
                set(gca,'FontSize',handles.fontsize)
                set(get(gca,'XLabel'),'FontSize',handles.fontsize)
                set(get(gca,'YLabel'),'FontSize',handles.fontsize)

            case 3 % 95 % isosurface
                % zcol>2 so do a 3d plot
                % make a 3d histogram
                [xim,xn]=samp3imw([x' y' z']',handles.samples(1,isamps),handles.res);
                % plot 
                handles.iso=plot3dconf(xn(1,:),xn(2,:),xn(3,:),xim,0.95,5,1);
%               set(handles.iso,'color',handles.colors)
                % sort out the axis labels etc
                xlabel(handles.names(handles.xcol))
                ylabel(handles.names(handles.ycol))
                zlabel(handles.names(handles.zcol))
                set(gca,'FontSize',handles.fontsize)
                set(get(gca,'XLabel'),'FontSize',handles.fontsize)
                set(get(gca,'YLabel'),'FontSize',handles.fontsize)

        end
        

    end
end

axis tight; box on; % can't think when we don't want this?

%status=close(hwait);

% Sort out the color bar for matlab releases >14
if(1)
if (str2num(version('-release'))>=14)
    % attempt to sort out colorbar width and fit colorbar text on
    kids=get(gcf,'children');
    if (strcmp(get(kids(1),'tag'),'Colorbar'))
      tmp=get(kids(1),'position');
      set(kids(1),'position',[tmp(1) tmp(2) 0.05 tmp(4)]);
      tmp1=get(kids(1),'outerposition');
      lhscbar=1-tmp1(3);
      set(kids(1),'outerposition',[lhscbar tmp1(2) tmp1(3) tmp1(4)]);
      tmp=get(kids(2),'outerposition');
      set(kids(2),'outerposition',[tmp(1) tmp(2) (lhscbar-tmp(1)) tmp(4)]);
    end
end
end

