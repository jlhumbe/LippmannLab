function [xout,yout,iout,jout,zout,XSdist,XSoffset]=get_xsect(X,Y,Z,nXsect)
% Josh Humberston 2018 UNH/USACE
%  
% Based on the input grid parameters (X,Y,Z - for which Z must be a matrix), this code allows the user to
% define beginning and end points of a desired number of cross sections
% (nXsect) 
%
% Input : X (grid x values in meters) - example Eastings
%         Y (grid y values in meters) - example Northings
%         Z (grid z values) - example depth/bathymetry/intensity/etc.
%         nXsect (number of Cross Sections to make) - example 3
%
% Output : xout     - cross section x locations
%          yout     - cross section y locations
%          iout     - cross section i index
%          jout     - cross section j index
%          zvalue   - cross section values ; example: depth
%          XSdist   - distance of each zvalue relative to start of cross
%                       section
%          XSoffset - the distance the grid nodes used are offset from the
%                        exact cross-section line (a sort of error measurements)
%
%
% ********** Output Use Examples **********
% % Ex 1) plots z values of cross section #1 at distances from starting point:
%          plot(XSdist{1},zout{1}) % 
%
% % Ex 2) plots offsets of cross section #3 based on from starting point:
%          plot(XSdist{3},XSoffset{3}) % 
%
% % Ex 3) If Z_3d is a 3-D matrix, create a time stack along cross section #1:
%          [xout,yout,iout,jout,zout,XSdist,XSoffset]=get_xsect(X,Y,Z_3d(:,:,1),1);
%               % and pick x-section...
%          [~,~,layers] = size(Z_3d);
%          xi=iout{1}; yi=iout{1};
%          timestack = ones(layers,length(xi)).*NaN;
%          for id=1:length(iout{1});
%               timestack(:,id) = squeeze(Z_3d(xi(id),yi(id),:));
%          end
%          pcolor(XSdist{1},1:layers,timestack);shading flat
%          xlabel('X-Section Distance')

%% ***************\\\\\\\\\\ CODE //////////*************** %%
% if x/y aren't matricies, make them so
if isvector(X) && isvector(Y)
    [X,Y] = meshgrid(X,Y);
end

% Find grid spacing
dx = nanmean(nanmean(diff(X,1,2)));
dy = nanmean(nanmean(diff(Y,1,1)));
dxy = sqrt(dx^2 + dy^2);

% Plot Figure
figure;
pcolor(X,Y,Z)
shading('interp')
hold on

% Preallocate space
x_XSends = ones(nXsect,2).*NaN;
y_XSends = ones(nXsect,2).*NaN;

% Run through each cross section
for id = 1:nXsect
    for xsid=1:2 % Allow user to pick points
        [x_XSends(id,xsid),y_XSends(id,xsid)] = ginput(1);
        if xsid==1 % plot first point only
            plot(x_XSends(id,xsid),y_XSends(id,xsid),'.--k','markersize',20,'linewidth',1.5);
        else % plot both points with a dashed line between
            plot(x_XSends(id,:),y_XSends(id,:),'.r','markersize',20);
            plot(x_XSends(id,:),y_XSends(id,:),'--k','linewidth',0.5);
        end
    end
    
    % determine angle of line drawn
    difang = atan2(y_XSends(id,2)-y_XSends(id,1),x_XSends(id,2)-x_XSends(id,1));
    % determine incredments at which to find points along this line. 
    % This method used is a bit clunky... but works. It simply uses increments based 
    % on half the diagonal length of each grid node to increase sampling
    % percision by a factor of 2.
    xinc = cos(difang)*(dxy/2);
    yinc = sin(difang)*(dxy/2);
    x_XSraw = x_XSends(id,1):xinc:x_XSends(id,2);
    y_XSraw = y_XSends(id,1):yinc:y_XSends(id,2);
    
    % Go through increments and find corresponding groid points and values
    % that are closest
    for idfill=1:length(x_XSraw)
        xdiff = X-x_XSraw(idfill);
        ydiff = Y-y_XSraw(idfill);
        totdiff = sqrt(xdiff.^2+ydiff.^2);
        [iraw(idfill,1), jraw(idfill,1)] = find(abs(totdiff)==nanmin(nanmin(abs(totdiff))));
        zraw(idfill,1) = Z(iraw(idfill),jraw(idfill));
    end
    
    % Due to spacing, there may be repeats, so remove those
    [~,ia,~] = unique([iraw jraw],'rows');
    zout{id} = zraw(ia);
    iout{id} = iraw(ia);
    jout{id} = jraw(ia);

    % Get x/y values associated with indicies
    for idfill=1:length(iout{id})
        xraw(idfill) = X(iraw(ia(idfill)),jraw(ia(idfill)));
        yraw(idfill) = Y(iraw(ia(idfill)),jraw(ia(idfill)));
    end
    
    % put values in cell array for output
    xout{id}=xraw;
    yout{id}=yraw;
    
    % plot dots on the grid nodes used for cross section so user can see if
    % its whack or not
    plot(xout{id},yout{id},'.r')
    XSdist{id} = sqrt((xout{id}-x_XSraw(1)).^2 + (yout{id}-y_XSraw(1)).^2);
    XSoffset{id} = sqrt((xout{id}-x_XSraw(ia)).^2 + (yout{id}-y_XSraw(ia)).^2);
    clear xraw yraw iraw jraw
end
end
