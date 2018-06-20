function [xout,yout,iout,jout,zvalue]=get_station_xyz_ij(X,Y,Z,npts,varargin)
%
% NOTE: if choosing a "same length segment case" (Cases B or D) the units
% of the segment are in the same units as X and Y. Be careful about
% latitude/longitude vs UTM!
%
% Input : X (grid x values in meters) - example Eastings
%         Y (grid y values in meters) - example Northings
%         Z (grid z values) - example depth/bathymetry
%       npts (# of points ) - example 10
%
% Variable Inputs
%
% Case A: If you want the segments to have variable length : no more
%         inputs; length(varargin)=0
%         [xout,yout,iout,jout,zvalue] = get_station_xyz_ij(X,Y,Z,npts)
%
% Case B: If you want the segments to have SAME length; length(varargin)=1;
%         Example : seglen = 100 (m)
%         [xout,yout,iout,jout,zvalue] = get_station_xyz_ij(X,Y,Z,npts,100)
%
% Case C: If you want the segments to have variable length, but use previous
%         data points; length(varargin)=2; Note: can be greater than 1
%         point
%         Example : xprev (previous x location(s))
%                   yprev (previous y location(s))
%         [xout,yout,iout,jout,zvalue] = get_station_xyz_ij(X,Y,Z,npts,xprev,yprev)
%
% Case D: If you want the segments to have SAME length, but use previous
%         data points; length(varargin)=3; Note: can be greater than 1
%         point
%         Example : seglen = 100 (m)
%                   xprev (previous x location(s))
%                   yprev (previous y location(s))
%
%         [xout,yout,iout,jout,zvalue] = get_station_xyz_ij(X,Y,Z,npts,xprev,yprev,seglen)
%
% Output : xout  - station x locations
%          yout  - station y locations
%          iout  - station i index
%          jout  - station j index
%          zvalue - station z value ; example: station depth

% Check to see if coordinates were put in as arrays, and if so, change to
% matrix
if isvector(X) && isvector(Y)
    [X,Y] = meshgrid(X,Y);
end

%%%% Set up Differnt Cases based on input %%%%
if nargin == 4 % CASE A
    xloc = []; yloc = [];
    
elseif nargin == 5 % CASE B : Same segment length
    xloc = []; yloc = [];
    seglen = varargin{1};
    
elseif nargin > 5 % CASE C : Use previous points ; variable segment length
    xloc = varargin{1}; % previous x points (xloc = xprev)
    yloc = varargin{2}; % previous y points (yloc = yprev)
    
    if nargin == 7 % CASE D : Use previous points, same segment length
        seglen = varargin{3};
    end
end

% Plot Figure
Ftemp = figure;
pcolor(X,Y,Z)
shading('interp')
hold on

% Plot previous points
plot(xloc,yloc,'ro','markersize',8)
if nargin>5
    % Plot last point
    plot(xloc(end),yloc(end),'ro','markersize',8,'markerfacecolor','r')
end

% Pick your points
xnew = ones(1,npts)*nan;
ynew = ones(1,npts)*nan;

for id = 1:npts
    [xnew(id),ynew(id)] = ginput(1);
    h = plot(xnew(id),ynew(id),'.b','MarkerSize',20);
    
    if nargin > 4
        
        if id ==1 % Your first point (id==1)
            
            if nargin == 5
                xvs = xnew(1); % already picked
                yvs = ynew(1); % already picked
                continue
            elseif nargin > 5
                xvs = xloc(1,end); % from previous
                yvs = yloc(1,end); % from previous
            end
            
        else % if not starter (id>1)
            xvs = xnew(id-1); % Starting Point X
            yvs = ynew(id-1); % Starting Point Y
        end
        
        
        xdiff = xnew(id) - xvs;
        ydiff = ynew(id) - yvs;
        totdiff = sqrt((xdiff^2) + (ydiff^2));
        
        if nargin == 5 || nargin == 7 % SAME LENGTH SEGMENTS
            if totdiff < seglen
                % Get rid of bad point
                % Pick a new point
                while totdiff < seglen
                    [xnew(id), ynew(id)] = ginput(1);
                    delete(h)
                    h = plot(xnew(id),ynew(id),'.b','MarkerSize',20);
                    xdiff = xnew(id) - xvs;
                    ydiff = ynew(id) - yvs;
                    totdiff = sqrt((xdiff^2) + (ydiff^2));
                end
                difang = atan2(ydiff,xdiff);
                disp(num2str(rad2deg(difang)))
                xnew(id) = xvs + (cos(difang)*seglen);
                ynew(id) = yvs + (sin(difang)*seglen);
                delete(h)
                h = plot(xnew(id),ynew(id),'.r','MarkerSize',20);
                
            elseif totdiff > seglen
                delete(h)
                difang = atan2(ydiff,xdiff);
                disp(num2str(rad2deg(difang)))
                xnew(id) = xvs + (cos(difang)*seglen);
                ynew(id) = yvs + (sin(difang)*seglen);
                h = plot(xnew(id),ynew(id),'.r','MarkerSize',20);
            end
        end
    end
    
end

xout = [xloc xnew];
yout = [yloc ynew];

% Should be short, but preallocate some space anyways...
iout = ones(length(xout),1).*NaN;
jout=iout;
zvalue = iout;

% Fill indicies and Z values arrays
for i=1:length(xout)
    xdiff = X-xout(i);
    ydiff = Y-yout(i);
    totdiff = sqrt(xdiff.^2+ydiff.^2);
    [iout(i), jout(i)] = find(abs(totdiff)==nanmin(nanmin(abs(totdiff))));
    zvalue(i) = Z(iout(i),jout(i));
end

close(Ftemp);
end
