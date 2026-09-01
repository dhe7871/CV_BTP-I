function colorvfield3(X,Y,Z,U,V,W,varargin)
%
% COLORVFIELD3 Colored 3D Vector Field Plotter.
%
% COLORVFIELD3(X,Y,Z,U,V,W) plots colored vectors with components (u,v,w) 
%   at the points (x,y,z). The vectors are colored using the jet color map 
%   (with the smallest vectors colored blue the largest colored red) and 
%   are divided into 32 discrete color levels.
%
% COLORVFIELD3(X,Y,Z,U,V,W,NUMBER_COLOR_LEVELS) divides the color map into 
%   NUMBER_COLOR_LEVELS dicrete color levels.
%
% Example:
%	X = linspace(-10,10,32);
%	Y = linspace(0, 0, 32);
%	Z = linspace(0, 0, 32);
%	U = linspace(-1,-5,32);
%	V = linspace(1,5,32);
%	W = linspace(0, 0, 32);
%	colorvfield3(X,Y,Z,U,V,W,16);	% plot 16 color levels
%	title('3D Vectors Plotted with Color Map (colorvfield3.m)')
%	axis([-12.25,11.75,-0.25,4.75,-0.25,0.25])
%	xlabel('X')
%	ylabel('Y')
%	zlabel('Z')
%

% Written by M. B. Sullivan, George Mason University (mbsullivan@gmail.com)
% 
% Revision 1.0. 
% Released 26 November 2007.
%

% error checking
if nargin < 6       % too few arguments
    disp('COLORVFIELD3 Error: Too few input arguments.');
    return
elseif nargin > 7   % too many arguments
    disp('COLORVFIELD3 Error: Too many input arguments.');
    return
end
                    % improperly sized input vectors
if length(X) ~=  length(Y) || length(X) ~= length(Z)...
        || length(X) ~= length(U) || length(X) ~= length(V)...
        || length(X) ~= length(W) || length(Y) ~= length(Z)...
        || length(Y) ~= length(U) || length(Y) ~= length(V)...
        || length(Y) ~= length(W) || length(Z) ~= length(U)...
        || length(Z) ~= length(V) || length(Z) ~= length(W)...
        || length(U) ~= length(V) || length(U) ~= length(W)...
        || length(V) ~= length(W) 
    disp('COLORVFIELD Error: X, Y, Z, U, V and W are not the same length!');
    return
end

% define constants
if  isempty(varargin) || length(varargin) == 0    % number of color levels not given 
    NUMBER_COLOR_LEVELS = 32;   % default value
else
    NUMBER_COLOR_LEVELS = varargin{1}(:);   % set number of color levels
end
ARRAY_LENGTH = length(X);   % length of X,Y,Z,U,V and W

% form color map
vcolormap = jet(NUMBER_COLOR_LEVELS);    % blue is smallest, red largest

% create color level bounds
vnorm = zeros(1,ARRAY_LENGTH);  % preallocate for speed
for i=1:ARRAY_LENGTH
    vnorm(i) = norm([U(i),V(i),W(i)]);
end
vmin = min(vnorm);  % smallest vector
vmax = max(vnorm);  % largest vector
vbounds = linspace(vmin, vmax, NUMBER_COLOR_LEVELS+1);    % color bounds

% plot each point in the appropriate color
for j=1:ARRAY_LENGTH
    % figure out what color level is appropriate
    color_level = NUMBER_COLOR_LEVELS;   % start at max index
    for k=1:NUMBER_COLOR_LEVELS
        if vnorm(j) >= vbounds(k)
            color_level = k;
        else
            break
        end
    end
    % color_level now holds the approprate color level
        
    % plot the vector in its approprate color
    quiver3(X(j),Y(j),Z(j),U(j),V(j),W(j),'Color',vcolormap(color_level,:))
    hold on
end
