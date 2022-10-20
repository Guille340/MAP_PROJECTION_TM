function [N,E] = tmDirect_Redfearn(lat,lon,lat0,lon0,N0,E0,k0,ellipsoidName)

%  [N,E] = tmDirect_Redfearn(lat,lon,lat0,lon0,N0,E0,k0,ellipsoidName)
%
%  DESCRIPTION: Calculates the projected coordinates (Easting and Northing)
%  of a point from its geographic coordinates (latitude and longitude)
%  using the Transverse Mercator Projection (TMP). The function applies the 
%  Redfearn series to calculate the projected coordinates. Apart from the
%  datum (latitude, longitude and reference ellipsoid) it is necessary to 
%  specify the parameters for the TMP (grid origin, false Easting/Northing 
%  and scale factor).
%
%  INPUT VARIABLES
%  - lat: latitude of the point [deg]
%  - lon: longitude of the point [deg]
%  - lat0: latitude of the grid origin [deg]
%  - lon0: longitude of the grid origin [deg]
%  - N0: false Northing [m]. In UTM projection N0 = 10,000,000 m for
%    South hemisphere and N0 = 0 m for North hemisphere.
%  - E0: false Easting [m]. In UTM projection E0 = 500,000 m
%  - k0: scale factor. In UTM projection k0 = 0.9996
%  - ellipsoidName: name of reference ellipsoid. Check refEllip.m help for 
%    string options
%
%  OUTPUT VARIABLES
%  - N: Northing [m]
%  - E: Easting [m]
%
%  INTERNALLY CALLED FUNCTIONS
%  - refEllip
%  - meridianarc
%
%  CONSIDERATIONS & LIMITATIONS
%  - Error < 1 mm for a longitude difference lambda < 3 deg.
%  - The OSGB datum is used for the calculation of the meridian arc, based
%    on the example provided in the Wikipedia article about the Redfearn
%    series.
%
%  REFERENCES
%  - http://en.wikipedia.org/wiki/Transverse_Mercator:_Redfearn_series 
%  - http://en.wikipedia.org/wiki/Transverse_Mercator_projection

%  VERSION HISTORY
%  =============== 
%  VERSION 1.0.0, 09 Jan 2020
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
% _______________________________

latRadians = lat*pi/180; % latitude of point [rad]
lonRadians = lon*pi/180; % longitude of point [rad]
lon0Radians = lon0*pi/180; % longitude of grid origin [rad]
lambda = (lonRadians - lon0Radians); % longitude difference [rad]

% Ellipsoid Parameters
[a,b,f,~,~] = refEllip(ellipsoidName); % parameters of the ref. ellipsoid
ec = sqrt(f*(2 - f)); % eccentricity
M = meridianarc(lat,[a b],'OSGB'); % meridian arc at point [m]
M0 = meridianarc(lat0,[a b],'OSGB'); % meridian arc at grid origin [m]
v = a/sqrt(1 - (ec*sin(latRadians))^2); % radius of curvature in prime vertical                 
p = v^3*(1-ec^2)/a^2;                  ... at footpoint latitude
beta = v/p;

% Pre-Stored Constants
tanLat = tan(latRadians);
tanLatPow2 = tanLat * tanLat;
tanLatPow4 = tanLatPow2 * tanLatPow2;
tanLatPow6 = tanLatPow4 * tanLatPow2;
lambdaCos = lambda*cos(latRadians);
betaPow2 = beta * beta;
betaPow3 = betaPow2 * beta;

% Direct Formula (Constants)
W3 = beta - tanLatPow2;
W4 = 4*betaPow2 + beta - tanLatPow2;
W5 = 4*betaPow3*(1 - 6*tanLatPow2) + betaPow2*(1 + 8*tanLatPow2) ...
   - 2*beta*tanLatPow2 + tanLatPow4;
W6 = 8*betaPow2^2*(11 - 24*tanLatPow2) - 28*betaPow3*(1 - 6*tanLatPow2) ...
   + betaPow2*(1 - 32*tanLatPow2) - 2*beta*tanLatPow2 + tanLatPow4;
W7 = 61 - 479*tanLatPow2 + 179*tanLatPow4 - tanLatPow6;
W8 = 1385 - 3111*tanLatPow2 + 543*tanLatPow4 - tanLatPow6;

% Direct Formula
y = k0*M + k0*v*tanLat*(lambdaCos^2/2 + lambdaCos^4*W4/24 ...
  + lambdaCos^6*W6/720 + lambdaCos^8*W8/40320); % true Northing (+/-)
x = k0*v*(lambdaCos + lambdaCos^3*W3/6 + lambdaCos^5*W5/120 ...
  + lambdaCos^7*W7/5040); % true Easting (+/-)

% Grid
N = y + N0 - k0*M0; % corrected Northing (+, ref. to grid origin)
E = x + E0; % corrected Easting (+)

