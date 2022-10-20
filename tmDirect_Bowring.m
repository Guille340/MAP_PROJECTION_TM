
function [N,E] = tmDirect_Bowring(lat,lon,lat0,lon0,N0,E0,k0,ellipsoidName)

%  [N,E] = tmDirect_Bowring(lat,lon,lat0,lon0,N0,E0,k0,ellipsoidName)
%
%  DESCRIPTION: Calculates the projected coordinates (Northing and Easting)
%  of a point from its geographic coordinates (latitude and longitude)
%  using the Transverse Mercator Projection (TMP). The function applies the 
%  Bowring series to calculate the projected coordinates. Apart from the
%  datum (latitude, longitude and reference ellipsoid) it is necessary to 
%  specify the parameters for the TMP (grid origin, false Easting/Northing 
%  and scale factor).
%
%  INPUT VARIABLES
%  - lat: latitude of the point [deg]
%  - lon: longitude of the point [deg]
%  - lat0: latitude of the grid origin [deg]
%  - lon: longitude of the grid origin [deg]
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
%  REFERENCES
%  - http://en.wikipedia.org/wiki/Transverse_Mercator:_Bowring_series
%  - http://en.wikipedia.org/wiki/Transverse_Mercator_projection

%  VERSION HISTORY
%  ===============
%  VERSION 1.0.1, 09 Jan 2020
%  - Variables renamed.
%  - Updated comments and help
%
%  VERSION 1.0.0, 09 Jul 2014
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
% ______________________________

latRadians = lat*pi/180; % latitude of point [rad]
lonRadians = lon*pi/180; % longitude of point [rad]
lon0Radians = lon0*pi/180; % longitude of grid origin [rad]

% Ellipsoid Parameters
[a,b,f,~,~] = refEllip(ellipsoidName); % parameters of the ref. ellipsoid
ec = sqrt(f*(2-f)); % eccentricity
ec2 = ec/(1-f); % second eccentricity
M = meridianarc(lat,[a b],'Bowring'); % meridian arc at point [m]
M0 = meridianarc(lat0,[a b],'Bowring'); % meridian arc at grid origin [m]

% Pre-Stored Constants
cosLat = cos(latRadians); % pre-stored constant
cosLatPow2 = cosLat^2; % pre-stored constant
sinLat = sin(latRadians); % pre-stored constant
lambda = (lonRadians - lon0Radians); % longitude difference [rad]
lambdaPow2 = lambda^2; % pre-stored constant

% Direct Formula (Parameters)
v = a*sqrt((1 + ec2)/(1 + ec2*cosLatPow2)); % radius of curvature in prime vertical 
z = (ec2*lambda^3*cosLat^5)/6; 
theta2 = atan( (2*sinLat*cosLat*sin(lambda/2)^2) / ...
        (sinLat^2 + cosLatPow2*cos(lambda)) );
    
% Direct Formula
y = k0*(M + v*theta2 + z*v*lambda*sinLat/4*(9 + 4*ec2*cosLatPow2 ...
    - 11*lambdaPow2 + 20*lambdaPow2*cosLatPow2)); % true Northing (+/-)
x = k0*v*(atanh(cosLat*sin(lambda)) ...
    + z*(1+lambdaPow2/10*(36*cosLatPow2-29))); % true Easting (+/-)

% Corrected Eastings & Northings
N = y + N0 - k0*M0; % corrected Northing (+, ref. grid origin)
E = x + E0; % corrected Easting (+)

