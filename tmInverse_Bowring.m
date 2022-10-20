
function [lat,lon] = tmInverse_Bowring(N,E,lat0,lon0,N0,E0,k0,ellipsoidName)

%  [lat,lon] = tmInverse_Bowring(N,E,lat0,lon0,N0,E0,k0,ellipsoidName)
%
%  DESCRIPTION: Calculates the geographic coordinates (latitude and 
%  longitude) of a point from its projected coordinates (Easting and
%  Northing) on a Transverse Mercator Projection (TMP). The function 
%  applies the Bowring series to calculate the geographic coordinates. 
%  The following inputs are required: projected coordinates (Northing and
%  Easting), TMP parameters (grid origin, false Easting/Northing, scale
%  factor) and the reference ellipsoid.
%
%  INPUT VARIABLES
%  - N: Northing [m]
%  - E: Easting [m]
%  - lat0: latitude of the grid origin [deg]
%  - lon0: longitude of the grid origin [deg]
%  - N0: false Northing [m]. In UTM projection N0 = 10,000,000 m for
%    South hemisphere and N0 = 0 m for North hemisphere.
%  - E0: false Easting [m]. In UTM projection E0 = 500,000 m
%  - k0: scale factor. In UTM projection k0 = 0.9996
%  - ellipsoidName: reference ellipsoid. Check refEllip.m help for string 
%    options.
%
%  OUTPUT VARIABLES
%  - lat: latitude of the point [deg]
%  - lon: longitude of the point [deg]
%
%  INTERNALLY CALLED FUNCTIONS
%  - refEllip
%  - meridianarc
%  - meridianarcInverse
%
%  CONSIDERATIONS & LIMITATIONS
%  - Error < 1 mm for a longitude difference lambda < 3 deg.
%
%  REFERENCES
%  - http://en.wikipedia.org/wiki/Transverse_Mercator:_Bowring_series
%  - http://en.wikipedia.org/wiki/Transverse_Mercator_projection
%

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

lon0Radians = lon0*pi/180; % longitude of grid origin [rad]

% Ellipsoid Parameters
[a,b,f,~,~] = refEllip(ellipsoidName); % parameters of the ref. ellipsoid
ec = sqrt(f*(2-f)); % eccentricity
ec2 = ec/(1-f); % second eccentricity
M0 = meridianarc(lat0,[a b]); % meridian arc at grid origin [m]

% Easting & Northing
y = N - N0 + k0*M0; % true Northing (+/-, ref. to equator)
x = E - E0; % true Easting (+/-)

% Inverse Formula (Parameters & Pre-Stored Constants)
Mf = y/k0; % distance from the equator to the footpoint [m]
latf = meridianarcInverse(Mf,[a b],'Bowring'); % footpoint latitude [deg]          
latfRadians = latf*pi/180; % footpoint latitude [rad]
cosLatf = cos(latfRadians); % pre-stored constant
cosLatfPow2 = cosLatf^2; % pre-stored constant
tanLatf = tan(latfRadians); % pre-stored constant
vf = a*sqrt((1 + ec2)/(1 + ec2*cosLatfPow2)); % radius of curvature in prime                                   
q = x/(k0*vf);                            ... vertical at footpoint latitude
theta4 = atan(sinh(q)/cosLatf);
theta5 = atan(tanLatf*cos(theta4));

% Inverse Formula
latRadians = (1 + ec2*cosLatfPow2)*(theta5 ...
    - ec2/24*q^4*tanLatf*(9 - 10*cosLatfPow2)) ...
    - ec2*cosLatfPow2*latfRadians; % latitude of point [rad]
lonRadians = theta4 ...
    - ec2/60*q^3*cosLatf*(10 - q^2*(4/cosLatfPow2 + cosLatfPow2)) ...
    + lon0Radians; % longitude of point [rad]

% Latitudes & Longitudes
lat = latRadians*180/pi; % latitude of point [deg]
lon = lonRadians*180/pi; % longitude of point [deg]

