function [lat,lon] = tmInverse_Redfearn(N,E,lat0,lon0,N0,E0,k0,ellipsoidName)
                   
% [lat,lon] = tmInverse_Redfearn(N,E,lat0,lon0,E0,N0,k0,ellipsoidName)
%
%  DESCRIPTION: Calculates the geographic coordinates (latitude and longitude) 
%  of a point from its projected coordinates (Easting and Northing) on a 
%  Transverse Mercator Projection (TMP). The function applies the Redfearn 
%  series to calculate the geographic coordinates. The following inputs are 
%  required: projected coordinates (Northing and Easting), TMP parameters (grid 
%  origin, false Easting/Northing, scale factor) and the reference ellipsoid.
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
%  - The OSGB iterative method is used for the calculation of the footpoint 
%    latitude, based on the example provided in the Wikipedia article about 
%    the Redfearn series.
%
%  REFERENCES
%  - http://en.wikipedia.org/wiki/Transverse_Mercator:_Redfearn_series 
%  - http://en.wikipedia.org/wiki/Transverse_Mercator_projection

%  VERSION HISTORY
%  ===============
%  VERSION 1.0: 09 Jan 2020
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
% _______________________________
 
lon0Radians = lon0*pi/180; % longitude of grid origin [rad]

% Ellipsoid Parameters
[a,b,f,~,~] = refEllip(ellipsoidName); % parameters of the ref. ellipsoid
ec = sqrt(f*(2 - f)); % eccentricity
M0 = meridianarc(lat0,[a b],'OSGB'); % meridian arc at grid origin [m]

% Easting & Northing
y = N - N0 + k0*M0; % true Northing (+/-, ref. to equator)
x = E - E0; % true Easting (+/-)

% Inverse Formula (Parameters)
Mf = y/k0; % distance from the equator to the footpoint [m]
latf = meridianarcInverse(Mf,[a b],'OSGB'); % footpoint latitude [deg]     
latfRadians = latf*pi/180; % footpoint latitude [rad]                                                   
vf = a/sqrt(1 - (ec*sin(latfRadians))^2); % radius of curvature in prime
                                  ... vertical at footpoint latitude [m]
pf = vf^3*(1-ec^2)/a^2; % radius of curvature in the meridian [m]               
betaf = vf/pf;

% Pre-Stored Constants
cosLatf = cos(latfRadians);
tanLatf = tan(latfRadians);
tanLatfPow2 = tanLatf * tanLatf;
tanLatfPow4 = tanLatfPow2 * tanLatfPow2;
tanLatfPow6 = tanLatfPow4 * tanLatfPow2;
betafPow2 = betaf * betaf;
betafPow3 = betafPow2 * betaf;
k0vf = k0*vf;

% Inverse Formula (Constants)
V3 = betaf + 2*tanLatfPow2;
V5 = 4*betafPow3*(1 - 6*tanLatf^2) - betafPow2*(9 - 68*tanLatfPow2) ...
   - 72*betaf*tanLatfPow2 - 24*tanLatfPow4;
V7 = 61 + 662*tanLatfPow2 + 1320*tanLatfPow4 + 720*tanLatfPow6;
U4 = 4*betafPow2 - 9*betaf*(1-tanLatfPow2) - 12*tanLatfPow2;
U6 = 8*betafPow2^2*(11 - 24*tanLatfPow2) ...
   - 12*betafPow3*(21 - 71*tanLatfPow2) ...
   + 15*betafPow2*(15 - 98*tanLatfPow2 + 15*tanLatfPow4) ...
   + 180*betaf*(5*tanLatfPow2 - 3*tanLatfPow4) + 360*tanLatfPow4;
U8 = - 1385 - 3633*tanLatfPow2 - 4095*tanLatfPow4 - 1575*tanLatfPow6;

% Inverse Formula
latRadians = latfRadians ...
   - x^2*betaf*tanLatf/k0vf^2*(1/2 + x^2*U4/(24*k0vf^2) ...
   + x^4*U6/(720*k0vf^4) + x^6*U8/(40320*k0vf^6)); % latitude of point [rad]
lonRadians = x/(cosLatf*k0vf) * (1 - x^2*V3/(6*k0vf^2) ...
   - x^4*V5/(120*k0vf^4) - x^6*V7/(5040*k0vf^6)) ...
   + lon0Radians; % longitude of point [rad]
        
% Geographic Coordinates
lat = latRadians*180/pi; % latitude of point [deg]
lon = lonRadians*180/pi; % longitude of point [deg]


