function lat = meridianarcInverse(M,ellip,varargin)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  lat = meridianarcInverse(M,ellip,varargin)
%
%  DESCRIPTION: The meridian arc is the distance between two points on the earth 
%  with the same longitude. This function calculates the latitude of a point by 
%  the meridian arc defined between this point and the equator.
%
%  INPUT VARIABLES
%  - M: Meridian arc between the point and the equator [m]
%  - ellip: parameters for the reference ellipsoid [m]. Vector of the form
%    [semi-major axis , semi-minor axis] 
%  - formula (varargin{1}): formulation for the estimation of the latitude of a 
%    point knowing its meridian arc to the equator (Inverse Meridian Arc Formula). 
%    Two methods can be selected:
%    ¬ 'Bowring'  -> Bowring series (default method)
%    ¬ 'Snyder' -> Series used in the USGS Map Projections working manual.
%    ¬ 'UTM' -> Iterative method with paramters from UTM transform.
%    ¬ 'OSGB' -> Iterative method with parameters from OSGB transform.
%
%  OUTPUT VARIABLES
%  - lat: latitude of the point [deg]
%
%  INTERNALLY CALLED FUNCTIONS
%  - None
%
%  REFERENCES
%  - http://en.wikipedia.org/wiki/Transverse_Mercator:_Bowring_series
%  - http://en.wikipedia.org/wiki/Transverse_Mercator:_Redfearn_series
%  - Snyder, John P. (1987). Map Projections – A Working Manual. U.S. Geological 
%    Survey Professional Paper 1395. United States Government Printing Office, 
%    Washington, D.C.

%  VERSION HISTORY
%  ===============
%  VERSION 2.0.1: 13 Jan 2020
%  - Updated Bowring routine to accept negative values of meridian arc.
%
%  VERSION 2.0: 07 Jan 2020
%  - Added iterative calculation of the latitude for datums UTM and OSGB. 
%
%  VERSION 1.0: 08 Jul 2014
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

switch nargin
    case {0 1}
        error('Not enough input arguments')
    case 2
        formula = 'UTM';     
    case 3
        formula = varargin{1};
    otherwise
        error('Too many input arguments')
end

% Parameters
Mmax = 10001965.72931272; % maximum meridian distance [m]
a = ellip(1); % semi-major axis of the ellipsoid [m]
b = ellip(2); % semi-minor axis of the ellipsoid [m]
n = (a-b)/(a+b); % third flattening

% Pre-Stored Constants
n2 = n^2; % pre-store constant
n3 = n^3; % pre-store constant
n4 = n2*n2; % pre-store constant
n5 = n2*n3; % pre-store constant

switch formula
    case 'Bowring'
        if abs(M) > Mmax
            warning(['Meridian distance of %0.3f [m] is higher than the maximum '...
                'allowed. Maximum value of %0.3f [m] will be assumed'],M,Mmax)
            M = Mmax;
        end

        theta = M*(1+n)/(a*(1+n2/8)^2);
        p = 1 - 33/20*n*cos(2*theta);
        q = 33/20*n*sin(2*theta);
        Z = 5/4*(1 - 9/16*n2)*(p +q*1i)^(8/33);
        psi = theta + imag(Z);

        latRadians = atan(tan(psi)*(1+n)/(1-n));

    case 'Snyder'
        B0 = b*(1 + n + 5/4*(n2 + n3));
        Mp = pi/2*B0; % Meridian distance from the equator to the pole
        mu = pi/2*(M/Mp);
        
        D2 = 3/2*n - 27/32*n3;
        D4 = 21/16*n2 - 55/32*n4;
        D6 = 151/96*n3;
        D8 = 1097/512*n4;
        
        latRadians = mu + D2*sin(2*mu) + D4*sin(4*mu) + D6*sin(6*mu) + D8*sin(8*mu);
        
    case {'UTM','OSGB'}
        % Calculation of Parameter B0
        if strcmp(formula,'UTM')
            B0 = a*(1 - n + 5/4*(n2 - n3) + 81/64*(n4 - n5)); 
        else % 'OSGB'
            B0 = b*(1 + n + 5/4*(n2 + n3));
        end
        
        % Iterative Calculation of the Meridian Distance
        latRadians = M/B0;
        latRadians0 = latRadians;
        lat0 = latRadians0 * 180/pi;
        M0 = meridianarc(lat0,[a b],formula);
        tol = 1e-5; % tolerance [m]
        while abs(M-M0) > tol
            latRadians = latRadians0 + (M - M0)/B0;
            latRadians0 = latRadians;
            lat0 = latRadians0 * 180/pi;
            M0 = meridianarc(lat0,[a b],formula);
        end   
end
        
lat = latRadians*180/pi;

end