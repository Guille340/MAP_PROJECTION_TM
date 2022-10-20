function M = meridianarc(lat,ellip,varargin)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  M = meridianarc(lat,ellip,varargin)
%
%  DESCRIPTION: The meridian arc is the distance in the vertical direction  between
%  two points on the earth which share the same longitude. This function calculates 
%  the meridian arc between a point of latitude 'lat' and the equator (latitude=0).
%
%  INPUT VARIABLES
%  - lat: latitude of the point [deg]
%  - ellip: parameters for the reference ellipsoid [m]. Vector of the form
%    [semi-major axis , semi-minor axis] 
%  - formula (varargin{1}): formulation for the estimation of the meridian arc. 
%    There are 4 options to choose from. All of them use third flattening 'n' 
%    (Series in Third Flattening 'n'):
%    ¬ 'Bessel'  -> Bessel formulation [1837]
%    ¬ 'Helmert' -> Helmert formulation [1880] (default). Is an expanded and 
%                   simplified version of Bessel formula.
%    ¬ 'UTM'     -> Hinks formulation [1927]. Fully expanded form of Bessel series 
%                   used by the U.S. Defense Mapping Agency. Uses definition of UTM. 
%                   Slow convergence
%    ¬ 'OSGB'    -> Similar to UTM. Fully expanded series. Adopted by the Ordnance 
%                   Survey of Great Britain. Slow convergence.
%    ¬ 'Bowring' -> Bowring formulation (1983).
%
%  OUTPUT VARIABLES
%  - M: meridian arc [m]
%
%  INTERNALLY CALLED FUNCTIONS
%  - None
%
%  OTHER FORMULAS (Series in Eccentricity 'e')
%
%   M=a*[(1-ec^2/4-3*ec^4/64-5*ec^6/256)*lat
%   ... - (3*ec^2/8+3*ec^4/32+45*ec^6/1024)*sin(2*lat)
%   ... + (15*ec^4/256+45*ec^6/1024)*sin(4*lat) 
%   ... - (35*ec^6/3072)*sin(6*lat)];
% 
%  REFERENCES
%  - http://en.wikipedia.org/wiki/Meridian_arc#cite_note-13
%  - Bessel, F. W (1837). "Bestimmung der Axen des elliptischen 
%    Rotationssphäroids, welches den vorhandenen Messungen von 
%    Meridianbögen der Erde am meisten entspricht" [Estimation of the axes
%    of the ellipsoid through measurements of the meridian arc]. 
%    Astronomische Nachrichten, 14 (333):333–346

%  VERSION HISTORY
%  ===============
%
%  VERSION 1.0.1: 07 Jan 2020
%  - Renamed lat as latRadians where the latitude variable is in radians.
%  - Replaced the default formula of 'Helmert' by 'UTM'.
%
%  VERSION 1.0.0: 8 Jul 2014
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%**************************************************************************

latRadians = lat*pi/180; % latitude [rad]
a = ellip(1); % semi-major axis of the ellipsoid [m]
b = ellip(2); % semi-minor axis of the ellipsoid [m]
n = (a-b)/(a+b); % third flattening

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

n2 = n^2; % pre-store constant
n3 = n^3; % pre-store constant
n4 = n2*n2; % pre-store constant
n5 = n2*n3; % pre-store constant

switch formula
    case 'Bessel'
        D0 = 1 + 9/4*n2 + 225/64*n4;
        D2 = 3/2*n + 45/16*n3 + 525/128*n5;
        D4 = 15/16*n2 + 105/64*n4;
        D6 = 35/48*n3 + 315/256*n5;
        
        M = a*(1 - n)*(1-n2)*(D0*latRadians - D2*sin(2*latRadians) + D4*sin(4*latRadians) - D6*sin(6*latRadians));
        
        % --------------------------------------------------------------------------
        %  ALTERNATIVE BESSEL FORMULATION (see medianarc.m)
        %  D0= 1 + 9/4*n2 + 225/64*n4;
        %  K1=3/2*n - 9/16*n3;
        %  K2=15/16*n2 - 15/32*n4;
        %  K3=35/48*n3;
        %  K4=315/512*n4;
        %        
        %  r=a*(1 - n)*(1 - n2)*D0;
        %  mu=lat - K1*sin(2*lat) + K2*sin(4*lat) - K3*sin(6*lat) + K4*sin(8*lat);
        %  M=r*mu;
        % --------------------------------------------------------------------------
        
    case 'Helmert'
        H0 = 1 + n2/4 + n4/64;
        H2 = 3/2*(n - n3/8);
        H4 = 15/16*(n2 - n4/4);
        H6 = 35/48*n3;
        H8 = 315/512*n4;
        
        M = a/(1+n)*(H0*latRadians - H2*sin(2*latRadians) + H4*sin(4*latRadians) - H6*sin(6*latRadians) + H8*sin(8*latRadians));
        
    case 'UTM'
        B0 =  1 - n + 5/4*(n2 - n3) + 81/64*(n4 - n5);
        B2 = -3/2*(n - n2 + 7/8*(n3 - n4) + 55/64*n5);
        B4 =  15/16*(n2 - n3 + 3/4*(n4 - n5));
        B6 = -35/48*(n3 - n4 + 11/16*n5);
        B8 =  315/512*(n4 - n5);
        
        M = a*(B0*latRadians + B2*sin(2*latRadians) + B4*sin(4*latRadians) + B6*sin(6*latRadians) + B8*sin(8*latRadians));
        
    case 'OSGB'
        B0 =  b*(1 + n + 5/4*(n2 + n3));
        B2 = -b*(3/2*(n+ n2) + 21/16*n3);
        B4 =  b*15/16*(n2 + n3);
        B6 = -b*35/48*n3;
        
        M = B0*latRadians + B2*sin(2*latRadians) + B4*sin(4*latRadians) + B6*sin(6*latRadians);
        
    case 'Bowring'
        psi = atan((1-n)/(1+n)*tan(latRadians));
        p = 1 - 3/4*n*cos(2*psi);
        q = 3/4*n*sin(2*psi);
        Z = (1-3/8*n2)*(p + q*1i)^(2/3);
        theta = psi - imag(Z);
        
        M = a*theta/(1+n)*(1+n2/8)^2;    
end  
        
end

        


