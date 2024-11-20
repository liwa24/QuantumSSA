function mjd2000 = epoch2mjd2000(year,fracYear)
% Description
%   Converter from epoch format to MJD2000
%
% Prototype
%   mjd2000 = epoch2mjd2000(year,fracYear)
%
% Input
%   year - year number
%   fracYear - days of the year [num of days]
%
% Output
%   mjd2000 - MJD2000 format of the date [num of days]
%
% Functions used
%   fracday2hms
%   month
%   fracday2hms
%   date2mjd2000
%
% Author and Version
%	Ver. 1 - W. Litteri - 01-2024


day = floor(fracYear);
fracDay = fracYear - day;
m = month(day);
[h,min,s] = fracday2hms(fracDay);

mjd2000 = date2mjd2000([year,m,day,h,min,s]);