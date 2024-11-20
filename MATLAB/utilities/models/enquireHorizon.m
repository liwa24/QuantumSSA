function [carVec, kepVec, t0, tvec, catalogued_sat] = enquireHorizon(TLE_input, start_time, stop_time, time_step)
% DESCRIPTION
% Tool for the propagation a spacecraft state (Cartesian and Keplerian) using
% the SGP4/SDP4 model of the JPL Horizons on-line ephemeris system API,
% privided an initial input in text format containing the TLE of the
% specified set of spacecraft. 
% 
% PROTOTYPE
%   [carVec, kepVec, t0, tvec, catalogued_sat] = enquireHorizon(TLE_input, start_time, stop_time, time_step)
% INPUT
%   TLE_input   [file_name.txt] text file containing the list of TLE for the n_sat
%       satellites of interest. The format is the following :
%       Header (Line 1)
%       All TLEs (Lines 2 to end, without header)    
%   start_time  [string] start time of propagation in format "YYYY-MM-DD" or "YYYY-MM-DD-hh-mm-ss"
%   stop_time   [string] end time of propagation in format "YYYY-MM-DD" or "YYYY-MM-DD-hh-mm-ss"
%   time_step   [string] time step of the propagation "XX d,m,h" (XX days, or minutes, or hours)
%
% OUTPUT
%   carVec [6,n_t,n_sat] Cartesian state vector by columns in Geocentric Inertial
%       coordinates, for each of the n_sat satellites and n_t time instants
%       [X,Y,Z,Vx,Vy,Vz] [km..., km/s...]
%   kepVec [6,n_t,n_sat] Keplerian state vector by columns in Geocentric Inertial
%       coordinates, for each of the n_sat satellites and n_t time instants
%       [a,e,i,Om,om,th] [km,#,rad]
%   t0              [1,1] start instant of the propagation in Julian date []
%   tvec            [n_t,1] instant of propagation from t0 (t-t0) in average days (1d=24h) []
%   catalogued_sat  [n_sat,1] list of the propagated satellites []
%
% DEPENDENCIES
% -
% 
% NOTE
% this function can be enlarged including all the functionalities of the
% Horizon system API ...
% more info on the Horizon tool : https://ssd.jpl.nasa.gov/horizons/manual.html
% more info on the Horizon API : https://ssd-api.jpl.nasa.gov/doc/horizons.html
%
% AUTHOR AND VERSION
%	Ver. 1 - W. Litteri - 01-2024


[TLE_input_encoded, catalogued_sat] = URLencodeTLE(TLE_input);
n_sat = length(catalogued_sat);
for k = 1:n_sat
%url_car = "2011-10-01 2011-10-14 60%20m %27SC-1%0A1%2087820U%2011053A%20%20%2011273.79990913%20%20.00099611%20%2000000-0%20%2064461-3%200%20%209991%0A2%2087820%20042.7843%20189.7738%200014383%20039.8647%20002.5266%2015.74868665%20%20%20196%0A1%2087820U%2011053A%20%20%2011273.86983630%20-.00085102%20+00000-0%20-55758-3%200%20%209998%0A2%2087820%20042.7804%20189.3478%200014258%20040.7498%20038.5752%2015.74826749000204%27";
url_car = strcat("https://ssd.jpl.nasa.gov/api/horizons.api?format=text&COMMAND=%27TLE%27&OBJ_DATA=%27NO%27&MAKE_EPHEM=%27YES%27&EPHEM_TYPE=%27VECTOR%27&CENTER=%27500@399%27&START_TIME=%27", ...
    start_time, "%27&STOP_TIME=%27", stop_time, "%27&STEP_SIZE=%27", regexprep(time_step, ' ', '%20'), "%27&TLE=", TLE_input_encoded(k));
url_kep = strcat("https://ssd.jpl.nasa.gov/api/horizons.api?format=text&COMMAND=%27TLE%27&OBJ_DATA=%27NO%27&MAKE_EPHEM=%27YES%27&EPHEM_TYPE=%27ELEMENTS%27&CENTER=%27500@399%27&START_TIME=%27", ...
    start_time, "%27&STOP_TIME=%27", stop_time, "%27&STEP_SIZE=%27", regexprep(time_step, ' ', '%20'), "%27&TLE=", TLE_input_encoded(k));

data_car = webread(url_car);
data_kep = webread(url_kep);
% eliminating the header of the URL page 
data_car = data_car((strfind(data_car, "$$SOE") + 6):(strfind(data_car, "$$EOE") -2)); 
data_kep = data_kep((strfind(data_kep, "$$SOE") + 6):(strfind(data_kep, "$$EOE") -2)); 

% create the matrix correspondent to the satellite

data_car_split = splitlines(data_car);
data_kep_split = splitlines(data_kep);


if k == 1 %initialize the data structures for n_sat spacecraft, assuming the same time instant for everyone
    n_lines_car = 4;
    n_lines_kep = 5;
    n_t = size(data_car_split,1)/n_lines_car;

    carVec = zeros(6,n_t, n_sat);
    kepVec = zeros(6,n_t, n_sat);
    tvec = zeros(n_t,1);
    t0 = str2double(data_car(1:17));
end

for j = 1:n_t

    if k == 1 %the time axis is the same for all the satellites
        tvec(j) = sscanf(data_car_split{(j-1)*n_lines_car + 1}, "%f") - t0; %instant in julian days = 24h
    end

    XYZ = sscanf(data_car_split{(j-1)*n_lines_car + 2}, " X =%f Y =%f Z =%f", [3,1]);
    VXYZ = sscanf(data_car_split{(j-1)*n_lines_car + 3}, " VX=%f VY=%f VZ=%f", [3,1]);
    carVec(:,j, k) = [XYZ;VXYZ];

    L1kep = sscanf(data_kep_split{(j-1)*n_lines_kep + 2}, " EC=%f QR=%f IN=%f", [3,1]);
    L2kep = sscanf(data_kep_split{(j-1)*n_lines_kep + 3}, " OM=%f W =%f Tp=%f", [3,1]);
    L3kep = sscanf(data_kep_split{(j-1)*n_lines_kep + 4}, " N =%f MA=%f TA=%f", [3,1]);
    L4kep = sscanf(data_kep_split{(j-1)*n_lines_kep + 5}, " A =%f AD=%f PR=%f", [3,1]);
    kepVec(:,j, k) = [L4kep(1); L1kep(1) ;deg2rad(L1kep(3)); deg2rad(L2kep(1:2)); deg2rad(L3kep(3))];
end

end






end 

function [TLE_input_encoded, catalogued_sat] = URLencodeTLE(TLE_input)
% this function allows an URL encode of the input TLE data for a correct
% interaction with the Horizon API.
% Spaces are encoded with %20 and new lines with %0A
%TLE_output_encoded = "";
fd = fopen(TLE_input);
A0 = fgetl(fd);
A1 = fgetl(fd);
A2 = fgetl(fd);

catalogued_sat = [];
j = 0;
while ischar(A2) 
    j = j + 1;
    satnum = str2double(A1(3:7));
    %assert(chksum(A1), 'Checksum failure on line 1')
    %assert(chksum(A2), 'Checksum failure on line 2')
    if ~ismember(satnum,catalogued_sat)
        %fprintf('New satellite: %d\n', satnum);
        catalogued_sat = [catalogued_sat; satnum];
        if j == 1
            TLE_input_encoded = strcat("%27SAT", num2str(satnum), "%0A", regexprep(A1, ' ', '%20'), ...
                "%0A", regexprep(A2, ' ', '%20'));
        else
            TLE_input_encoded = [TLE_input_encoded; strcat("%27SAT", num2str(satnum), "%0A", regexprep(A1, ' ', '%20'), ...
                "%0A", regexprep(A2, ' ', '%20'))];
        end
    


    else 
        ind = find(satnum == catalogued_sat); 
        TLE_input_encoded(ind) = strcat(TLE_input_encoded(ind),"%0A", regexprep(A1, ' ', '%20'), ...
            "%0A", regexprep(A2, ' ', '%20'));

    end
    A1 = fgetl(fd);
    A2 = fgetl(fd);
end
  fclose(fd);
  %close the strings
for k=1:length(catalogued_sat)
    TLE_input_encoded(k) = strcat(TLE_input_encoded(k), "%27");
end


end