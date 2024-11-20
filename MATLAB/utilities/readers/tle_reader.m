function out_data = tle_reader(file)
% DESCRIPTION
% Read satellite ephemeris data from a NORAD two-line element (TLE) text
% file and sort useful data in a struct, for each of the listed
% satellites.
%
% PROTOTYPE
%   out_data = tle_reader(file)
%
% INPUT
%   file [file_name.txt] text file containing the list of TLE for the n_sat
%       satellites of interest. The format is the following :
%       Header (Line 1)
%       All TLEs (Lines 2 to end, without header)    
%
% OUTPUT
% out_data {n_sat,1} cell array of structs relative to each satellite, for
% which the TLE have been provided. Information are on semi-major axis
% (km), eccentricity (e), Inclination, RAAN, arg. of perigee, Mean anomaly,
% and true anomaly (rad), and the epoch (in MJD2000). Each entry is a
% vector of [n_t,1], where n_t is the number of epochs provided in the TLE
% files for each of the n_sat spacecraft.
%
% DEPENDENCIES
%
% NOTES
%   All the other information can also be added in a subsequent version
%
% AUTHOR AND VERSION
%	Ver. 1 - W. Litteri - 01-2024, adapted from 
%       @ Brett Pantalone North Carolina State University

  fd = fopen(file,'r');
  if fd < 0, fd = fopen([file '.txt'],'r'); end
  assert(fd > 0,['Can''t open file ' file ' for reading.'])
  j = 0;
  A0 = fgetl(fd);
  A1 = fgetl(fd);
  A2 = fgetl(fd);
  catalogued_sat = [];
  out_data = {};
  while sum(size(A2)) > 0 
    j = j + 1;
    satnum = str2double(A1(3:7));
    assert(chksum(A1), 'Checksum failure on line 1')
    assert(chksum(A2), 'Checksum failure on line 2')
    if ~ismember(satnum,catalogued_sat)
        fprintf('New satellite: %d\n', satnum);
        catalogued_sat = [catalogued_sat; satnum];

        
        [Year, days, Incl, Omega, e, w, M, n, a, theta] = getInfo(A1,A2);
        

        tle_decod = struct('name', satnum, 'epoch', epoch2mjd2000(Year,days), 'semi_major_axis', a, 'eccentricity', e, 'inclination', Incl, 'right_ascension', Omega, ...
             'arg_of_perigee', w, 'true_anomaly', theta,'mean_anomaly', M, 'mean_motion', n);
        out_data{length(catalogued_sat),1} = tle_decod;

    else 
        ind = find(satnum == catalogued_sat); 
        sat_struct = out_data{ind};

        [Year, days, Incl, Omega, e, w, M, n, a, theta] = getInfo(A1,A2);
        sat_struct.epoch = [sat_struct.epoch; epoch2mjd2000(Year,days)]; 
        sat_struct.semi_major_axis = [sat_struct.semi_major_axis; a]; 
        sat_struct.eccentricity = [sat_struct.eccentricity; e];
        sat_struct.inclination = [sat_struct.inclination; Incl];
        sat_struct.right_ascension = [sat_struct.right_ascension; Omega];
        sat_struct.arg_of_perigee = [sat_struct.arg_of_perigee; w];
        sat_struct.mean_anomaly = [sat_struct.mean_anomaly; M];
        sat_struct.mean_motion = [sat_struct.mean_motion; n];
        sat_struct.true_anomaly = [sat_struct.true_anomaly; theta];

        out_data{ind} = sat_struct;


    end
    %A0 = fgetl(fd);
    A1 = fgetl(fd);
    A2 = fgetl(fd);
  end
  fclose(fd);
  fprintf('Processed: %d TLE entries on %d satellites. \n', j, length(catalogued_sat));
end

%%
% Checksum (Modulo 10)
% Letters, blanks, periods, plus signs = 0; minus signs = 1
function result = chksum(str)
  result = false; c = 0;
  
  for k = 1:68
    if str(k) > '0' && str(k) <= '9'
      c = c + str(k) - 48;
    elseif str(k) == '-'
      c = c + 1;
    end
  end
  if mod(c,10) == str(69) - 48
    result = true;
  end
  
end

function [Year, days, Incl, Omega, e, w, M, n, a, theta] = getInfo(A1,A2)
        
        Year = str2double(['20' A1(19:20)]);
        days = str2double(A1(21:32));
        Incl = str2double(A2(9:16))*pi/180;
        Omega = str2double(A2(18:25))*pi/180;
        e = str2double(['0.' A2(27:33)]);
        w = str2double(A2(35:42))*pi/180;
        M = str2double(A2(44:51))*pi/180;
        n = str2double(A2(53:63));
        T = 86400/n;
        a = ((T/(2*pi))^2*astroConstants(13))^(1/3);

        theta = kepler_solver('from_M', e, 1e-6, astroConstants(13), M); 

end



