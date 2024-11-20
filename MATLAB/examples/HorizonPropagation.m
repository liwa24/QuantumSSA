% EXAMPLE on how to generate a database of spacecraft ephemerides using the
% Horizon API, and a TLE input file. Note that the online tool only works
% with a few TLE per spacecraft, and at the date of the example, only 2 TLE
% were passed to tthe API. More TLE would cause a server down. Howewer the
% propagation is consistent for all the 8 spacecraft, part of the Starlink
% constellation

clc
clear all 

TLE_input = "new_starlink.txt";
start_time = "2024-01-21"; 
stop_time = "2024-01-31";
time_step = "1 m"; %one minute

% Propagation with Horizon API 
[carVec, kepVec, t0, tvec, catalogued_sat] = enquireHorizon(TLE_input, start_time, stop_time, time_step);

% plot the results for the 8 satellites
Terra_3D; hold on; 
for j= 1:length(catalogued_sat)
    plot3(carVec(1,:, j), carVec(2,:, j), carVec(3,:, j)); 
end
legend(["";num2str(catalogued_sat)])

%Write a csv file 
n_t = length(tvec);
newCar = zeros(n_t + 1, 6*length(catalogued_sat) + 1);
newKep = zeros(n_t + 1, 6*length(catalogued_sat) + 1);

newCar(2:end,1) = tvec;
newKep(2:end,1) = tvec;

for k = 1:length(catalogued_sat)
    newCar(2:end,(k-1)*6 + (2:7)) = transpose(carVec(:,:,k));
    newKep(2:end,(k-1)*6 + (2:7)) = transpose(kepVec(:,:,k));

end

writematrix(newCar, "files\output\cartesian_states.csv")
writematrix(newKep, "files\output\keplerian_states.csv")