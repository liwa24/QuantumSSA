clc 
clear
F = findall(0,'type','figure','tag','TMWWaitbar');
output_name = "em_orbits_try.h5";
delete(F);
close all 

% load the data to work on 
load ("files\PERIODIC ORBITS\EM\EM_IC_ARRAY.mat");
% class of the orbit (1) | Initial state (6) | Jacobi constant (1) | Period
% (1)| stability index (1)
% all the quantities are normalized depending on the type of problem. EM
% means Earth-moon.
[mu, LU,TU,VU,LPs] = constants_3BP("EM");
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13); 

% strategy to follow: each orbit is propagated for Np periods (the period
% is indicated inside the Ics array  (9th position) with N points in total

Np = 5; 
N = 1500*Np; 

% batch propagation estimation 
% just propagate a couple of orbits to estimate the time required for the
% propagation of the entire batch. 
Norb = size(out_EM,1);
nmax = 1000;
fprintf('Esimating the required time propagating for %d orbits \n', nmax)
tic
for i= 1:nmax
    j = randi(Norb);
    [t, X] = ode113(@(t, X) eq_motion_CR3BP(t, X, mu), linspace(0, Np*out_EM(j,9),N), out_EM(j,2:7), options);

end 
t_est = toc/nmax; %average propagation time per orbit 
fprintf('Estimated propagation time: %f hours \n', t_est*Norb/3600)
fprintf('Memory space required for propagation: %f GB \n', N*6*8/(1e9)*Norb + N*Norb*8/1e9)

if strcmp(input('Do you want to continue? Y/N \n', "s"),'N')
    fprintf('ciao')
    return 
end

tic
%T = zeros(N,Norb); % time nodes for each orbit
%out_orbits = zeros(N,6,Norb);
file_name = strcat("files\PERIODIC ORBITS\",output_name);
database_location = "/files/PERIODIC ORBITS";
h5create(file_name,database_location,[N 7 Norb], 'ChunkSize', [N,7,32],'Deflate',6)

start = 1;

f = waitbar(0,'1','Name','Batch propagation',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

setappdata(f,'canceling',0);
not_comp = [];
%%

for i=3755:Norb
    if getappdata(f,'canceling')
        delete(f)
        fprintf('Operation aborted by user after %f min. \n', toc/60);
        return
    end
    
    % Update waitbar and message
    waitbar(i/Norb,f,sprintf('Progress: %.2f %%', (i/Norb)*100), 'interpreter', 'none')
    % propagation 
    [t, X] = ode113(@(t, X) eq_motion_CR3BP(t, X, mu), linspace(0, Np*out_EM(i,9),N), out_EM(i,2:7), options);
    % store in output array 
    if length(t) < N 
        t = [t;zeros(N-length(t),1)];
        X = [X;zeros(N-size(X,1),6)];
        not_comp = [not_comp;i];
    end
    h5write(file_name,database_location,t,[1,1,i],[N,1,1]);
    h5write(file_name,database_location,X,[1,2,i],[N,6,1]);
    
end

t_simu = toc;
waitbar(1,f,'End')
fprintf('Operation completed. Elapsed time:  %f hours\n', t_simu/3600);
delete(f)

data = h5read(file_name,"/files/PERIODIC ORBITS");
