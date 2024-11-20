clc 
clear
run('C:\Users\walth\OneDrive - University of Strathclyde\Documents\ThreeBodyProblem\CR3BP\halo_continuation.m');
load('C:\Users\walth\OneDrive - University of Strathclyde\Documents\ORBITGTP\files\PERIODIC ORBITS\EM\ICs\EM_IC_ARRAY.mat');

categories = {'BN'; 'BS'; 'DN'; 'DPO'; 'DRO'; 'DS'; 'L1\_A';'L1\_HN';'L1\_HS';'L1\_L';'L1\_V';'L2\_A';'L2\_HN';'L2\_HS';'L2\_L';'L2\_V';'L3\_A';'L3\_HN';'L3\_HS';'L3\_L';'L3\_V';'L4\_A';'L4\_LP';'L4\_SP';'L4\_V';'L5\_A';'L5\_LP';'L5\_SP';'L5\_V';'LPOE';'LPOW'};
figure; 
% types of orbits
histogram(out_EM(:,1),1:32)
title('Class of orbits in the database. Earth-Moon system. Total = 36071')
xticks(1.5:31.5');
xticklabels(categories); 
ylabel('Number of orbits')

% for i=1:length(categories)
%     gcf; 
%     text(i+0.5,600,categories(i),'Rotation',90, 'Interpreter','latex')
% end

TU = 382981/(24*3600); %TU in days
LU = 389703; %km
Mr = 1737.1; %km

% Lagrangian Points positions (in LU)

L1 = [0.83691513;	0.00000000;	 0.00000000];
L2 = [1.15568217;	0.00000000;	 0.00000000];
L3 = [-1.00506265;	0.00000000;	 0.00000000];
L4 = [0.48784941;	0.86602540;	 0.00000000];
L5 = [0.48784941;   -0.86602540; 0.00000000];

means = zeros(length(categories),3);
stds = zeros(length(categories),3);
maxs = zeros(length(categories),3);
mins = zeros(length(categories),3);
for i = 1:length(categories)
    means(i,:) = mean(out_EM(out_EM(:,1)==i, 8:10),1);
    stds(i,:) = std(out_EM(out_EM(:,1)==i, 8:10),1);
    maxs(i,:) = max(out_EM(out_EM(:,1)==i, 8:10),[],1);
    mins(i,:) = min(out_EM(out_EM(:,1)==i, 8:10),[],1);
end
means(:,2)= means(:,2)*TU;
stds(:,2)= stds(:,2)*TU;
maxs(:,2)= maxs(:,2)*TU;
mins(:,2)= mins(:,2)*TU;

figure; 
bar(1:31, means(:,1)); hold on
errorbar(1:31, means(:,1),stds(:,1), 'b*')
plot(1:31, mins(:,1), 'r*');
plot(1:31, maxs(:,1), 'g*');
yline(max(maxs(:,1)), 'g--');
yline(min(mins(:,1)), 'r--');
legend('Mean', 'std', 'min', 'max')
xticks(1:31');
xticklabels(categories); 
ylabel('Jacobi Constant [$LU^{2}/TU^{2}$]');
title('Distribution of the Jacobi Constant')

figure; 
bar(1:31, means(:,2)); hold on
errorbar(1:31, means(:,2),stds(:,2), 'b*')
plot(1:31, mins(:,2), 'r*');
plot(1:31, maxs(:,2), 'g*');
yline(max(maxs(:,2)), 'g--');
yline(min(mins(:,2)), 'r--');
legend('Mean', 'std', 'min', 'max')
xticks(1:31');
xticklabels(categories); 
ylabel('Period [days]');
title('Distribution of the Period')

figure; 
bar(1:31, means(:,3)); hold on
errorbar(1:31, means(:,3),stds(:,3), 'b*')
plot(1:31, mins(:,3), 'r*');
plot(1:31, maxs(:,3), 'g*');
yline(max(maxs(:,3)), 'g--');
yline(min(mins(:,3)), 'r--');
legend('Mean', 'std', 'min', 'max')
xticks(1:31');
xticklabels(categories); 
ylabel('Stability Index []');
title('Distribution of the Stability Index')


options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13); 
num_orb = 5;
figure; hold on;

Xshow = zeros(1000,6,num_orb*31);
k=1;
for i = 1:31
    min_i = find(out_EM(:,1)==i,1);
    max_i = find(out_EM(:,1)==i,1, 'last');
    plots = floor(linspace(min_i, max_i,num_orb));
    for j=plots
    [~, X] = ode113(@(t, X) eq_motion(t, X, mu), linspace(0, out_EM(j,9),1000), out_EM(j,2:7), options);
    Xshow(:,:,k) = X;
    k=k+1;
    plot3(X(:,1), X(:,2),X(:,3), 'Color',colors{i});
    end    
end
axis equal
plot3(L1(1),L1(2),L1(3));
plot3(L2(1),L2(2),L2(3));
plot3(L3(1),L3(2),L3(3));
plot3(L4(1),L4(2),L4(3));
plot3(L5(1),L5(2),L5(3));



