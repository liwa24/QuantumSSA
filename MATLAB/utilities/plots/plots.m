function [orbit_dep, orbit_arr, orbit_transfer] = plots(mj2000_departure_window, mj2000_arrival_window, ...
    dv_matrix, TOF, trans_opt, par_opt, trans_opt_constr, par_opt_constr, levels, tag_P1, tag_P2)

mu_sun = astroConstants(4);
AU = astroConstants(2);
dep = trans_opt(1,:);
arr = trans_opt(2,:);
dep_constr = trans_opt_constr(1,:);
arr_constr = trans_opt_constr(2,:);
dv_opt = par_opt(1);  
t_opt = par_opt(2);
dv_opt_constr = par_opt_constr(1);  
t_opt_constr = par_opt_constr(2);

switch tag_P1
    case 1 
        str1 = ' Mercury';
    case 2
        str1 = ' Venus';
    case 3
        str1 = ' Earth';
    case 4
        str1 = ' Mars';
    case 5
        str1 = ' Jupiter';
    case 6
        str1 = ' Saturn';
    case 7
        str1 = ' Uranus';
    case 8   
       str1= ' Neptune';
    case 9
       str1 = ' Pluto';
                  
end

switch tag_P2
    case 1 
        str2 = ' Mercury';
    case 2
        str2 = ' Venus';
    case 3
        str2 = ' Earth';
    case 4
        str2 = ' Mars';
    case 5
        str2 = ' Jupiter';
    case 6 
        str2 = ' Saturn';
    case 7
        str2 = ' Uranus';
    case 8   
       str2= ' Neptune';
    case 9
       str2 = ' Pluto';
                  
end
datesx = zeros(length(mj2000_departure_window), 6);
datesy = datesx;
for i= 1:length(mj2000_departure_window)
datesx(i,:) = mjd20002date(mj2000_departure_window(i));
datesy(i,:) = mjd20002date(mj2000_arrival_window(i));
end
datesx= datenum(datesx);
datesy= datenum(datesy);
%dv_matrix(dv_matrix>levels(end)) = NaN;
figure
contour(datesx', datesy', dv_matrix, levels, 'LineColor', 'flat');

h1 = gcf;
hold on
xtickangle(45);
ytickangle(45);
datetick('x', 'yyyy mmm dd', 'keeplimits', 'keepticks');
datetick('y', 'yyyy mmm dd', 'keeplimits', 'keepticks');

TOF(TOF<0) = NaN;
contour(datesx', datesy', TOF, 'LineColor', 'black', 'showtext', 'on');
caxis([levels(1) levels(end)]);
% h4 = gcf;
h3 = colorbar;
h3.Limits = [levels(1) levels(end)];
h2 = plot(datenum(dep), datenum(arr), 'bo','MarkerFaceColor', [0 0 1], 'Markersize', 5);
h4 = plot(datenum(dep_constr), datenum(arr_constr), 'bo','MarkerFaceColor', [1 0 0], 'Markersize', 5);
%axes
%h4 = legend('delta_v transfer', 'Time of flight','delta_v opt',  'Location', 'southWest');
title(['porkchop plot of a transfer between' str1 ' and' str2]);
xlabel('Departure Date'); ylabel('Arrival Date');

grid on

%% static plot - non constrained orbit
[orbit_dep,~] = uplanet(date2mjd2000(dep), tag_P1);
[orbit_arr, ~] = uplanet(date2mjd2000(arr), tag_P2);
[r_dep, ~] = kep2car(orbit_dep, mu_sun);
[r_arr, ~] = kep2car(orbit_arr, mu_sun);
[~,~,~,~,VI,~,~,THETA] = lambertMR(r_dep,r_arr,t_opt*24*3600,mu_sun,0,0,0,1);
orbit_transfer = car2kep(r_dep, VI', mu_sun);

[Xi,Yi,Zi] = plotOrbit(orbit_dep(1:5), mu_sun, 0, 2*pi, 0.01);
[Xf, Yf, Zf] = plotOrbit(orbit_arr(1:5), mu_sun, 0, 2*pi, 0.01);
[Xt, Yt, Zt] = plotOrbit(orbit_transfer(1:5), mu_sun, orbit_transfer(6), orbit_transfer(6) + THETA, 0.01);
figure
hold on
plot3(Xi/AU, Yi/AU, Zi/AU,  'b--');
plot3(Xf/AU, Yf/AU, Zf/AU,  'r--');
plot3(Xt/AU, Yt/AU, Zt/AU,  'g-');
plot3(0,0,0, 'yo', 'MarkerFaceColor', [1 1 0], 'Markersize', 5)
axis equal
xlabel('x [AU]');
ylabel('y [AU]');
title(['non - constrained transfer orbit between' str1 ' and' str2]);
