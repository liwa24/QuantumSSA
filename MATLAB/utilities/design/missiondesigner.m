function [mj2000_departure_window, mj2000_arrival_window, dv_matrix, TOF, trans_opt, par_opt, trans_opt_constr, par_opt_constr,  error_matrix] = missiondesigner (dep_window, arr_window, tag_P1, tag_P2, discr_interval, v_inf)

mu_sun = astroConstants(4);
% mu_P1 = astroConstants(10 + tag_P1);  %valid for planets!
% mu_p2 = astroConstants(10 + tag_P2);

%% conversion to julian dates from gregorian calendar
dim_dep = floor((date2mjd2000(dep_window(2,:)) - date2mjd2000(dep_window(1,:)))/discr_interval);
dim_arr = floor((date2mjd2000(arr_window(2,:)) - date2mjd2000(arr_window(1,:)))/discr_interval);
dim = max(dim_dep, dim_arr);  %dovendo fare differenze temporali mi servono vettori di stesse dimensioni
mj2000_departure_window = linspace(date2mjd2000(dep_window(1,:)), date2mjd2000(dep_window(2,:)), dim);
mj2000_arrival_window = linspace(date2mjd2000(arr_window(1,:)), date2mjd2000(arr_window(2,:)), dim);
dv_matrix = zeros(dim, dim);
error_matrix = dv_matrix;
TOF = dv_matrix;
if nargin == 6
dv_matrix_constrained = TOF;
end
for i=1:dim %partenza
    for j = 1:dim %arrivo
        TOFij = (mj2000_arrival_window(j)- mj2000_departure_window(i))*24*3600; %[s]
        TOF(j,i) = TOFij/(24*3600);
        [kep_dep_i , ~] = uplanet(mj2000_departure_window(i), tag_P1);
        [kep_arr_j , ~] = uplanet(mj2000_arrival_window(j), tag_P2);
        [RIi,VIi] = kep2car(kep_dep_i , mu_sun);
        [RFj,VFj] = kep2car(kep_arr_j , mu_sun);
        
        [~,~,~,error_matrix(j,i),VIt,VFt,~,~] = lambertMR(RIi,RFj,TOFij,mu_sun,0,0,0,1);
        dv_matrix(j,i) = norm(VIt'-VIi) + norm(VFj - VFt');
        if nargin == 6
            dv_matrix_constrained(j,i) = dv_matrix(j,i);
                if norm(VIt'-VIi) > v_inf
                    dv_matrix_constrained(j,i) = NaN;
                end
        end
    end
end

[dvvec,j] = min(dv_matrix);
[dv_opt,i] = min(dvvec);
j= j(i);
dep = mjd20002date(mj2000_departure_window(i));
arr = mjd20002date(mj2000_arrival_window(j));
t_opt = TOF(j,i);

trans_opt = [dep;arr];
par_opt = [dv_opt,  t_opt];

if nargin == 6
    [dvvec,j] = min(dv_matrix_constrained);
    [dv_opt_constr,i] = min(dvvec);
    j= j(i);
    dep_constr = mjd20002date(mj2000_departure_window(i));
    arr_constr = mjd20002date(mj2000_arrival_window(j));
    t_opt_constr = TOF(j,i);
    trans_opt_constr = [dep_constr;arr_constr];
    par_opt_constr = [dv_opt_constr,  t_opt_constr];
else
    trans_opt_constr = [];
    par_opt_constr = [];
end
