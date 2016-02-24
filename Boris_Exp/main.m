close all
clc
clear

% parameters of the learning algorithm
sigma = 0.2;
cutoff = 0;

% reading the data file
fileID = fopen('data/data_green_4.txt','r');
% fileID = fopen('data/data_white_4.txt','r');
% fileID = fopen('data/data_desk_4_3.txt','r');
formatSpec = '%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f';
size_data = [25 inf];
data = fscanf(fileID,formatSpec,size_data);
fclose(fileID);

% extracting the variables from the data structure
time = data(1,:)';
P = data(2:4,:);
V = data(5:7,:);
w = data(8:10,:);
Rot = data(11:19,:);
F = data(20:25,:);

% contact points on the foot
m = 5;
r(:,1) = [0.05;0.05;0];
r(:,2) = [-0.05;0.05;0];
r(:,3) = [0.05;-0.05;0];
r(:,4) = [-0.05;-0.05;0];

% calculating P and V for the points on the foot
for i = 1:size(data,2)
    Ri = [Rot(1,i), Rot(4,i), Rot(7,i);
        Rot(2,i), Rot(5,i), Rot(8,i);
        Rot(3,i), Rot(6,i), Rot(9,i)];
    for j = 1:(m-1)
        P(3*j+1:3*j+3,i) = P(1:3,i) + Ri*r(:,j);
        omegai = skew(w(1:3,i));
        V(3*j+1:3*j+3,i) = V(1:3,i) + omegai*Ri*r(:,j);
    end
end

% learning part
st_no = 1500;
end_no = size(P,2);
window = 250;
selection = linspace(3,3*m,m);  %selection matrix to choose z positions

F_predict(:,1:st_no) = F(:,1:st_no);
f_predict(1:3,1:st_no) = F(1:3,1:st_no);
for i = st_no:(end_no-1)
    st_ind = max(1, i-window);
    
    % estimating Fx
    posdiff = diff(P(selection-2,st_ind:i)');
    veldiff = diff(V(selection-2,st_ind:i)');
    Diff_X_Xd = [posdiff, veldiff];
    Diff_Fx = diff(F(1,st_ind:i))';        
    lwrmodel_fx = LWRModel(Diff_X_Xd, Diff_Fx, sigma, cutoff);
    delta_x = [ (P(selection-2,i+1) - P(selection-2,i))',...
        (V(selection-2,i+1) - V(selection-2,i))'];
    F_predict(1,i+1) = F(1,i) + lwrmodel_fx.predict(delta_x);
    
    % estimating Fy
    posdiff = diff(P(selection-1,st_ind:i)');
    veldiff = diff(V(selection-1,st_ind:i)');
    Diff_Y_Yd = [posdiff, veldiff];
    Diff_Fy = diff(F(2,st_ind:i))';        
    lwrmodel_fy = LWRModel(Diff_Y_Yd, Diff_Fy, sigma, cutoff);
    delta_y = [ (P(selection-1,i+1) - P(selection-1,i))',...
        (V(selection-1,i+1) - V(selection-1,i))'];
    F_predict(2,i+1) = F(2,i) + lwrmodel_fy.predict(delta_y);
    
    % estimating Fz
    posdiff = diff(P(selection,st_ind:i)');
    veldiff = diff(V(selection,st_ind:i)');
    Diff_Z_Zd = [posdiff, veldiff];
    Diff_Fz = diff(F(3,st_ind:i))';        
    lwrmodel_fz = LWRModel(Diff_Z_Zd, Diff_Fz, sigma, cutoff);
    delta_z = [ (P(selection,i+1) - P(selection,i))',...
        (V(selection,i+1) - V(selection,i))'];
    F_predict(3,i+1) = F(3,i) + lwrmodel_fz.predict(delta_z);
    
    Diff_XYZ_ZYZD = [Diff_X_Xd, Diff_Y_Yd, Diff_Z_Zd];
    delta_xyz = [delta_x, delta_y, delta_z];
    
    % estimating Taux
    Diff_Tx = diff(F(4,st_ind:i))';        
    lwrmodel_tx = LWRModel(Diff_XYZ_ZYZD, Diff_Tx, sigma, cutoff);
    F_predict(4,i+1) = F(4,i) + lwrmodel_tx.predict(delta_xyz);
    
    % estimating Tauy
    Diff_Ty = diff(F(5,st_ind:i))';        
    lwrmodel_ty = LWRModel(Diff_XYZ_ZYZD, Diff_Ty, sigma, cutoff);
    F_predict(5,i+1) = F(5,i) + lwrmodel_ty.predict(delta_xyz);
    
    % estimating Tauz
    Diff_Tz = diff(F(6,st_ind:i))';        
    lwrmodel_tz = LWRModel(Diff_XYZ_ZYZD, Diff_Tz, sigma, cutoff);
    F_predict(6,i+1) = F(6,i) + lwrmodel_tz.predict(delta_xyz);
    
    % estimating Fx, Fy and Fz with one point only
    % estimating Fx
    Diff_X_Xd = [diff(P(1,st_ind:i)'), diff(V(1,st_ind:i)')];
    Diff_Fx = diff(F(1,st_ind:i))';        
    lwrmodel_fx1 = LWRModel(Diff_X_Xd, Diff_Fx, sigma, cutoff);
    delta_x = [ P(1,i+1) - P(1,i), V(1,i+1) - V(1,i)];
    f_predict(1,i+1) = F(1,i) + lwrmodel_fx1.predict(delta_x);
    
    % estimating Fy
    Diff_Y_Yd = [diff(P(2,st_ind:i)'), diff(V(2,st_ind:i)')];
    Diff_Fy = diff(F(2,st_ind:i))';        
    lwrmodel_fy1 = LWRModel(Diff_Y_Yd, Diff_Fy, sigma, cutoff);
    delta_y = [ P(2,i+1) - P(2,i), V(2,i+1) - V(2,i)];
    f_predict(2,i+1) = F(2,i) + lwrmodel_fy1.predict(delta_y);
    
    % estimating Fz
    Diff_Z_Zd = [diff(P(3,st_ind:i)'), diff(V(3,st_ind:i)')];
    Diff_Fz = diff(F(3,st_ind:i))';        
    lwrmodel_fz1 = LWRModel(Diff_Z_Zd, Diff_Fz, sigma, cutoff);
    delta_z = [ P(3,i+1) - P(3,i), V(3,i+1) - V(3,i)];
    f_predict(3,i+1) = F(3,i) + lwrmodel_fz1.predict(delta_z);
end


figure
for i = 1:3
    subplot(3,1,i);
    plot(time(1:end_no), F_predict(i,1:end_no), 'r');
    hold on
    plot(time(1:end_no), f_predict(i,1:end_no), 'k');
    plot(time(1:end_no), F(i,1:end_no),'b');
    legend('predicted', 'predicted1', 'actual');
    ylabel('force (N)');
end

figure
for i = 1:3
    subplot(3,1,i);
    plot(time(1:end_no), F_predict(i+3,1:end_no), 'r');
    hold on
    plot(time(1:end_no), F(i+3,1:end_no),'b');
    legend('predicted', 'actual');
    ylabel('torque (N.m)');
end
