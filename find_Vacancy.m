% By Philippe Baron
% Algorithm to locate a vacant particle in a colloidal crystal
% ** It is important to note that this algorithm assumes the input of a
%    colloidal crystal with global 6-fold symmetry (n > 2)

clear all  
close all 

%% Load Simulation Information
path = '';

filename = strcat(path,'run.txt');
run = textread(filename,'%s','delimiter','/n');
Matrix = dlmread(strcat(path,'mc_xyz_1.txt')); % Simulation coordinates and angles
Neighbors = dlmread(strcat(path,'VoronoiTess.txt')); % Particle neighbors determined by voronoi tessellation

Np = str2double(run(7)); % number of particles
r_x = str2double(run(49))/1000; % x-radius
r_y = str2double(run(51))/1000; % y-radius
r_z = r_y; 
k = r_x/r_y; % aspect ratio
z_coord = 1.0666667;
rad = [r_x r_y];

width = str2double(run(92))*r_y; % system width
height = 50; % system height
n_super = str2double(run(47)); % super-ellipse curvature parameter

Mat = [Matrix(:,3) Matrix(:,4)]*(min(rad));
angles = [Matrix(:,6) Matrix(:,7)];
box_x = str2double(run(56))*(min(rad)); 
box_y = str2double(run(59))*(min(rad));

Nstep = ((str2double(run(9))-str2double(run(25)))/str2double(run(17))); % number of simulation frames


SHOWFIGURE = 'off';
videoname = strcat(path,'video','.avi');
writer = VideoWriter(videoname);
interval = 10;
writer.FrameRate = 10;
open(writer);

frame = Nstep; % selects last simulation frame

% selects relevant coordinate information
M = Mat((frame-1)*Np+1:frame*Np,:);
neighbors = Neighbors((frame-1)*Np+1:frame*Np,:);

% set minimum interaction radius
r_max = 3*r_y;

%% DELETE FAR AWAY NEIGHBOR PARTICLES
for i = 1:Np
    neig = neighbors(i,3:end);
    adjust = find(neig);
    neig = neig(adjust);
    x1 = M(i,1);
    y1 = M(i,2);
    zeta = angles((frame-1)*Np+i,2);
    chi = angles((frame-1)*Np+i,1);
    % Find rotated particle angle from quaternion values zeta & chi
    azimuthal = atan2(2*(chi*zeta),(1-2*(zeta^2)));
    % Inverse rotation matrix
    R = [cos(-azimuthal), -sin(-azimuthal); sin(-azimuthal), cos(-azimuthal)];
    
    for j = 1:length(neig)
        x2 = M(neig(j),1);
        y2 = M(neig(j),2);
        r2 = [x2-x1 y2-y1]';
        r2 = R*r2;
        r2(1) = r2(1)/k;
        r2 = inv(R)*r2;
        % Checks if neighbor particle is outside interaction radius
        if norm(r2) > r_max
            neighbors(i,2+j)=0;
        end
    end
end
 
%% FIND VACANCY BORDERS
vac = [];
vac2 = [];
density = Np/(box_x*box_y/k); 
for i = 1:Np
    dr = 0.25;
    dth = (2*pi)/(r_max/dr);
    
    % Call to radial distribution function 'get_grt' defined below
    g_r = get_grt(i,M,angles,dr,dth,k,r_max,density,Np,frame);
    gr = sum(g_r');
    gt = sum(g_r); % Sum of the radial distrubution values over angular bins
    ind = find(gt>0);
    
    % If there are less than six non-zero values in 'gt' that particle is 
    % a border particle of some kind
    if length(ind) ~= 6
        vac = [vac 1];
    else
        vac = [vac 0];
    end
    
    % If there are exactly five non-zero values in 'gt' that particle is
    % potentially bordering a vacany
    if length(ind) == 5
        vac2 = [vac2 1];
    else
        vac2 = [vac2 0];
    end  
end

% If a particle with exactly five close packed particles around it has
% neigbors that also have exactly five close packed particles around them,
% then those neighbors are surrounding a vacancy

C5 = zeros(1,length(vac2)); % pseudo order parameter of vacancy proximity
for i = 1:Np
    neig = neighbors(i,3:end);
    adjust = find(neig);
    neig = neig(adjust);
    count = 0;
    for x = 1:length(neig)
        if vac2(i)==1 && vac2(neig(x))==1
            C5(i) = C5(i)+1;
        else
            C5(i) = C5(i);
        end
    end
end

%% APPROXIMATE VACANCIES

V = find(C5==max(C5) | C5==(max(C5)-1));
density = Np/(box_x*box_y/k);
Vxy = [];
% Examines only particles which are definitely bordering a vacancy
for i = 1:length(V)
    r_max = 3*r_y;
    dr = 0.05;
    dth = (2*pi)/(r_max/dr);
    g_r = get_grt(V(i),M,angles,dr,dth,k,r_max,density,Np,frame);
    gt = sum(g_r);
    idx = find(gt>0)*dth; % find angles at which neighbors are distributed
    peaks = zeros(1,6);
    % elliptical particle geometry favors a 6-fold symmetry so a one is
    % placed in the corresponding index of 'peaks' for each particle
    % present at the expected 60 degree increments
    for j = 1:length(idx)
        if round(idx(j))==6
            peaks(1) = 1;
        else
            peaks(round(idx(j))+1) = 1;
        end
    end
    % the index in 'peaks' which is 0 gives the direction in which the
    % vacancy is located
    gap = find(peaks==0);
    if length(gap)>1
        gap = gap(1);
    end
    angle = (gap-1)*(pi/3);
    x = M(V(i),1) + 2*r_x*cos(angle);
    y = M(V(i),2) + 2*r_y*sin(angle);
    Vxy = [Vxy [x y]'];
end

%% FIND VACANT POINTS

% Since each particle neighboring a vacancy predicts its location, this
% segment of the code averages all of the predictions of particles
% neighboring the same vacancy into one single point prediction

tol = 1;
a = length(Vxy(1,:));
y = 0;
for i = 1:a
    if y == 1
        i = 1;
    elseif y==2
        break;
    end
    count = 0;
    for j = 1:a
        d = sqrt((Vxy(1,j)-Vxy(1,i))^2+(Vxy(2,j)-Vxy(2,i))^2);
        if i~=j & d<tol
            Vxy(:,j) = (Vxy(:,i)+Vxy(:,j))/2;
            count = count + 1;
        end
    end
    if count > 0 
        Vxy(:,i) = [];
        y = 1;
        a = a-1;
    else
        y = 2;
    end
end
Nv = length(Vxy(1,:)); % Number of vacancies

%% PLOT PARTICLE SYSTEM
if Nv < 10
    scatter(Vxy(1,:),Vxy(2,:),15,[0.25 0.25 1],'filled');
end
for j = 1:Np
[mesh_x,mesh_y,mesh_z] = MeshSuperEllipse(r_x,r_y,r_z,n_super,40);
mesh_x = mesh_x*min(rad);
mesh_y = mesh_y*min(rad);
clear mesh_j
    zeta = angles((frame-1)*Np+j,2);
    eta = 0;
    chi = angles((frame-1)*Np+j,1);
    xi = 0;
    At = zeros(3,3);
    % Creation of quaternion rotation matrix to rotate the superellipse
    % mesh
    At(1,1) = -zeta.^2 + eta.^2 - xi.^2 + chi.^2;
    At(1,2) = 2*(zeta.*chi - xi.*eta);
    At(1,3) = 2*(eta.*zeta + xi.*chi);

    At(2,1) = -2*(xi.*eta + zeta.*chi);
    At(2,2) = -zeta.^2 - eta.^2 + xi.^2 + chi.^2;
    At(2,3) = 2*(eta.*chi - xi.*zeta);

    At(3,1) = 2*(eta.*zeta - xi.*chi);
    At(3,2) = -2*(xi.*zeta + eta.*chi);
    At(3,3) = zeta.^2 - eta.^2 - xi.^2 + chi.^2;

    At_i = inv(At);

    mesh_j = [At_i(1,1)*mesh_x' + At_i(1,2)*mesh_y' + At_i(1,3)*mesh_z',  ...
            At_i(2,1)*mesh_x' + At_i(2,2)*mesh_y' + At_i(2,3)*mesh_z', ...
            At_i(3,1)*mesh_x' + At_i(3,2)*mesh_y' + At_i(3,3)*mesh_z'];

    mesh_j(:,1) = mesh_j(:,1)+M(j,1);    
    mesh_j(:,2) = mesh_j(:,2)+M(j,2);   
    
    % Particles bordering vacancies will be colored shades of green and
    % crystal border particles will be red while interior particles not near a
    % vacancy will be blue
    
    if vac(j) == 1 && C5(j) == 0
        patch(mesh_j(:,1),mesh_j(:,2),[1 0.5 0.5])
    elseif C5(j) ~= 0
        patch(mesh_j(:,1),mesh_j(:,2),[0.5 C5(j)/max(C5) 0.5])
    else
        patch(mesh_j(:,1),mesh_j(:,2),[0.5 0.5 1])
    end
    hold on

    plot(mesh_j(:,1),mesh_j(:,2),'-k')
    text(M(j,1),M(j,2),num2str(j));
    hold off
end
ylim([-25 25]);
axis equal

%% FUNCTIONS
function grt = get_grt(j,temp,angles,dr,dth,k,r_max,density,N,i)
    % Function to calculate radial and angular distribution of particles 
    % around a central particle 
    grt = zeros(round(r_max/dr),(2*pi)/dth);
    p_temp(:,1) = temp(:,1)-temp(j,1);
    p_temp(:,2) = temp(:,2)-temp(j,2);
    zeta = angles((i-1)*N+j,2); 
    chi = angles((i-1)*N+j,1);
    azimuthal = atan2(2*(chi*zeta),(1-2*(zeta^2)));
    if abs(round(azimuthal,1)) == round(pi,1)
        azimuthal = 0;
    end
    R = [cos(-azimuthal), -sin(-azimuthal); sin(-azimuthal), cos(-azimuthal)];
    p_temp = (R*(p_temp'))';
    p_temp(:,1) = p_temp(:,1)./k;
    for a = dr:dr:r_max
        for b = dth:dth:(2*pi)
            count = 0;
            for z = 1:N
                x = p_temp(z,1);
                y = p_temp(z,2);
                d = sqrt(x^2+y^2);
                theta = atan2(y,x);
                area = (0.5*dth)*((a)^2-(a-dr)^2);
                if theta < 0
                    theta = (2*pi)+theta;
                end
                if z~=j && theta < b && theta >= (b-dth) && d < a && d >= (a-dr)
                    count = count + 1;
                end
            end
            grt(floor(a/dr),floor(b/dth)) = grt(floor(a/dr),floor(b/dth))+(count/area)/density;
        end
    end 
end
        
        
    
    






