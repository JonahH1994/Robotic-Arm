%% Robotic Arm:

% Number of links:
p.N = 4 ;

% Order of variables: theta_i, dtheta_i , i = 1 : N

% Define link directions:
p.i = [1;0;0] ; % er_i of each link
p.j = [0;1;0] ; % e0_i of each link
p.k = [0;0;1] ; % ez_i of each link

% Define link properties:
p.t = 0.1 ; % radius of the rod in m
p.A = pi * p.t^2 ; % Cross sectional area of rod in m

mass = 1 ; % mass of each link in Kg
p.m = mass * ones(p.N,1) ;
%p.m(1) = 0 ;

p.len = 1.5 ; % length of each link in meters
p.l = p.len * ones(p.N,1) ;
%p.l(1) = 2 ;

l_com = p.len - 0.5 ; % Position from one end of link to COM:
p.G = l_com * ones(p.N,1) ;

% Define initial inertias assuming the links are all oriented along the
% y-axis:
I1 = 1/12*p.m.*(3*p.t^2+p.len.^2) ;
I2 = 1/2*p.m*p.t^2 + p.m.*((p.len/2).^2-p.G.^2) ;
I3 = I1 ;
p.I = zeros(3*p.N,3) ;

for i = 1 : p.N
   
    I = p.I((3*i-2):3*i,:) ;
    I(1,1) = I1(i) ;
    I(2,2) = I2(i) ;
    I(3,3) = I3(i) ;
    p.I(3*i-2:3*i,:) = I ;
    
end

% DEFINE INITIAL CONDITIONS %%%%%%%%%%:

% initial angle:
th0 = zeros(3*p.N,1) ;
th0(:) = (0:3*p.N-1)' * pi/20 + pi/2 ;
w0 = zeros(3*p.N,1) ;

% Create initial rotation matrices:
R0 = zeros(3,3*p.N) ;
for  i = 1 : p.N
   
    % rotation matrix about z:
    Rr0z = [ cos(th0(3*i-2)), -sin(th0(3*i-2)), 0 ; ...
             sin(th0(3*i-2)),  cos(th0(3*i-2)), 0 ;...
                           0,                0, 1 ] ;
    
    % Rotation matrix about y:
    Rr0y = [ cos(th0(3*i-1)), 0,-sin(th0(3*i-1)) ; ...
                           0, 1,               0 ; ...
             sin(th0(3*i-1)), 0, cos(th0(3*i-1)) ] ;
        
    % Rotation matrix about x:
    Rr0x = [ 1,             0,              0 ; ...
             0, cos(th0(3*i)),  sin(th0(3*i)) ; ...
             0,-sin(th0(3*i)),  cos(th0(3*i)) ] ;
                    
    R0(:,3*i-2:3*i) = Rr0z * Rr0y * Rr0x ;
    
end

% Setup ode
tSpan = linspace(0,50,200) ; %[0 10] ;
opts = odeset( 'RelTol', 1e-10, 'AbsTol', 1e-10 ) ;
z0 = [ th0; w0; R0(:) ] ;

% integrate:
[t,xState] = ode45( @dyn, tSpan, z0, opts, p ) ;

% Simulate each link:
r = calculateState(xState,p) ;
[sz,~] = size(xState) ;

% figure
% plot3( r{1}(:,1), r{1}(:,2), r{1}(:,3) ) ;
% 
% figure
% plot3( r{2}(:,1), r{2}(:,2), r{2}(:,3) ) ;

figure 
hold on;
grid on ;

view(3) ;
axis([-p.len p.len -p.len p.len -p.len p.len]*(p.N+0.25) ) ;

handll = gobjects(1,p.N) ;

for i = 1 : p.N
    handll(i) = plot3( [ r{i}(1,1) r{i+1}(1,1)], [r{i}(1,2) r{i+1}(1,2)], ...
        [r{i}(1,3) r{i+1}(1,3)], '-o', 'LineWidth', 4, ...
        'MarkerFaceColor',[.49 1 .63], 'MarkerSize',5 ) ;
end

for i = 2 : sz
   
    pause(0.01) ;
    
    for j = 1 : p.N
        handll(j).XData = [ r{j}(i,1) r{j+1}(i,1)] ;
        handll(j).YData = [ r{j}(i,2) r{j+1}(i,2)] ;
        handll(j).ZData = [ r{j}(i,3) r{j+1}(i,3)] ;
    end
    
    drawnow ;
    
end