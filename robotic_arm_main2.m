%% Robotic Arm:

% Number of links:
p.N = 6 ;
p.Xtarg = rand(3,1) ; %[3.678;2.148;3.542] ; %[2;2;4] ;
p.Xtarg = p.Xtarg / norm(p.Xtarg)*p.N/2 ;
p.theTarg = zeros(3,1) ;

% Use jacobian to find target:
p.jacob = 1 ;

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
p.g = 9.8 ;
%p.m(1) = 0 ;

p.len = 1 ; % length of each link in meters
p.Xtarg = p.Xtarg * p.len ;
p.l = p.len * ones(p.N,1) ;
p.l(1) = 4 ;

l = 2 ; %max(p.l) ;

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
th0 = zeros(p.N,1) ;
th0(:) = (0:p.N-1)' * pi/20 ; % + pi/2;
w0 = zeros(p.N,1) ;
p.n = zeros(3,p.N) ;
%w0(1) = 3 ;

% Create initial rotation matrices:
R0 = zeros(3,3*p.N) ;
for  i = 1 : p.N
    
    if i == 1 
        p.nor(:,i) = p.k ;
    else
        p.nor(:,i) = p.j ;
    end
   
%     if i == 1
%         
%         p.nor(:,i) = p.k ;
%         Rr = [ cos(th0(i)),  sin(th0(i)), 0 ; ...
%               -sin(th0(i)),  cos(th0(i)), 0 ; ...
%                          0,            0, 1 ] ;
%     else
%         p.nor(:,i) = p.j ;
%         Rr = [ cos(th0(i)), 0, sin(th0(i)) ; ...
%                          0, 1,           0 ; ...
%               -sin(th0(i)), 0, cos(th0(i)) ] ;
%         
%     end

    Rr = rotMat(p.nor(:,i),th0(i)) ;
    
    R0(:,3*i-2:3*i) = Rr ;
    
%     % rotation matrix about z:
%     Rr0z = [ cos(th0(3*i-2)), -sin(th0(3*i-2)), 0 ; ...
%              sin(th0(3*i-2)),  cos(th0(3*i-2)), 0 ;...
%                            0,                0, 1 ] ;
%     
%     % Rotation matrix about y:
%     Rr0y = [ cos(th0(3*i-1)), 0,-sin(th0(3*i-1)) ; ...
%                            0, 1,               0 ; ...
%              sin(th0(3*i-1)), 0, cos(th0(3*i-1)) ] ;
%         
%     % Rotation matrix about x:
%     Rr0x = [ 1,             0,              0 ; ...
%              0, cos(th0(3*i)),  sin(th0(3*i)) ; ...
%              0,-sin(th0(3*i)),  cos(th0(3*i)) ] ;
%                     
%     R0(:,3*i-2:3*i) = Rr0z * Rr0y * Rr0x ;
    
end

% Setup ode
dt = 1/60 ;
tf = 15 ;
tSpan = linspace(0,tf,tf/dt) ; %[0 10] ;
opts = odeset( 'RelTol', 1e-10, 'AbsTol', 1e-10 ) ;
z0 = [ th0; w0; R0(:) ] ;
%disp(z0(1:3)) ;

%p.Xtarg = z0(1:3) + [0;cos(pi/4);cos(pi/4)]* 1 ;
%p.Xtarg = [1.34;-1.2;2.05] ; %[2;-2;1] ;

% integrate:
% [t,xState] = ode45( @dyn, tSpan, z0, opts, p ) ;

xState = zeros(length(tSpan),length(z0)) ;
xState(1,:) = z0' ;

nd = length(tSpan) ;
p.tol = 0.05 ;

disp(nd)
res = 0 ;
rann = 0 ;
while rann == 0
    
    for i = 2 : length(tSpan) 

        if (res == 1)
            
           %res = 0 ;
           th0 = rand(p.N,1) ;
           th0 = th0/max(th0)*pi ;
           z0 = [ th0; w0; R0(:) ] ;
           xState = zeros(length(tSpan),length(z0)) ;
           xState(1,:) = z0' ;
           res = 0 ;
           break ;

        end
        [dz, res] = dyn2(tSpan,xState(i-1,:)',p) ;
        %dz = dyn1(tSpan,xState(i-1,:)',p) ;
        xState(i,:) = dz'*(tSpan(i)-tSpan(i-1)) + xState(i-1,:) ;

        [r,~,~,~,~,~,~] = findKinematics(xState(i,1:p.N),xState(i,p.N+1:2*p.N),zeros(p.N),p) ;
        if abs(r(:,end)-p.Xtarg) < p.tol

            nd = i ;
            res = 2 ;
            break ;

        elseif i >= length(tSpan)-1
            res = 2;
        end

    end
    if (res == 2)
        rann = 1;
    end
end
disp(nd) ;
xState = xState(1:nd,:) ;
tSpan = tSpan(1:nd) ;
t = tSpan ;

% Simulate each link:
r = calculateState(xState,p) ;
[sz,~] = size(xState) ;

% figure
% plot(t,r{6}(:,1)) ;
% hold on ;
% plot(t,ones(size(r{6}(:,1)))*p.Xtarg(1)) ;
% hold off ;

% figure
% plot3( r{1}(:,1), r{1}(:,2), r{1}(:,3) ) ;
% 
% figure
% plot3( r{2}(:,1), r{2}(:,2), r{2}(:,3) ) ;

figure 
hold on;
grid on ;

view(3) ;
a = (p.N-1)*p.l(2) ;
h = sum(p.l) ;
axis([-a a -a a 0 h]) ;
%axis([-l l -l l -l l]*(p.N+0.25) ) ;
xlabel( 'X Position (m)' ) ;
ylabel( 'Y Position (m)' ) ;
zlabel( 'Z Positoin (m)' ) ;

handll = gobjects(1,p.N) ;

for i = 1 : p.N
    handll(i) = plot3( [ r{i}(1,1) r{i+1}(1,1)], [r{i}(1,2) r{i+1}(1,2)], ...
        [r{i}(1,3) r{i+1}(1,3)], '-o', 'LineWidth', 4, ...
        'MarkerFaceColor',[.49 1 .63], 'MarkerSize',5 ) ;
end
plot3( r{end}(1,1), r{end}(1,2), r{end}(1,3), '.k' ) ;
plot3( p.Xtarg(1), p.Xtarg(2), p.Xtarg(3), '*r', 'MarkerSize', 10 ) ;

for i = 2 : sz
   
    pause(0.05) ;
    
    for j = 1 : p.N
        handll(j).XData = [ r{j}(i,1) r{j+1}(i,1)] ;
        handll(j).YData = [ r{j}(i,2) r{j+1}(i,2)] ;
        handll(j).ZData = [ r{j}(i,3) r{j+1}(i,3)] ;
    end
    
    plot3( r{end}(i,1), r{end}(i,2), r{end}(i,3), '.k' ) ;
    drawnow ;
    
end
