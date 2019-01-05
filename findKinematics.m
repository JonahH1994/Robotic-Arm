%% Kinematic Calculator:
% Given theta, dtheta, ddtheta

function [r,v,a,w,alpha,R,n] = findKinematics(the,dthe,ddthe,p)

r = zeros(3,p.N) ;
v = zeros(3,p.N) ;
a = zeros(3,p.N) ;
w = zeros(3,p.N) ;
R = zeros(3,3*p.N) ;
alpha = zeros(3,p.N) ;
n = zeros(3,p.N) ;

rG = r ;
vG = v ;
aG = a ;

for i = 1 : p.N
   
   	if i == 1
    
        R(:,3*i-2:3*i) = rotmat(p.nor(:,i),the(i)) ; 
        n(:,i) = R(:,3*i-2:3*i) * p.nor(:,i) ;
        r(:,i) = R(:,3*i-2:3*i)*p.k*p.l(i) ;
        rG(:,i) = R(:,3*i-2:3*i)*p.k*p.G(i) ;
        w(:,i) = R(:,3*i-2:3*i)*p.nor(:,i)*dthe(i) ;
        %v(:,i) = cross(w(:,i),r(:,i)) ;
        alpha(:,i) = ddthe(i)*R(:,3*i-2:3*i)*p.nor(:,i) ;
        %a(:,i) = cross(alpha(:,i),r(:,i)) ;
        %a_krelj_F = cross(omega_j, v_krelj_F) + cross(alpha_j, r_krelj_F);
        %as_F(:,k) = as_F(:,j) + a_krelj_F;
    
    else
        
        nj = R(:,3*(i-1)-2:3*(i-1)) * p.nor(:,i) ;
        Rji = rotmat(nj,the(i)) ;
        R(:,3*i-2:3*i) = Rji * R(:,3*(i-1)-2:3*(i-1)) ;
        
        n(:,i) = nj ; %R(:,3*i-2:3*i) * p.nor(:,i) ;
        
        r_j_i = R(:,3*i-2:3*i) * p.l(i) * p.i ;
        rG_j_i = R(:,3*i-2:3*i)* p.G(i) * p.i ;
        r(:,i) = r(:,i-1) + r_j_i ;
        rG(:,i) = r(:,i-1) + rG_j_i ; 
        
        w_j_i = dthe(i) * nj ;
        w(:,i) = w(:,i-1) + w_j_i ;
        
        v_j_i =  cross(w(:,i),r_j_i) ;
        vG_j_i = cross(w(:,i),rG_j_i) ;
        v(:,i) = v(:,i-1) + v_j_i ;
        vG(:,i) = v(:,i-1) + vG_j_i ;
        
        alp_j_i = ddthe(i) * nj ;
        alpha(:,i) = alpha(:,i-1) + alp_j_i + cross(alpha(:,i-1),alp_j_i) ;
        
        a_j_i = cross(w(:,i), v_j_i) + cross(alpha(:,i), r_j_i) ;
        aG_j_i =cross(w(:,i), vG_j_i) + cross(alpha(:,i), rG_j_i) ;
        a(:,i) = a(:,i-1) + a_j_i ;
        aG(:,i) = a(:,i-1) + aG_j_i ;
        
    end
    
end

end

function R = rotmat(n, theta)
    n = n/norm(n);
    R = n*n' + cos(theta)*(eye(3)-n*n') + sin(theta)*skew(n);
end

function S = skew(v)

S = [  0   -v(3)  v(2) ; ...
      v(3)  0    -v(1) ; ...
     -v(2)  v(1)    0  ];

end