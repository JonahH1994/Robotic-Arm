function torques = invDynamics(the, dthe, ddthe, p)
   
    Fg = zeros(3, p.N); Fg(3, :) = -p.m*p.g;
    F_hand = [0; 0; 0];
    F  = [Fg F_hand];
    
    M_hand = [0; 0; 0];
    
    [r, rG, ~, ~, ~, aG, R, w, alphas, ns] = kinematics(the, dthe, ddthe, p);
    
    torques = zeros(p.N, 1);
    for j = 1:p.N
        Mg = [0; 0; 0];
        Hdot = 0;
        for q = j:p.N
            rGqj = rG(:,q) - rG(:,j);
            mq        = p.m(q);
            aGq     = aG(:,q);
            Rq        = reshape(R(:,q), [3,3]);
            alphaq    = alphas(:,q);
            wq      = w(:,q);
            
            Iq        = Rq' * p.I(:,3*q-2:3*q) * Rq;
            
            % angular momentum
            Hdot = Hdot + cross(rGqj, mq*aGq) + Iq * alphaq + ...
                cross(wq, Iq * wq);
            
            % gravity terms
            Mg = Mg + cross(rGqj, F(:,q));
        end
        
        % hand terms
        rhandj = r(:,end) - r(:,j);
        M_hand = M_hand + cross(rhandj, F_hand);
        
        torques(j) = dot(Hdot - Mg - M_hand, ns(:,j));
    end

end