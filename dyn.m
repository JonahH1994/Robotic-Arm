function dz = dyn(t,z,p)
    
    % Define derivative state:
    dw = zeros(3*p.N,1);
    dR = zeros(3,3*p.N) ;
    w = z(3*p.N+1:3*(2*p.N)) ;
    
    for i = 1 : p.N
       
        if sum(sum(p.I((3*i-2):3*i,:) == zeros(3,3))) ~= 9
            %w_i = z(p.N*i+1) ;
            %I = p.I(3*i-2:3*i,:) ;
            %dz(i) = -inv(I)*cross(w_i,I*w_i) ;
            if i == 1
                dw(i) = 0 ;
            elseif ( mod(i,2) == 1 )
                dw(i) = exp(-0.5*i*t) ;
            else
                dw(i) = cos(t*i + i*pi/10 ) ;
            end
            dR(:,3*i-2:3*i) = skeww(w(3*i-2:3*i)) * reshape(z(3*(2*p.N) + ...
                1+9*(i-1): 3*(2*p.N) + 9*i),[3,3]) ;
        else
            dR(:,3*i-2:3*i) = skeww(w(3*i-2:3*i)) * reshape(z(3*(2*p.N) + ...
                1+9*(i-1): 3*(2*p.N) + 9*i),[3,3]) ;
        end
        
    end
    
    if (p.N == 7)
        dw = zeros(3*p.N,1);
        dw(1:3:3*p.N) = [ 0; 3*exp(t)/(1+exp(t))^2; cos(3*t); cos(1.5*t); sin(3*t); ...
            sin(1.5*t); exp(0.05*t) ] ;
    end
    
    dz = [ w; dw; dR(:)] ;

end

function S = skeww(w)

    S = [ 0, -w(3), w(2); w(3), 0, -w(1); -w(2), w(1), 0 ] ;

end