function [dz,res] = dyn(t,z,p)
    
    res=0;
    % Define derivative state:
    dw = zeros(p.N,1);
    dR = zeros(3,p.N*3) ;
    w = z(p.N+1:(2*p.N)) ;
    
    if  (p.jacob == 1)
        
        J = zeros(6,p.N) ;
        ddthe = zeros(p.N,1) ;
        
%         [r,~,~,~,~,~,n] = findKinematics(z(1:p.N),z(p.N+1:2*p.N),zeros(p.N),p) ;
%         for i = 1 : p.N
%            
%             if i == p.N
%                 J(1:3,1) = cross(n(:,i),r(:,i)) ;
%                 J(4:end,1) = n(:,i) ;
%             else
%                 J(1:3,i+1) = cross(n(:,i),r(:,end)-r(:,i)) ;
%                 J(4:end,i+1) = n(:,i) ;
%             end
%             
%         end
        
        for i = 1 : p.N
           
            dthe = zeros(p.N,1) ;
            dthe(i) = 1 ;
            
            %r = calculateState([z(1:p.N);dthe;z(2*p.N+1:end)]',p) ;
            [~,v,~,w,~,~,~] = findKinematics(z(1:p.N),dthe,zeros(p.N),p) ;
            %J(:,i) = [r{end}(1,4:6)' ; r{end}(1,7:9)' ] ;
            J(:,i) = [v(:,end);w(:,end)] ;
            
        end
        
        Jp = pinv(J) ; %((J'*J)\J') ;
%         if (p.N == 6)
%             rc = rcond(Jp) ;
%         else
%             rc = 1 ;
%         end
%         if rc <= 10e-22
%            res = 1 ; 
%         end
        [r,~,~,~,~,~,~] = findKinematics(z(1:p.N),z(p.N+1:2*p.N),zeros(p.N),p) ;
        %vDes = 0.2*(p.Xtarg-r{end}(1,1:3)')/norm(p.Xtarg-r{end}(1,1:3)') ;
        %vDes = 0.75*(p.Xtarg-r(:,end))/norm(p.Xtarg-r(:,end)) ;
        vDes = 0.5*(p.Xtarg-r(:,end)) ;
        dq = Jp * ([vDes;zeros(3,1)]) ;
        %ddthe = zeros(p.N,1) ; %0.1*(dq-z(p.N+1:2*p.N)) ;
        %dq = pinv(J(1:3,:)) * vDes ;
        %dz = [dq;zeros(size(dq));dR(:)];
       
%         R = zeros(3,p.N*3) ;
%         O = zeros(3,p.N) ;
%         J = zeros(6,p.N) ;
%         
%         for i = 1 : p.N
%            
%             if i == 1
%                 
%                 R(:,(3*i-2):3*i) = rotMat(z(1),3) ;
%                 O(:,i) = R(:,(3*i-2):3*i)*p.l(i)*p.k ;
%                 
%             else
%                 
%                 R(:,3*i-2:3*i) = R(:,(3*(i-1)-2):3*(i-1)) * rotMat(z(i),2) ;
%                 O(:,i) = R(:,(3*i-2):3*i)*p.l(i) * p.i + O(:,i-1) ;
%                 
%             end
%             
%         end
%         
%         for i= 1 : p.N
%            
%             if i == 1
%             
%                 zi = R(:,end-2:end)*p.k ; 
%                 J(:,1) = [cross(zi,O(:,end));zi] ;
%                 
%             else
%                 
%                 zi = R(:,(3*(i-1)-2):3*(i-1)) * p.k ;
%                 J(:,i) = [cross(zi,(O(:,end)-O(:,i-1)));zi] ;
%                 
%             end
%             
%         end
%         
%         Jp = ((J'*J)\J') ;
%         dq = Jp * ([p.Xtarg;p.theTarg]-[O(:,end);z(p.N);0;z(1)]) ;
%         disp(norm(p.Xtarg-O(:,end))) ;
%         %disp(O(:,end))
        
        %dz = [z(p.N+1:2*p.N);dq;dR(:)] ;
        dz = [dq;ddthe;dR(:)] ;
        %dz = [-dq+z(p.N+1:2*p.N);zeros(size(dq));dR(:)] ;
        %dz = [dq;-dq+z(p.N+1:2*p.N);dR(:)] ;
        
    else

        for i = 1 : p.N

            dw(i) = -p.g/p.l(i) * sin(z(i)) ;

    %         if sum(sum(p.I((3*i-2):3*i,:) == zeros(3,3))) ~= 9
    %             %w_i = z(p.N*i+1) ;
    %             %I = p.I(3*i-2:3*i,:) ;
    %             %dz(i) = -inv(I)*cross(w_i,I*w_i) ;
    %             if i == 1
    %                 dw(i) = 0 ; %sin(t*i + i*pi/10) ; %0 ;
    %             elseif ( mod(i,2) == 1 )
    %                 dw(i) = exp(-0.5*i*t) ;
    %             else
    %                 dw(i) = cos(t*i + i*pi/10 ) ;
    %             end
    %             %dR(:,3*i-2:3*i) = skeww(w(3*i-2:3*i)) * reshape(z(3*(2*p.N) + ...
    %             %    1+9*(i-1): 3*(2*p.N) + 9*i),[3,3]) ;
    %             dR(:,3*i-2:3*i) = zeros(3,3) ;
    %         else
    %             %dR(:,3*i-2:3*i) = skeww(w(3*i-2:3*i)) * reshape(z(3*(2*p.N) + ...
    %             %    1+9*(i-1): 3*(2*p.N) + 9*i),[3,3]) ;
    %             dR(:,3*i-2:3*i) = zeros(3,3) ;
    %         end

        end

    %     if (p.N == 7)
    %         %dw = zeros(p.N,1);
    %         dw = [ 0; 3*exp(t)/(1+exp(t))^2; cos(3*t); cos(1.5*t); sin(3*t); ...
    %             sin(1.5*t); exp(0.05*t) ] ;
    %     end

        dz = [ w; dw; dR(:)] ;
    end

end

function S = skeww(w)

    S = [ 0, -w(3), w(2); w(3), 0, -w(1); -w(2), w(1), 0 ] ;

end

function R = rotMat(ang,ax)

if ax == 1 
    R = eye(3) ;
elseif ax == 2 
    R = [ cos(ang), 0, sin(ang) ; ...
                    0, 1,           0 ; ...
         -sin(ang), 0, cos(ang) ] ;
else
    R = [ cos(ang),  sin(ang), 0 ; ...
         -sin(ang),  cos(ang), 0 ; ...
                 0,         0, 1 ] ;
end

end