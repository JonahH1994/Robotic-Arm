function r = calculateState(xState,p)

[sz,~] = size(xState) ;
r = cell(1,p.N+1) ;

for i = 1 : p.N
    r{i} = zeros(sz,3) ;
end

R = zeros(3,3) ;

for j = 1 : sz 
    for i = 2 : p.N+1
        if i == 2
            R = reshape(xState(j,3*(2*p.N) + 1+9*(i-2): 3*(2*p.N) + 9*(i-1)),[3,3]) ;
            r{i}(j,:) = p.l(i-1)*R*p.i ;
        else
            R = R * reshape(xState(j,3*(2*p.N) + 1+9*(i-2): 3*(2*p.N) ...
                + 9*(i-1)),[3,3]) ;
            
            r{i}(j,:) = r{i-1}(j,:) + (p.l(i-1)* R * p.i)' ;
%             for l = flip(2 : i)
%                 if l ~= i
%                     r{i}(j,:) = r{i}(j,:) + r{l}(j,:) ;
%                 else
%                     r{i}(j,:) = p.l(i-1) * R * p.i ;
%                 end
%             end
        end
    end
end

end