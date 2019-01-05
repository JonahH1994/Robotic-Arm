function R = rotMat(n, theta)
    n = n/norm(n);
    R = n*n' + cos(theta)*(eye(3)-n*n') + sin(theta)*skew(n);
end