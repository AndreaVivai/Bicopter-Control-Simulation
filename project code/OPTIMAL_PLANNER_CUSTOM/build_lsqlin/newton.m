function [Deltax, lambda] = newton(G, H, A)

    M = [H,  A'; 
         A,  zeros(size(A, 1))];

    b = [-G; 
          zeros(size(A, 1), 1)];

    solution = M \ b;

    Deltax = solution(1:length(G));
    lambda = solution(length(G) + 1:end);

end

