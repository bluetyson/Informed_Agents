function lG = average_path_length(A)
% Calculate the average path length in a graph.
% (see http://en.wikipedia.org/wiki/Average_path_length for a definition)
%
% INPUT
% A: [n n]: adjacency matrix
%
% OUTPUT
% lG: [1]: average path length
%
% This function makes use of the Parallel Computing Toolbox, i.e. the outer loop
% can be divided among several workers.

n = size(A, 1);
lGs = zeros(n, 1);
parfor i = 1:n
    active = i;
    todo = ones(n,1);
    todo(i) = 0;
    dist = 0;
    total = 0;
    while ~isempty(active)
        dist = dist + 1;
        active = find(sum(A(:,active), 2) .* todo);
        todo(active) = 0;
        total = total + size(active, 1) * dist;
    end
    lGs(i) = total;
end
lG = sum(lGs) / (n * (n - 1));

end % average_path_length(...)
