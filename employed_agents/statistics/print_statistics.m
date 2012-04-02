function print_statistics(A)
% Print statics for a undirected, loop-free graph.
%
% INPUT
% A: [n n]: adjacency matrix

n = size(A, 1);
fprintf('number of nodes = %d\n', n);

m = nnz(A);
assert(mod(m, 2) == 0);
fprintf('number of (undirected) edges = %d\n', m/2);

k = m/n;
fprintf('average node degree = %.4f\n', k);

maxDeg = max(full(sum(A)));
fprintf('maximum node degree = %d\n', maxDeg);

avgPathLength = average_path_length(A);
fprintf('average path length = %.4f\n', avgPathLength);

clusterCoeff = global_clustering_coefficient(A);
fprintf('global clustering coefficient = %.4f\n', clusterCoeff);

end % print_statistics(...)
