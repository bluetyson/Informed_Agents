function [neighbor, nNeighbor]= convert_graph(A)
% Convert an adjacency matrix into a concatenated vector listing all the
% neighbors and another vector containing the number of neighbors per node.
%
% INPUT
% A: [n n]
%  adjacency matrix
%
% OUTPUT
% neighbor: [m 1]: concatenated list of neighbors, m=nnz(A)
% nNeighbor: [n 1]: node i has neighbor(i) neighbors
n = size(A, 1);
nNeighbor = full(sum(A, 2));
[neighbor, ~] = ind2sub(n, full(find(A)));
end % convert_graph(...)
