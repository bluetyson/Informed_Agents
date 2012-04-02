function B = largest_connected_component(A)
% Returns the adjacency matrix of the largest connected component of an
% undirected graph.
%
% INPUT
% A: [n,n]: symmetric adjacency matrix with entries 0 resp. 1
%
% OUTPUT
% B: [m,m]: Adjacency matrix of the largest connected component of the graph
%  represented by the matrix A.

n = size(A, 1);
% vector to record visited nodes: node i was visited iff visited(i)==1
visited = zeros(n, 1);
todo = ones(n,1);
indexLargestComponent = 0;
sizeLargestComponent = 0;
for i = 1:n
    % if node i wasn't visited yet, start a new component search at node i
    if ~visited(i)
        % vector with indices of currently active nodes
        active = i;
        todo(i) = 0;
        count = 1;
        visited(i) = i;
        while ~isempty(active)
            % Set active nodes to all not visited nodes reachable from ...
            % currently active nodes.
            active = find(sum(A(:, active), 2) .* todo);
            count = count + length(active);
            todo(active) = 0;
            visited(active) = i;
        end
        if count>sizeLargestComponent
            sizeLargestComponent = count;
            indexLargestComponent = i;
        end
    end
end

% extract indices of nodes belonging to the largest connected component
largestComponent = find(visited==indexLargestComponent);
B = A(largestComponent, largestComponent);

end % largest_connected_component(...)
