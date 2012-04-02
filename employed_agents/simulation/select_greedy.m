function employedAgents = select_greedy(A, n0)
% Select greedily (as described in the report) the n0 agents corresponding to
% the nodes with the largest degrees.
%
% INPUT
% A: [n n]: adjacency representation of a graph describing the connections
%  between agents, i.e. A(i, j) is 1 iff there is an edge between node i and j.
%  A is assumed to be symmetric.
% n0: [1]: number of employed agents to be selected, n0 <= n
%
% OUTPUT
% employedAgents: [1 n0]: indices of the selected employed agents

[~, J] = sort(sum(A));
employedAgents = J(end-n0+1:end);

end % select_greedy(...)
