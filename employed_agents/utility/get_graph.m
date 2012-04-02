function [A, name, stream, state] = get_graph(graphID, graphSize, ...
 generateNew, verbose)
% Function for loading, generating and storing graphs. If the requested graph
% does not exist, the function generates a new one and stores it to a mat-file.
% For proper operation, the folders 'generatedGraphs' and 'facebook100' have to
% be within the Matlab path.
% 
% INPUT
% graphID: [1]: identifier for the graph type which should be loaded/generated:
%  1   <--> Random graph
%  2   <--> Scale free graph
%  3   <--> Small world graph
%  4-6 <--> Real world graph
% graphSize: [1]: determines size of graph:
%  0 <-> small graph, i.e. ~1000 nodes
%  1 <-> large graph, i.e. ~35000 nodes
% generateNew: [1]: flag to generate a new graph, i.e. the graph is newly
%  generated even if it has been stored before
% verbose: [1]: enable showing what is going on
%
% OUTPUT
% A: [n,n] (sparse): adjacency matrix of the loaded/generated graph. If the
%  graph is not connected, its largest component is returned.
% name: string: describing the returned graph
% stream: [1 RandStream]: the random stream used for generating the graph (as
%  returned by RandStream.getDefaultStream).
% state: the state of stream before the graph is generated

% Cell array for graph filenames
generatedGraphsFilenames = {{'Random_small', 'ScaleFree_small', ...
 'SmallWorld_small'}; {'Random_large', 'ScaleFree_large', 'SmallWorld_large'}};
% Cell array for real world graph filenames
realWorldGraphsFilenames = {{'Simmons81', 'Reed98', 'Caltech36'}; ...
 {'Penn94', 'UF21', 'Texas84'}};

% Parameters for graph generation. The number of nodes n, edges m and the
% average nodal degree k have been chosen such that they best match the
% corresponding characteristics of the real world graphs.
n = [1000 35000];
m = [40000 2900000];
k = [40 76];
k2 = ceil(k/2);

% Cell array with functions for graph generation.
generationFunctions = {{@(n) random_graph(n, m(1)/((n-1)*(n-1))), ...
 @(n) scale_free(n, k2(1)+1, k2(1)), @(n) small_world(n, k(1),0.25)} ...
 {@(n) random_graph(n, m(2)/((n-1)*(n-1))), @(n) scale_free(n, k2(2)+1, ...
 k2(2)), @(n) small_world(n, k(2), 0.25)}};

% determine random stream and its state
stream = RandStream.getDefaultStream;
state = stream.State;

name = '';
A = [];
if graphID <= 3
    % Determine filename with cell array and check if file exists.
    name = generatedGraphsFilenames{graphSize+1}{graphID};
    filename = [name '.mat'];
    if ~generateNew && exist(filename, 'file') == 2
        if verbose
            fprintf('Loading graph from: %s \n', filename);
        end
        load(filename, 'A', 'stream', 'state');
    else % New graph to be generated or file does not exist.  
        if graphSize>1
            error 'Invalid value for graphSize'; 
        end;
        if verbose
            fprintf(['New graph with %d nodes is generated and saved to ' ...
             '%s\n'], n(graphSize+1), filename);
        end
        A = generationFunctions{graphSize+1}{graphID}(n(graphSize+1));
        oldSize = size(A, 1);
        A = largest_connected_component(A);
        if verbose && size(A, 1) ~= oldSize
            fprintf(['dropped %d nodes from graph of size %d resulting in ' ...
             'a graph of size %d\n'], oldSize-size(A, 1), oldSize, size(A, 1));
        end
        save(['generatedGraphs/' filename], 'A', 'stream', 'state'); 
    end
elseif graphID <= 6 % real world graph
    name = realWorldGraphsFilenames{graphSize+1}{graphID-3};
    filename = [name '.mat'];
    if verbose
        fprintf('Loading real world graph: %s\n', filename);
    end
    load(filename, 'A');
    oldSize = size(A, 1);
    A = largest_connected_component(A);
    if verbose && size(A, 1) ~= oldSize
        fprintf(['dropped %d nodes from graph of size %d resulting in a ' ...
         'graph of size %d\n'], oldSize-size(A, 1), oldSize, size(A, 1));
    end
else
    error 'Wrong graphID';
end

end % get_graph(...)
