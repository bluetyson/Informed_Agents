function [opinion, consensusTime, t, intermediateOpinion, employedAgent, ...
 filename] = simulate_and_store(A, m, maxT, dt, n0, mu0, u0, ...
 consensusFraction, epsilon, employedAgent, getIntermediateOpinion, ...
 breakOnConsensus, action, graphName)
% Wrapper for simulate(...). Perform a simulation of the "employed" agent model
% on a given graph and save input and output as well as the random stream and
% its state for later processing / reproducing the received results.
% The file is written to the directory 'simulationResult', for proper operation,
% this folder has to be within the Matlab path.
%
% INPUT
% Same as input to simulate(...). One additional argument has to be provided:
% graphName: string: name of the graph represented by A
%
% OUTPUT
% Same as output of simulate(...). One additional variable is returned:
% filename: string: the name of the file that input and output have been saved
%  to. The name of the output file is a concatenation of graphName and the
%  parameters n0, mu0 and u0 being used. To make sure that things work properly
%  using the Parallel Computing Toolbox, the filename also includes the worker
%  id.

muStr = 'var';
if size(mu0, 1) == 1, muStr = num2str(mu0); end;
uStr = 'var';
if size(u0, 1) == 1, uStr = num2str(u0); end;

worker = '';
w = getCurrentTask();
if ~isempty(w)
    worker = ['_' num2str(get(w, 'ID'))];
end

% determine random stream and its state
stream = RandStream.getDefaultStream;
state = stream.State;

[opinion, consensusTime, t, intermediateOpinion, employedAgent] = ...
 simulate(A, m, maxT, dt, n0, mu0, u0, consensusFraction, epsilon, ...
 employedAgent, getIntermediateOpinion, breakOnConsensus, action);

for i=1:1024
    filename = sprintf('simulationResult/%s_n0=%d_mu=%s_u=%s_[%d%s].mat', ...
     graphName, n0, muStr, uStr, i, worker);
    if exist(filename, 'file') ~= 2
        save(filename, 'opinion', 'consensusTime', 't', ...
         'intermediateOpinion', 'employedAgent', 'A', 'm', 'maxT', 'dt', ...
         'mu0', 'u0', 'consensusFraction', 'epsilon', 'employedAgent', ...
         'getIntermediateOpinion', 'breakOnConsensus', 'action', ...
         'graphName', 'stream', 'state');
        return
    end
end

fprintf('Error writing parameters and output to file\n');

end % simulate_and_store(...)
