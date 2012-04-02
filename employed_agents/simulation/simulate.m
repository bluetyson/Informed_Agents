function [opinion, consensusTime, t, intermediateOpinion, employedAgent] = ...
 simulate(A, m, maxT, dt, n0, mu0, u0, consensusFraction, epsilon, ...
 employedAgent, getIntermediateOpinion, breakOnConsensus, action)
% Perform a simulation of the "employed" agent model on a given graph.
% n: [1]: number of agents
%
% INPUT
% A: [n n]: adjacency representation of a graph describing the connections
%  between agents, i.e. A(i, j) is 1 iff there is an edge between node i and j.
%  A is assumed to be symmetric.
% m: [1]: number of simulations done in parallel
% maxT: [1]: maximal number of time steps simulated
% dt: [1]: number of time steps done in a single call to vectorized_simulate 
% n0: [1]: number of employed agents, 0 <= n0 <= n.
% mu0, u0: [1] / [m 1]: model constants. Either a single constant is provided
%  and used for every simulation or a vector containing (different) constants
%  providing the u and mu-values for the different simulations
% consensusFraction: [1]: fraction of agents with positive opinion needed for
%  consensus to be reached
% employedAgent: empty / [n0 m]: if employedAgent is empty, for every
%  simulation the employed agents are selected uniformly at random. Otherwise
%  the i-th row should be a subset of size n0 of the numbers 1 to n
%  corresponding to the employed agent in simulation i
% epsilon: [1]: if average absolute change over dt simulation steps is smaller
%  than epsilion, i.e. deltaOpinion / dt < epsilon, in all m simulation then
%  stop the simulation. Use epsilon = -1 for no early break.
% getIntermediateOpinion: [1]: flag to get the intermediate opinions of all m
%  simulations after every dt time steps
% breakOnConsensus: [1]: flag to break when consensus is reached in every
%  simulation
% action: string: eval(action) is executed after each simulation step. This
%  parameter can be used to individually print or animate things without
%  changing the whole function. Note that eval(action) possibly has side-effects
%  that might change the function's behavior.
%
% OUTPUT
% opinion: [m n]: opinion(i, :) gives the final opinion of the agents in
%  simulation i
% consensusTime: [m 1]: consensusTime(i) gives the time consensus is reached in
%  simulation i. If no consensus is reached, consensusTime(i) = Inf.
% t: [1]: number of time steps simulated
% intermediateOpinion: {1 nIntermediateSteps}: if getIntermediateOpinion is set,
%  intermediateOpinion{i} holds a matrix of size [m n] representing the opinion
%  of the agents before the i-th iteration of the while-loop
% employedAgent: [n0 m]: employedAgent(:, i) gives the agents selected as
%  employed agent in simulation i

% get data structure for generation of meetings
[neighbor, nNeighbor] = convert_graph(A);
n = size(A, 1);

% initial opinion
opinion = 2*rand(m, n)-1; % U(-1, 1)

% simulation parameters
if size(mu0, 1) == 1, mu = repmat(mu0, m, 1); else mu = mu0; end
if size(u0, 1) == 1, u = repmat(u0, m, 1); else u = u0; end

% agent parameters
selfOpinionFactor = ones(m, n);
otherOpinionFactor = zeros(m, n);
opinionBias = zeros(m, n);
if isempty(employedAgent)
    employedAgent = zeros(n0, m);
    for i=1:m
        employedAgent(:, i) = randsample(n, n0);
    end
end
nPositiveAgent = zeros(m, 1);
normalAgentFlag = ones(m, n);

% select special agents and modify their parameter
for i = 1:m
    opinion(i, employedAgent(:, i)) = ones(1, n0);
    normalAgentFlag(i, employedAgent(:, i)) = zeros(1, n0);
    nPositiveAgent(i) = sum(opinion(i, :)>0);
    opinionBias(i, employedAgent(:, i)) = u(i);
    selfOpinionFactor(i, employedAgent(:, i)) = zeros(1, n0);
    otherOpinionFactor(i, employedAgent(:, i)) = ones(1, n0);
end

t = 0;
consensusTime = inf(m, 1);
if getIntermediateOpinion ~= 0
    intermediateOpinion = cell(1, ceil(maxT/dt));
else
    intermediateOpinion = {};
end
whileStep = 0;

while t < maxT
    whileStep = whileStep + 1;
    if getIntermediateOpinion ~= 0
        intermediateOpinion{whileStep} = opinion;
    end
    % generate m * dt pairs of meeting agents
    [meetingA, meetingB] = generate_meetings(neighbor, nNeighbor, m, dt);
    % perform simulation
    [opinion, stepConsensusTime, nPositiveAgent, deltaOpinion] = ...
     vectorized_simulate(meetingA, meetingB, opinion, mu, u, ...
     selfOpinionFactor, otherOpinionFactor, opinionBias, normalAgentFlag, ...
     nPositiveAgent, consensusFraction);
    t = t + dt;
    noConsensusSimulation = find(stepConsensusTime==dt+1);
    stepConsensusTime(noConsensusSimulation) = inf;
    consensusTime = min(consensusTime, stepConsensusTime+t-dt);
    eval(action);
    if (breakOnConsensus && isempty(noConsensusSimulation)) ...
         || max(deltaOpinion)/dt<epsilon
        break;
    end
end

if getIntermediateOpinion ~= 0
    intermediateOpinion = intermediateOpinion(1:whileStep);
end

end % simulate(...)
