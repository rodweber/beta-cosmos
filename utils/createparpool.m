function [activeParpool] = createparpool(maxNumWorkers)
%CREATEPARPOOL Summary of this function goes here
% Detailed explanation goes here.

% Check for existing parallel pool.
currentParpool = gcp('nocreate'); % If no pool, do not create new one.
if isempty(currentParpool)
  poolsize = 0;
else
  poolsize = currentParpool.NumWorkers;
end

% Return pool with the required number of workers.
if poolsize == maxNumWorkers
  % Current pool already has the required number of workers.
  activeParpool = currentParpool;
else
  delete(currentParpool); % Delete current pool, if it exists.
  % Create a cluster object according to the specified profile.
  activeCluster = parcluster('local');
  activeCluster.NumWorkers = maxNumWorkers;
  activeParpool = parpool(activeCluster);
end

end % Function createparpool()
