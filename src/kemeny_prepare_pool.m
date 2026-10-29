function [pool, worker_count] = kemeny_prepare_pool(k, workers)
%KEMENY_PREPARE_POOL Prepare a local process pool outside method timings.
%   WORKERS=[] selects min(k, MATLAB's hardware and configured pool limits).
%   A positive integer WORKERS overrides the automatic count exactly.
%   Saved profiles and operating-system CPU affinity are unchanged.

if nargin < 2
    workers = [];
end
validateattributes(k, {'numeric'}, {'real', 'scalar', 'finite', ...
    'integer', 'positive'}, mfilename, 'k');
if ~isempty(workers)
    validateattributes(workers, {'numeric'}, {'real', 'scalar', 'finite', ...
        'integer', 'positive'}, mfilename, 'workers');
end
if ~license('test', 'Distrib_Computing_Toolbox')
    error('kemeny:ParallelToolboxRequired', ...
        'The five-method comparison requires Parallel Computing Toolbox.');
end

if verLessThan('matlab', '9.13') % R2021b/R2022a use the local profile.
    cluster = parcluster('local');
else
    cluster = parcluster('Processes');
end
job_storage = getenv('KEMENY_JOB_STORAGE');
if ~isempty(job_storage)
    if ~isfolder(job_storage), mkdir(job_storage); end
    cluster.JobStorageLocation = job_storage;
end
cluster.NumThreads = 1;
if isempty(workers)
    % Respect both hardware defaults and saved profile limits.
    hardware_defaults = parallel.cluster.Local();
    default_workers = min(cluster.NumWorkers, hardware_defaults.NumWorkers);
    if isprop(cluster, 'PreferredPoolNumWorkers')
        default_workers = min(default_workers, cluster.PreferredPoolNumWorkers);
    end
    worker_count = min(double(k), default_workers);
else
    worker_count = double(workers);
    cluster.NumWorkers = max(cluster.NumWorkers, worker_count);
end

pool = gcp('nocreate');
if ~isempty(pool)
    same_configuration = isa(pool, 'parallel.ProcessPool') && ...
        isa(pool.Cluster, 'parallel.cluster.Local') && ...
        pool.NumWorkers == worker_count && pool.Cluster.NumThreads == 1;
    if ~same_configuration
        if ~isempty(pool.FevalQueue.RunningFutures) || ...
                ~isempty(pool.FevalQueue.QueuedFutures)
            error('kemeny:PoolBusy', ...
                ['The existing parallel pool has running or queued futures. ' ...
                 'Finish those tasks before changing its worker count or type.']);
        end
        evalc('delete(pool);');
        pool = [];
    end
end
if isempty(pool)
    evalc('pool = parpool(cluster, worker_count);');
end

% Warm the pool before algorithm timing.
warmup = zeros(worker_count, 1);
parfor q = 1:worker_count
    warmup(q) = sqrt(q); %#ok<PFBNS>
end
end
