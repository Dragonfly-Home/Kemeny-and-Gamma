function varargout = run_kemeny(varargin)
%RUN_KEMENY Independent repeats and multiple partition groups per method.
% Columns: dataset, Direct repeats, Bini-Rec repeats, Bini-DirRec repeats,
%          Du [k repeats] pairs, New-full-rhs [k repeats] pairs.
% [[4,1],[6,2]] is a flat row in MATLAB; adjacent numbers form one pair.
% Counts are MEASURED calls; each enabled method/k group warms once first.
% Count 0 skips only that group, including its warmup.
TEST_CASES = {
    'as-22july06', 1, 1, 1, [[4, 1], [6, 1], [8, 1], [10, 1], [12, 1], [14, 1]], [[4, 1], [6, 1], [8, 1], [10, 1], [12, 1], [14, 1]];
    'big', 1, 1, 1, [[4, 1], [6, 1], [8, 1], [10, 1], [12, 1], [14, 1]], [[4, 1], [6, 1], [8, 1], [10, 1], [12, 1], [14, 1]];
    'ca-CondMat', 1, 1, 1, [[4, 1], [6, 1], [8, 1], [10, 1], [12, 1], [14, 1]], [[4, 1], [6, 1], [8, 1], [10, 1], [12, 1], [14, 1]];
    'ca-HepPh', 1, 1, 1, [[4, 1], [6, 1], [8, 1], [10, 1], [12, 1], [14, 1]], [[4, 1], [6, 1], [8, 1], [10, 1], [12, 1], [14, 1]];
    'cond-mat-2003', 1, 1, 1, [[4, 1], [6, 1], [8, 1], [10, 1], [12, 1], [14, 1]], [[4, 1], [6, 1], [8, 1], [10, 1], [12, 1], [14, 1]];
    'cs4', 1, 1, 1, [[4, 1], [6, 1], [8, 1], [10, 1], [12, 1], [14, 1]], [[4, 1], [6, 1], [8, 1], [10, 1], [12, 1], [14, 1]];
    'de2010', 1, 1, 1, [[4, 1], [6, 1], [8, 1], [10, 1], [12, 1], [14, 1]], [[4, 1], [6, 1], [8, 1], [10, 1], [12, 1], [14, 1]];
    'gre_1107', 1, 1, 1, [[4, 1], [6, 1], [8, 1], [10, 1], [12, 1], [14, 1]], [[4, 1], [6, 1], [8, 1], [10, 1], [12, 1], [14, 1]];
    'hi2010', 1, 1, 1, [[4, 1], [6, 1], [8, 1], [10, 1], [12, 1], [14, 1]], [[4, 1], [6, 1], [8, 1], [10, 1], [12, 1], [14, 1]];
    'vsp_c-30_data_data', 1, 1, 1, [[4, 1], [6, 1], [8, 1], [10, 1], [12, 1], [14, 1]], [[4, 1], [6, 1], [8, 1], [10, 1], [12, 1], [14, 1]];
    'worms20_10NN', 1, 1, 1, [[4, 1], [6, 1], [8, 1], [10, 1], [12, 1], [14, 1]], [[4, 1], [6, 1], [8, 1], [10, 1], [12, 1], [14, 1]];
};
WORKERS = [];             % HPC: min(allocated CPUs-1, largest enabled k).
SC_BATCH_SIZE = Inf;       % New stochastic-complement RHS batch width.
FEEDBACK_BATCH_SIZE = 64;  % New feedback RHS batch width.

ip = inputParser;
addParameter(ip, 'Tests', TEST_CASES);
addParameter(ip, 'Workers', WORKERS, @(x) isempty(x) || integer_scalar(x, 1));
addParameter(ip, 'SCBatchSize', SC_BATCH_SIZE, @batch_width);
addParameter(ip, 'BatchSize', FEEDBACK_BATCH_SIZE, @(x) integer_scalar(x, 1));
addParameter(ip, 'OutputDir', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(ip, 'PlanOnly', false, @(x) islogical(x) && isscalar(x));
addParameter(ip, 'MemoryLimitGB', Inf, @batch_width);
parse(ip, varargin{:});
opt = ip.Results;
tests = normalize_tests(opt.Tests);
plan = make_plan(tests);
configuration = struct('Tests', {tests}, 'Plan', plan, ...
    'RequestedWorkers', opt.Workers, 'SCBatchSize', opt.SCBatchSize, ...
    'BatchSize', opt.BatchSize, 'WarmupsPerEnabledGroup', 1, ...
    'Reference', 'mean of successful measured Direct results in the same TEST_CASES row; warmup excluded', ...
    'Timer', 'caller tic/toc around one algorithm call; warmup excluded');
if opt.PlanOnly
    disp(plan);
    if nargout > 0, varargout{1} = plan; end
    if nargout > 1, varargout{2} = configuration; end
    return
end

root = fileparts(mfilename('fullpath'));
old_path = path;
path_cleanup = onCleanup(@() path(old_path)); %#ok<NASGU>
addpath(fullfile(root, 'third_party', 'bini'), '-begin');
addpath(fullfile(root, 'third_party', 'du'), '-begin');
addpath(fullfile(root, 'src'), '-begin');
assert(~isempty(which('dissect')), 'KemenyRunner:MissingDissect', 'MATLAB dissect is required.');
[workers, allocated] = choose_workers(plan, opt.Workers);
memory_gib = memory_limit(opt.MemoryLimitGB, allocated);
configuration.Workers = workers;
configuration.AllocatedCPUs = allocated;
configuration.MemoryLimitGiB = memory_gib;

output = char(opt.OutputDir);
if isempty(output)
    parent = fullfile(root, 'runs');
    if ~isfolder(parent), mkdir(parent); end
    output = tempname(parent);
end
if ~isfolder(output), mkdir(output); end
reserved = {'results.csv', 'results.mat', 'summary.csv', 'configuration.mat', ...
    'results-partial.tsv', 'warmup-results.tsv', 'timing-matlab.tsv'};
for j = 1:numel(reserved)
    assert(~isfile(fullfile(output, reserved{j})), 'KemenyRunner:OutputExists', ...
        'Output already exists: %s. Use a new OutputDir.', fullfile(output, reserved{j}));
end
configuration.OutputDir = output;
save(fullfile(output, 'configuration.mat'), 'configuration');
writetable(plan, fullfile(output, 'plan.csv'));
stage('reset', output);
launch = str2double(getenv('KEMENY_MATLAB_LAUNCH_EPOCH'));
if isfinite(launch)
    now_epoch = double(javaMethod('currentTimeMillis', 'java.lang.System')) / 1000;
    stage('record', 'MATLAB_start_to_runner', now_epoch-launch);
end
fprintf('RUNNER_SCHEMA: multi-partitions-v2\nMATLAB: %s\nResults: %s\n', version, output);
fprintf('BLAS: %s\nLAPACK: %s\n',version('-blas'),version('-lapack'));
metadata = struct('MATLAB',version,'Computer',computer,'BLAS',version('-blas'), ...
    'LAPACK',version('-lapack'),'JobID',getenv('SLURM_JOB_ID'), ...
    'Account',getenv('SLURM_JOB_ACCOUNT'),'Partition',getenv('SLURM_JOB_PARTITION'), ...
    'AllocatedCPUs',allocated,'Workers',workers,'WorkerThreads',1);
save(fullfile(output,'environment.mat'),'metadata');
fprintf('RESOURCES: cpus=%d workers=%d memory_limit_GiB=%g\n', allocated, workers, memory_gib);
fprintf('COUNTS: measured=%d warmup=%d; each enabled method/k group warms once before its own repeats.\n', ...
    sum(plan.RequestedRepeats), nnz(plan.RequestedRepeats));
fprintf('ClientCPUSeconds excludes worker processes. TOTAL timings include their child stages.\n');
for j = 1:height(plan)
    fprintf('PLAN: G%d T%d %s %s k=%g measured=%d warmup=%d\n', ...
        plan.PlanIndex(j), plan.TestIndex(j), plan.Matrix{j}, plan.Method{j}, plan.k(j), ...
        plan.RequestedRepeats(j), plan.Warmups(j));
end
old_threads = maxNumCompThreads;
thread_cleanup = onCleanup(@() maxNumCompThreads(old_threads)); %#ok<NASGU>
pool = [];
rows = cell(0, 17);
warm_rows = cell(0, 17);
batch_clock = stage('start', 'batch_TOTAL');
try
    files = cell(size(tests, 1), 1);
    for q = 1:size(tests, 1)
        selected = plan.TestIndex == q & plan.RequestedRepeats > 0;
        if any(selected), files{q} = resolve_matrix(tests{q, 1}, root); end
    end
    for q = 1:size(tests, 1)
        selected = find(plan.TestIndex == q);
        enabled = selected(plan.RequestedRepeats(selected) > 0);
        if isempty(enabled)
            fprintf('CASE_SKIPPED: T%d %s; all counts are zero.\n', q, tests{q, 1});
            continue
        end
        name = [tests{q, 1} '.mat'];
        stage('context', {q, name, NaN, 0, 0});
        fprintf('CASE_BEGIN: T%d %s\n', q, name);
        set_threads('matrix_prepare', allocated);
        token = stage('start', 'matrix_prepare_TOTAL');
        [P, pi_vec] = prepare_matrix(files{q});
        stage('end', 'matrix_prepare_TOTAL', token);
        n = size(P, 1); nz = nnz(P);
        fprintf('MATRIX: %s N=%d Nonzeros=%d\n', name, n, nz);
        for pos = enabled'
            method = plan.Method{pos}; k = plan.k(pos);
            if isfinite(k)
                assert(k <= n, 'KemenyRunner:Partitions', ...
                    '%s k=%d exceeds N=%d for %s.', method, k, n, name);
            end
            if strcmp(method, 'Direct')
                dense_gib = 8 * double(n)^2 / 2^30;
                fprintf('DIRECT_MEMORY: one_dense_array_GiB=%.3f allocation_GiB=%g; temporaries need extra memory.\n', dense_gib, memory_gib);
                assert(dense_gib < memory_gib, 'KemenyRunner:DirectMemory', ...
                    ['%s Direct needs %.3f GiB for ONE dense array, exceeding the %.3f GiB limit. ', ...
                     'Set its Direct repeat count to 0 or explicitly allocate sufficient memory. ', ...
                     'No algorithm was silently skipped.'], name, dense_gib, memory_gib);
                if 3*dense_gib > memory_gib
                    warning('KemenyRunner:MemoryHeadroom', ...
                        'Direct dense array %.3f GiB leaves limited room for inverse/temporary arrays in %.3f GiB. Peak memory is not known.', dense_gib, memory_gib);
                end
            end
        end
        reference = NaN;
        new_options = struct('UseParallel', true, 'BatchSize', opt.BatchSize, ...
            'SCBatchSize', opt.SCBatchSize);
        for pos = selected'
            method = plan.Method{pos}; k = plan.k(pos); count = plan.RequestedRepeats(pos);
            if count == 0
                fprintf('GROUP_SKIPPED: G%d T%d %s k=%g; requested repeats=0; no warmup.\n', pos, q, method, k);
                continue
            end
            if isfinite(k)
                stage('context', {q, name, k, 0, pos});
                set_threads('pool_setup', allocated);
                if isempty(pool)
                    pool = open_pool(workers, output);
                    pool_cleanup = onCleanup(@() close_pool(pool)); %#ok<NASGU>
                else
                    verify_pool(pool, workers);
                    fprintf('POOL_REUSED: workers=%d method=%s k=%d\n', workers, method, k);
                end
                fprintf('PARALLEL_BLOCK_LIMIT: %s k=%d; at most %d workers execute blocks.\n', method, k, min(k, workers));
            end
            switch method
                case 'Direct', call = @() kemenydirect(P);
                case 'Bini-Rec', call = @() recursivekemeny(P, pi_vec);
                case 'Bini-DirRec', call = @() recursivekemenydirect(P, pi_vec);
                case 'Du', call = @() kemeny_k_parallel(P, k, pi_vec, struct('solver', 'lu'));
                case 'New-full-rhs', call = @() kemeny_feedback(P, k, pi_vec, new_options);
            end
            for r = 0:count
                if r == 0, phase = 'warmup'; else, phase = 'measured'; end
                if isfinite(k), verify_pool(pool, workers); end
                threads = set_threads(method, allocated);
                stage('context', {q, name, k, r, pos});
                fprintf('CALL_BEGIN: G%d T%d %s %s k=%g repeat=%d/%d\n', pos, q, method, phase, k, r, count);
                [value, seconds, timing, status, message] = timed_call(call, [phase '.' method]);
                relative = NaN;
                if isfinite(reference) && strcmp(status, 'OK')
                    relative = abs(value-reference) / max(abs(reference), eps);
                end
                row = {q, name, n, nz, method, k, r, value, seconds, relative, ...
                    threads, timing.cpu_seconds, timing.mean_cpus, status, message, pos, reference};
                checkpoint(output, phase, row);
                if r > 0
                    rows(end+1, :) = row; %#ok<AGROW>
                else
                    warm_rows(end+1, :) = row; %#ok<AGROW>
                end
                stage('save');
                fprintf('CALL_DONE: G%d T%d %s %s k=%g %d/%d Status=%s K=%.17g Seconds=%.6f\n', ...
                    pos, q, method, phase, k, r, count, status, value, seconds);
                if ~strcmp(status, 'OK') && r == 0
                    error('KemenyRunner:WarmupFailed', '%s warmup failed: %s', method, message);
                end
            end
            if strcmp(method, 'Direct')
                % Use all measured Direct results; exclude warmup.
                [rows, warm_rows, reference, valid_direct] = finalize_reference(rows, warm_rows, q);
                refresh_checkpoints(output, rows, warm_rows);
                fprintf('REFERENCE: T%d DirectMean=%.17g successful_measured=%d requested=%d; warmup excluded.\n', ...
                    q, reference, valid_direct, count);
            end
        end
        clear P pi_vec call
    end
    stage('context', {0, '', NaN, 0, 0});
    stage('end', 'batch_TOTAL', batch_clock);
    close_pool(pool);
    [results, summary] = save_results(output, rows, warm_rows, plan, configuration);
    stage('save');
    assert(~any(strcmp(results.Status, 'FAILED')), 'KemenyRunner:MethodFailed', ...
        'Some measured calls failed; partial results and summary have been saved.');
    fprintf('KEMENY_BENCHMARK_SUCCESS: measured=%d, skipped_groups=%d\n', ...
        height(results), nnz(plan.RequestedRepeats == 0));
catch ME
    fprintf(2, 'KEMENY_BENCHMARK_FAILED: %s: %s\n', ME.identifier, ME.message);
    [results, summary] = save_results(output, rows, warm_rows, plan, configuration); %#ok<ASGLU>
    stage('save');
    rethrow(ME);
end
if nargout > 0, varargout{1} = results; end
if nargout > 1, varargout{2} = configuration; end
if nargout > 2, varargout{3} = summary; end
end

function close_pool(pool)
if isempty(pool) || ~isvalid(pool), return; end
try
    token = stage('start','pool.shutdown');
    delete(pool);
    stage('end','pool.shutdown',token);
    stage('save');
catch cleanup_error
    fprintf(2,'POOL_CLEANUP: %s\n',cleanup_error.message);
end
end

function tf = integer_scalar(x, lower)
tf = isnumeric(x) && isreal(x) && isscalar(x) && isfinite(x) && x >= lower && x == fix(x);
end

function tf = batch_width(x)
tf = (isnumeric(x) && isreal(x) && isscalar(x) && isinf(x) && x > 0) || integer_scalar(x, 1);
end

function tests = normalize_tests(tests)
assert(iscell(tests) && ismatrix(tests) && ~isempty(tests) && size(tests,2) == 6, ...
    'KemenyRunner:Tests', 'Each row must be {dataset, DirectCount, BiniRecCount, BiniDirRecCount, [DuK DuCount], [NewK NewCount]}.');
for q = 1:size(tests,1)
    name = tests{q,1};
    if isstring(name) && isscalar(name) && ~ismissing(name), name = char(name); end
    assert(ischar(name) && isrow(name) && ~isempty(strtrim(name)), ...
        'KemenyRunner:Tests', 'Row %d requires a dataset name.', q);
    name = strtrim(name);
    if endsWith(name, '.mat', 'IgnoreCase', true), name = name(1:end-4); end
    assert(~isempty(name) && isempty(regexp(name, '[\\/\x00-\x1F]', 'once')), ...
        'KemenyRunner:Tests', 'Row %d: use a dataset basename, without a directory.', q);
    tests{q,1} = name;
    for c = 2:4
        assert(integer_scalar(tests{q,c},0), 'KemenyRunner:Tests', ...
            'Row %d column %d must be a nonnegative integer repeat count.', q, c);
        tests{q,c} = double(tests{q,c});
    end
    for c = 5:6
        pairs = tests{q,c};
        assert(isnumeric(pairs) && isreal(pairs) && ~isempty(pairs) && ismatrix(pairs), ...
            'KemenyRunner:Tests', 'Row %d column %d requires numeric [k repeats] pairs.', q, c);
        if isrow(pairs) && mod(numel(pairs),2) == 0
            pairs = reshape(pairs,2,[])';
        else
            assert(size(pairs,2) == 2, 'KemenyRunner:Tests', ...
                'Row %d column %d must contain an even-length row or an N-by-2 matrix of [k repeats].', q, c);
        end
        for g = 1:size(pairs,1)
            assert(integer_scalar(pairs(g,1),2) && integer_scalar(pairs(g,2),0), ...
                'KemenyRunner:Tests', 'Row %d column %d pair %d requires integer k>=2 and repeats>=0.', q, c, g);
        end
        tests{q,c} = double(pairs);
    end
end
end

function plan = make_plan(tests)
methods = {'Direct','Bini-Rec','Bini-DirRec','Du','New-full-rhs'};
rows = cell(0,7);
for q = 1:size(tests,1)
    for m = 1:5
        if m <= 3
            pairs = [NaN tests{q,m+1}];
        else
            pairs = tests{q,m+1};
        end
        for g = 1:size(pairs,1)
            id = size(rows,1)+1;
            rows(end+1,:) = {id, q, [tests{q,1} '.mat'], methods{m}, ...
                pairs(g,1), pairs(g,2), double(pairs(g,2)>0)}; %#ok<AGROW>
        end
    end
end
plan = cell2table(rows,'VariableNames', ...
    {'PlanIndex','TestIndex','Matrix','Method','k','RequestedRepeats','Warmups'});
end

function [workers, allocated] = choose_workers(plan, requested)
allocation = getenv('SLURM_CPUS_PER_TASK');
if isempty(allocation)
    allocated = maxNumCompThreads;
else
    allocated = str2double(allocation);
    assert(integer_scalar(allocated,1), 'KemenyRunner:CPUs', 'Invalid SLURM_CPUS_PER_TASK.');
end
ks = plan.k(plan.RequestedRepeats > 0 & isfinite(plan.k));
if isempty(ks), workers = 0; return; end
assert(allocated >= 2, 'KemenyRunner:CPUs', 'Parallel methods require at least two allocated CPUs.');
if isempty(requested)
    workers = min(allocated-1, max(ks));
else
    workers = requested;
end
assert(workers < allocated, 'KemenyRunner:CPUs', ...
    'Workers=%d requires at least %d allocated CPUs; got %d.', workers, workers+1, allocated);
end

function limit = memory_limit(requested, allocated)
limit = requested;
per_node = str2double(getenv('SLURM_MEM_PER_NODE'));
per_cpu = str2double(getenv('SLURM_MEM_PER_CPU'));
if isfinite(per_node) && per_node > 0
    limit = min(limit, per_node/1024);
elseif isfinite(per_cpu) && per_cpu > 0
    limit = min(limit, per_cpu*allocated/1024);
end
% Local runs also guard against a known physical memory limit.
if ismac
    [rc, bytes] = system('sysctl -n hw.memsize');
    if rc == 0 && isfinite(str2double(bytes)), limit = min(limit, str2double(bytes)/2^30); end
elseif isunix && isfile('/proc/meminfo')
    hit = regexp(fileread('/proc/meminfo'), 'MemTotal:\s+(\d+)\s+kB', 'tokens', 'once');
    if ~isempty(hit), limit = min(limit, str2double(hit{1})/2^20); end
end
end

function threads = set_threads(method, allocated)
if any(strcmp(method, {'Du','New-full-rhs','pool_setup'})), target = 1; else, target = allocated; end
maxNumCompThreads(target);
threads = maxNumCompThreads;
assert(threads == target, 'KemenyRunner:Threads', 'Failed to set client thread budget.');
fprintf('THREAD_BUDGET: method=%s client_max_threads=%d worker_threads=1\n', method, threads);
end

function pool = open_pool(workers, output)
assert(license('test','Distrib_Computing_Toolbox'), 'KemenyRunner:Toolbox', 'Parallel Computing Toolbox is required.');
assert(isempty(gcp('nocreate')), 'KemenyRunner:ExistingPool', ...
    'An existing local pool is open. Close it explicitly before this benchmark.');
token = stage('start','pool.prepare_TOTAL');
if verLessThan('matlab','9.13'), cluster = parcluster('local'); else, cluster = parcluster('Processes'); end
cluster.NumThreads = 1;
cluster.NumWorkers = max(cluster.NumWorkers,workers);
storage = fullfile(output,'parallel_jobs');
if ~isfolder(storage), mkdir(storage); end
cluster.JobStorageLocation = storage;
create = stage('start','pool.create');
pool = parpool(cluster,workers,'IdleTimeout',Inf);
stage('end','pool.create',create);
fprintf('POOL_CREATED: workers=%d worker_threads=1 IdleTimeout=Inf\n',workers);
warm = stage('start','pool.warmup_and_verify');
try
spmd
    worker_threads = maxNumCompThreads;
end
actual = [worker_threads{:}];
fprintf('WORKER_THREADS_ACTUAL: %s\n',mat2str(actual));
assert(all(actual == 1),'KemenyRunner:WorkerThreads','Workers must have NumThreads=1.');
catch ME
    delete(pool);
    rethrow(ME);
end
stage('end','pool.warmup_and_verify',warm);
stage('end','pool.prepare_TOTAL',token);
end

function verify_pool(pool, workers)
current = gcp('nocreate');
assert(~isempty(current) && isequal(current,pool) && ...
    pool.NumWorkers == workers && pool.Cluster.NumThreads == 1, ...
    'KemenyRunner:PoolLost','The warmed process pool was lost or changed.');
end

function [value, seconds, timing, status, message] = timed_call(call, name)
token = stage('start',name);
status = 'OK'; message = ''; value = NaN; seconds = NaN;
try
    clock = tic;
    value = call();
    seconds = toc(clock);
catch ME
    status = 'FAILED'; message = [ME.identifier ': ' ME.message];
end
timing = stage('end',name,token);
if strcmp(status,'OK') && ~(isnumeric(value) && isreal(value) && isscalar(value) && isfinite(value))
    status = 'FAILED'; message = 'Nonfinite, nonreal or nonscalar Kemeny result.';
end
if strcmp(status,'FAILED'), value = NaN; seconds = NaN; end
end

function names = result_columns
names = {'TestIndex','Matrix','N','Nonzeros','Method','k','Repeat','Kemeny','Seconds', ...
    'RelativeError','ClientMaxThreads','ClientCPUSeconds','MeanClientCPUs','Status','ErrorMessage', ...
    'PlanIndex','DirectMean'};
end

function checkpoint(output, phase, row)
if strcmp(phase,'warmup'), file = 'warmup-results.tsv'; else, file = 'results-partial.tsv'; end
file = fullfile(output,file);
new_file = ~isfile(file);
fid = fopen(file,'a','n','UTF-8');
assert(fid ~= -1,'KemenyRunner:Checkpoint','Cannot open %s.',file);
closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
if new_file, fprintf(fid,'%s\n',strjoin(result_columns,sprintf('\t'))); end
write_row(fid,row);
check_write(fid,file);
end

function write_row(fid,row)
for j = [2 5 14 15]
    row{j} = regexprep(row{j}, '[\t\r\n]', ' ');
end
fprintf(fid,['%d\t%s\t%d\t%d\t%s\t%g\t%d\t%.17g\t%.9f\t%.17g', ...
    '\t%d\t%.6f\t%.6f\t%s\t%s\t%d\t%.17g\n'],row{:});
end

function check_write(fid,file)
[message, code] = ferror(fid);
assert(code == 0,'KemenyRunner:Checkpoint','Failed writing %s: %s',file,message);
end

function [rows, warm_rows, reference, count] = finalize_reference(rows, warm_rows, q)
reference = NaN; count = 0;
if ~isempty(rows)
    selected = cell2mat(rows(:,1)) == q & strcmp(rows(:,5),'Direct') & ...
        strcmp(rows(:,14),'OK') & isfinite(cell2mat(rows(:,8)));
    values = cell2mat(rows(selected,8));
    count = numel(values);
    if count > 0, reference = mean(values); end
end
rows = apply_reference(rows,q,reference);
warm_rows = apply_reference(warm_rows,q,reference);
end

function rows = apply_reference(rows,q,reference)
if isempty(rows), return; end
for i = find(cell2mat(rows(:,1)) == q)'
    rows{i,17} = reference;
    rows{i,10} = NaN;
    if isfinite(reference) && strcmp(rows{i,14},'OK')
        rows{i,10} = abs(rows{i,8}-reference) / max(abs(reference),eps);
    end
end
end

function refresh_checkpoints(output,rows,warm_rows)
% Finish writing a replacement before replacing the previous checkpoint.
files = {'results-partial.tsv','warmup-results.tsv'};
histories = {rows,warm_rows};
for j = 1:2
    if isempty(histories{j}), continue; end
    temporary = tempname(output);
    write_history(temporary,histories{j});
    [ok,message] = movefile(temporary,fullfile(output,files{j}),'f');
    assert(ok,'KemenyRunner:Checkpoint','Cannot finalize checkpoint: %s',message);
end
end

function write_history(file,rows)
fid = fopen(file,'w','n','UTF-8');
assert(fid ~= -1,'KemenyRunner:Checkpoint','Cannot open %s.',file);
closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',strjoin(result_columns,sprintf('\t')));
for i = 1:size(rows,1), write_row(fid,rows(i,:)); end
check_write(fid,file);
end

function [results, summary] = save_results(output, rows, warm_rows, plan, configuration)
% Preserve measured Direct references and incomplete groups after failure.
for q = unique(plan.TestIndex)'
    [rows,warm_rows] = finalize_reference(rows,warm_rows,q);
end
refresh_checkpoints(output,rows,warm_rows);
% An empty result table still has stable numeric/text column types.
prototype = {0,'',0,0,'',NaN,0,NaN,NaN,NaN,0,0,0,'','',0,NaN};
results = cell2table([prototype; rows],'VariableNames',result_columns);
results(1,:) = [];
stats_rows = cell(height(plan),17);
for j = 1:height(plan)
    sub = results(results.PlanIndex == plan.PlanIndex(j),:);
    good = strcmp(sub.Status,'OK'); t = sub.Seconds(good);
    expected = plan.RequestedRepeats(j);
    failed = nnz(~good); missing = max(0,expected-height(sub));
    if expected == 0
        status = 'SKIPPED';
    elseif failed > 0
        status = 'FAILED';
    elseif missing > 0
        status = 'INCOMPLETE';
    else
        status = 'COMPLETE';
    end
    stats = NaN(1,6);
    if ~isempty(t)
        stats(1:4) = [median(t),mean(t),min(t),max(t)];
        if numel(t)>1, stats(5)=std(t,0); end
        errors = sub.RelativeError(good & isfinite(sub.RelativeError));
        if ~isempty(errors), stats(6)=max(errors); end
    end
    direct = results(results.TestIndex == plan.TestIndex(j) & strcmp(results.Method,'Direct') & strcmp(results.Status,'OK'),:);
    reference = NaN;
    if ~isempty(direct), reference = mean(direct.Kemeny); end
    stats_rows(j,:) = [{plan.TestIndex(j),plan.Matrix{j},plan.Method{j},plan.k(j), ...
        expected,numel(t),failed,missing,status},num2cell(stats),{plan.PlanIndex(j),reference}];
end
summary = cell2table(stats_rows,'VariableNames',{'TestIndex','Matrix','Method','k', ...
    'RequestedRepeats','ValidRepeats','FailedRepeats','MissingRepeats','Status', ...
    'MedianSeconds','MeanSeconds','MinSeconds','MaxSeconds','StdSeconds','MaxRelativeError','PlanIndex','DirectMean'});
writetable(results,fullfile(output,'results.csv'));
writetable(summary,fullfile(output,'summary.csv'));
save(fullfile(output,'results.mat'),'results','summary','configuration');
end

function token = stage(action, name, value)
persistent origin context rows output
token = [];
switch action
    case 'reset'
        origin = tic; context = {0,'',NaN,0,0}; rows = cell(0,10); output = name;
    case 'context'
        context = name;
    case 'start'
        fprintf('[BEGIN +%.3fs G%d T%d R%d] %s\n',toc(origin),context{5},context{1},context{4},name);
        token = struct('wall',tic,'cpu',cputime);
    case 'end'
        seconds = toc(value.wall); cpu = cputime-value.cpu;
        token = struct('seconds',seconds,'cpu_seconds',cpu,'mean_cpus',cpu/max(seconds,eps));
        rows(end+1,:) = [context,{name,seconds,toc(origin),cpu,token.mean_cpus}];
        fprintf('[TIME +%.3fs G%d T%d R%d] %s %.6f s; client_cpu=%.6f s\n', ...
            toc(origin),context{5},context{1},context{4},name,seconds,cpu);
    case 'record'
        rows(end+1,:) = [context,{name,value,toc(origin),NaN,NaN}];
        fprintf('[TIME] %s %.6f s\n',name,value);
    case 'save'
        report = cell2table(rows,'VariableNames',{'TestIndex','Matrix','k','Repeat','PlanIndex','Stage', ...
            'Seconds','RunnerElapsedSeconds','ClientCPUSeconds','MeanClientCPUs'});
        writetable(report,fullfile(output,'timing-matlab.tsv'),'FileType','text','Delimiter','\t');
end
end

function file = resolve_matrix(name, root)
name = [name '.mat'];
if isfile(fullfile(root, name))
    file = fullfile(root, name);
elseif isfile(fullfile(root, 'matrix', name))
    file = fullfile(root, 'matrix', name);
else
    error('KemenyRunner:FileNotFound', 'Matrix file not found: %s', name);
end
end

function [P, pi_vec, meta] = prepare_matrix(file)
token = stage('start','matrix.load_mat');
data = load(file);
stage('end','matrix.load_mat',token);
token = stage('start','matrix.graph_preprocess');
if isfield(data, 'Problem') && isstruct(data.Problem) && isfield(data.Problem, 'A')
    % Match Du's input pipeline; choose exactly one largest component.
    A = spones(sparse(double(data.Problem.A)));
    validate_adjacency(A);
    original_n = size(A, 1);
    if ~issymmetric(A), A = max(A, A'); end
    G = graph(A);
    [labels, sizes] = conncomp(G);
    [~, largest] = max(sizes);
    G = subgraph(G, labels == largest);
    A = adjacency(G);
    kind = 'Problem.A -> unweighted undirected largest component';
elseif isfield(data, 'A')
    A = sparse(double(data.A));
    validate_adjacency(A);
    original_n = size(A, 1);
    kind = 'A -> original weights, diagonal and direction preserved';
else
    error('KemenyRunner:MatFormat', '%s must contain Problem.A or A.', file);
end
n = size(A, 1);
if n < 2
    error('KemenyRunner:TooSmall', 'At least two retained vertices are required.');
end
if issymmetric(A)
    labels = conncomp(graph(A));
else
    labels = conncomp(digraph(A), 'Type', 'strong');
end
if max(labels) ~= 1
    error('KemenyRunner:Reducible', ...
        'Raw A is not irreducible. Its weights and vertex set were preserved; supply a connected dataset.');
end
stage('end','matrix.graph_preprocess',token);
token = stage('start','matrix.dissect');
p = dissect(A);
A = A(p, p);
stage('end','matrix.dissect',token);
token = stage('start','matrix.transition_and_stationary');
degree = full(sum(A, 2));
if any(degree <= 0)
    error('KemenyRunner:ZeroDegree', 'The retained matrix contains a zero-degree row.');
end
P = spdiags(degree, 0, n, n) \ A;
if issymmetric(A)
    pi_vec = degree / sum(degree);
    pi_method = 'normalized weighted degree';
else
    % Replace one redundant stationarity equation with a sparse anchoring
    % equation, then normalize. This also works for periodic chains.
    M = speye(n) - P';
    M(1, :) = 0;
    M(1, 1) = 1;
    rhs = sparse(1, 1, 1, n, 1);
    pi_vec = full(M \ rhs);
    pi_vec = pi_vec / sum(pi_vec);
    pi_method = 'sparse stationary solve';
end
if any(~isfinite(pi_vec)) || any(pi_vec <= 0) || ...
        norm(P' * pi_vec - pi_vec, 1) > 1e-10
    error('KemenyRunner:Stationarity', 'Stationary distribution failed positivity or residual checks.');
end
meta = struct('kind', kind, 'original_n', original_n, 'pi_method', pi_method);
stage('end','matrix.transition_and_stationary',token);
end

function validate_adjacency(A)
if ~ismatrix(A) || size(A, 1) ~= size(A, 2) || isempty(A) || ~isreal(A) || ...
        any(~isfinite(nonzeros(A))) || any(nonzeros(A) < 0)
    error('KemenyRunner:Adjacency', 'The input must be a nonempty square real nonnegative finite matrix.');
end
end
