import subprocess, os, numpy as np, re, pickle, time

def addArguments(parser, combinedJobs=False):
    parser.add_argument('--nodes', type=int)
    parser.add_argument('--chainsPerNode', type=int)
    parser.add_argument('--coresPerNode', type=int)
    parser.add_argument('--mem_per_node', type=int, help="Memory in MB per node")
    parser.add_argument('--walltime')
    if combinedJobs:
        parser.add_argument('--combineOneJobName',
                            help='run all one after another, under one job submission (good for many fast operations)')
        parser.add_argument('--runsPerJob', type=int, default=int(os.environ.get('COSMOMC_runsPerJob', '1')),
                            help='submit multiple mpi runs at once from each job script (e.g. to get more than one run per node)')

    parser.add_argument('--job_template', help="template file for the job submission script")
    parser.add_argument('--program', help='actual program to run (default: ./cosmomc)')
    parser.add_argument('--queue', help='name of queue to submit to')
    parser.add_argument('--qsub', help='option to change qsub command to something else')

    parser.add_argument('--dryrun', action='store_true')
    parser.add_argument('--no_sub', action='store_true')


def replacePlaceholders(txt, vals):
    txt = txt.replace('\r', '')
    for name, value in vals.iteritems():
        txt = txt.replace('##' + name + '##', str(value))
    return txt

def extractValue(template, name):
    match = re.search('##' + name + ':(.*)##', template)
    if match: return match.group(1).strip()
    return None

def getDefaulted(key_name, default=None, tp=str, template=None, ext_env=None, **kwargs):
    val = kwargs.get(key_name)
    if val is None and template is not None:
            val = extractValue(template, 'DEFAULT_' + key_name)
    if val is None: val = os.environ.get('COSMOMC_' + key_name, None)
    if val is None and ext_env: val = os.environ.get('ext_env', None)
    if val is None: val = default
    if val is None: return None
    return tp(val)

def checkArguments(**kwargs):
    submitJob(None, None, msg=True, **kwargs)

class jobSettings():

    def __init__(self, jobName, msg=False, **kwargs):
        self.jobName = jobName
        self.job_template = getDefaulted('job_template', 'job_script', **kwargs)
        with open(self.job_template, 'r') as f:
            template = f.read()

        self.coresPerNode = getDefaulted('coresPerNode', 16, tp=int, template=template, **kwargs)
        self.chainsPerNode = getDefaulted('chainsPerNode', 4, tp=int, template=template, **kwargs)
        self.nodes = getDefaulted('nodes', 1, tp=int, template=template, **kwargs)
        self.nchains = self.nodes * self.chainsPerNode

        self.runsPerJob = getDefaulted('runsPerJob', 1, tp=int, template=template, **kwargs)
        # also defaulted at input so should be set here unless called programmatically

        self.omp = self.coresPerNode / (self.chainsPerNode * self.runsPerJob)
        if self.omp != np.floor(self.omp): raise Exception('Chains must each have equal number of cores')
        if msg:
            print 'Job parameters: %i cosmomc runs of %i chains on %i nodes, each node with %i MPI chains, each chain using %i OpenMP cores (%i cores per node)' % (
             self.runsPerJob, self.nchains, self.nodes, self.chainsPerNode, self.omp, self.coresPerNode)

        self.mem_per_node = getDefaulted('mem_per_node', 63900, tp=int, template=template, **kwargs)
        self.walltime = getDefaulted('walltime', '24:00:00', template=template, **kwargs)
        self.program = getDefaulted('program', './cosmomc', template=template, **kwargs)
        self.queue = getDefaulted('queue', '', template=template, **kwargs)
        self.gridEngine = getDefaulted('GridEngine', 'PBS', template=template, **kwargs)
        self.qsub = getDefaulted('qsub', ('qsub', 'msub')[self.gridEngine == 'MOAB'], template=template, **kwargs)
        self.runCommand = extractValue(template, 'RUN')


class jobIndex():
    """
     Stores the mappings between job Ids, jobNames
    """
    def __init__(self):
        self.jobSettings = dict()
        self.jobNames = dict()
        self.rootNames = dict()
        self.jobSequence = []

def loadJobIndex(batchPath, must_exist=False):
    if batchPath is None: batchPath = './scripts/'
    fileName = os.path.join(batchPath, 'jobIndex.pyobj')
    if os.path.exists(fileName):
        with open(fileName, 'rb') as inp:
            return pickle.load(inp)
    else:
        if not must_exist: return jobIndex()
        return None

def saveJobIndex(obj, batchPath=None):
    if batchPath is None: batchPath = './scripts/'
    with open(os.path.join(batchPath, 'jobIndex.pyobj'), 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)


def addJobIndex(batchPath, jobName, j):
    if batchPath is None: batchPath = './scripts/'
    index = loadJobIndex(batchPath)
    index.jobSettings[j.jobId] = j
    index.jobNames[j.jobName] = j.jobId
    for name in j.names:
        index.rootNames[name] = j.jobId
    index.jobSequence.append(j.jobId)
    saveJobIndex(index, batchPath)

def submitJob(jobName, paramFiles, sequential=False, msg=False, **kwargs):

    j = jobSettings(jobName, msg, **kwargs)
    if kwargs.get('dryrun', False) or paramFiles is None: return

    if isinstance(paramFiles, basestring): paramFiles = [paramFiles]
    paramFiles = [ini.replace('.ini', '') for ini in paramFiles]

    j.runsPerJob = (len(paramFiles), 1)[sequential]
    # adjust omp for the actual number (may not be equal to input runsPerJob because of non-integer multiple)
    j.omp = int(j.coresPerNode / (j.chainsPerNode * j.runsPerJob))

    j.path = os.getcwd()
    j.onerun = (0, 1)[len(paramFiles) == 1 or sequential]
    # mem = mem_per_node * numnodes
    vals = dict()
    vals['JOBNAME'] = jobName
    vals['OMP'] = j.omp
    vals['MEM_MB'] = j.mem_per_node
    vals['WALLTIME'] = j.walltime
    vals['NUMNODES'] = j.nodes
    vals['NUMRUNS'] = j.runsPerJob
    vals['NUMMPI'] = j.nchains
    vals['CHAINSPERNODE'] = j.chainsPerNode
    vals['PPN'] = j.chainsPerNode * j.runsPerJob * j.omp
    vals['MPIPERNODE'] = j.chainsPerNode * j.runsPerJob
    vals['NUMTASKS'] = j.nchains * j.runsPerJob
    vals['ROOTDIR'] = j.path
    vals['ONERUN'] = j.onerun
    vals['PROGRAM'] = j.program
    vals['QUEUE'] = j.queue

    j.names = [ os.path.basename(param) for param in paramFiles]

    commands = []
    for param, name in zip(paramFiles, j.names):
        ini = param + '.ini'
        if j.runCommand is not None:
            vals['INI'] = ini
            vals['INIBASE'] = name
            command = replacePlaceholders(j.runCommand, vals) + (' &', '')[sequential]
        else:
            command = ('time mpirun -np %i %s %s > ./scripts/%s.log 2>&1 %s' %
                         (j.nchains, j.program, ini, os.path.basename(ini), ('&', '')[sequential]))
        commands.append(command)

    vals['COMMAND'] = "\n".join(commands)
    with open(j.job_template, 'r') as f:
        template = f.read()
        script = replacePlaceholders(template, vals)
        scriptRoot = './scripts/' + jobName
        scriptName = scriptRoot + '_subscript'
        open(scriptName, 'w').write(script)
        if len(paramFiles) > 1:
            open(scriptRoot + '.batch', 'w').write("\n".join(paramFiles))
        if not kwargs.get('no_sub', False):
            res = subprocess.check_output(j.qsub + ' ' + scriptName, shell=True).strip()
            if not res: print 'No qsub output'
            else:
                j.paramFiles = paramFiles
                j.jobId = res
                j.subTime = time.time()
                open(scriptRoot + '.sub', 'w').write(res)
                addJobIndex(kwargs.get('batchPath'), jobName, j)


def queue_job_details(batchPath=None, running=True, queued=True, warnNotBatch=True):
    """
    Return: list of jobIds, list of jobNames, list of list names
    """
    index = loadJobIndex(batchPath)
    if not index:
        print 'No existing job index found'
        return []
    res = subprocess.check_output('showq -u $USER', shell=True).strip()
    res = res.split("\n")
    names = []
    jobNames = []
    ids = []
    infos = []
    for line in res[2:]:
        if ' ' + os.environ.get('USER') + ' ' in line and (queued and not ' Running ' in line or running and ' Running ' in line):
            items = line.split()
            jobId = items[0].split('.')
            if jobId[0].upper() == 'TOTAL': continue
            if len(jobId) == 1 or jobId[0].isdigit(): jobId = jobId[0]
            else: jobId = jobId[1]
            j = index.jobSettings.get(jobId)
            if j is None:
                if warnNotBatch: print '...Job ' + jobId + ' not in this batch, skipping'
                continue

            names += [j.names]
            jobNames += [j.jobName]
            ids += [jobId]
            infos += [line]
    return ids, jobNames, names, infos

def queue_job_names(batchPath=None, running=False, queued=True):
    lists = queue_job_details(batchPath, running, queued)[2]
    names = []
    for nameset in lists:
        names += nameset
    return names
