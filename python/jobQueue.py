import subprocess, os, numpy as np

def addArguments(parser, combinedJobs=False):
    parser.add_argument('--nodes', type=int, default=int(os.environ.get('COSMOMC_nodes', '1')))
    parser.add_argument('--chainsPerNode', type=int, default=int(os.environ.get('COSMOMC_chainsPerNode', '4')))
    parser.add_argument('--coresPerNode', type=int, default=int(os.environ.get('COSMOMC_coresPerNode', '16')))
    parser.add_argument('--mem_per_node', type=int, default=int(os.environ.get('COSMOMC_walltime', '63900')))
    if combinedJobs:
        parser.add_argument('--combineOneJobName', default=None,
                            help='run all one after another, under one job submission (good for many fast operations)')
        parser.add_argument('--runsPerJob', type=int, default=1,
                            help='submit multiple mpi runs at once from each job script (e.g. to get more than one run per node)')

    parser.add_argument('--walltime', default=os.environ.get('COSMOMC_walltime', '24:00:00'))
    parser.add_argument('--job_template', default=os.environ.get('COSMOMC_job_script', 'job_script'),
                        help="template file for the job submission script")
    parser.add_argument('--program', default='./cosmomc')
    parser.add_argument('--dryrun', action='store_true')
    parser.add_argument('--no_sub', action='store_true')

def getArgsOmp(args, msg=True, runsPerJob=1):
    omp = args.coresPerNode / (args.chainsPerNode * runsPerJob)
    if omp != np.floor(omp): raise Exception('Chains must each have equal number of cores')
    if msg:
        print 'Job parameters: %i cosmomc runs on %i nodes, each with %i MPI chains, each chain using %i OpenMP cores (%i cores per node)' % (runsPerJob, args.nodes,
            args.chainsPerNode, omp, args.coresPerNode)
    return omp


def queued_jobs():
    res = subprocess.check_output('qstat -u $USER', shell=True)
    res = res.split("\n")
    names = []
    for line in res[4:]:
        if 'Q 00:00' in line:
            items = line.split()
            jobid = items[0].split('.')[0]
            output = subprocess.check_output('qstat -f ' + str(jobid) , shell=True).split('\n')
            pars = []
            current = ''
            for L in output:
                if '=' in L:
                    if len(current) > 0:
                        pars.append(current)
                    current = L.strip()
                else: current += L.strip()
            if len(current) > 0: pars.append(current)
            props = dict()
            for L in pars[1:]:
                (key, val) = L.split('=', 1)
                props[key.strip()] = val.strip()
            names.append(props['Job_Name'])
    return names


def queued_jobs_PBS():  # what used to use on Darwin
    res = subprocess.check_output('qstat -u $USER', shell=True)
    res = res.split("\n")
    names = []
    for line in res[4:]:
        if 'master' in line:
            items = line.split()
            jobid = items[0].split('.')[0]
            output = subprocess.check_output('qstat -f ' + str(jobid) , shell=True).split('\n')
            pars = []
            current = ''
            for L in output:
                if '=' in L:
                    if len(current) > 0:
                        pars.append(current)
                    current = L.strip()
                else: current += L.strip()
            if len(current) > 0: pars.append(current)
            props = dict()
            for L in pars[1:]:
                (key, val) = L.split('=', 1)
                props[key.strip()] = val.strip()
            names.append(props['Job_Name'])
    return names

def replacePlaceholders(txt, vals):
    txt = txt.replace('\r', '')
    for name, value in vals.iteritems():
        txt = txt.replace('##' + name + '##', str(value))
    return txt


def submitJob(jobName, paramFiles, sequential=False, **kwargs):
    if isinstance(paramFiles, basestring): paramFiles = [paramFiles]
    job_template = kwargs.get('job_template', 'job_script')
    nodes = int(kwargs.get('nodes', 1))
    omp = int(kwargs.get('omp', 4))
    chainsPerNode = int(kwargs.get('chainsPerNode', 1))
    mem_per_node = int(kwargs.get('mem_per_node', 63900))
    walltime = kwargs.get('walltime', '24:00:00')
    qsub = kwargs.get('qsub', 'qsub')
    program = kwargs.get('program', './cosmomc')

    ppn = chainsPerNode * omp
    nchains = nodes * chainsPerNode
    runsPerJob = (len(paramFiles), 1)[sequential]
    # mem = mem_per_node * numnodes
    vals = dict()
    vals['JOBNAME'] = jobName
    vals['OMP'] = omp
    vals['MEM_MB'] = mem_per_node
    vals['WALLTIME'] = walltime
    vals['NUMNODES'] = nodes
    vals['PPN'] = ppn * runsPerJob
    vals['NUMMPI'] = nchains
    vals['NUMTASKS'] = nchains * runsPerJob
    vals['ROOTDIR'] = os.getcwd()
    vals['ONERUN'] = (0, 1)[len(paramFiles) == 1 or sequential]

    commands = []
    for param in paramFiles:
        ini = param
        if ini[-4:] != '.ini': ini += '.ini'
        commands.append('time mpirun -ppn %i -np %i %s %s > ./scripts/%s.log 2>&1 %s' %
                         (chainsPerNode * runsPerJob, nchains, program, ini, os.path.basename(ini), ('&', '')[sequential]))
    vals['COMMAND'] = "\n".join(commands)
    script = replacePlaceholders(open(job_template, 'r').read(), vals)
    scriptName = './scripts/' + jobName + '_subscript'
    open(scriptName, 'w').write(script)
    if len(paramFiles) > 1:
        open('./scripts/' + jobName + '.batch', 'w').write("\n".join(paramFiles))
    if not kwargs.get('no_sub', False): os.system(qsub + ' ' + scriptName)

