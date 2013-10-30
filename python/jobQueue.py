import subprocess, os

def queued_jobs():
    res = subprocess.check_output('qstat -u $USER', shell=True)
    res = res.split("\n")
    names=[]
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
    for name, value in vals.iteritems():
        txt = txt.replace('##'+name+'##',str(value))
    return txt

def submitJob(jobName, paramFiles, pbs_template='job_script', numnodes=1, omp = 4, chainsPerNode=1,mem_per_node=63900,walltime='24:00:00',
               qsub='qsub' ):
    ppn = chainsPerNode * omp
    nchains = numnodes * chainsPerNode
    mem = mem_per_node * numnodes 
    vals = dict()
    vals['JOBNAME'] = jobName
    vals['OMP'] = omp
    vals['MEM_MB']=mem
    vals['WALLTIME']=walltime
    vals['NUMNODES']=numnodes
    vals['PPN']=ppn   
    vals['ROOTDIR'] = os.getcwd()
    commands=[]
    if isinstance(paramFiles, basestring): paramFiles=[paramFiles]
    for param in paramFiles:
        ini=param
        if ini[-4:] != '.ini': ini+= '.ini'
        commands.append('time mpirun -ppn %i -np %i ./cosmomc %s > ./scripts/%s.log 2>&1' % (ppn, nchains, ini, jobName ))
    vals['COMMAND']="\n".join(commands)
    script = replacePlaceholders(open(pbs_template, 'r').read(),vals)
    scriptName = './scripts/' +jobName+'_subscript'
    open(scriptName, 'w').write(script)
    if len(paramFiles)>1:
        open('./scripts/'+jobName+'.batch', 'w').write("\n".join(paramFiles))
    os.system(qsub+' '+scriptName)
    
