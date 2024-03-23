
import sys
import subprocess
import shlex
import time
import os
import shutil

def build_all_atom(dir,pulchar_path):
    path = os.path.join(dir,'chain_split/')
    path = os.path.abspath(path)
    pulchar_path = os.path.abspath(pulchar_path)
    process_list = []
    n_job_per_node = 30
    filelist = os.listdir(path)
    for f in filelist:
        prefix = f.split('.')[0]
        
        if 'rebuilt' not in f and 'pdb' in f and not os.path.exists(os.path.join(path , prefix + '.rebuilt.pdb')):
            command = pulchar_path + ' {} -c '.format(f)
            args = shlex.split(command)
            # print(command)
            # print(path+p)
            with open(os.path.join(path,'{}.log'.format(prefix)), 'w') as log:
                if len(process_list) < n_job_per_node:
                    process_list.append(subprocess.Popen(args, cwd = path, stdout=log))
                else:
                    have_finished = False
                    while True:
                        for i in range(len(process_list)):
                            if process_list[i].poll() is not None:
                                process_list[i] = subprocess.Popen(args, cwd = path, stdout=log)
                                have_finished = True
                                break
                        if have_finished:
                            break
                        time.sleep(0.5)
            # print(f)
    for p in process_list:
        p.wait()

def split_chain(dir,file):
    #path = '/home/chens/projects/AutoEM/data/outputs'
    filelist = os.listdir(dir)
    cspath = os.path.join(dir,'chain_split/')
    if os.path.exists(cspath):
        shutil.rmtree(cspath)
    os.makedirs(cspath)

    pdbfile=os.path.join(dir,file)
    if not os.path.exists(pdbfile):
        return

    prefix = pdbfile.split('/')[-1].split('.')[0]
    lastcid = ''
    lastrid = ''
    newlines = []
    if not os.path.exists(cspath):
        os.makedirs(cspath)

    with open(pdbfile, 'r') as aafile:
        lines = aafile.readlines()
        for l in lines:
            #print(l)
            if not l.startswith('ATOM'):
                #newlines.append(l)
                continue
            cid = l[21]
            rid = int(l[22:26])
            if lastcid == '':
                lastcid = cid

            if lastrid == '':
                lastrid = rid
            
            if lastcid != cid or rid-lastrid not in [0,1]:
                if len(newlines)>3:
                    rid_name=f'{lastrid//1000%10}{lastrid//100%10}{lastrid//10%10}{lastrid%10}'
                    with open(os.path.join(cspath,f'{prefix}_{lastcid}_{rid_name}.pdb'), 'w') as newfile:
                        for nl in newlines:
                            newfile.writelines(nl)
                newlines=[]
                lastcid = cid
            newlines.append(l)
            lastrid = rid

        
        if len(newlines)>3:
            rid_name=f'{lastrid//1000%10}{lastrid//100%10}{lastrid//10%10}{lastrid%10}'
            with open(os.path.join(cspath,f'{prefix}_{lastcid}_{rid_name}.pdb'), 'w') as newfile:
                    for nl in newlines:
                        newfile.writelines(nl)


def run_pulchra(dir,file,pulchra_path):
    split_chain(dir,file)
    build_all_atom(dir,pulchra_path)
