from glob import glob
import os
from subprocess import run, DEVNULL
import time

def merge_chain(cspath,output_path,protocol):
    files = os.listdir(cspath)
    
    file_list = []
    for f in files:
        if '.rebuilt.pdb' in f:
            if (f.split('.rebuilt.pdb')[0]+'_real_space_refined_000.pdb') not in files:
                if len(open(os.path.join(cspath,f)).readlines())>60:
                    file_list.append(f)
        elif 'real_space_refined_000.pdb' in f:
            file_list.append(f)
    file_list.sort()
    aid = 1
    
    with open(output_path,'w') as acf:
        for f in file_list:
            cid = f.split(protocol)[-1].split('_')[1]
            with open(cspath + '/'+ f, 'r') as pf:
                lines = pf.readlines()
                for l in lines:
                    if l.startswith('ATOM') and 'nan' not in l:
                        if len(l)<70:
                            acf.write(l[:4] + str(aid).rjust(7, ' ') + l[11:21] + cid + l[22:54])
                            acf.write(f'  1.00  0.00           {l[13]}\n')
                        else:
                            acf.write(l[:4] + str(aid).rjust(7, ' ') + l[11:21] + cid + l[22:])
                        aid+=1

def run_refine(dir,em_path,resol,output_path,phenix_act,protocol):
    em_path=os.path.abspath(em_path)
    file_path=os.path.join(dir,'chain_split')
    for file in glob(os.path.join(file_path,'*.rebuilt.pdb')):
        name=file.split('/')[-1]
        cmd=f'bash refine.sh {file_path} "phenix.real_space_refine {name} {em_path} resolution={resol}" {phenix_act}'
        task_num=60
    
        while True:
            try:
                results=os.popen('ps -aux|grep real_space_refine.py|wc -l').readlines()
                task_num=int(results[0].strip())-1
            except:
                task_num=60
            if task_num<=32:
                time.sleep(2)
                break
            time.sleep(10)
        run([cmd], stdout=DEVNULL, stderr=DEVNULL, shell=True)
        # print(cmd)
    while True:
        try:
            results=os.popen('ps -aux|grep real_space_refine.py|wc -l').readlines()
            task_num=int(results[0].strip())-1
            print('task_num:',task_num)
        except:
            task_num=1
        if task_num<=1:
            break
        else:
            time.sleep(10)
    merge_chain(file_path,output_path,protocol)
