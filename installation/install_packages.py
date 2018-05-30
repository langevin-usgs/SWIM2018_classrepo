import os
import platform
import sys
import shutil

try:
    # For Python 3.0 and later
    from urllib.request import urlopen
except ImportError:
    # Fall back to Python 2's urllib2
    from urllib2 import urlopen
import platform
from subprocess import Popen, PIPE, STDOUT

pythonpath = os.path.join('..', 'miniconda3', 'python')
scriptpath = os.path.join('..', 'miniconda3', 'Scripts')
condacommand = os.path.join(scriptpath, 'conda')
pipcommand = os.path.join(scriptpath, 'pip')

exeloc = {'Windows': 'python',
          'Darwin': 'python'}


# simple function to print a message to STDOUT
def printmsg(msg):
    print(msg)
    return


def run_and_print(cmds):
    for cmd in cmds:
        print(' {}'.format(cmd))
        cmd_list = cmd.split()
        p = Popen(cmd_list, stdout=PIPE, stderr=STDOUT)
        while True:
            line = p.stdout.readline()
            c = line.decode('utf-8')
            if c != '':
                c = c.rstrip('\r\n')
                print('{}'.format(c))
            else:
                break
    return


def root_install():
    pip_list = []
    conda_list = ['jupyter',
                  'scipy',
                  'pyshp',
                  'nose',
                  'pandas',
                  'flopy']

    # prepare the pip installs to run in a single command (after activating env)
    cmds = ['{} config --add channels conda-forge'.format(condacommand)]
    cmds.append('{} config --set ssl_verify false'.format(condacommand))
    cmds.append('{} update conda -y'.format(condacommand))
    cmds.append('{} update -y --all'.format(condacommand))
    for c in conda_list:
        cmds.append('{} install {} -y'.format(condacommand, c))
    cmds.append('{} info'.format(condacommand))
    # add pip installs
    cmds.append('{} -m pip install --upgrade pip'.format(pythonpath))
    for p in pip_list:
        cmds.append('{} install {}'.format(pipcommand, p))
    
    run_and_print(cmds)
    
    print('\nRunning tests of installed python packages in root installation...')
    print('    using python installed in "{}"'.format(pythonpath))
    cmds = [pythonpath + ' -m nose -v test_root_install.py']
    run_and_print(cmds)    

    return

    
if __name__ == "__main__":

    install_root = True
    
    # install packages to root environment
    if install_root:
        root_install()
