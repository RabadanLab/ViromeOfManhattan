#!/usr/bin/env python

"""
    Miscellaneous helper functions (shell-related, SGE-related,
    bioinformatics-related) to be imported as needed
    ~~~~~~
"""

import sys
import subprocess
import os
import time
import ConfigParser
from distutils import spawn

# -------------------------------------

### shell-related functions

# -------------------------------------

def mkdirp(mydir):
    """The loose equivalent of mkdir -p"""

    try:
        os.mkdir(mydir)
    except:
        pass

# -------------------------------------

def check_dependencies(required_programs):
    """Check for errors, check dependencies"""

    # required_programs - define required programs in list, e.g., ["java", "samtools", "tabix"]

    for j in required_programs:
        if not spawn.find_executable(j):
            print("[ERROR] Can't find " + j + ". Please add it to your PATH")
            sys.exit(1)

# -------------------------------------

def check_path(myfile):
    """Check for the existence of a file"""

    # expanduser handle tilde
    for f in myfile.split(','):
        if not os.path.isfile(os.path.expanduser(f)):
            print("[ERROR] Can't find the file " + f)
            sys.exit(1)

# -------------------------------------

def check_null(i):
    """Check that value is not None """

    # expects a tuple like this:
    # ('input', args.input)

    if i[1] is None:
        print("[ERROR] " + i[0] + " must not be null")
        sys.exit(1)

# -------------------------------------

def check_file_exists_and_nonzero(myfile):
    """Check for the existence and nonzero-ness of a file"""

    # loop through comma-delimited list of files
    for i in myfile.split(','):
        # if (os.path.isfile(i)):
        if os.path.isfile(os.path.expanduser(i)):
            if os.path.getsize(os.path.expanduser(i)) == 0:
                print(i + " is empty. Exiting")
                sys.exit(1)
        else:
            print("Can't find " + i + ". Exiting")
            sys.exit(1)

# -------------------------------------

def run_cmd(cmd, bool_verbose, bool_getstdout):
    """Run a system (i.e., shell) command"""

    # if verbose, print command
    if bool_verbose:
        print("[command] " + cmd)

    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc.wait() 
    # print return code
    # print(proc.returncode) 
    # print stdout stderr tuple
    # proc.communicate()

    (stdout, stderr) = proc.communicate()

    # if error, print it
    if stderr:
        print("[stderror] " + stderr),

    # return stdout
    if bool_getstdout: 
        return stdout.rstrip()
    else:
        return "0" # note: this must return a str

# -------------------------------------

def run_long_cmd(cmd, bool_verbose, myfile):
    """Run a system command that prints a lot of stuff to stdout or stderr"""

    # https://thraxil.org/users/anders/posts/2008/03/13/Subprocess-Hanging-PIPE-is-your-enemy/

    # "A closer inspection of the subprocess.Popen docs revealed a warning:
    # "Note: The data read is buffered in memory, so do not use this method if the data size is large or unlimited." ...
    # Apparently, that warning actually means "Note: if there's any chance that the data read will be more than a couple pages, 
    # this will deadlock your code." What was happening was that the memory buffer was a fixed size. 
    # When it filled up, it simply stopped letting the child process write to it. 
    # The child would then sit and patiently wait to be able to write the rest of its output.
    # Luckily the solution is fairly simple. Instead of setting stdout and stderr to PIPE, 
    # they need to be given proper file (or unix pipe) objects that will accept a reasonable amount of data."

    # if verbose, print command
    if bool_verbose:
        print("[command] " + cmd)

    with open(myfile, 'w') as f:
        proc = subprocess.Popen(cmd, shell=True, stdout=f, stderr=f)
        proc.wait() 

# -------------------------------------

### SGE-related functions

# -------------------------------------

def getjid(x):
    """Parse out and return SGE job id from string"""
    # Note: the SGE string must look like this: 'Your job 8379811 ("test") has been submitted'
    return x.split('Your job ')[1].split()[0]

# -------------------------------------

### bioinformatics-related functions

# -------------------------------------

def fastajoinlines(infile, outfile, idbase):
    """Put all the sequence portions of a fasta file onto a single line"""

    # a counter for entries
    counter = 1
    # a flag
    flag = 0

    with open(outfile, 'w') as f:
        with open(infile, 'r') as g:
            for line in g:
                line = line.rstrip()
                # if id line
                if line[0] == '>':
                    if flag:
                        f.write('\n')
                        counter += 1
                    f.write('>' + idbase + '_' + str(counter) + '\n')
                    flag = 1
                # if sequence line
                else:
                    f.write(line)
        f.write('\n')

    # return the number of entries in the fasta file
    return counter

# -------------------------------------

def fastasplit(infile, filename, cutoff):
    """Split a fasta file to produce a new file per entry such that seq length > cutoff (assume fastajoinlines)"""

    # a counter for entries
    counter = 0
    id = ''

    with open(infile, 'r') as g:
        for line in g:
            line = line.rstrip()
            if line[0] == '>':
                id = line
            elif len(line) > cutoff:
                counter += 1
                with open(filename + '_' + str(counter) + '.fasta', 'w') as f:
                    f.write(id + '\n')
                    f.write(line + '\n')

    # return counter for contigs above threshold length
    return counter

# -------------------------------------

### Misc

# -------------------------------------

def ConfigSectionMap(Config, section):
    """Parse a configuration file and return dict"""
    # https://wiki.python.org/moin/ConfigParserExamples

    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1

# -------------------------------------

def mytimer(myfunc):
    """Decorator for timing a function"""
    # http://stackoverflow.com/questions/5478351/python-time-measure-function

    def mynewfunc(*args, **kwargs):
        starttime = time.time()
        myfunc(*args, **kwargs)
        print('[step delta t] {} sec'.format(int(time.time() - starttime)))

    return mynewfunc

# -------------------------------------

if __name__ == "__main__":

    pass
