# -*- coding: utf-8 -*-

"""
@package
@brief
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2015
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

# Standard library imports
# Third party imports
# Local imports

from subprocess import Popen, PIPE


#~~~~~~~COMMAND LINE UTILITIES~~~~~~~#

def run_command(cmd, stdin=None, ret_stderr=False, ret_stdout=True):
    """
    Run a command line in the default shell and return the standard output
    @param  cmd A command line string formated as a string
    @param  stdinput    Facultative parameters to redirect an object to the standard input
    @param  ret_stderr  If True the standard error output will be returned
    @param  ret_stdout  If True the standard output will be returned
    @note If ret_stderr and ret_stdout are True a tuple will be returned and if both are False
    None will be returned
    @return If no standard error return the standard output as a string
    @exception  OSError Raise if a message is return on the standard error output
    @exception  (ValueError,OSError) May be raise by Popen
    """

    # Execute the command line in the default shell
    if stdin:
        proc = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        stdout, stderr = proc.communicate(input=stdin)
    else:
        proc = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        stdout, stderr = proc.communicate()

    if proc.returncode == 1:
        msg = "An error occured during execution of following command :\n"
        msg += "COMMAND : {}\n".format(cmd)
        msg += "STDERR : {}\n".format(stderr)
        raise Exception (msg)

    # Else return data according to user choices is returned
    if ret_stdout and ret_stderr:
        return stdout, stderr
    elif ret_stdout:
        return stdout
    elif ret_stderr:
        return stderr
    else:
        return None

def make_cmd_str(prog_name, opt_dict={}, opt_list=[]):
    """
    Create a Unix like command line string from a
    @param prog_name Name (if added to the system path) or path of the programm
    @param opt_dict Dictionnary of option arguments such as "-t 5". The option flag have to
    be the key (without "-") and the the option value in the dictionnary value. If no value is
    requested after the option flag "None" had to be asigned to the value field.
    @param opt_list List of simple command line arguments
    @exemple make_cmd_str("bwa", {"b":None, t":6, "i":"../idx/seq.fa"}, ["../read1", "../read2"])
    """

    # Start the string by the name of the program
    cmd = "{} ".format(prog_name)

    # Add options arguments from opt_dict
    if opt_dict:
        for key, value in opt_dict.items():
            if value:
                cmd += "-{} {} ".format(key, value)
            else:
                cmd += "-{} ".format(key)

    # Add arguments from opt_list
    if opt_list:
        for value in opt_list:
            cmd += "{} ".format(value)

    return cmd
