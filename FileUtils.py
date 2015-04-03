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
from os import access, R_OK, remove, path, mkdir
from gzip import open as gopen
from shutil import copy as shutilCopy
from shutil import Error as shutilError
from sys import stdout
from time import time
from io import open as iopen
from zlib import decompressobj as zlib_decompressobj
from zlib import MAX_WBITS as zlib_MAX_WBITS

# Local imports

#~~~~~~~ PREDICATES ~~~~~~~#

def is_readable_file (fp):
    """ Verify the readability of a file or list of file """
    if not access(fp, R_OK):
        raise IOError ("{} is not a valid file".format(fp))

def is_gziped (fp):
    """ Return True if the file is Gziped else False """
    return fp[-2:].lower() == "gz"

#~~~~~~~ PATH MANIPULATION ~~~~~~~#

def file_basename (path):
    """Return the basename of a file without folder location and extension """
    return path.rpartition('/')[2].partition('.')[0]

def file_extension (path):
    """ Return The extension of a file in lowercase """
    return path.rpartition(".")[2].lower()

def file_name (path):
    """ Return The complete name of a file with the extension but without folder location """
    return path.rpartition("/")[2]

def dir_name (path):
    """ Return the complete path where is located the file without the file name """
    return path.rpartition("/")[0].rpartition("/")[2]

def rm_blank (name, replace=""):
    """ Replace blank spaces in a name by a given character (default = remove)
    Blanks at extremities are always removed and nor replaced """
    return replace.join(name.split())

#~~~~~~~ FILE MANIPULATION ~~~~~~~#

def copyFile(src, dest):
    """
    Copy a single file to a destination file or folder (with error handling/reporting)
    @param src Source file path
    @param dest Path of the folder where to copy the source file
    """
    try:
        shutilCopy(src, dest)
    # eg. src and dest are the same file
    except shutilError as :
        print('Error: %s' % e)
    # eg. source or destination doesn't exist
    except IOError as E:
        print('Error: %s' % E.strerror)

def gzip_file (in_path, out_path=None):
    """
    @param in_path Path of the input uncompressed file
    @param out_path Path of the output compressed file (facultative)
    @exception  OSError Can be raise by open
    """
    # Generate a automatic name if none is given
    if not out_path:
        out_path = in_path +".gz"

    # Try to initialize handle for
    try:
        in_handle = open(in_path, "rb")
        out_handle = gopen(out_path, "wb")
        # Write input file in output file
        print ("Compressing {}".format(in_path))
        out_handle.write (in_handle.read())
        # Close both files
        in_handle.close()
        out_handle.close()
        return path.abspath(out_path)

    except IOError as E:
        print(E)
        if path.isfile (out_path):
            try:
                remove (out_path)
            except OSError:
                print "Can't remove {}".format(out_path)

def gunzip_file (in_path, out_path=None):
    """
    @param in_path Path of the input compressed file
    @param out_path Path of the output uncompressed file (facultative)
    @exception  OSError Can be raise by open
    """
    # Generate a automatic name without .gz extension if none is given
    if not out_path:
        out_path = in_path[0:-3]

    try:
        # Try to initialize handle for
        in_handle = gzip.GzipFile(in_path, 'rb')
        out_handle = open(out_path, "wb")
        # Write input file in output file
        print ("Uncompressing {}".format(in_path))
        out_handle.write (in_handle.read())
        # Close both files
        out_handle.close()
        in_handle.close()
        return path.abspath(out_path)

    except IOError as E:
        print(E)
        if path.isfile (out_path):
            try:
                remove (out_path)
            except OSError:
                print "Can't remove {}".format(out_path)

class gunzip_iterator(object):
    """
    This class decompress gziped files much more quickly than the gzip package from python standard
    library. It was adapted from the code corner blog m.wolf@code-corner.de
    """

    def __init__(self, buffer_size=1024*1024*8):
        self.dobj = zlib_decompressobj(16+zlib_MAX_WBITS) #16+zlib.MAX_WBITS -> zlib can decompress gzip
        self.decomp = []
        self.lines = []
        self.buffer_size = buffer_size

    def open(self, filename):
        self.fhwnd = iopen(filename, "rb")
        self.eof = False

    def close(self):
        self.fhwnd.close()
        self.dobj.flush()
        self.decomp = []

    def decompress(self):
        raw = self.fhwnd.read(self.buffer_size)
        if not raw:
            self.eof = True
            self.decomp.insert(0, self.dobj.flush())

        else:
            self.decomp.insert(0, self.dobj.decompress(raw))

    def readline(self):

        out_str = []

        while True:
            if len(self.lines) > 0:
                return self.lines.pop() + "\n"

            elif len(self.decomp) > 0:
                out = self.decomp.pop()
                arr = out.split("\n")

                if len(arr) == 1:
                    out_str.append(arr[0])

                else:
                    self.decomp.append(arr.pop())
                    arr.reverse()
                    out_str.append(arr.pop())
                    self.lines.extend(arr)

                    out_str.append("\n")
                    return "".join(out_str)

            else:
                if self.eof: break
                self.decompress()

        if len(out_str) > 0:
            return "".join(out_str)


def expand_file (infile, outdir="./"):
    """
    expand file in outdir if the file are gzipped
    Else the file won't be modified and won't be moved from it's current directory
    @param infile Path to a file eventually gzipped
    @param outdir Path of the directory in which to uncompress or copy the files
    @return A path to an uncompressed file
    """
    assert path.isfile(infile), "{} is not a valid file".format(infile)

    # Extract if gziped
    if is_gziped(infile):
        return fgunzip (in_path=infile, out_path=path.join(outdir,file_name(infile)[:-3]))

    # else just return the original file path
    else:
        return infile

def merge_files (inpath_list, outpath="out", compress_output=True, bufsize = 100000):
    """
    Merge a list of text file (gzip or not) in a single file taht can be compress or not
    @param input_list List of files to merge
    @param outpath Destination file
    @param compress_output Gzip the output file. Slower if true
    @param bufline Size of the output file write buffer in line (positive integer)
    @return path of the output merged file
    """
    stime = time()
    # Creating and storing a file for writting output
    outpath = path.abspath(outpath)+".gz" if compress_output else path.abspath(outpath)
    openout = gopen if compress_output else open

    with openout(outpath, "wb") as out_handle:
        # Iterate over files in the input list
        for inpath in inpath_list:

            # Open according to the compression
            openin = gopen if is_gziped(inpath) else open
            with openin (inpath, "rb") as in_handle:
                stdout.write("\t+ {}  ".format(file_name(inpath)))
                stdout.flush()

                # Init a line counter and a text buffer
                lineno = 0
                linebuf = ""

                # Store line in the buffer until the line size is full then flush in out_handle
                for line in in_handle:
                    lineno += 1
                    linebuf += line
                    if lineno % bufsize == 0:
                        out_handle.write(linebuf)
                        linebuf = ""
                    if lineno % 1000000 == 0:
                        stdout.write("*")
                        stdout.flush()

                # Flush the remaining lines in the buffer
                stdout.write("*\n")
                stdout.flush()
                out_handle.write(linebuf)

    print ("{} files merged in {}s\n".format (len(inpath_list), round(time()-stime,3)))
    return outpath

#~~~~~~~ DIRECTORY MANIPULATION ~~~~~~~#

def mkdir(fp):
    """
    Create a directory at the indicated path\n
    Reproduce the ability of UNIX "mkdir -p" command
    (ie if the path already exits no exception will be raised).
    @param  fp path name where the folder should be created
    @exception  OSError Can be raise by os.mkdir
    """
    if path.exists(fp) and path.isdir(fp):
        #print ("'{}' already exist in the current directory".format(fp))
        return fp
    else:
        #print ("Creating '{}' in the current directory".format(fp))
        mkdir(fp)
        return fp
