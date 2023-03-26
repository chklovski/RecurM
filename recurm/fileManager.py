import os
import errno
import sys
import logging
import shutil
import tempfile
import subprocess


def remove_intermediates(list_of_files):
    for entry in list_of_files:
        os.remove(entry)

def list_assembly_folder(assembly_folder, ext):
    assembly_files = []
    if assembly_folder is not None:
        all_files = os.listdir(assembly_folder)
        for f in all_files:
            if f.endswith(ext):
                assembly = os.path.join(assembly_folder, f)
                if os.stat(assembly).st_size == 0:
                    logging.warning("Skipping assembly {} as it has a size of 0 bytes.".format(f))
                else:
                    assembly_files.append(assembly)

    if not assembly_files:
        logging.error("No assemblies found. Check the extension (-x) used to identify assemblies.")
        sys.exit(1)

    return sorted(assembly_files)

def bash_sort_file(file, column_loc, outdir, threads, numerical):
    if not numerical:
        numerical = ''
    else:
        numerical = 'n'
    try:
        sortedfile = '{}/TMP.hashes'.format(outdir)
        cmd = "sort -t$'\\t' -k{} --parallel {} -{}r {} > {} && mv {} {}" \
            .format(column_loc, threads, numerical, file, sortedfile, sortedfile, file)

        logging.debug(cmd)
        subprocess.call(cmd, shell=True)
        logging.debug('Finished sorting')
    except Exception as e:
        logging.error('An error occured while sorting file {}: {}'.format(file, e))
        sys.exit(1)

def check_empty_dir(input_dir, overwrite=False):
    """Check the the specified directory is empty and create it if necessary."""
    if not os.path.exists(input_dir):
        make_sure_path_exists(input_dir)
    else:
        # check if directory is empty
        files = os.listdir(input_dir)
        if len(files) != 0:
            if overwrite:
                for root, dirs, files in os.walk(input_dir):
                    for f in files:
                        os.unlink(os.path.join(root, f))
                    for d in dirs:
                        shutil.rmtree(os.path.join(root, d))
            else:
                logging.error('Output directory must be empty: ' + input_dir + ' Use --force if you wish to overwrite '
                                                                               'existing directory. \n')
                sys.exit(1)

def check_if_file_exists(inputFile):
    """Check if file exists."""
    if not os.path.exists(inputFile):
        logging.error('File does not exist: ' + inputFile + '\n')
        sys.exit(1)


def make_sure_path_exists(path):
    """Create directory if it does not exist."""
    if not path:
        return

    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            logging.error('Specified path does not exist: ' + path + '\n')
            sys.exit(1)