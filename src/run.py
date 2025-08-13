#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# Python standard library
from __future__ import print_function
from shutil import copytree
import os, re, json, sys, subprocess

# Local imports
from utils import (git_commit_hash,
    join_jsons,
    fatal,
    which,
    exists,
    err)

from . import version as __version__


def init(repo_path, output_path, links=[], required=['workflow', 'resources', 'config'], sym_link = True):
    """Initialize the output directory. If user provides a output
    directory path that already exists on the filesystem as a file
    (small chance of happening but possible), a OSError is raised. If the
    output directory PATH already EXISTS, it will not try to create the directory.
    @param repo_path <str>:
        Path to installation source code and its templates
    @param output_path <str>:
        Pipeline output path, created if it does not exist
    @param links list[<str>]:
        List of files to symlink into output_path
    @param required list[<str>]:
        List of folder to copy over into output_path
    """
    if not exists(output_path):
        # Pipeline output directory
        # does not exist on filesystem
        os.makedirs(output_path)

    elif exists(output_path) and os.path.isfile(output_path):
        # Provided Path for pipeline
        # output directory exists as file
        raise OSError("""\n\tFatal: Failed to create provided pipeline output directory!
        User provided --output PATH already exists on the filesystem as a file.
        Please run {} again with a different --output PATH.
        """.format(sys.argv[0])
        )

    # Copy over templates are other required resources
    copy_safe(source = repo_path, target = output_path, resources = required)

    # Create renamed symlinks for each rawdata
    # file provided as input to the pipeline
    inputs = sym_safe(input_data = links, target = output_path, link = sym_link)

    return inputs


def copy_safe(source, target, resources = []):
    """Private function: Given a list paths it will recursively copy each to the
    target location. If a target path already exists, it will NOT over-write the
    existing paths data.
    @param resources <list[str]>:
        List of paths to copy over to target location
    @params source <str>:
        Add a prefix PATH to each resource
    @param target <str>:
        Target path to copy templates and required resources
    """

    for resource in resources:
        destination = os.path.join(target, resource)
        if not exists(destination):
            # Required resources do not exist
            copytree(os.path.join(source, resource), destination)


def sym_safe(input_data, target, link):
    """Creates re-named symlinks for each FastQ file or cellranger output folder provided
    as input. If a symlink already exists, it will not try to create a new symlink.
    If relative source PATH is provided, it will be converted to an absolute PATH.
    It is currently forcing a link to be created for cellranger output folders, even if the provided
    link parameter is Fals
    @param input_data <list[<str>]>:
        List of input files to symlink to target location
    @param target <str>:
        Target path to copy templates and required resources
    @return input_fastqs list[<str>]:
        List of renamed input FastQs
    """
    input_fastqs = [] # store renamed fastq file names
    for file in input_data:
        if os.path.isdir(file): #Checking if provided file is a directory. If so, assumes it is a cellranger outs folder
            if os.path.exists(os.path.join(file, 'outs')):
                #filename = os.path.join(os.path.basename(os.path.dirname(file)), os.path.basename(file))
                filename = os.path.basename(file)
                link = True
            else:
                raise NameError("""\n\tFatal: Provided input '{}' does not match expected format!
                Cannot determine if existing folder is a cellranger output folder. 
                Please check the folder name and structure before trying again.
                Here is example of expected cellranger output folder structure:
                  input: sampleName     structure: sampleName/outs
                """.format(file, sys.argv[0])
                ) 
        else:
            filename = os.path.basename(file)
        try:
            if not link:
                renamed = rename(filename)
                renamed = os.path.join(target, renamed)
            else:
                renamed = os.path.join(target, filename)
        except NameError as e:
            if not link:
                # Don't care about creating the symlinks
                err('Warning: Skipping over provided input {} file.'.format(filename))
                continue # goto next file
            else:
                raise e

        input_fastqs.append(renamed)
        print(filename, file, renamed)

        if not exists(renamed) and link:
            # Create a symlink if it does not already exist
            # Follow source symlinks to resolve any binding issues
            os.symlink(os.path.abspath(os.path.realpath(file)), renamed, target_is_directory=True)

    return input_fastqs


def rename(filename):
    """Dynamically renames FastQ file to have one of the following extensions: *.R1.fastq.gz, *.R2.fastq.gz
    To automatically rename the fastq files, a few assumptions are made. If the extension of the
    FastQ file cannot be infered, an exception is raised telling the user to fix the filename
    of the fastq files.
    @param filename <str>:
        Original name of file to be renamed
    @return filename <str>:
        A renamed FastQ filename
    """
    # Covers common extensions from SF, SRA, EBI, TCGA, and external sequencing providers
    # key = regex to match string and value = how it will be renamed
    extensions = {
        # Matches: _R[12]_fastq.gz, _R[12].fastq.gz, _R[12]_fq.gz, etc.
        ".R1.f(ast)?q.gz$": ".R1.fastq.gz",
        ".R2.f(ast)?q.gz$": ".R2.fastq.gz",
        # Matches: _R[12]_001_fastq_gz, _R[12].001.fastq.gz, _R[12]_001.fq.gz, etc.
        # Capture lane information as named group
        ".R1.(?P<lane>...).f(ast)?q.gz$": ".R1.fastq.gz",
        ".R2.(?P<lane>...).f(ast)?q.gz$": ".R2.fastq.gz",
        # Matches: _[12].fastq.gz, _[12].fq.gz, _[12]_fastq_gz, etc.
        "_1.f(ast)?q.gz$": ".R1.fastq.gz",
        "_2.f(ast)?q.gz$": ".R2.fastq.gz",
        f"{os.path.sep}outs": ""
    }

    if (filename.endswith('.R1.fastq.gz') or
        filename.endswith('.R2.fastq.gz')):
        # Filename is already in the correct format
        return filename

    converted = False
    for regex, new_ext in extensions.items():
        matched = re.search(regex, filename)
        if matched:
            # regex matches with a pattern in extensions
            converted = True
            filename = re.sub(regex, new_ext, filename)
            break # only rename once

    if not converted:
        raise NameError("""\n\tFatal: Failed to rename provided input '{}'!
        Cannot determine the extension of the user provided input file.
        Please rename the file list above before trying again.
        Here is example of acceptable input file extensions:
          sampleName.R1.fastq.gz      sampleName.R2.fastq.gz
          sampleName_R1_001.fastq.gz  sampleName_R2_001.fastq.gz
          sampleName_1.fastq.gz       sampleName_2.fastq.gz
        Please also check that your input files are gzipped?
        If they are not, please gzip them before proceeding again.
        """.format(filename, sys.argv[0])
        )

    return filename


def setup(sub_args, ifiles, repo_path, output_path):
    """Setup the pipeline for execution and creates config file from templates
    @param sub_args <parser.parse_args() object>:
        Parsed arguments for run sub-command
    @param repo_path <str>:
        Path to installation or source code and its templates
    @param output_path <str>:
        Pipeline output path, created if it does not exist
    @return config <dict>:
         Config dictionary containing metadata to run the pipeline
    """
    # Check for mixed inputs,
    # inputs which are a mixture
    # of FastQ and BAM files
    mixed_inputs(ifiles)

    # Check if inputs are folders
    folder_inputs(ifiles)

    # Resolves PATH to reference file
    # template or a user generated
    # reference genome built via build
    # subcommand
    genome_config = os.path.join(repo_path,'config','genome.json')
    if sub_args.genome.endswith('.json'):
        # Provided a custom reference genome generated by build pipline
        genome_config = os.path.abspath(sub_args.genome)
        sub_args.genome = list(json.load(open(sub_args.genome))['references'].keys())[0]

    required = {
        # Base configuration file
        "base": os.path.join(repo_path,'config','config.json'),
        # Template for project-level information
        "project": os.path.join(repo_path,'config','containers.json'),
        # Template for genomic reference files
        # User provided argument --genome is used to select the template
        "genome": genome_config,
        # Template for tool information
        "tools": os.path.join(repo_path,'config', 'modules.json'),
    }

    # Create the global or master config
    # file for pipeline, config.json
    config = join_jsons(required.values()) # uses templates in config/*.json
    config['project'] = {}
    config = add_user_information(config)
    config = add_rawdata_information(sub_args, config, ifiles)

    # Resolves if an image needs to be pulled
    # from an OCI registry or a local SIF exists
    config = image_cache(sub_args, config, repo_path)

    # Add other runtime info for debugging
    config['project']['version'] = __version__
    config['project']['workpath'] = os.path.abspath(sub_args.output)
    git_hash = git_commit_hash(repo_path)
    config['project']['git_commit_hash'] = git_hash   # Add latest git commit hash
    config['project']['pipeline_path'] = repo_path    # Add path to installation

    # Add all cli options for data provenance
    for opt, v in vars(sub_args).items():
        if opt == 'func':
            # Pass over sub command's handler
            continue
        elif not isinstance(v, (list, dict)):
            # CLI value can be converted to a string
            v = str(v)
        config['options'][opt] = v

    # Run a check for conditional parameters with all metadata full compiled
    check_conditional_parameters(config)

    # Save config to output directory
    with open(os.path.join(output_path, 'config.json'), 'w') as fh:
        json.dump(config, fh, indent = 4, sort_keys = True)

    return config


def unpacked(nested_dict):
    """Generator to recursively retrieves all values in a nested dictionary.
    @param nested_dict dict[<any>]:
        Nested dictionary to unpack
    @yields value in dictionary
    """
    # Iterate over all values of
    # given dictionary
    for value in nested_dict.values():
        # Check if value is of dict type
        if isinstance(value, dict):
            # If value is dict then iterate
            # over all its values recursively
            for v in unpacked(value):
                yield v
        else:
            # If value is not dict type
            # then yield the value
            yield value


def get_fastq_screen_paths(fastq_screen_confs, match = 'DATABASE', file_index = -1):
    """Parses fastq_screen.conf files to get the paths of each fastq_screen database.
    This path contains bowtie2 indices for reference genome to screen against.
    The paths are added as singularity bind points.
    @param fastq_screen_confs list[<str>]:
        Name of fastq_screen config files to parse
    @param match <string>:
        Keyword to indicate a line match [default: 'DATABASE']
    @param file_index <int>:
        Index of line line containing the fastq_screen database path
    @return list[<str>]:
        Returns a list of fastq_screen database paths
    """
    databases = []
    for file in fastq_screen_confs:
        with open(file, 'r') as fh:
            for line in fh:
                if line.startswith(match):
                        db_path = line.strip().split()[file_index]
                        databases.append(db_path)
    return databases


def resolve_additional_bind_paths(search_paths):
    """Finds additional singularity bind paths from a list of random paths. Paths are
    indexed with a compostite key containing the first two directories of an absolute
    file path to avoid issues related to shared names across the /gpfs shared network
    filesystem. For each indexed list of file paths, a common path is found. Assumes
    that the paths provided are absolute paths, the build sub command creates reference
    files with absolute filenames.
    @param search_paths list[<str>]:
        List of absolute file paths to find common bind paths from
    @return common_paths list[<str>]:
        Returns a list of common shared file paths to create additional singularity bind paths
    """
    common_paths = []
    indexed_paths = {}

    for ref in search_paths:
        # Skip over resources with remote URI and
        # skip over strings that are not file PATHS as
        # build command creates absolute resource PATHS
        if ref.lower().startswith('sftp://') or \
        ref.lower().startswith('s3://') or \
        ref.lower().startswith('gs://') or \
        not ref.lower().startswith(os.sep):
            continue

        # Break up path into directory tokens
        path_list = os.path.abspath(ref).split(os.sep)
        try: # Create composite index from first two directories
            # Avoids issues created by shared /gpfs/ PATHS
            index = path_list[1:3]
            index = tuple(index)
        except IndexError:
            index = path_list[1] # ref startswith /
        if index not in indexed_paths:
            indexed_paths[index] = []
        # Create an INDEX to find common PATHS for each root
        # child directory like /scratch or /data. This prevents
        # issues when trying to find the common path betweeen
        # these two different directories (resolves to /)
        indexed_paths[index].append(str(os.sep).join(path_list))

    for index, paths in indexed_paths.items():
        # Find common paths for each path index
        p = os.path.dirname(os.path.commonprefix(paths))
        if p == os.sep:
            # Aviods adding / to bind list when
            # given /tmp or /scratch as input
            p = os.path.commonprefix(paths)
        common_paths.append(p)

    return list(set(common_paths))


def bind(sub_args, config):
    """Resolves bindpaths for singularity/docker images.
    @param sub_args <parser.parse_args() object>:
        Parsed arguments for run sub-command
    @param configfile dict[<any>]:
        Config dictionary generated by setup command.
    @return bindpaths list[<str>]:
        List of singularity/docker bind paths
    """
    bindpaths = []
    for value in unpacked(config):
        if not isinstance(value, str):
            continue
        if exists(value):
            if os.path.isfile(value):
                value = os.path.dirname(value)
            if value not in bindpaths:
                bindpaths.append(value)

    # Bind input file paths, working
    # directory, and other reference
    # genome paths
    rawdata_bind_paths = [os.path.realpath(p) for p in config['project']['datapath'].split(',')]
    working_directory =  os.path.realpath(config['project']['workpath'])
    genome_bind_paths = resolve_additional_bind_paths(bindpaths)
    bindpaths = [working_directory] + rawdata_bind_paths +  genome_bind_paths
    bindpaths = list(set([p for p in bindpaths if p != os.sep]))

    return bindpaths


def mixed_inputs(ifiles):
    """Check if a user has provided a set of input files which contain a
    mixture of FastQ and BAM files. The pipeline does not support processing
    a mix of FastQ and BAM files.
    @params ifiles list[<str>]:
        List containing pipeline input files (renamed symlinks)
    """
    bam_files, fq_files = [], []
    fastqs = False
    bams = False
    for file in ifiles:
        if file.endswith('.R1.fastq.gz') or file.endswith('.R2.fastq.gz'):
            fastqs = True
            fq_files.append(file)
        elif file.endswith('.bam'):
            bams = True
            bam_files.append(file)

    if fastqs and bams:
        # User provided a mix of FastQs and BAMs
        raise TypeError("""\n\tFatal: Detected a mixture of --input data types.
            A mixture of BAM and FastQ files were provided; however, the pipeline
            does NOT support processing a mixture of input FastQ and BAM files.
            Input FastQ Files:
                {}
            Input BAM Files:
                {}
            Please do not run the pipeline with a mixture of FastQ and BAM files.
            This feature is currently not supported within '{}', and it is not
            recommended to process samples in this way either. If this is a priority
            for your project, please run the set of FastQ and BAM files separately
            (in two separate output directories). If you feel like this functionality
            should exist, feel free to open an issue on Github.
            """.format(" ".join(fq_files), " ".join(bam_files), sys.argv[0])
        )

def folder_inputs(ifiles):
    """Check if a user has provided directories as input. 
    @params ifiles list[<str>]:
        List containing pipeline input files (renamed symlinks)
    """
    folder_files, file_files = [], []
    folders = False
    files = False
    for file in ifiles:
        if os.path.isdir(file):
            folders = True
            folder_files.append(file)
        else:
            files = True
            file_files.append(file)

    if folders and files:
        # User provided a mix of folders and files
        raise TypeError("""\n\tFatal: Detected a mixture of --input data types.
            A mixture of folders and files were provided; however, the pipeline
            does NOT support processing a mixture of input FastQ files and 
            cellranger outputs.
            Input Folders:
                {}
            Input Files:
                {}
            Please do not run the pipeline with a mixture of files and folders.
            This feature is currently not supported within '{}'. If you feel like 
            this functionality should exist, feel free to open an issue on Github.
            """.format(" ".join(folder_files), " ".join(file_files), sys.argv[0])
        )
    return(folders)

def add_user_information(config):
    """Adds username and user's home directory to config.
    @params config <dict>:
        Config dictionary containing metadata to run pipeline
    @return config <dict>:
         Updated config dictionary containing user information (username and home directory)
    """
    # Get PATH to user's home directory
    # Method is portable across unix-like
    # OS and Windows
    home = os.path.expanduser("~")

    # Get username from home directory PATH
    username = os.path.split(home)[-1]

    # Update config with home directory and
    # username
    config['project']['userhome'] = home
    config['project']['username'] = username

    return config


def add_sample_metadata(input_files, config, group=None):
    """Adds sample metadata such as sample basename, label, and group information.
    If sample sheet is provided, it will default to using information in that file.
    If no sample sheet is provided, it will only add sample basenames and labels.
    @params input_files list[<str>]:
        List containing pipeline input fastq files
    @params config <dict>:
        Config dictionary containing metadata to run pipeline
    @params group <str>:
        Sample sheet containing basename, group, and label for each sample
    @return config <dict>:
        Updated config with basenames, labels, and groups (if provided)
    """
    import re

    # TODO: Add functionality for basecase
    # when user has samplesheet
    added = []
    config['samples'] = []
    for file in input_files:
        # Split sample name on file extension
        sample = re.split('\.R[12]\.fastq\.gz', os.path.basename(file))[0]
        if sample not in added:
            # Only add PE sample information once
            added.append(sample)
            config['samples'].append(sample)

    return config


def parse_libraries(config, libraries_file, delimeter = ','):
    """Adds sample information from the libraries
    file. The libraries file is a CSV file containing
    information about each library. It contains each
    sample's name, flowcell, demultiplexed name, and
    library type. The relationship between samples
    provided to the --input option and samples listed
    in the libraries file is 1:many. This is because
    you can have a set of FastQ file for different
    single-cell applications, like Gene Expression and
    Antibody Capture, etc.
    @params config <dict>:
        Config dictionary containing metadata to run pipeline
    @params libraries_file <string>:
        10x libraries file containing information about each library
    @return config <dict>:
         Updated config dictionary containing library file information
         (Name, Flowcell)
    """
    def _require(fields, d, lib):
        """Private function that checks to see if all required fields
        are provided in the libraries file. If nan item in fields does
        not  exist in d, then the user forget to add this required field.
        """
        missing = []
        for f in fields:
            try:
                i = d[f]
            except KeyError:
                missing.append(f)
                pass
        if missing:
            fatal(
                "Error: Missing required fields in --libraries {} file!\n \
                └── Please add information for the following field(s): {}".format(
                    lib,
                    ','.join([f.title() for f in missing])
                )
            )

        return

    config['libraries'] = {}
    # Get file extension to determine
    # the appropriate file delimeter
    extension = os.path.splitext(libraries_file)[-1].lower()
    if extension in ['.tsv', '.txt', '.text', '.tab']:
        # file is tab seperated
        delimeter = '\t'
    # Find index of file dynamically,
    # makes it so the order of the
    # columns does not matter
    indices = {}

    with open(libraries_file) as fh:
        try:
            header = next(fh).strip().split(delimeter)
        except StopIteration:
            fatal(
                'Error: --libraries {} cannot be empty!\n \
            └── Please ensure the file is not empty before proceeding again.'.format(libraries_file)
            )
        for i in range(len(header)):
            colname = header[i].strip().lower()
            indices[colname] = i
        _require(['name', 'flowcell', 'sample', 'type'], indices, libraries_file)
        for line in fh:
            linelist = line.strip().split(delimeter)
            # Get indices of fields and
            # parse value in file
            i_name = indices['name']
            name = linelist[i_name]
            i_flowcell = indices['flowcell']
            flowcell = linelist[i_flowcell]
            if name not in config['libraries']:
                config['libraries'][name] = []
            config['libraries'][name].append(flowcell)

    return config


def check_reference_file(reference_file, flag, delimeter = ','):
    """Check reference information from the features
    or cmo reference file. The reference file is a CSV
    file containing information about each feature or
    hashtag/cmo. It contains each feature's id, name,
    read that the barcode is located on, pattern to
    specify the location of the barcode in the read,
    barcode sequence, and feature type.
    @params reference_file <string>:
        10x reference file containing information about each feature
    @params flag <string>:
        Config flag that was used to provide the reference file

    """
    def _require(fields, d, lib, flag):
        """Private function that checks to see if all required fields
        are provided in the reference file. If nan item in fields does
        not  exist in d, then the user forget to add this required field.
        """
        missing = []
        for f in fields:
            try:
                i = d[f]
            except KeyError:
                missing.append(f)
                pass
        if missing:
            fatal(
                f"Error: Missing required fields in --{flag} {{}} file!\n \
                └── Please add information for the following field(s): {{}}".format(
                    lib,
                    ','.join([f.lower() for f in missing])
                )
            )

        return

    # Get file extension to determine
    # the appropriate file delimeter
    extension = os.path.splitext(reference_file)[-1].lower()
    if extension in ['.tsv', '.txt', '.text', '.tab']:
        # file is tab seperated
        delimeter = '\t'
    # Find index of file dynamically,
    # makes it so the order of the
    # columns does not matter
    indices = {}
    with open(reference_file) as fh:
        try:
            header = next(fh).strip().split(delimeter)
        except StopIteration:
            fatal(
                f'Error: --{flag} {{}} cannot be empty!\n \
            └── Please ensure the file is not empty before proceeding again.'.format(reference_file)
            )
        for i in range(len(header)):
            colname = header[i].strip().lower()
            indices[colname] = i
        _require(['id', 'name', 'read', 'pattern', 'sequence', 'feature_type'], indices, reference_file, flag)



def check_rename_file(config, rename_file, delimeter = ','):
    """Check sample information from the rename file.
    The rename file is a CSV file containing information
    about each sample. It contains each sample's name
    and its associated demultiplexed (FastQ) name. The
    relationship between samples provided to the --input
    option and samples listed in the rename file is 1:many.
    This is because sets of FastQ files with different
    names that are from the same sample. It contains each
    unique set of FASTQ name and sample name.
    @params config <dict>:
        Config dictionary containing metadata to run pipeline
    @params rename_file <string>:
        rename file containing information about each sample

    """
    def _require(fields, d, lib):
        """Private function that checks to see if all required fields
        are provided in the reference file. If nan item in fields does
        not  exist in d, then the user forget to add this required field.
        """
        missing = []
        for f in fields:
            try:
                i = d[f]
            except KeyError:
                missing.append(f)
                pass
        if missing:
            fatal(
                f"Error: Missing required fields in --rename {{}} file!\n \
                └── Please add information for the following field(s): {{}}".format(
                    lib,
                    ','.join([f.lower() for f in missing])
                )
            )

        return

    # Get file extension to determine
    # the appropriate file delimeter
    extension = os.path.splitext(rename_file)[-1].lower()
    if extension in ['.tsv', '.txt', '.text', '.tab']:
        # file is tab seperated
        delimeter = '\t'
    # Find index of file dynamically,
    # makes it so the order of the
    # columns does not matter
    indices = {}
    with open(rename_file) as fh:
        try:
            header = next(fh).strip().split(delimeter)
        except StopIteration:
            fatal(
                f'Error: --rename {{}} cannot be empty!\n \
            └── Please ensure the file is not empty before proceeding again.'.format(rename_file)
            )
        for i in range(len(header)):
            colname = header[i].strip().lower()
            indices[colname] = i
        _require(['fastq', 'name'], indices, rename_file)

def check_forcecells_file(config, forcecells_file, delimeter = ','):
    """Check sample information from the force cells file.
    The force cells file is a CSV file containing information
    about the samples that required the force cells flag to be
    used during analysis. It contains the sample's name
    and its associated number of cells. A third column is included
    if the samples are part of the multiplexed library where the
    hashtags will be called by the Cell Ranger analysis.
    @params config <dict>:
        Config dictionary containing metadata to run pipeline
    @params forcecells_file <string>:
        force cells file containing information about each sample
    @params flag <string>:
        Config flag that was used to provide the rename file

    """
    def _require(fields, d, lib):
        """Private function that checks to see if all required fields
        are provided in the reference file. If nan item in fields does
        not  exist in d, then the user forget to add this required field.
        """
        missing = []
        for f in fields:
            try:
                i = d[f]
            except KeyError:
                missing.append(f)
                pass
        if missing:
            fatal(
                f"Error: Missing required fields in --rename {{}} file!\n \
                └── Please add information for the following field(s): {{}}".format(
                    lib,
                    ','.join([f.lower() for f in missing])
                )
            )
        return

    # Get file extension to determine
    # the appropriate file delimeter
    extension = os.path.splitext(forcecells_file)[-1].lower()
    if extension in ['.tsv', '.txt', '.text', '.tab']:
        # file is tab seperated
        delimeter = '\t'
    # Find index of file dynamically,
    # makes it so the order of the
    # columns does not matter
    indices = {}
    with open(forcecells_file) as fh:
        try:
            header = next(fh).strip().split(delimeter)
        except StopIteration:
            fatal(
                f'Error: --forcecells {{}} cannot be empty!\n \
            └── Please ensure the file is not empty before proceeding again.'.format(forcecells_file)
            )
        for i in range(len(header)):
            colname = header[i].strip().lower()
            indices[colname] = i
        _require(['sample', 'cells'], indices, forcecells_file)

def finalcheck(config, flag, delimeter=','):
    """Check the contents of the rename or libraries
    file against input. This function checks to see if
    the input files are not used in the rename/libraries
    file and prints a warning if that occurs. If either
    file lists a FASTQ file or path that is not included
    in the input and throws an error if that is detected.
    @params config <dict>:
        Config dictionary containing metadata to run pipeline
    @params flag <string>:
        Config flag that was used to provide the input file
    """
    filename = config['options'][flag]

    extension = os.path.splitext(filename)[-1].lower()
    if extension in ['.tsv', '.txt', '.text', '.tab']:
        # file is tab seperated
        delimeter = '\t'

    # Find index of file dynamically,
    # makes it so the order of the
    # columns does not matter
    indices = {}

    # Dictionary holding unique contents from files to use for comparisons
    contents = {}
    name_type = {}
    with open(filename) as fh:
        try:
            header = next(fh).strip().split(delimeter)
        except StopIteration:
            fatal(
                f'Error: --rename {{}} cannot be empty!\n \
            └── Please ensure the file is not empty before proceeding again.'.format(filename)
            )
        for i in range(len(header)):
            colname = header[i].strip().lower()
            indices[colname] = i
        for line in fh:
            linelist = line.strip().split(delimeter)
            for i in indices:
                values = contents.get(i, set())
                values.add(linelist[indices[i]])
                contents[i] = values
            value = name_type.get(linelist[indices['name']], set())
            value.add(linelist[indices['type']])
            name_type[linelist[indices['name']]] = value

    # Compiles the sample names and fastq paths from the input (config)
    samples  = set([re.sub("_S[0-9]+_L00[0-9]", "", i) for i in config['samples']])
    fastq_paths = set([os.path.dirname(i) for i in config['options']['input']])

    for index_name in indices:
        comparison = contents[index_name]

        #Check the FASTQ names against the sample (fastq) names provided in the input files
        if index_name in ['sample', 'fastq']:
            if samples != comparison:
                if len(samples-comparison) > 0:
                    print(f"\nWarning: Some FASTQs will be skipped! \nWarning: --{{}} {{}} does not contain values for all provided FASTQ files.\n \
            └── Please note that no sample names have been provided for FASTQ files with the following id(s): {{}} \n \
            These FASTQ files will be skipped when running the pipeline.".format(flag, filename, ','.join(samples-comparison)))
                if len(comparison-samples) > 0:
                    fatal(
                        f'\nError: --{{}} {{}} contains values in FASTQ column that is not in the provided FASTQ files!\n \
            └── Please note that the followed listed FASTQ names are not found in the input files: {{}} '.format(flag, filename, ','.join(comparison-samples))
                    )

        if index_name == 'flowcell':
            # Check to see which values in file are not found in fastq_paths
            missing_file = set([i for i in comparison if sum([i in fastq_path for fastq_path in fastq_paths]) == 0])

            # Check to see which fastq_paths from input are not found in the flowcell values in file
            missing_path = set([fastq_path for fastq_path in fastq_paths if sum([i in fastq_path for i in comparison]) == 0])

            if len(missing_path) > 0:
                print(f"\nWarning: Some FASTQs will be skipped! \nWarning: --{{}} {{}} does not contain values for all provided FASTQ paths.\n \
            └── Please note that no samples contain flowcells that are on the following path(s): \n \
            {{}} \n \
            Any FASTQ files in these paths will be skipped when running the pipeline.".format(flag, filename, ','.join(missing_path)))
            if len(missing_file) > 0:
                fatal(
                    f'\nError: --{{}} {{}} contains values in FASTQ column that is not in the provided FASTQ files!\n \
            └── Please note that the followed listed FASTQ names are not found in the input files: {{}} '.format(flag, filename, ','.join(missing_file))
                )

    names = [name for name in name_type if (len(name_type[name]) <= 1)]
    if len(names) > 0:
        print(f"\nWarning: Some samples only have one feature type associated with them! \nWarning: --{{}} {{}} only contains one feature type for some of the samples.\n \
            └── Please note that only one feature type was provided for the following sample(s): {{}} \n \
            If this is correct, these samples do not need to be run using cellranger multi.".format(flag, filename, ','.join(names)))


def check_conditional_parameters(config):
    """Check the compiled config fictionary to ensure
    that any parameters that are only required for
    certain pipelines are provided if that is the version
    of the pipeline that is being run.
    @params config <dict>:
        Config dictionary containing metadata to run pipeline
    """
    errorMessage = []
    input_folders = folder_inputs(config['options']['input'])
    #Check if cellranger version is provided when required
    if config['options']['pipeline'] in ['gex', 'cite', 'multi'] and config['options']['cellranger'] == '':
        errorMessage += [
            "Error: Version of cellranger to use is required for {} pipeline\n \
            └── Please use the --cellranger flag to select one of the available versions: {}".format(
                config['options']['pipeline'],
                ', '.join(['7.1.0', '7.2.0', '8.0.0', '9.0.0'])
            )
        ]

    #Check if libraries file is provided when required
    if config['options']['pipeline'] in ['cite', 'multi', 'multiome'] and config['options']['libraries'] == 'None' and not input_folders:
        errorMessage += [
            "Error: Libraries file is required for {} pipeline\n \
            └── Please use the --libraries flag to provide the CSV file with the columns: {}".format(
                config['options']['pipeline'],
                ','.join(['Name', 'Flowcell', 'Sample', 'Type'])
            )
        ]

    #Check if features file is provided when required
    if config['options']['pipeline'] in ['cite'] and config['options']['features'] == 'None' and not input_folders:
        errorMessage += [
            "Error: Features file is required for {} pipeline\n \
            └── Please use the --features flag to provide the CSV file with the columns: {}".format(
                config['options']['pipeline'],
                ','.join(['id', 'name', 'sequence', 'feature_type', 'read', 'pattern'])
            )
        ]

    #Check reference
    if config['options']['genome'] in ['hg2024', 'mm2024']:
        if config['options']['pipeline'] in []:
            errorMessage += [
                "Error: The {} reference is not available for the {} pipeline\n \
                └── Please use the --genome flag to select one of the available references: {}".format(
                    config['options']['genome'], config['options']['pipeline'],
                    ', '.join(['hg38', 'mm10'])
                )
            ]


    if len(errorMessage) > 0:
        errorMessage += ["\nAdditional information about flags can be found via cell-seek run --help or at https://openomics.github.io/cell-seek/usage/run/"]
        fatal('\n'.join(errorMessage))


def add_rawdata_information(sub_args, config, ifiles):
    """Adds information about rawdata provided to pipeline.
    Determines whether the dataset is paired-end or single-end and finds the set of all
    rawdata directories (needed for -B option when running singularity). If a user provides
    paired-end data, checks to see if both mates (R1 and R2) are present for each sample.
    @param sub_args <parser.parse_args() object>:
        Parsed arguments for run sub-command
    @params ifiles list[<str>]:
        List containing pipeline input files (renamed symlinks)
    @params config <dict>:
        Config dictionary containing metadata to run pipeline
    @return config <dict>:
         Updated config dictionary containing user information (username and home directory)
    """

    # Determine whether dataset is paired-end
    # or single-end
    # Updates config['project']['nends'] where
    # 1 = single-end, 2 = paired-end, -1 = bams
    convert = {1: 'single-end', 2: 'paired-end', -1: 'bam', -2: 'cellranger'}
    nends = get_nends(ifiles)  # Checks PE data for both mates (R1 and R2)
    config['project']['nends'] = nends
    config['project']['filetype'] = convert[nends]

    # Finds the set of rawdata directories to bind
    rawdata_paths = get_rawdata_bind_paths(input_files = sub_args.input)
    config['project']['datapath'] = ','.join(rawdata_paths)

    # Add each sample's basename
    # from the list of supplied
    # FastQ files provided via
    # the --input option
    config = add_sample_metadata(input_files = ifiles, config = config)
    # Add the list of samples
    # provided in the libraries
    # file, i.e. libraries.csv
    if sub_args.libraries != None:
        libraries = sub_args.libraries
        config = parse_libraries(config = config, libraries_file = libraries)

    # Check to see if submitted
    # features reference file has
    # required columns
    if sub_args.features != None:
        reference = sub_args.features
        check_reference_file(reference_file = reference, flag = "features")

    # Check to see if submitted
    # HTO reference file has
    # required columns
    if sub_args.cmo_reference != None:
        reference = sub_args.cmo_reference
        check_reference_file(reference_file = reference, flag = "cmo_reference")

    if sub_args.rename != None:
        rename = sub_args.rename
        check_rename_file(rename_file = rename, config=config)

    if sub_args.forcecells != None:
        forcecells = sub_args.forcecells
        check_forcecells_file(forcecells_file = forcecells, config=config)

    return config


def image_cache(sub_args, config, repo_path):
    """Adds Docker Image URIs, or SIF paths to config if singularity cache option is provided.
    If singularity cache option is provided and a local SIF does not exist, a warning is
    displayed and the image will be pulled from URI in 'config/containers.json'.
    @param sub_args <parser.parse_args() object>:
        Parsed arguments for run sub-command
    @params config <file>:
        Docker Image config file
    @param repo_path <str>:
        Path to installation or source code and its templates
    @return config <dict>:
         Updated config dictionary containing user information (username and home directory)
    """
    images = os.path.join(repo_path, 'config','containers.json')

    # Read in config for docker image uris
    with open(images, 'r') as fh:
        data = json.load(fh)
    # Check if local sif exists
    for image, uri in data['images'].items():
        if sub_args.sif_cache:
            sif = os.path.join(sub_args.sif_cache, '{}.sif'.format(os.path.basename(uri).replace(':', '_')))
            if not exists(sif):
                # If local sif does not exist on in cache,
                # print warning and default to pulling from
                # URI in config/containers.json
                print('Warning: Local image "{}" does not exist in singularity cache'.format(sif), file=sys.stderr)
            else:
                # Change pointer to image from Registry URI
                # to local SIF
                data['images'][image] = sif

    config.update(data)

    return config


def get_nends(ifiles):
    """Determines whether the dataset is paired-end or single-end.
    If paired-end data, checks to see if both mates (R1 and R2) are present for each sample.
    If single-end, nends is set to 1. Else if paired-end, nends is set to 2.
    @params ifiles list[<str>]:
        List containing pipeline input files (renamed symlinks)
    @return nends_status <int>:
         Integer reflecting nends status: 1 = se, 2 = pe, -1 = bams
    """
    # Determine if dataset contains paired-end data
    paired_end = False
    bam_files = False
    cellranger = False
    nends_status = 1

    for file in ifiles:
        if file.endswith('.bam'):
            bam_files = True
            nends_status = -1
            break
        elif file.endswith('.R2.fastq.gz'):
            paired_end = True
            nends_status = 2
            break # dataset is paired-end
        elif os.path.isdir(file):
            cellranger = True
            nends_status = -2

    # Check to see if both mates (R1 and R2)
    # are present paired-end data
    if paired_end:
        nends = {} # keep count of R1 and R2 for each sample
        for file in ifiles:
            # Split sample name on file extension
            sample = re.split('\.R[12]\.fastq\.gz', os.path.basename(file))[0]
            if sample not in nends:
                nends[sample] = 0

            nends[sample] += 1

        # Check if samples contain both read mates
        missing_mates = [sample for sample, count in nends.items() if count == 1]
        if missing_mates:
            # Missing an R1 or R2 for a provided input sample
            raise NameError("""\n\tFatal: Detected pair-end data but user failed to provide
                both mates (R1 and R2) for the following samples:\n\t\t{}\n
                Please check that the basename for each sample is consistent across mates.
                Here is an example of a consistent basename across mates:
                  consistent_basename.R1.fastq.gz
                  consistent_basename.R2.fastq.gz

                Please do not run the pipeline with a mixture of single-end and paired-end
                samples. This feature is currently not supported within {}, and it is
                not recommended either. If this is a priority for your project, please run
                paired-end samples and single-end samples separately (in two separate output
                directories). If you feel like this functionality should exist, feel free to
                open an issue on Github.
                """.format(missing_mates, sys.argv[0])
            )
    elif not bam_files and not cellranger:
        # Provided only single-end data
        # not supported or recommended
        raise TypeError("""\n\tFatal: Single-end data detected.
            {} does not support single-end data. Calling variants from single-end
            data is not recommended either. If you feel like this functionality should
            exist, feel free to open an issue on Github.
            """.format(sys.argv[0])
        )

    return nends_status


def get_rawdata_bind_paths(input_files):
    """
    Gets rawdata bind paths of user provided fastq files.
    @params input_files list[<str>]:
        List containing user-provided input fastq files
    @return bindpaths <set>:
        Set of rawdata bind paths
    """
    bindpaths = []
    for file in input_files:
        # Get directory of input file
        rawdata_src_path = os.path.dirname(os.path.abspath(os.path.realpath(file)))
        if rawdata_src_path not in bindpaths:
            bindpaths.append(rawdata_src_path)

    return bindpaths


def dryrun(outdir, config='config.json', snakefile=os.path.join('workflow', 'Snakefile')):
    """Dryruns the pipeline to ensure there are no errors prior to runnning.
    @param outdir <str>:
        Pipeline output PATH
    @return dryrun_output <str>:
        Byte string representation of dryrun command
    """
    try:
        dryrun_output = subprocess.check_output([
            'snakemake', '-npr',
            '-s', str(snakefile),
            '--use-singularity',
            '--rerun-incomplete',
            '--cores', str(1),
            '--configfile={}'.format(config)
        ], cwd = outdir,
        stderr=subprocess.STDOUT)
    except OSError as e:
        # Catch: OSError: [Errno 2] No such file or directory
        #  Occurs when command returns a non-zero exit-code
        if e.errno == 2 and not which('snakemake'):
            # Failure caused because snakemake is NOT in $PATH
            err('\n\x1b[6;37;41mError: Are snakemake AND singularity in your $PATH?\x1b[0m')
            fatal('\x1b[6;37;41mPlease check before proceeding again!\x1b[0m')
        else:
            # Failure caused by unknown cause, raise error
            raise e
    except subprocess.CalledProcessError as e:
        print(e, e.output.decode("utf-8"))
        raise(e)

    return dryrun_output


def runner(mode, outdir, alt_cache, logger, additional_bind_paths = None,
    threads=2,  jobname='pl:master', submission_script='run.sh',
    tmp_dir = '/lscratch/$SLURM_JOBID/'):
    """Runs the pipeline via selected executor: local, slurm, uge.
    If 'local' is selected, the pipeline is executed locally on a compute node/instance.
    If 'slurm' is selected, jobs will be submited to the cluster using SLURM job scheduler.
    If 'uge' is selected, jobs will be submited to the cluster using UGE job scheduler.
    Support for additional job schedulers (i.e. PBS, SGE, LSF) may be added in the future.
    @param outdir <str>:
        Pipeline output PATH
    @param mode <str>:
        Execution method or mode:
            local runs serially a compute instance without submitting to the cluster.
            slurm will submit jobs to the cluster using the SLURM job scheduler.
            uge will submit jobs to the cluster using the UGE job scheduler
    @param additional_bind_paths <str>:
        Additional paths to bind to container filesystem (i.e. input file paths)
    @param alt_cache <str>:
        Alternative singularity cache location
    @param logger <file-handle>:
        An open file handle for writing
    @param threads <str>:
        Number of threads to use for local execution method
    @param jobname <str>:
        Name of the master job
    @param submission_script <str>:
        Path to pipeline submission script, see src/run.sh
    @return masterjob <subprocess.Popen() object>:
    """
    # Add additional singularity bind PATHs
    # to mount the local filesystem to the
    # containers filesystem, NOTE: these
    # PATHs must be an absolute PATHs
    outdir = os.path.abspath(outdir)
    # Add any default PATHs to bind to
    # the container's filesystem, like
    # tmp directories, /lscratch
    bindpaths = "{},{}".format(outdir, os.path.dirname(tmp_dir.rstrip('/')))
    # Set ENV variable 'SINGULARITY_CACHEDIR'
    # to output directory
    my_env = {}; my_env.update(os.environ)
    cache = os.path.join(outdir, ".singularity")
    my_env['SINGULARITY_CACHEDIR'] = cache
    if alt_cache:
        # Override the pipeline's default
        # cache location
        my_env['SINGULARITY_CACHEDIR'] = alt_cache
        cache = alt_cache

    if additional_bind_paths:
        # Add Bind PATHs for rawdata directories
        bindpaths = "{},{}".format(additional_bind_paths,bindpaths)

    if not exists(os.path.join(outdir, 'logfiles')):
        # Create directory for logfiles
        os.makedirs(os.path.join(outdir, 'logfiles'))

    # Create .singularity directory for
    # installations of snakemake without
    # setuid which creates a sandbox in
    # the SINGULARITY_CACHEDIR
    if not exists(cache):
        # Create directory for sandbox
        # and image layers
        os.makedirs(cache)

    # Run on compute node or instance
    # without submitting jobs to a scheduler
    if mode == 'local':
        # Run pipeline's main process
        masterjob = subprocess.Popen([
                'snakemake', '-pr', '--rerun-incomplete',
                '--use-singularity',
                '--singularity-args', "'-B {}'".format(bindpaths),
                '--cores', str(threads),
                '--configfile=config.json'
            ], cwd = outdir, stderr=subprocess.STDOUT, stdout=logger, env=my_env)

    # Submitting jobs to cluster via a job scheduler
    elif mode in ['slurm', 'uge']:
        # Run pipeline's main process
        masterjob = subprocess.Popen([
                str(os.path.join(outdir, 'resources', str(submission_script))), mode,
                '-j', jobname, '-b', str(bindpaths),
                '-o', str(outdir), '-c', str(cache),
                '-t', "'{}'".format(tmp_dir)
            ], cwd = outdir, stderr=subprocess.STDOUT, stdout=logger, env=my_env)

    return masterjob
