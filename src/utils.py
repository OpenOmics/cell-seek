#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# Python standard library
from __future__ import print_function
from shutil import copytree
import os, sys, hashlib, re, socket
import subprocess, json


def md5sum(filename, first_block_only = False, blocksize = 65536):
    """Gets md5checksum of a file in memory-safe manner.
    The file is read in blocks/chunks defined by the blocksize parameter. This is 
    a safer option to reading the entire file into memory if the file is very large.
    @param filename <str>:
        Input file on local filesystem to find md5 checksum
    @param first_block_only <bool>:
        Calculate md5 checksum of the first block/chunk only
    @param blocksize <int>:
        Blocksize of reading N chunks of data to reduce memory profile
    @return hasher.hexdigest() <str>:
        MD5 checksum of the file's contents
    """
    hasher = hashlib.md5()
    with open(filename, 'rb') as fh:
        buf = fh.read(blocksize)
        if first_block_only:
            # Calculate MD5 of first block or chunck of file.
            # This is a useful heuristic for when potentially 
            # calculating an MD5 checksum of thousand or 
            # millions of file.
            hasher.update(buf)
            return hasher.hexdigest()
        while len(buf) > 0:
            # Calculate MD5 checksum of entire file
            hasher.update(buf)
            buf = fh.read(blocksize)

    return hasher.hexdigest()


def permissions(parser, path, *args, **kwargs):
    """Checks permissions using os.access() to see the user is authorized to access
    a file/directory. Checks for existence, readability, writability and executability via:
    os.F_OK (tests existence), os.R_OK (tests read), os.W_OK (tests write), os.X_OK (tests exec).
    @param parser <argparse.ArgumentParser() object>:
        Argparse parser object
    @param path <str>:
        Name of path to check
    @return path <str>:
        Returns abs path if it exists and permissions are correct
    """
    if not exists(path):
        parser.error("Path '{}' does not exists! Failed to provide vaild input.".format(path))
    if not os.access(path, *args, **kwargs):
        parser.error("Path '{}' exists, but cannot read path due to permissions!".format(path))

    return os.path.abspath(path)


def standard_input(parser, path, *args, **kwargs):
    """Checks for standard input when provided or permissions using permissions().
    @param parser <argparse.ArgumentParser() object>:
        Argparse parser object
    @param path <str>:
        Name of path to check
    @return path <str>:
        If path exists and user can read from location
    """
    # Checks for standard input
    if not sys.stdin.isatty():
        # Standard input provided, set path as an
        # empty string to prevent searching of '-'
        path = ''
        return path

    # Checks for positional arguments as paths
    path = permissions(parser, path, *args, **kwargs)

    return path


def exists(testpath):
    """Checks if file exists on the local filesystem.
    @param parser <argparse.ArgumentParser() object>:
        argparse parser object
    @param testpath <str>:
        Name of file/directory to check
    @return does_exist <boolean>:
        True when file/directory exists, False when file/directory does not exist
    """
    does_exist = True
    if not os.path.exists(testpath):
        does_exist = False # File or directory does not exist on the filesystem

    return does_exist


def ln(files, outdir):
    """Creates symlinks for files to an output directory.
    @param files list[<str>]:
        List of filenames
    @param outdir <str>:
        Destination or output directory to create symlinks
    """
    # Create symlinks for each file in the output directory
    for file in files:
        ln = os.path.join(outdir, os.path.basename(file))
        if not exists(ln):
                os.symlink(os.path.abspath(os.path.realpath(file)), ln)


def which(cmd, path=None):
    """Checks if an executable is in $PATH
    @param cmd <str>:
        Name of executable to check
    @param path <list>:
        Optional list of PATHs to check [default: $PATH]
    @return <boolean>:
        True if exe in PATH, False if not in PATH
    """
    if path is None:
        path = os.environ["PATH"].split(os.pathsep)

    for prefix in path:
        filename = os.path.join(prefix, cmd)
        executable = os.access(filename, os.X_OK)
        is_not_directory = os.path.isfile(filename)
        if executable and is_not_directory:
            return True
    return False


def err(*message, **kwargs):
    """Prints any provided args to standard error.
    kwargs can be provided to modify print functions 
    behavior.
    @param message <any>:
        Values printed to standard error
    @params kwargs <print()>
        Key words to modify print function behavior
    """
    print(*message, file=sys.stderr, **kwargs)



def fatal(*message, **kwargs):
    """Prints any provided args to standard error
    and exits with an exit code of 1.
    @param message <any>:
        Values printed to standard error
    @params kwargs <print()>
        Key words to modify print function behavior
    """
    err(*message, **kwargs)
    sys.exit(1)


def require(cmds, suggestions, path=None):
    """Enforces an executable is in $PATH
    @param cmds list[<str>]:
        List of executable names to check
    @param suggestions list[<str>]:
        Name of module to suggest loading for a given index
        in param cmd.
    @param path list[<str>]]:
        Optional list of PATHs to check [default: $PATH]
    """
    error = False
    for i in range(len(cmds)):
        available = which(cmds[i])
        if not available:
            c = Colors
            error = True
            err("""\n{}{}Fatal: {} is not in $PATH and is required during runtime!{}
            └── Possible solution: please 'module load {}' and run again!""".format(
                c.bg_red, c.white, cmds[i], c.end, suggestions[i])
            )

    if error: fatal()

    return


def is_skyline():
    """Determines if running on a Skyline interactive node (ai-hpcnXXX,
    where XXX is any integer), mirroring the equivalent is_skyline() check
    in src/run.sh. Deliberately false on the head/submit node
    (ai-hpcsubmitN) and any other host, even if SLURM's ClusterName also
    resolves to skyline-prod there.
    @return <boolean>:
        True if on a Skyline interactive node, False otherwise
    """
    host = socket.gethostname().split(".")[0]
    if not re.match(r"^ai-hpcn[0-9]+$", host):
        return False

    try:
        out = subprocess.run(
            ["scontrol", "show", "config"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True, timeout=5,
        ).stdout
    except Exception:
        return False

    m = re.search(r"^ClusterName\s*=\s*(\S+)", out, re.MULTILINE)
    return bool(m) and m.group(1) == "skyline-prod"


def default_tmp_dir():
    """Returns a cluster-appropriate default for --tmp-dir. Skyline's
    scratch space is node-local and already exported via $TMPDIR (e.g.
    /data/scratch/<user>/<hostname>), unlike Biowulf, which requires
    building a path from $SLURM_JOBID (/lscratch/$SLURM_JOBID/). Both
    variables are deliberately left unexpanded here so they resolve
    correctly inside whichever node a given job actually lands on,
    rather than being baked in at CLI-parse time on the submitting node.
    @return <str>:
        Default --tmp-dir value for this cluster
    """
    if is_skyline():
        return "$TMPDIR"
    return "/lscratch/$SLURM_JOBID/"


def container_pull_required(config):
    """Determines whether singularity still needs to pull a Cell Ranger image
    from ghcr.io, or whether image_cache() already resolved the selected
    pipeline/version to a local SIF from --sif-cache. Fails safe: if the
    image for the selected pipeline/version can't be resolved in config,
    this assumes a pull (and GHCR auth) is required.
    @param config <dict>:
        Pipeline config dictionary, after setup() (and image_cache()) have run
    @return <boolean>:
        True if a docker:// pull is still required, False if a local SIF
        was already found in --sif-cache
    """
    pipeline = config.get("options", {}).get("pipeline", "")
    version = config.get("options", {}).get("cellranger", "")
    tool = {
        "atac": "cellranger-atac",
        "multiome": "cellranger-arc",
    }.get(pipeline, "cellranger")

    uri = config.get(tool, {}).get(version)
    if uri is None:
        # Could not resolve the image for this pipeline/version,
        # fail safe and assume a pull (and auth) is required.
        return True

    if uri.startswith("docker://"):
        return True

    if pipeline in ("gex", "cite", "multi"):
        r_uri = config.get("images", {}).get("r_base")
        return r_uri is None or r_uri.startswith("docker://")

    return False


def require_github_auth():
    """Enforces that GitHub Container Registry (ghcr.io) credentials are present
    in the environment. Every pipeline branch pulls a Cell Ranger image from a
    private ghcr.io/openomics/cell-seek package via singularity. This checks for
    the GITHUB_TOKEN and GITHUB_ACTOR (or GITHUB_USER) environment variables,
    which configure_singularity_github_auth() uses to log singularity into
    ghcr.io, prior to runtime.
    """
    c = Colors
    token = os.environ.get("GITHUB_TOKEN", "").strip()
    actor = (os.environ.get("GITHUB_ACTOR") or os.environ.get("GITHUB_USER") or "").strip()

    missing = []
    if not token:
        missing.append("GITHUB_TOKEN")
    if not actor:
        missing.append("GITHUB_ACTOR or GITHUB_USER")

    if missing:
        fatal("""\n{}{}Fatal: missing required GitHub Container Registry authentication!{}
        This pipeline pulls Cell Ranger images from a private ghcr.io/openomics/cell-seek
        package, so the following environment variable(s) must be set: {}
        └── Possible solution: export GITHUB_TOKEN=<personal access token with read:packages scope>
            and export GITHUB_ACTOR=<your github username>, then run again!
        Need a token? See: https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens#creating-a-personal-access-token-classic""".format(
            c.bg_red, c.white, c.end, ", ".join(missing))
        )


def require_singularity_registry_login_support():
    """Enforces a minimum Singularity/SingularityCE version that supports the
    `singularity registry login` subcommand, used to authenticate to ghcr.io.
    That subcommand did not exist prior to SingularityCE 4.0.0: OCI/Docker
    registry auth was previously performed with `singularity remote login`
    (added in Singularity 3.7.0), until the registry-related commands were
    split out of `remote` into their own dedicated `registry` command group
    in SingularityCE 4.0.0. Installs older than 4.0.0 would otherwise fail
    with a raw CalledProcessError from an unrecognized subcommand; fail with
    a clear, actionable message instead.
    @calls fatal() if singularity is older than 4.0.0, or its version cannot
        be determined
    """
    c = Colors
    try:
        out = subprocess.run(
            ["singularity", "--version"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True, timeout=5,
        ).stdout
    except Exception as e:
        out = ""

    m = re.search(r"(\d+)\.(\d+)\.(\d+)", out)
    if not m or int(m.group(1)) < 4:
        fatal("""\n{}{}Fatal: installed singularity version is too old for GitHub Container Registry authentication!{}
        `singularity registry login` (used to authenticate to ghcr.io) requires
        SingularityCE >= 4.0.0; detected version: {}
        └── Possible solution: upgrade to SingularityCE >= 4.0.0, or provide a
            pre-pulled/cached SIF via --sif-cache so no ghcr.io pull is needed.""".format(
            c.bg_red, c.white, c.end, out.strip() or "unknown")
        )


def configure_singularity_github_auth():
    """Logs into ghcr.io with GitHub Container Registry credentials, scoped
    to that registry only. SINGULARITY_DOCKER_USERNAME/PASSWORD are deliberately
    NOT used here: those env vars are applied by singularity to every docker://
    pull regardless of registry, which breaks anonymous pulls of this pipeline's
    other, public DockerHub images (e.g. docker://skchronicles/chicyte). Using
    `singularity registry login` instead stores credentials keyed to ghcr.io
    only (in the user's default registry config), leaving DockerHub pulls
    anonymous. Call require_github_auth() first to ensure the GITHUB_*
    variables are actually present.
    @calls fatal() if the login to ghcr.io fails, or if singularity is too old
        to support `singularity registry login` (see
        require_singularity_registry_login_support())
    """
    token = os.environ.get("GITHUB_TOKEN", "")
    actor = os.environ.get("GITHUB_ACTOR") or os.environ.get("GITHUB_USER") or ""
    if not token or not actor:
        # Nothing to configure; require_github_auth() should have
        # already enforced these are set before this is called.
        return

    require_singularity_registry_login_support()

    try:
        subprocess.run(
            [
                "singularity", "registry", "login",
                "--username", actor,
                "--password-stdin",
                "docker://ghcr.io",
            ],
            input=token,
            universal_newlines=True,
            check=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
    except subprocess.CalledProcessError as e:
        c = Colors
        fatal("""\n{}{}Fatal: failed to login to ghcr.io with the provided GitHub credentials!{}
        {}""".format(c.bg_red, c.white, c.end, e.stderr or e.stdout or str(e)))


def safe_copy(source, target, resources = []):
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


def git_commit_hash(repo_path):
    """Gets the git commit hash of the repo.
    @param repo_path <str>:
        Path to git repo
    @return githash <str>:
        Latest git commit hash
    """
    try:
        githash = subprocess.check_output(['git', 'rev-parse', 'HEAD'], stderr=subprocess.STDOUT, cwd = repo_path).strip().decode('utf-8')
        # Typecast to fix python3 TypeError (Object of type bytes is not JSON serializable)
        # subprocess.check_output() returns a byte string
        githash = str(githash)
    except Exception as e:
        # Github releases are missing the .git directory,
        # meaning you cannot get a commit hash, set the 
        # commit hash to indicate its from a GH release
        githash = 'github_release'
    return githash


def join_jsons(templates):
    """Joins multiple JSON files to into one data structure
    Used to join multiple template JSON files to create a global config dictionary.
    @params templates <list[str]>:
        List of template JSON files to join together
    @return aggregated <dict>:
        Dictionary containing the contents of all the input JSON files
    """
    # Get absolute PATH to templates in git repo
    repo_path = os.path.dirname(os.path.abspath(__file__))
    aggregated = {}

    for file in templates:
        with open(os.path.join(repo_path, file), 'r') as fh:
            aggregated.update(json.load(fh))

    return aggregated


def check_cache(parser, cache, *args, **kwargs):
    """Check if provided SINGULARITY_CACHE is valid. Singularity caches cannot be
    shared across users (and must be owned by the user). Singularity strictly enforces
    0700 user permission on on the cache directory and will return a non-zero exitcode.
    @param parser <argparse.ArgumentParser() object>:
        Argparse parser object
    @param cache <str>:
        Singularity cache directory
    @return cache <str>:
        If singularity cache dir is valid
    """
    if not exists(cache):
        # Cache directory does not exist on filesystem
        os.makedirs(cache)
    elif os.path.isfile(cache):
        # Cache directory exists as file, raise error
        parser.error("""\n\t\x1b[6;37;41mFatal: Failed to provided a valid singularity cache!\x1b[0m
        The provided --singularity-cache already exists on the filesystem as a file.
        Please run {} again with a different --singularity-cache location.
        """.format(sys.argv[0]))
    elif os.path.isdir(cache):
        # Provide cache exists as directory
        # Check that the user owns the child cache directory
        # May revert to os.getuid() if user id is not sufficent
        if exists(os.path.join(cache, 'cache')) and os.stat(os.path.join(cache, 'cache')).st_uid != os.getuid():
                # User does NOT own the cache directory, raise error
                parser.error("""\n\t\x1b[6;37;41mFatal: Failed to provided a valid singularity cache!\x1b[0m
                The provided --singularity-cache already exists on the filesystem with a different owner.
                Singularity strictly enforces that the cache directory is not shared across users.
                Please run {} again with a different --singularity-cache location.
                """.format(sys.argv[0]))

    return cache


def unpacked(nested_dict):
    """Generator to recursively retrieves all values in a nested dictionary.
    @param nested_dict dict[<any>]:
        Nested dictionary to unpack
    @yields value in dictionary 
    """
    # Iterate over all values of given dictionary
    for value in nested_dict.values():
        # Check if value is of dict type
        if isinstance(value, dict):
            # If value is dict then iterate over
            # all its values recursively
            for v in unpacked(value):
                yield v
        else:
            # If value is not dict type then
            # yield the value
            yield value


class Colors():
    """Class encoding for ANSI escape sequeces for styling terminal text.
    Any string that is formatting with these styles must be terminated with
    the escape sequence, i.e. `Colors.end`.
    """
    # Escape sequence
    end = '\33[0m'
    # Formatting options
    bold   = '\33[1m'
    italic = '\33[3m'
    url    = '\33[4m'
    blink  = '\33[5m'
    higlighted = '\33[7m'
    # Text Colors
    black  = '\33[30m'
    red    = '\33[31m'
    green  = '\33[32m'
    yellow = '\33[33m'
    blue   = '\33[34m'
    pink  = '\33[35m'
    cyan  = '\33[96m'
    white = '\33[37m'
    # Background fill colors
    bg_black  = '\33[40m'
    bg_red    = '\33[41m'
    bg_green  = '\33[42m'
    bg_yellow = '\33[43m'
    bg_blue   = '\33[44m'
    bg_pink  = '\33[45m'
    bg_cyan  = '\33[46m'
    bg_white = '\33[47m'


def hashed(l):
    """Returns an MD5 checksum for a list of strings. The list is sorted to
    ensure deterministic results prior to generating the MD5 checksum. This 
    function can be used to generate a batch id from a list of input files.
    It is worth noting that path should be removed prior to calculating the 
    checksum/hash.
    @Input:
        l list[<str>]: List of strings to hash
    @Output:
        h <str>: MD5 checksum of the sorted list of strings
    Example:
        $ echo -e '1\n2\n3' > tmp
        $ md5sum tmp
        # c0710d6b4f15dfa88f600b0e6b624077  tmp
        hashed([1,2,3])   # returns c0710d6b4f15dfa88f600b0e6b624077
    """
    # Sort list to ensure deterministic results
    l = sorted(l)
    # Convert everything to strings
    l = [str(s) for s in l]
    # Calculate an MD5 checksum of results
    h = hashlib.md5()
    # encode method ensure cross-compatiability 
    # across python2 and python3 
    h.update("{}\n".format("\n".join(l)).encode())
    h = h.hexdigest()

    return h


if __name__ == '__main__':
    # Calculate MD5 checksum of entire file 
    print('{}  {}'.format(md5sum(sys.argv[0]), sys.argv[0]))
    # Calcualte MD5 cehcksum of 512 byte chunck of file,
    # which is similar to following unix command: 
    # dd if=utils.py bs=512 count=1 2>/dev/null | md5sum 
    print('{}  {}'.format(md5sum(sys.argv[0], first_block_only = True, blocksize = 512), sys.argv[0]))
