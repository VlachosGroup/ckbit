#!/bin/bash
#
# pystan-install.sh
#
# Download the latest PyStan source, switch to the CVODES branch, and build and install into
# the current Python environment.
#

BUILD_ROOT="./"
VERBOSE=""
NO_CLONE=0

function doCmd() {
    local   output error_str forced_rc rc cmd_rc verbose no_exit
    local   no_output_capture success_str
    
    forced_rc=0
    verbose=0
    no_exit=0
    no_output_capture=0
    while [ $# -gt 0 ]; do
        case "$1" in
            --no-output-capture|-O)
                no_output_capture=1
                ;;
            --no-exit|-E)
                no_exit=1
                ;;
            --success-string|-s)
                shift
                if [ $# -eq 0 ]; then
                    echo "ERROR:  no argument to --success-string/-s in doCmd()"
                    exit 22
                fi
                success_str="$1"
                ;;
            --error-string|-e)
                shift
                if [ $# -eq 0 ]; then
                    echo "ERROR:  no argument to --error-string/-e in doCmd()"
                    exit 22
                fi
                error_str="$1"
                ;;
            --rc|-r)
                shift
                if [ $# -eq 0 ]; then
                    echo "ERROR:  no argument to --rc/-r in doCmd()"
                    exit 22
                fi
                forced_rc=1
                rc="$1"
                ;;
            --verbose|-v)
                verbose=1
                ;;
            *)
                break
                ;;
        esac
        shift
    done
    if [ $# -gt 0 ]; then
        if [ $no_output_capture -eq 0 ]; then
            output="$("$@" 2>&1)"
            cmd_rc=$?
        else
            "$@"
            cmd_rc=$?
        fi
        if [ $cmd_rc -ne 0 ]; then
            if [ -n "$error_str" ]; then
                echo "ERROR:  $error_str"
                if [ $verbose -ne 0 ]; then
                    printf "        %s" "$output"
                fi
            fi
            if [ $no_exit -eq 0 ]; then
                if [ $forced_rc -ne 0 ]; then
                    exit $rc
                fi
                exit $cmd_rc
            fi
            if [ $forced_rc -ne 0 ]; then
                return $rc
            fi
            return $cmd_rc
        fi
    fi
    if [ -n "$success_str" ]; then
        echo "$success_str"
    fi
    return 0
}

function usage() {
    cat <<EOT
usage:

    $1 {--verbose|-v} {--help|-h} {--build-root|-p <path>}
        {--no-clone|-C}
EOT
    exit 0
}

while [ $# -gt 0 ]; do
    case "$1" in
        --help|-h)
            usage "$0"
            ;;
        --verbose|-v)
            VERBOSE="-v"
            ;;
        --no-clone|-C)
            NO_CLONE=1
            ;;
        --build-root|-p)
            shift
            if [ $# -eq 0 ]; then
                echo "ERROR:  no argument provided with --build-root/-p"
                exit 22
            fi
            BUILD_ROOT="$1"
            ;;
        *)
            break
            ;;
    esac
    shift
done

if [ ! -d "$BUILD_ROOT" ]; then
    echo "ERROR:  build root directory does not exist:  $BUILD_ROOT"
    exit 1
fi
BUILD_DIR="${BUILD_ROOT}/pystan-src"
CURRENT_DIR="$(pwd)"

if [ $NO_CLONE -eq 0 ]; then
    doCmd $VERBOSE -e "failed to clone PyStan github repository" -s "OK - cloned git repository" \
            git clone --recursive https://github.com/stan-dev/pystan "$BUILD_DIR"
fi
pushd "$BUILD_DIR" >/dev/null 2>&1
if [ $? -ne 0 ]; then
    echo "ERROR:  unable to change to PyStan source directory: $BUILD_DIR"
    exit 1
fi
doCmd $VERBOSE -e "unable to checkout CVODES branch" -s "OK - checked-out cvodes branch" \
        git checkout cvodes
doCmd $VERBOSE -e "unable to update submodules in repository" -s "OK - updated submodules in repository"  \
        git submodule update --recursive
doCmd $VERBOSE -e "failed while cleaning unneeded components from source repository" -s "OK - removed unneeded components from repository"  \
        python -c "\"import os, shutil; [shutil.rmtree(r'\\?\{}'.format(os.path.abspath(dirname)), ignore_errors=True) for dirname in [dirname for dirname, *_ in os.walk('pystan/stan') if any(dirname.endswith(ends) for ends in ['doc', 'test'])]]\""
doCmd $VERBOSE -e "failed while building PyStan+CVODES module" -s "OK - PyStan+CVODES module built"  \
        python setup.py build
doCmd $VERBOSE -e "failed while installing PyStan+CVODES module" -s "OK - PyStan+CVODES module installed"  \
        python setup.py install
popd >/dev/null 2>&1
