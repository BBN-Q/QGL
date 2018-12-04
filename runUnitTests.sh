#!/bin/bash
echo ""
echo "#---------- $0 (QGL) start..."
echo ""
#----- Simple script to capture unittest preparation and invocation logic.

pwd

export szGit_LFS_Path=`which git-lfs`
echo "szGit_LFS_Path: [${szGit_LFS_Path}]"


# Instantiate the git lfs cached data files...
export CMD="git lfs pull"
echo ""
if [ -z ${szGit_LFS_Path} ]; then
    echo ""
    echo "#--------------------"
    echo "#   W A R N I N G   -- git-lfs extension library unavailable;"
    echo "#                      -- QGL h5 data file processing may fail."
    echo "#                         << [${szGit_LFS_Path}] << 'which git-lfs'"
    echo "#--------------------"
else
    echo "git-lfs available;  retrieving/instantiating cached QGL h5 data files via ~:"
    echo $CMD
    echo ""
    $CMD
fi
echo ""


# Don't forget to set the BBN_MEAS_FILE reference.
export BBN_MEAS_FILE=tests/test_measure.yml

echo ""
echo "#-----Unit-Testing general QGL (BBN_MEAS_FILE=${BBN_MEAS_FILE}) via ~:"
# -f option fails fast
#export CMD="python -m unittest discover . -v -f"
#export CMD="python -m unittest discover . -v"
# Trimm'ed down (non LFS unitest calls) syntax as follows:
export CMD="python -m unittest  tests/test_A*.py tests/test_C*.py tests/test_Scheduler.py tests/test_config.py tests/test_pulse_types.py -v"
echo $CMD
echo ""
$CMD

if [ -z ${szGit_LFS_Path} ]; then
    echo ""
    echo "#-----   W A R N I N G:"
    echo "    Fast-failing unittest modules test_QGL.py, and test_Sequences.py,"
    echo "    currently, due to git lfs data file dependencies;  in-short, the"
    echo "    docker continuumio/miniconda load appears to omit the necessary"
    echo "    git-lfs library installation."
    echo ""
    echo "    Without the git-lfs library extension, h5 data calls from these"
    echo "    modules typically error out with ~:"
    echo "    OSError: Unable to open file (file signature not found)"
    echo ""
    echo "    If/when the git-lfs library becomes available in the docker load, "
    echo "    remove the \"-f\"from the invocation, below:"
fi

echo ""
echo "#----- Testing LFS data file dependent QGL modules via ~:"
# Careful -- in this format (tests.moduleName) DON'T cite the .py suffix
# (it will render an odd error regarding missing 'py' attribute)
export CMD="python -m unittest tests.test_QGL tests.test_Sequences -v -f"
echo $CMD
echo ""
$CMD

echo ""
echo "#---------- $0 (QGL) stop."
echo ""

