#!source
SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
export PYTHONPATH=$(cat ${SCRIPTPATH}/path/PYTHONPATH.bk)
export R_LIBS=$(cat ${SCRIPTPATH}/path/R_LIBS.bk)
echo "PYTHONPATH and R_LIBS was restored"
