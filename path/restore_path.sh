#!source
#SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
export PYTHONPATH=$(cat ~/.PYTHONPATH.bk)
export R_LIBS=$(cat ~/.R_LIBS.bk)
echo "PYTHONPATH and R_LIBS was restored"
