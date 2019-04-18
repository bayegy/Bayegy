#!source
#SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
echo $PYTHONPATH > ~/.PYTHONPATH.bk
echo $R_LIBS > ~/.R_LIBS.bk
export PYTHONPATH=
export R_LIBS=
echo "PYTHONPATH and R_LIBS was cleaned"
