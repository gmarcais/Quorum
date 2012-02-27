#! /bin/sh

set -x # Show what is being done
set -e # Exit in case of failure


# Paths
JF=jellyfish
COMB=combine_jf_dbs
QUORUM=quorum

# Parameters
JF_SIZE=
NUM_THREADS=16
MIN_Q_CHAR=33
MIN_QUALITY=5
WINDOW=10
MAX_ERR_PER_WINDOW=3

$JF count -t $NUM_THREADS -m 24 -s $JF_SIZE -C -r -p 126 -o pe_trim --quality-start $MIN_Q_CHAR --min-quality $MIN_QUALITY "$@"
[ -f pe_trim_1 ] && (echo >&2 "Increase JF_SIZE parameter"; false)
$JF count -t $NUM_THREADS -m 24 -s $JF_SIZE -C -r -p 126 -o pe_all "$@"
[ -f pe_all_1 ] && (echo >&2 "Increase JF_SIZE parameter"; false)

$COMB -m 1 pe_trim_0 pe_all_0 -o combined
$QUORUM -d combined_0 -c 2 -C -m 1 -s 1 -g 1 -a 3 -t $NUM_THREADS -w $WINDOW -e $MAX_ERR_PER_WINDOW "$@"
