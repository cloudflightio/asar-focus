(
SCRIPTDIR=$(dirname $(readlink -f ${BASH_SOURCE[0]}))
cd "$SCRIPTDIR"
docker build .. -f Dockerfile -t asar-focus \
    --build-arg USER_ID=$(id -u) \
    --build-arg GROUP_ID=$(id -g) \
    --progress=plain 2>&1
) | tee build.log
