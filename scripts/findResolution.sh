#/bin/bash
xvfbStr=$(ps -ef | grep screen | grep Xvfb)
[[ $xvfbStr =~ ([0-9]+)x([0-9]+)x([0-9]+) ]] && array=(${BASH_REMATCH//x/ })
echo  ${array[0]} ${array[1]}  