#/bin/bash
widget=$1
category=$2

[ -z $widget ] && echo "no widget given" && exit 1
[ -z $category ] && echo "no category given" && exit 1

echo "rm -rf  /widgets/$category/$widget"
rm -rf /widgets/$category/$widget
echo "rm -f /biodepot/$category/OW$widget.py"
rm -f /biodepot/$category/OW$widget.py
