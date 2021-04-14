#/bin/bash
widget=$1
category=$2

[ -z $widget ] && echo "no widget given" && exit 1
[ -z $category ] && echo "no category given" && exit 1

echo "cp -r $widget /widgets/$category/."
cp -r $widget /widgets/$category/.
barewidget=${widget##*/}
echo "ln -sf ../../widgets/$category/$barewidget/$barewidget.py /biodepot/$category/OW$barewidget.py"
ln -sf ../../widgets/$category/$barewidget/$barewidget.py /biodepot/$category/OW$barewidget.py
