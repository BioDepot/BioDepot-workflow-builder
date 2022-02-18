#!/bin/bash

orderFile=drawers.txt

mkdir -p /biodepot
setup=/biodepot/setup.py
createLinks () {
	workflow=/biodepot/$1
	widgets=($(cd $workflow/widgets && find $1 -mindepth 1 -maxdepth 1 -type d ! -name icon))
	for widget in "${widgets[@]}"; do
		widgetName=$(basename $widget)
		cd $workflow && ln -sf widgets/$widget/$widgetName.py OW$widgetName.py
	done
}
createInit () {
	initFile=/biodepot/$1/__init__.py
	iconPath=($(cd /biodepot/$1/widgets/$1 && ls icon/*))
	echo "import sysconfig" > $initFile
	echo 'ICON = ''"'"$iconPath"'"' >> $initFile
	echo 'BACKGROUND = ''"'light-purple'"'	>> $initFile
}
processLine () {
 #expect 1 or 2 whitespace delimited parameters
 eval "parms=($1)"
 dirName="${parms[0]}"
 drawName="${parms[1]}" 
 if [ -d "/biodepot.temp/$dirName" ]; then
	dirPath="/biodepot.temp/$dirName"
	widgetPath="/widgets.temp/$dirName"
	cp -r  $dirPath /biodepot/$dirName || exit 1
	cp -r  $widgetPath /widgets/$dirName || exit 1
 elif [ -d "/extras/$dirName" ]; then
	dirPath="/extras/$dirName"
	cp -r  $dirPath /biodepot/$dirName || exit 1
	cp -r  /biodepot/$dirName/widgets/$dirName/icon  /biodepot/$dirName/.
	createLinks $dirName
	createInit $dirName
 fi

 [ -n "$drawName" ]  || drawName="${parms[0]}" >> $setup
 echo "setup(" >> $setup
 echo '    name=''"'"$drawName"'",' >> $setup
 echo '    packages=''["'"$dirName"'"],' >> $setup
 echo '    package_data={''"'"$dirName"'": [''"icons/*.svg''"]},' >> $setup
 echo '    entry_points={''"'orange.widgets'": ''"'"$drawName"' = '"$dirName"'"}' >> $setup
 echo ')' >> $setup
 pip3 install -e /biodepot
}
echo "from setuptools import setup" > $setup

while IFS= read -r LINE; do
    processLine "$LINE"
done < $orderFile
