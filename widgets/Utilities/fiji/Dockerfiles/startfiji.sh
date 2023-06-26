#!/bin/bash
set -e
trap 'catch $?' SIGINT SIGTERM EXIT

catch(){
 errCode=$?
 [ -n "$(jobs -pr)" ] && kill $(jobs -pr)
 [ -n "$errCode" ] && exit $errcode
 echo "errcode is $errCode"
 exit 0
}
if [ -n "$installfiji" ]; then
    dest="$installfiji/Fiji.app"
	if [ -d $dest ]; then
	   if [ -n "$overwrite" ]; then
	    rm -r "$dest"
	    echo "installing Fiji.app to $dest"
		cp -r /usr/local/bin/Fiji.app $dest
	   else
	    echo "$dest exists - will not overwrite unless the overwrite option is checked"
	   fi 
	else
	  echo "installing Fiji.app to $dest"
	  cp -r /usr/local/bin/Fiji.app $dest
	fi
fi
if [ -n "$fijidir" ]; then
	cmdString="$fijidir/ImageJ-linux64 "
else
    echo "Running read-only copy of ImageJ/fiji"
    cmdString="/usr/local/bin/Fiji.app/ImageJ-linux64"
fi
if [ -n "$updatesites" ]; then
  #remove "," between each element 
  sites=($(echo $updatesites | sed  's/\"\,\"/ /g' | sed  's/\"//g' | sed 's/\[//g' | sed 's/\]//g'))
  for site in "${sites[@]}"; do
    updateparam=$(echo $site | sed 's/\,/ /g')
    echo "$cmdString --update add-update-site $updateparam"
    exec "$cmdString" --update add-update-site $updateparam &
    wait
  done
fi
if [ -n "$updatefiji" ]; then
	echo "Updating FIJI..."
	exec "$cmdString" --update update &
	wait
fi

if [ -n "$macro" ]; then
	if [ -n "$quitimmediately" ]; then
		timestamp=$(date '+%Y_%m_%d__%H_%M_%S');
		echo 'runMacro(''"'"$macro"'",''"'"$param"'")' > "/tmp/macro.$timestamp"
		echo 'eval("script", "System.exit(0);");' >> "/tmp/macro.$timestamp"
		macro="/tmp/macro.$timestamp"
		cat "$macro"
	fi
	echo "$cmdString  $@ -macro $macro $param "
	$cmdString $@ -macro $macro $param
	exit "$?"
elif [ -n "$script" ]; then
	eval "$cmdString $@ -run $macro $param "
else
	exec "$cmdString" "$@"
fi
