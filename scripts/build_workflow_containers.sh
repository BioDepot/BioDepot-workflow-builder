#!/bin/bash

docker_tag_exists() {
    curl --silent -f -lSL https://index.docker.io/v1/repositories/$1/tags/$2 > /dev/null
}

docker_build(){
 dockerDir=$workflow/widgets/$workflow/$widget/Dockerfiles
 currentDir=$(pwd)
 if test -f "$dockerDir/build.sh"; then
	if [ -z "$SCRIPT_ONLY" ]; then
		echo "Building widget container with script $dockerDir/build.sh $1"
		cd $dockerDir && build.sh $1  
		cd $currentDir
	else
	   scriptLine="cd $dockerDir && build.sh $1; cd $currentDir"
	fi
 else
	if [ -z "$SCRIPT_ONLY" ]; then
		echo "Building widget container from provided Dockerfile"
		cd $dockerDir && docker build -t $1 . 
		cd $currentDir
	else
	   scriptLine="cd $dockerDir && docker build -t $1 . && echo cd $currentDir"		
	fi
 fi
}

getContainer(){
 imageName=$(jq -r '.docker_image_name' $workflow/widgets/$workflow/$widget/$widget.json)
 tag=$(jq -r '.docker_image_tag' $workflow/widgets/$workflow/$widget/$widget.json) 
 containerString=$imageName:$tag
 if [[ -v "containerStringSeen[$containerString]" ]]; then
	return
 fi
 containerStringSeen["$containerString"]=$widget
 noLocalContainer=""
 $(docker inspect --type=image $containerString &> /dev/null) || noLocalContainer=True
 if [ -z "$SCRIPT_ONLY" ] && [ -z "$noLocalContainer" ]  && [ -z "$FORCE_BUILD" ]; then
	echo "$containerString exists"
	return
 fi
 if [ -z "$SCRIPT_ONLY" ]; then
	attempts=0
	while [ -z "$USE_DOCKERFILE" ] && [ $attempts -lt $maxAttempts ]; do
		attempts=$((attempts+1))
		echo "Attempt $attempts to pull $containerString"
		docker pull $containerString &> /dev/null && echo "Successfully pulled $containerString" && return || echo "pull failed"
	done
	docker_build $containerString
 else
	docker_build $containerString
	if [ -z "$USE_DOCKERFILE" ]; then
	   if [ -z "$scriptLine" ]; then
			echo 'echo "Pulling '"$containerString"'"'
			echo "docker pull $containerString &> /dev/null "
	   else 
	        echo  'echo "Pulling '"$containerString"' if possible - otherwise build from Dockerfile"'
			echo "docker pull $containerString &> /dev/null || ( $scriptLine )"
	   fi
	else
	   if [ -z "$scriptLine" ]; then
			echo 'echo "missing Dockerfile - pulling '"$containerString"' instead"'
			echo "docker pull $containerString &> /dev/null "
	   else
			echo 'echo "Building container '"$containerString"' from Dockerfile"'
			echo "$scriptLine"
	   fi	 
	fi
 fi
}

declare -A containerStringSeen
workflow=$1
maxAttempts=$2

if [ -z "$maxAttempts" ]; then
	maxAttempts=3
fi
if [ -z "$SCRIPT_ONLY" ]; then 
	echo "Working on $workflow"
else
	echo "#!/bin/bash"
	echo "#build container script for $workflow"
fi
widgetList=($( cat $workflow/$workflow.ows | grep -o -P '(?<= name=").*?(?=")' | uniq ))
for widget in "${widgetList[@]}"; do
 scriptLine=""
 if [ -z "$SCRIPT_ONLY" ]; then
	echo "Checking $widget"
 fi
 getContainer
done
