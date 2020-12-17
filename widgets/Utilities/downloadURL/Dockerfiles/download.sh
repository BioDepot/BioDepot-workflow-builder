#!/bin/bash
decompress=1
function getFilename(){
	filename=""
	echo "finding filename for url $url"
	tempDir="$(mktemp -d /tmp/XXXXXXXXX)"
	#make a temporary directory without write permissions to force curl to quit after obtaining filename
	chmod -w $tempDir
	filename=$( (cd $tempDir; curl -JLO "$url" |& grep file | grep -o -P '(?<=file ).*(?=:)') )
	rm $tempDir -rf
	if [ -z "$filename" ]; then
	    filename="${url##*/}"
	fi
}

function decompString(){
    dcmd=""
    zipFlag=""
    if [ -n "$decompress" ]
    then
        case $1 in
        *.tar.bz2 | *.tbz2 )
            if [ -n "$concatenateFile" ]; then
                dcmd="| tar xjO >> $concatenateFile"
            else
                dcmd='| tar xj'
            fi
            return
            ;;
        *.tar.gz | *.tgz) 
            if [ -n "$concatenateFile" ]; then
                dcmd="| tar xzO >> $concatenateFile"
            else
                dcmd='| tar xz'
            fi
            return
            ;;  
        *.tar)
            if [ -n "$concatenateFile" ]; then
                dcmd="| tar xO >> $concatenateFile"
            else
                dcmd='| tar x'
            fi
            return
            ;;
        *.gz)
            if [ -n "$concatenateFile" ]; then
                dcmd="| gzip -d >> $concatenateFile"
            else
                outputname=$(basename "$1" .gz)
                dcmd="| gzip -d > $outputname"
            fi
            return
            ;;
        *.bz2)
            if [ -n "$concatenateFile" ]; then
                dcmd="| bzip2 -d >> $concatenateFile"
            else
                outputname=$(basename "$1" .bz2)
                dcmd="| bzip2 -d > $outputname"
            fi
            return
            ;;
        *.zip)
			zipFlag=1
            if [ -n "$concatenateFile" ]; then
                dcmd="&& unzip -p >> $concatenateFile"
            else
                dcmd="&& unzip -o '$filename' && rm '$filename'"
            fi
            return
            ;;
        esac
    fi
    if [ -n "$concatenateFile" ]; then
        dcmd=">> $concatenateFile"
    else
        dcmd="-o $filename"
    fi
    return
}
while [[ $# -gt 0 ]] ; do
    case $1 in
    --decompress)
        decompress=1
        shift 
        ;;
    --directory)
        mkdir -p $2
        cd $2
        shift
        shift 
        ;;
    --concatenateFile)
        concatenateFile=$2
        shift
        shift 
        ;;
    *)   
    urls+=("$1")
    shift # past argument
    ;;
    esac
done


function findGoogleFilename(){
 #check if it fits the
    if [[  $1 == *drive.google.com/file/d/* ]]; then
        fileID=$(echo "$1" | sed -n -e  's/.*drive\.google\.com\/file\/d\///p' | sed  's:/.*::')
    else  
        fileID=$(echo "$1" | sed -n -e 's/.*\?id\=//p')
        fileID=${fileID%%/*}
    fi
    echo fileID is ${fileID}
    filename=$(curl -s -L "$1" | sed -n -e 's/.*<meta property\="og\:title" content\="//p' | sed -n -e 's/">.*//p')
    echo filename is "$filename"
    bash -c "curl -c ./cookie 'https://drive.google.com/uc?export=download&id=${fileID}' &> /dev/null" 
    code=$(cat ./cookie | grep -o 'download_warning.*' | cut -f2)
    if [[ -n "$code" ]]; then
        echo code is "$code"
    else
		echo "no code found"
    fi

}
#empty the concatenateFile if it exists
#do it here instead of in parse loop because in case the directory change comes after the concatenate
if [ -n "$concatenateFile" ]; then
    bash -c "> $concatenateFile"
fi

#loop through the urls
status=0
for url in  "${urls[@]}" ; do
    #if it falls through all code - then there is an error.
	curlret=1
	zipFlag=""
    if [[ $url == *drive.google.com* ]]  
    then
        #find filename and fileID and keep cookie
        findGoogleFilename $url
        decompString "$filename"
        echo "google drive url is $url filename is $filename fileID is $fileID code is $code dcmd is $dcmd"
        if [[ -n "$filename" ]]; then    
            if [ -z $code ]; then
                rm ./cookie
                echo No problem with virus check no verification needed
                decompString "$filename"
                if [[ -n $zipFlag ]]; then
					echo "curl  -L 'https://docs.google.com/uc?export=download&id=${fileID}' -o '$filename' $dcmd"
					bash -c "curl -L 'https://docs.google.com/uc?export=download&id=$fileID' -o '$filename' $dcmd"                
                else
					echo "curl  -L 'https://docs.google.com/uc?export=download&id=${fileID}' $dcmd"
					bash -c "curl -L 'https://docs.google.com/uc?export=download&id=$fileID' $dcmd"
				fi
                curlret=$?
            else
                echo "Verification code to bypass virus scan is $code "
                decompString "$filename"
                if [[ -n $zipFlag ]]; then
					echo "curl  -Lb ./cookie 'https://drive.google.com/uc?export=download&confirm=${code}&id=$fileID' -o '$filename' $dcmd"
					bash -c "curl  -Lb ./cookie 'https://drive.google.com/uc?export=download&confirm=${code}&id=$fileID' -o '$filename' $dcmd"                
                else
					echo "curl  -Lb ./cookie 'https://drive.google.com/uc?export=download&confirm=${code}&id=$fileID' $dcmd"
					bash -c "curl -Lb ./cookie 'https://drive.google.com/uc?export=download&confirm=${code}&id=$fileID' $dcmd"
				fi
                curlret=$?
                rm ./cookie
            fi
        else
            echo "did not download $url - can't find filename - authentication may be required"
            curlret=1
        fi
    else
        echo "url $url is not from google drive"
        getFilename
        if [ -n "$decompress" ] 
        then
            
            echo filename is "$filename"
            decompString "$filename"
            
            if [[ -n $zipFlag ]]; then
				echo "curl -JLO $url $dcmd" 
				bash -c "curl -JLO $url $dcmd"
			else
				echo "curl $url $dcmd" 
				bash -c "curl -L $url $dcmd"
			fi
            curlret=$?     
        else
            if [ -n "$concatenateFile" ]; then
                echo  'curl -L $url >>' "$concatenateFile"
                curl url >> $concatenateFile
                curlret=$?
            else
                echo "curl -JLO $url"
                curl -JLO $url
                curlret=$?
            fi
        fi
    fi
    if [ $curlret -ne 0 ]; then
		status=1
	fi
    shift
done
exit $status
