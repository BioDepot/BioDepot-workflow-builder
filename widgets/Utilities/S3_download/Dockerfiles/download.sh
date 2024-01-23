#!/bin/bash

function checkFilename(){
    echo "check zero length"
    [[ -n "$filename" ]] || return 1
    echo "check too long"
    [[ $(echo "${#filename}") -lt 255 ]] || return 1
    echo "check character"
    [[ $filename =~ ^[0-9a-zA-Z._-]+$ ]] || return 1
    echo "check first char"
    [[ $(echo $filename | cut -c1-1) =~ ^[0-9a-zA-Z]+$ ]] || return 1
    return 0
}
function getFilename(){
    curlret=0
    unset filename
    echo "finding filename for url $url"
    tempDir="$(mktemp -d /tmp/XXXXXXXXX)"
    #make a temporary directory without write permissions to force curl to quit after obtaining filename
    chmod -w $tempDir
    filename=$(cd $tempDir; su user -c "wget --content-disposition $url |& grep denied | sed 's/.*denied //; s/:.*//'")
    checkFilename && return
    filename=$(su user -c "curl -JLO $url |& grep -m 1 Warning | sed 's/.* file //; s/:.*//'")
    checkFilename && return
    filename=$(basename "$url")
    checkFilename && return
    filename="${url##*/}"
}
function getRequestURL(){
    local content_url="https://drive.google.com/uc?export=download&id=$1"
    local html_content=$(curl -c ./cookies.txt -s -L "$content_url")
    # Extract values using grep and sed
    local id=$(echo "$html_content" | grep -o 'name="id" value="[^"]*' | sed 's/name="id" value="//')
    local export=$(echo "$html_content" | grep -o 'name="export" value="[^"]*' | sed 's/name="export" value="//')
    local confirm=$(echo "$html_content" | grep -o 'name="confirm" value="[^"]*' | sed 's/name="confirm" value="//')
    local uuid=$(echo "$html_content" | grep -o 'name="uuid" value="[^"]*' | sed 's/name="uuid" value="//')

    # Construct the request URL
    request_url="https://drive.usercontent.google.com/download?id=${id}&export=${export}&confirm=${confirm}&uuid=${uuid}"
    
}

function findGoogleFilename(){
    if [[  $1 == *drive.google.com/file/d/* ]]; then
        fileID=$(echo "$1" | sed -n -e  's/.*drive\.google\.com\/file\/d\///p' | sed  's:/.*::')
    else
        fileID=$(echo "$1" | sed -n -e 's/.*\?id\=//p')
        fileID=${fileID%%/*}
    fi
    echo "fileID is ${fileID}"
    filename=$(curl -s -L "$1" | sed -n -e 's/.*<meta property\="og\:title" content\="//p' | sed -n -e 's/">.*//p')
    echo curl -s -L "$1" 
    [ -z "$filename" ] && echo "Unable to find filename - cannot download from google" &&  exit 1
    #first see if the file can be downloaded directly by checking the header when downloading to a non-writable directory
    tempDir="$(mktemp -d /tmp/XXXXXXXXX)"
    #make a temporary directory without write permissions to force curl to quit after obtaining filename
    chmod -w $tempDir
    cookie=$(cd $tempDir; su user bash -c "curl -I -L 'https://docs.google.com/uc?export=download&id=$fileID'" | grep -o -P '(?<=set-cookie: ).*' | sed 's/;.*//')
}

function decompString(){
    unset dcmd
    unset zipFlag
    if [ -n "$decompress" ]; then
        case $1 in
        *.tar.bz2 | *.tbz2 )
            if [ -n "$concatenateFile" ]; then
                dcmd="| tar -xjOf - >> $concatenateFile"
            else
                dcmd='| tar -xjf -'
            fi
            return
            ;;
        *.tar.gz | *.tgz)
            if [ -n "$concatenateFile" ]; then
                dcmd="| tar -xzOf - >> $concatenateFile"
            else
                dcmd='| tar -xzf -'
            fi
            return
            ;;
        *.tar)
            if [ -n "$concatenateFile" ]; then
                dcmd="| tar -xOf - >> $concatenateFile"
            else
                dcmd='| tar -xf -'
            fi
            return
            ;;
        *.gz)
            if [ -n "$concatenateFile" ]; then
                dcmd="| gzip -d >> $concatenateFile"
            else
                local outputname=$(basename "$1" .gz)
                dcmd="| gzip -d > $outputname"
            fi
            return
            ;;
        *.bz2)
            if [ -n "$concatenateFile" ]; then
                dcmd="| bzip2 -d >> $concatenateFile"
            else
                local outputname=$(basename "$1" .bz2)
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
}

while [[ $# -gt 0 ]] ; do
    case $1 in
    --decompress)
        decompress=1
        ;;
    --directory)
        mkdir -p $2
        cd $2
        shift
        ;;
    --concatenateFile)
        concatenateFile=$2
        shift
        ;;
    --noClobber)
        noClobber=1
        ;;
    *)
        urls+=("$1")
        ;;
    esac
    shift
done

#empty the concatenateFile if it exists
#do it here instead of in parse loop because in case the directory change comes after the concatenate
if [ -n "$concatenateFile" ]; then
    bash -c "> $concatenateFile"
fi

#loop through the urls
status=0
for url in "${urls[@]}" ; do
    #if it falls through all code - then there is an error.
    curlret=1
    if [[ $url == *drive.google.com* ]]; then
        #find filename and fileID and keep cookie
        findGoogleFilename $url
        decompString "$filename"
        echo "google drive url is $url filename is $filename fileID is $fileID dcmd is $dcmd"
        if [[ -n "$filename" ]]; then
            if [ -n "$noClobber" ] && [ -f "$filename" ]; then
                echo "File $filename is already present, skipping download"
                shift
                continue
            fi
            if [ -z "$cookie" ]; then
                echo No problem with virus check no verification needed
                cmd="curl -L 'https://docs.google.com/uc?export=download&id=${fileID}' "
            else
                echo "We need to pass the virus check"
                getRequestURL $fileID
                echo "request url is $request_url"
                cmd="curl -Lb ./cookies.txt '${request_url}' "
            fi
            if [[ -n $zipFlag ]]; then
                cmd+="-o '$filename' "
            fi
            cmd+="$dcmd"
            echo "$cmd"
            bash -c "$cmd"
            curlret=$?
            rm -f ./cookie
        else
            echo "did not download $url - can't find filename - authentication may be required"
        fi
    else
        echo "url $url is not from google drive"
        getFilename
        if [ -n "$noClobber" ] && [ -f "$filename" ]; then
            echo "File $filename is already present, skipping download"
            shift
            continue
        fi
        echo "$filename"
        if [ -n "$decompress" ]; then
            # check for a log that should contain all extracted objects
            if [ -s "$filename.log" ] && [ -n "$noClobber" ]; then
                skipDownload=true
                while read f; do
                    if [ ! -e "$f" ]; then
                        # if we are here then one of the extracted objects does not exist
                        skipDownload=false
                        break
                    fi
                done < $filename.log
                if $skipDownload; then
                    shift
                    continue
                fi
            fi
            decompString "$filename"
            # make a temp directory to store file content then move to permanent location after
            tmpdir=$(mktemp -d -p $PWD)
            pushd $tmpdir > /dev/null
            if [[ -n $zipFlag ]]; then
                cmd="curl -JLO $url $dcmd"
            else
                cmd="curl -L $url $dcmd"
            fi
            echo "$cmd"
            bash -c "$cmd"
            curlret=$?
            if [ $curlret -eq 0 ]; then
                # store a log file to prevent script from downloading again
                find -not -name . > ../$filename.log
                # set dot glob to move hidden files and directories
                shopt -s dotglob
                mv * ../
                # unset dot glob to avoid trouble
                shopt -u dotglob
            fi
            popd > /dev/null
            rmdir $tmpdir
        else
            if [ -n "$concatenateFile" ]; then
                cmd="curl $url >> $concatenateFile"
            else
                cmd="curl -JLO $url"
            fi
            echo "$cmd"
            $cmd
            #curlret=$?
        fi
    fi   
    if [ $curlret -ne 0 ]; then
        status=1
    fi
    shift
done

exit $status
