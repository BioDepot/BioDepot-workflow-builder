#!/bin/bash
decompString(){
    dcmd=""
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
            ;;  
        *.tar)
            if [ -n "$concatenateFile" ]; then
                dcmd="| tar xO >> $concatenateFile"
            else
                dcmd='| tar x'
            fi
            ;;
        *.gz)
            if [ -n "$concatenateFile" ]; then
                dcmd="| gzip -d >> $concatenateFile"
            else
                dcmd="| gzip -d > ${1%.gz}"
            fi
            return
            ;;
        *.bz2)
            if [ -n "$concatenateFile" ]; then
                dcmd="| bzip2 -d >> $concatenateFile"
            else
                dcmd="| bzip2 -d > ${1%.gz}"
            fi
            return
            ;;
        *.zip)
            if [ -n "$concatenateFile" ]; then
                dcmd="&& unzip -p >> $concatenateFile"
            else
                dcmd="&& unzip -o' $filename '&& rm $filename"
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


function findFilename(){
 #check if it fits the
    if [[  $1 == *drive.google.com/file/d/* ]]; then
        fileID=$(echo "$1" | sed -n -e  's/.*drive\.google\.com\/file\/d\///p' | sed  's:/.*::')
    else  
        fileID=$(echo "$1" | sed -n -e 's/.*\?id\=//p')
        fileID=${fileID%%/*}
    fi
    echo fileID is ${fileID}
    filename=$(curl -s -L "$1" | sed -n -e 's/.*<meta property\="og\:title" content\="//p' | sed -n -e 's/">.*//p')
    curl -c ./cookie -r 0-0 -s -L "https://drive.google.com/uc?export=download&id=${fileID}" >& /dev/null
    code=$(cat ./cookie | grep -o 'download_warning.*' | cut -f2)
    echo filename is "$filename"
    if [[ -n "$code" ]]; then
        echo code is "$code"
    fi

}
#empty the concatenateFile if it exists
#do it here instead of in parse loop because in case the directory change comes after the concatenate
if [ -n "$concatenateFile" ]; then
    bash -c "> $concatenateFile"
fi

#loop through the urls
for url in  "${urls[@]}" ; do

    if [[ $url == *drive.google.com* ]]  
    then
        #find filename and fileID and keep cookie
        findFilename $url
        echo "google drive url is $url filename is $filename fileID is $fileID"
            #get a cookie and check if there is a verification code 
        if [[ -n "$filename" ]]; then    
            if [ -z $code ]; then
                rm ./cookie
                echo No problem with virus check no verification needed
                decompString "$filename"
                echo "curl  -L 'https://drive.google.com/uc?export=download&id=${fileID}' $dcmd"
                bash -c "curl  -L 'https://drive.google.com/uc?export=download&id=$fileID' $dcmd"
            else
                echo Verification code to bypass virus scan is $code 
                decompString "$filename"
                echo "curl  -Lb ./cookie 'https://drive.google.com/uc?export=download&confirm=${code}&id=$fileID' $dcmd"
                bash -c "curl -Lb ./cookie 'https://drive.google.com/uc?export=download&confirm=${code}&id=$fileID' $dcmd"
                rm ./cookie
            fi
        else
            echo "did not download $url - can't find filename - authentication may be required"
        fi
    else
        echo "url $url is not from google drive"
        if [ -n "$decompress" ] 
        then
            filename="${url##*/}"
            echo filename is "$filename"
            decompString "$filename"
            echo "curl $url $dcmd" 
            bash -c "curl $url $dcmd"     
        else
            if [ -n "$concatenateFile" ]; then
                echo  'curl $url >>' "$concatenateFile"
                bash -c "curl url >> $concatenateFile"
            else
                echo "curl -JLO $url"
                bash -c "curl -JLO $url"
            fi
        fi
    fi
    shift
done

