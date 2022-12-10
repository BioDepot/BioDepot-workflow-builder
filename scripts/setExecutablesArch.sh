#!/bin/bash

([ -z "$2" ] || [ $2 = "amd64" ] || [ $2 = "linux_amd64" ] || [ $2 = "x64" ]) && suffix="linux_x64"
[ -n "$2" ] && ([ $2 = "arm64" ] || [ $2 = "linux_arm64" ]) && suffix="linux_arm64"

readarray -d '' files < <(find $1 -name "*-$suffix" -print0)
for file in "${files[@]}"; do
    echo "making symbolic link for $file"
    path=${file%/*}
    $(ln -fs $file $path/${path##*/})
done