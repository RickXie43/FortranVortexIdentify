#!/bin/bash
export folder_path=$1
parameterspath=source
file_path=$(ls $folder_path | head -n 1)
first_line=$(head -n 1 $folder_path/$file_path)
regex=".* ([0-9]+) ([0-9]+) \".*"
if [[ $first_line =~ $regex ]]; then
	num1=${BASH_REMATCH[1]}
	num2=${BASH_REMATCH[2]}
fi

sed -i "s/#define NX_DEF [0-9]*/#define NX_DEF $num2/" $parameterspath/parameters.F90
sed -i "s/#define NY_DEF [0-9]*/#define NY_DEF $num1/" $parameterspath/parameters.F90

echo "替换完成"
