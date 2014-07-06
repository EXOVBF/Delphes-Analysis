#!/bin/sh

if [ ${1} = "all" ]; then
    names=$(ls cpp | tr " " "\n")
    
    for name_tmp in $names 
    do
	if [$name_tmp = "DelphesTreeReader.cpp" -o
		${name_tmp:${#name_tmp}} = "~" -o
		    ${name_tmp:${#name_tmp}} = "#"]; then
	    continue;
	fi
	echo compiling $name_tmp
	c++ -O2 -lm `root-config --cflags --glibs` -o ${name_tmp%".cpp"} cpp/$name_tmp 
    done

else
    exe_name=${1%".cpp"}
    exe_name=${exe_name#"cpp/"}
    c++ -O2 -lm `root-config --cflags --glibs` -o $exe_name ${1} 
fi