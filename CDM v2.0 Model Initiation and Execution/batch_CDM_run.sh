 #!/bin/bash         

n=3 # number of cores

echo $PWD
read -r -p "Is this the correct working directory? [y/N] " response
if [[ "$response" =~ ^([yY][eE][sS]|[yY])+$ ]]
then
#    echo "Printing files..."
	sub_dir="/model_iter*/"
	program_file_name="DUNE"
	parameter_file_name="params.par"
	
	prog_fullpath=$PWD$sub_dir$program_file_name
	param_fullpath=$PWD$sub_dir$parameter_file_name
#	echo $fullpath

	parallel -j $n -k --eta --link echo {}\; {} ::: $prog_fullpath ::: $param_fullpath
	
	echo "All processes complete"
else
    echo "You answered no. Exiting..."
fi
