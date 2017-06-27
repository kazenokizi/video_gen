# before using this code, you need:
# * pymol and avconv

# * get the name of matfile in data folder and the name of the 
#     protein object, change them in this code below
# * in png_make.pml, change pdb name, set_view, community_highlight_new.py directory
# * in community_highlight_new.py, change AA IDs and 2nd structure residues, image size
# * in plot_for_showcase.m, change image size

# if any modifications, search in codes files 'bazinga'

echo -n "Please input start time and end time: "
read t1 t2
echo -n "Do you need all shown in cartoon? "
read response

if [ -z $t1 ] && [ -z $t2 ]; then
    t1=1
    t2=100
fi
if [[ response == 'yes' ]]; then
    default=1
elif [[ response == 'no' ]]; then
    default=0
fi

# make folders
mkdir -p 0partitions 1stability 2diffpartitions 6combination

# time params
tstart=$t1
tend=$t2
# setup the directory
newpath='/home/heng/Documents/0work/2codes/' # bazinga
# newpath='/Users/ac/Documents/MATLAB/codes/'

#1) execute matlab code
cd 1stability
# Matlab script
echo "addpath(genpath('$newpath'));" > matlab_run.m
echo "load('../data/kinase_4AKE_whole_stability_pp.mat')" >> matlab_run.m
echo "plot_for_showcase(kinase_4AKE,2,$tstart,$tend)" >> matlab_run.m
echo "pymol_stability_output(kinase_4AKE,2)" >> matlab_run.m
echo "exit" >> matlab_run.m
# run the script
if [[ `uname -s` == 'Linux' ]]; then
	matlab -nodesktop -nodisplay -nosplash -r matlab_run
elif [[ `uname -s` == 'Darwin' ]]; then
	/Applications/MATLAB_R2014b.app/bin/matlab -nodesktop -nodisplay -nosplash -r matlab_run
fi

# get into the subdir and copy matfile to data
for f in *; do
    if [[ -d $f ]]; then
        cd $f
        for file in *; do
            matfile=`echo $(basename $file)`
        done
        cp * ../../data/
        cd ..
    fi
done

rm matlab_run.m
cd ../codes/
cp png_make.pml png_make_backup.pml
# need to copy for safety
if [[ `uname -s` == 'Linux' ]]; then
    sed -i "s/=0,/=$tstart,/g" png_make.pml 
    sed -i "s/=5,/=$tend,/g" png_make.pml 
    sed -i "s/matfile/$matfile/g" png_make.pml 
elif [[ `uname -s` == 'Darwin' ]]; then
    sed -i .bak "s/=0,/=$tstart,/g" png_make.pml 
    sed -i .bak "s/=5,/=$tend,/g" png_make.pml 
    sed -i .bak "s/matfile/$matfile/g" png_make.pml 
fi

#2) create partitions with different visualisation formats
cd ../2diffpartitions
mkdir -p 1cartoon 2sticks 3both 4alpha

find 1cartoon 2sticks 3both 4alpha -maxdepth 0 -exec cp ../codes/png_make.pml {} \;
find 1cartoon 2sticks 3both 4alpha -maxdepth 0 -exec cp ../data/* {} \;

if [[ default ]]; then
    cd 1cartoon
    pymol -c png_make.pml
    # change the file names
    cd subunit_png
    for FILE in `ls`; do mv $FILE `echo $FILE | sed -e 's:^0*::'`; done
    mv _subunit_*.png 0_subunit_*.png
    cd ../../
else
    k=0
    for f in *; do
        if [[ -d $f ]]; then
    	   cd $f
           k=$((k+1))
    	   # substitue the sig number
            if [[ `uname -s` == 'Linux' ]]; then
    	       sed -i "s/1)/$k)/g" png_make.pml 
            elif [[ `uname -s` == 'Darwin' ]]; then
                sed -i .bak "s/1)/$k)/g" png_make.pml
            fi
            pymol -c png_make.pml
            # change the file names
            cd subunit_png
            for FILE in `ls`; do mv $FILE `echo $FILE | sed -e 's:^0*::'`; done
            mv _subunit_*.png 0_subunit_*.png
            cd ../../
        fi
    done
fi

# copy images
cd ../1stability
x=`ls | wc -l`
numf=$((x-2))
cd ../2diffpartitions
cp 1cartoon/subunit_png/*.png ../0partitions                # bazinga
if [[ default == 0 ]]; then
    cp 2sticks/subunit_png/{0..49}_subunit_*.png ../0partitions 
    cp 3both/subunit_png/{10..29}_*.png ../0partitions 
    cp 4alpha/subunit_png/{30..49}_*.png ../0partitions 
fi

# combination
cd ../6combination
for (( i=0; i<=numf; i++ )); do ../codes/combine.sh $i; done

# call avidemux/avconv and create video
# avidemux
if [[ `which avconv` != '' ]]; then
    avconv -i %d.jpeg -r 25 -b 65536k video.avi
    mv video.avi ..
else
    echo "Need avconv! "
    exit 1
fi
cd ../codes
if [[ `uname -s` == 'Darwin' ]]; then
    rm png_make.pml.bak
fi
rm png_make.pml
mv png_make_backup.pml png_make.pml
cd ..
rm -rf 0partitions 1stability 2diffpartitions 6combination