cd ./
run /home/heng/Desktop/video_gen/codes/community_highlight_new.py
load 4ake_H.pdb

hide all
show cartoon

#############################################################
# Change here the position of the camera for the subunit view
#############################################################
set_view (\
     0.417308629,   -0.071697339,   -0.905933022,\
     0.437800139,    0.889436722,    0.131275997,\
     0.796357572,   -0.451400787,    0.402558148,\
    -0.000025660,    0.000025567, -198.414794922,\
    -2.271391630,   -3.156492949,  -16.358736038,\
   159.843780518,  236.987731934,  -20.000000000 )

#############################################################
part2png('matfile','subunit', start=0, end=5, sig=1)

quit
