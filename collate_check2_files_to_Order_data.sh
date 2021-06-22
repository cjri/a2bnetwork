#Script to prepare files for the step --calculation 12; see Instructions.txt
#Here $1 should be the name of a modified checkpoint file.  Suppose that the file Check2_0_0.out is created during --calculation 2.
#Then --calculation 10 would be used to convert this to Check2_0_0_alt.out.  This latter file would be inputted into this script as $1.
#Note here that the following script is an example.  The initial ^11 corresponds to the number of individuals in the transmission network.
#The number of files described in this script would need changing if this was not equal to 11; here the indices go up to CheckX_alt.out and Check2_X_X.out where X is 10, being 11 minus one.  More files would be needed in a larger network and fewer in a smaller network.
#Furthermore, the position of the -1 in each line may not apply in all circumstances.  This is set up, in this case, so that the second individual in the network is the root node, to which no other individual transmits.  The root node is ignored when the code goes on to read the Check2_X_Y_alt.out files; if it is not in the first or second position the -1 could be removed from the lines below.
#If you are not sure where the -1 should go, or if it is needed at all, have a look at your Check2_0_0_alt.out file.  The first number on each line is the number of indviduals, while subsequent numbers correspond to where each individual is infected from, or whether that individual is the root.  If the first or second individual is the root, that needs to be accounted for in the below.
#So in reality you might need to do quite a bit of rewriting the script below.  Sorry.
#Once all of these files have been created, they need to be moved to a directory with the name ../Order_data with respect to where you are running the code.
#If you are wondering why I am making you do all this to run the code, the motivation is that the Check2_0_0.out file can be very big i.e. a few Gb.  In running calculations, the code creates a short list of networks to look at, then looks through the Check2 file to see which of these have plausible orderings, reading in the orderings from those files.  Splitting the Check2 file as below saves a lot of file I/O, if not your time.
#You may also have noticed the split_likelihoods.sh script.  This essentially does the same thing for the previously-calculated likelihoods.  It ensures in a lot of the calculations done in large networks that a calculation is not repeated if the answer is known already.  This involves creating a ../Likelihood_data directory to go alongside ../Order_data...
head -n1 $1 > top.line
grep '^11 0 ' $1 > Check0_alt.out &
grep '^11 1 ' $1 > Check1_alt.out &
grep '^11 2 ' $1 > Check2_alt.out &
grep '^11 3 ' $1 > Check3_alt.out &
grep '^11 4 ' $1 > Check4_alt.out &
grep '^11 5 ' $1 > Check5_alt.out &
grep '^11 6 ' $1 > Check6_alt.out &
grep '^11 7 ' $1 > Check7_alt.out &
grep '^11 8 ' $1 > Check8_alt.out &
grep '^11 9 ' $1 > Check9_alt.out &
grep '^11 10' $1 > Check10_alt.out &
grep '^11 0 -1 0 ' Check0_alt.out > Check2_0_0_alt.out &
grep '^11 0 -1 1 ' Check0_alt.out > Check2_0_1_alt.out &
grep '^11 0 -1 2 ' Check0_alt.out > Check2_0_2_alt.out &
grep '^11 0 -1 3 ' Check0_alt.out > Check2_0_3_alt.out &
grep '^11 0 -1 4 ' Check0_alt.out > Check2_0_4_alt.out &
grep '^11 0 -1 5 ' Check0_alt.out > Check2_0_5_alt.out &
grep '^11 0 -1 6 ' Check0_alt.out > Check2_0_6_alt.out &
grep '^11 0 -1 7 ' Check0_alt.out > Check2_0_7_alt.out &
grep '^11 0 -1 8 ' Check0_alt.out > Check2_0_8_alt.out &
grep '^11 0 -1 9 ' Check0_alt.out > Check2_0_9_alt.out &
grep '^11 0 -1 10' Check0_alt.out > Check2_0_10_alt.out
grep '^11 1 -1 0 ' Check1_alt.out > Check2_1_0_alt.out &
grep '^11 1 -1 1 ' Check1_alt.out > Check2_1_1_alt.out &
grep '^11 1 -1 2 ' Check1_alt.out > Check2_1_2_alt.out &
grep '^11 1 -1 3 ' Check1_alt.out > Check2_1_3_alt.out &
grep '^11 1 -1 4 ' Check1_alt.out > Check2_1_4_alt.out &
grep '^11 1 -1 5 ' Check1_alt.out > Check2_1_5_alt.out &
grep '^11 1 -1 6 ' Check1_alt.out > Check2_1_6_alt.out &
grep '^11 1 -1 7 ' Check1_alt.out > Check2_1_7_alt.out &
grep '^11 1 -1 8 ' Check1_alt.out > Check2_1_8_alt.out &
grep '^11 1 -1 9 ' Check1_alt.out > Check2_1_9_alt.out &
grep '^11 1 -1 10' Check1_alt.out > Check2_1_10_alt.out
grep '^11 2 -1 0 ' Check2_alt.out > Check2_2_0_alt.out &
grep '^11 2 -1 1 ' Check2_alt.out > Check2_2_1_alt.out &
grep '^11 2 -1 2 ' Check2_alt.out > Check2_2_2_alt.out &
grep '^11 2 -1 3 ' Check2_alt.out > Check2_2_3_alt.out &
grep '^11 2 -1 4 ' Check2_alt.out > Check2_2_4_alt.out &
grep '^11 2 -1 5 ' Check2_alt.out > Check2_2_5_alt.out &
grep '^11 2 -1 6 ' Check2_alt.out > Check2_2_6_alt.out &
grep '^11 2 -1 7 ' Check2_alt.out > Check2_2_7_alt.out &
grep '^11 2 -1 8 ' Check2_alt.out > Check2_2_8_alt.out &
grep '^11 2 -1 9 ' Check2_alt.out > Check2_2_9_alt.out &
grep '^11 2 -1 10' Check2_alt.out > Check2_2_10_alt.out
grep '^11 3 -1 0 ' Check3_alt.out > Check2_3_0_alt.out &
grep '^11 3 -1 1 ' Check3_alt.out > Check2_3_1_alt.out &
grep '^11 3 -1 2 ' Check3_alt.out > Check2_3_2_alt.out &
grep '^11 3 -1 3 ' Check3_alt.out > Check2_3_3_alt.out &
grep '^11 3 -1 4 ' Check3_alt.out > Check2_3_4_alt.out &
grep '^11 3 -1 5 ' Check3_alt.out > Check2_3_5_alt.out &
grep '^11 3 -1 6 ' Check3_alt.out > Check2_3_6_alt.out &
grep '^11 3 -1 7 ' Check3_alt.out > Check2_3_7_alt.out &
grep '^11 3 -1 8 ' Check3_alt.out > Check2_3_8_alt.out &
grep '^11 3 -1 9 ' Check3_alt.out > Check2_3_9_alt.out &
grep '^11 3 -1 10' Check3_alt.out > Check2_3_10_alt.out
grep '^11 4 -1 0 ' Check4_alt.out > Check2_4_0_alt.out &
grep '^11 4 -1 1 ' Check4_alt.out > Check2_4_1_alt.out &
grep '^11 4 -1 2 ' Check4_alt.out > Check2_4_2_alt.out &
grep '^11 4 -1 3 ' Check4_alt.out > Check2_4_3_alt.out &
grep '^11 4 -1 4 ' Check4_alt.out > Check2_4_4_alt.out &
grep '^11 4 -1 5 ' Check4_alt.out > Check2_4_5_alt.out &
grep '^11 4 -1 6 ' Check4_alt.out > Check2_4_6_alt.out &
grep '^11 4 -1 7 ' Check4_alt.out > Check2_4_7_alt.out &
grep '^11 4 -1 8 ' Check4_alt.out > Check2_4_8_alt.out &
grep '^11 4 -1 9 ' Check4_alt.out > Check2_4_9_alt.out &
grep '^11 4 -1 10' Check4_alt.out > Check2_4_10_alt.out 
grep '^11 5 -1 0 ' Check5_alt.out > Check2_5_0_alt.out &
grep '^11 5 -1 1 ' Check5_alt.out > Check2_5_1_alt.out &
grep '^11 5 -1 2 ' Check5_alt.out > Check2_5_2_alt.out &
grep '^11 5 -1 3 ' Check5_alt.out > Check2_5_3_alt.out &
grep '^11 5 -1 4 ' Check5_alt.out > Check2_5_4_alt.out &
grep '^11 5 -1 5 ' Check5_alt.out > Check2_5_5_alt.out &
grep '^11 5 -1 6 ' Check5_alt.out > Check2_5_6_alt.out &
grep '^11 5 -1 7 ' Check5_alt.out > Check2_5_7_alt.out &
grep '^11 5 -1 8 ' Check5_alt.out > Check2_5_8_alt.out &
grep '^11 5 -1 9 ' Check5_alt.out > Check2_5_9_alt.out &
grep '^11 5 -1 10' Check5_alt.out > Check2_5_10_alt.out
grep '^11 6 -1 0 ' Check6_alt.out > Check2_6_0_alt.out &
grep '^11 6 -1 1 ' Check6_alt.out > Check2_6_1_alt.out &
grep '^11 6 -1 2 ' Check6_alt.out > Check2_6_2_alt.out &
grep '^11 6 -1 3 ' Check6_alt.out > Check2_6_3_alt.out &
grep '^11 6 -1 4 ' Check6_alt.out > Check2_6_4_alt.out &
grep '^11 6 -1 5 ' Check6_alt.out > Check2_6_5_alt.out &
grep '^11 6 -1 6 ' Check6_alt.out > Check2_6_6_alt.out &
grep '^11 6 -1 7 ' Check6_alt.out > Check2_6_7_alt.out &
grep '^11 6 -1 8 ' Check6_alt.out > Check2_6_8_alt.out &
grep '^11 6 -1 9 ' Check6_alt.out > Check2_6_9_alt.out &
grep '^11 6 -1 10' Check6_alt.out > Check2_6_10_alt.out
grep '^11 7 -1 0 ' Check7_alt.out > Check2_7_0_alt.out &
grep '^11 7 -1 1 ' Check7_alt.out > Check2_7_1_alt.out &
grep '^11 7 -1 2 ' Check7_alt.out > Check2_7_2_alt.out &
grep '^11 7 -1 3 ' Check7_alt.out > Check2_7_3_alt.out &
grep '^11 7 -1 4 ' Check7_alt.out > Check2_7_4_alt.out &
grep '^11 7 -1 5 ' Check7_alt.out > Check2_7_5_alt.out &
grep '^11 7 -1 6 ' Check7_alt.out > Check2_7_6_alt.out &
grep '^11 7 -1 7 ' Check7_alt.out > Check2_7_7_alt.out &
grep '^11 7 -1 8 ' Check7_alt.out > Check2_7_8_alt.out &
grep '^11 7 -1 9 ' Check7_alt.out > Check2_7_9_alt.out &
grep '^11 7 -1 10' Check7_alt.out > Check2_7_10_alt.out
grep '^11 8 -1 0 ' Check8_alt.out > Check2_8_0_alt.out &
grep '^11 8 -1 1 ' Check8_alt.out > Check2_8_1_alt.out &
grep '^11 8 -1 2 ' Check8_alt.out > Check2_8_2_alt.out &
grep '^11 8 -1 3 ' Check8_alt.out > Check2_8_3_alt.out &
grep '^11 8 -1 4 ' Check8_alt.out > Check2_8_4_alt.out &
grep '^11 8 -1 5 ' Check8_alt.out > Check2_8_5_alt.out &
grep '^11 8 -1 6 ' Check8_alt.out > Check2_8_6_alt.out &
grep '^11 8 -1 7 ' Check8_alt.out > Check2_8_7_alt.out &
grep '^11 8 -1 8 ' Check8_alt.out > Check2_8_8_alt.out &
grep '^11 8 -1 9 ' Check8_alt.out > Check2_8_9_alt.out &
grep '^11 8 -1 10' Check8_alt.out > Check2_8_10_alt.out
grep '^11 9 -1 0 ' Check9_alt.out > Check2_9_0_alt.out &
grep '^11 9 -1 1 ' Check9_alt.out > Check2_9_1_alt.out &
grep '^11 9 -1 2 ' Check9_alt.out > Check2_9_2_alt.out &
grep '^11 9 -1 3 ' Check9_alt.out > Check2_9_3_alt.out &
grep '^11 9 -1 4 ' Check9_alt.out > Check2_9_4_alt.out &
grep '^11 9 -1 5 ' Check9_alt.out > Check2_9_5_alt.out &
grep '^11 9 -1 6 ' Check9_alt.out > Check2_9_6_alt.out &
grep '^11 9 -1 7 ' Check9_alt.out > Check2_9_7_alt.out &
grep '^11 9 -1 8 ' Check9_alt.out > Check2_9_8_alt.out &
grep '^11 9 -1 9 ' Check9_alt.out > Check2_9_9_alt.out &
grep '^11 9 -1 10' Check9_alt.out > Check2_9_10_alt.out
grep '^11 10 -1 0 ' Check10_alt.out > Check2_10_0_alt.out &
grep '^11 10 -1 1 ' Check10_alt.out > Check2_10_1_alt.out &
grep '^11 10 -1 2 ' Check10_alt.out > Check2_10_2_alt.out &
grep '^11 10 -1 3 ' Check10_alt.out > Check2_10_3_alt.out &
grep '^11 10 -1 4 ' Check10_alt.out > Check2_10_4_alt.out &
grep '^11 10 -1 5 ' Check10_alt.out > Check2_10_5_alt.out &
grep '^11 10 -1 6 ' Check10_alt.out > Check2_10_6_alt.out &
grep '^11 10 -1 7 ' Check10_alt.out > Check2_10_7_alt.out &
grep '^11 10 -1 8 ' Check10_alt.out > Check2_10_8_alt.out &
grep '^11 10 -1 9 ' Check10_alt.out > Check2_10_9_alt.out &
grep '^11 10 -1 10' Check10_alt.out > Check2_10_10_alt.out
rm Check?_alt.out
rm Check10_alt.out
for c in Check2*; do
  cat top.line $c > temp
  mv temp $c
done

