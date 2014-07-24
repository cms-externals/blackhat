#!/usr/bin/perl -s

#option to use for primitive amplitude tests
#for this case the argument list are the names 
#of parent data files, not the actual executables
#have to be passed to print_parents executable

my $do_primitive_tests = 0;
if($primitive){
    $do_primitive_tests = $primitive;
}

my $do_long_primitive_tests = 0;
if($longprimitive){
    $do_long_primitive_tests = $longprimitive;
}


use strict;
use Term::ANSIColor;
use warnings;


    print "**********************************************\n" ;
print "Starting the tests ... \n" ;

my $nbr_tests = @ARGV;
my $nbr_failed_tests = 0;
my $nbr_improved_tests = 0;
my $nbr_worse_tests = 0;
my $parent_data_files = "../../datafiles/parents/";
my $target_data_files = "../../test/primitive_targets/";


foreach ( @ARGV ) {
    my $nbr_errors = 0;
    #target file name
    #my $target_file = substr $_ ,0 , index($_,'.dat');
    #$target_file=$target_data_files.$target_file.'_target.tex';
    if($do_primitive_tests){
       #printf( "primitive\n");
        printf( "check_parent_diagrams\t $_\t ");
        $nbr_errors = system ("./check_parent_diagrams -infile $parent_data_files$_ -dir $target_data_files -test 1> $_.log 2>&1 ") >> 8;
        #print ("./check_parent_diagrams -infile $parent_data_files$_ -dir $target_data_files -test 1> $_.log 2>&1 ");
    }
    elsif($do_long_primitive_tests){
        #printf( "do long primitive\n" );
         printf( "check_parent_diagrams\t -cutonly $_\t ");
        $nbr_errors = system ("./check_parent_diagrams -cutonly -infile $parent_data_files$_ -dir $target_data_files -test 1> $_.log 2>&1 ") >> 8;
        #print "./print_parent_diagrams $parent_data_files$_ $target_file 1> $_.log 2>&1 \n"
    }
    else {
        printf( 'Executing : %-25s ', $_);
        $nbr_errors = system ("./$_ 1> $_.log 2>&1 ") >> 8;
    }
    my $history_file = "$_.hist";
    my $old_best_error_number = 100000;
    my $best_error_number = -1 ;
    my $last_error_number = -1 ;
    
     if ( open( my $hist_in , "<", $history_file) ){
        my @lines = <$hist_in>;
    	(undef,undef,$last_error_number,$old_best_error_number)=split(' ',$lines[$#lines]);
         close( $hist_in);
         if ( $nbr_errors < $old_best_error_number){
    	     $best_error_number = $nbr_errors
     	} else {
             $best_error_number = $old_best_error_number
     	} 
     } else {
     	$best_error_number=$nbr_errors;
     };
     
     
     
     open( my $hist , ">>", $history_file) ;
     my $date=`date +"%F %T"`;
    chomp($date);
     my $histmsg = "$date $nbr_errors $best_error_number\n";
     print $hist $histmsg;

    if ($nbr_errors < $last_error_number ){
        $nbr_improved_tests++;
    }   
    if ($nbr_errors > $last_error_number && $last_error_number >= 0 ){
        $nbr_worse_tests++;
    }   
    if ($nbr_errors == 0 ){
	print  colored ("PASSED\n",'green') ;
    } else {
    	    print  colored ("FAILED\n",'red') ;
	   $nbr_failed_tests++; 
    	if ( $last_error_number >= 0 ) {
	    print " Nbr of errors: $nbr_errors. (Last: $last_error_number  Best: $old_best_error_number)\n";
	} else {
	    print "Nbr of errors: $nbr_errors. \n";
	}
    }
   
}


    print "**********************************************\n" ;
    if ( $nbr_failed_tests > 0){
    	printf "Number of failed tests: %1d out of %2d\n", $nbr_failed_tests, $nbr_tests;
    } else {
    	printf colored( "All %1d test(s) passed!\n",'bold green'), $nbr_tests;
    }
    if ( $nbr_improved_tests > 0 ){
        printf colored ("Number of tests improved: %1d. Well done!\n",'blue'), $nbr_improved_tests;
    }
    if ( $nbr_worse_tests > 0 ){
        printf colored ("Number of tests that are worse than last time: %1d. Caution!\n",'bold red'), $nbr_worse_tests;
    }
    
    print "**********************************************\n" ;
    exit ($nbr_failed_tests);
