#!/usr/bin/env perl
#-------------------------------------------------------------------------------
# Written by Tim Campbell - added to NUMA by Alex Reinecke
#-------------------------------------------------------------------------------
# PERL script to generate, for each input target object, a complete list of
# object and include prerequisites (dependencies).
#
# The source files in a directory (based on suffix) or listed in the Makefile
# SOURCEF, SOURCEC and SOURCEH lists are parsed for module and include
# prerequisites.  Since include targets do not have make commands associated
# with them a spanning tree of include prerequisites is generated for each
# object target.  Object targets with prerequisites are output to the
# dependency makefile fragment.
#
# By default, objects are prefixed with the name of the directory containing
# the object source or by the value of the OBJDIR variable defined in the local
# Makefile (if it is defined).  The option "-p NONE" will cause the objects to
# be output without any prefix.
#-------------------------------------------------------------------------------

use strict;
use Getopt::Std;
use vars qw( $opt_h $opt_v $opt_d $opt_s $opt_f $opt_o $opt_p $opt_l );

my $program = $0;
$program =~ s/.*\///g;

# usage
sub usage {
  my $c = shift;
  print STDERR "Usage: $program [options] [dirs]\n";
  print STDERR " -h            : help -- display usage\n";
  print STDERR " -v            : verbose -- prints progress to stderr\n";
  print STDERR " -d            : use Makefile to obtain directory lists\n";
  print STDERR " -s            : use Makefile to obtain source file lists\n";
  print STDERR " -f <pattern>  : use suffix match <pattern> to identify Fortran sources\n";
  print STDERR " -o <out_file> : output to <out_file>; default is ./depend.mk\n";
  print STDERR " -p <obj_dir>  : attach prefix <obj_dir>/ to objects (NONE = no prefix;\n";
  print STDERR "               : LIB = list objects as member of library)\n";
  print STDERR " -l            : also output VPATH settings to <out_file>\n";
  print STDERR " [dirs]        : directories for the recursive search; default is .\n";
  exit $c;
}

# process optional arguments
getopts( 'hvdsf:o:p:l' ) || &usage;
&usage(0) if $opt_h;


#-------------------------------------------------------------------------------
# global variables
#-------------------------------------------------------------------------------
# top-level directories
my @tdirs=(); my $d;
if ( defined($ARGV[0]) ) {
  while ( @ARGV ) {
    $d = realpath( shift @ARGV );
    if (! -d $d) { die "$program: ERROR cannot find $d\n"; }
    push(@tdirs, $d);
  }
}
else {
  push(@tdirs, getcwd());
}
# output makefile fragment
my $out_file;
if ( $opt_o ) {
  $out_file = $opt_o;
}
else {
  $out_file = join("/",getcwd(),"depend.mk");
}
# global hash of targets (keyed on directory)
my %tgts = ();
# global hash of source prerequisites (keyed on target)
my %src_prqs = ();
# global hash of include prerequisites (keyed on target)
my %inc_prqs = ();
# global hash of module prerequisites (keyed on target)
my %mod_prqs = ();
# global hash of module source objects (keyed on lowercase module name)
my %mod_sobj = ();
# pattern matching string of directories to skip
my $skip_dirs = qr/(?:^\.|config|build|prologues|run|bin|mod|obj|^lib$|^bak$|^save$)/;
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# main
#-------------------------------------------------------------------------------
# traverse directory trees and build hash of prerequisites
traverse_dir_trees( \&process_files, \@tdirs, $skip_dirs, $opt_d, $opt_v );
# reconcile module prerequisites for each object
reconcile_mod_prqs();
# reconcile include prerequisites for each object
reconcile_inc_prqs();
# output prerequisite lists
output_prqs( $out_file );
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# subroutine to process the source and include files in a directory and add
# to the prerequisite hashes
#-------------------------------------------------------------------------------
sub process_files {
  my ($d) = @_;
  my $sub = $program.": ".(caller(0))[3];
  my $suffix = qr/\.[^.]*/;
  my $curdir = qr/\$\(CURDIR\)/;
  my $t; my $f; my $ff; my $m; my @files=(); my @mods=(); my $odir; my $lib;
  print STDERR "\n\nBegin processing in directory: $d\n" if $opt_v;
  if ( -e "Makefile" ) {
    $odir = &get_make_var( "Makefile", "OBJDIR", $opt_v );
    if ( $odir and $odir =~ m/$curdir/ ) { $odir =~ s/$curdir/$d/; }
    if ( $opt_p ) {
      if ( $opt_p =~ m/LIB/ ) {
        $lib = &get_make_var( "Makefile", "LIB", $opt_v );
      }
    }
  }
#-process Fortran sources
  if ( $opt_s ) {
    @files = &get_make_var( "Makefile", "SOURCEF", $opt_v );
  }
  else {
    my $file_fsrc = qr/(?:\.(?:f|f90|F|F90)$)/; #suffixes for Fortran sources 
    if ( $opt_f ) { $file_fsrc = qr/(?:\.(?:$opt_f)$)/; }
    opendir( DIR, $d) or die "$sub: ERROR opening $d for reading: $!";
    @files = grep( /$file_fsrc/, readdir(DIR) );
    closedir( DIR );
  }
  foreach $f ( @files ) {
    print STDERR "Processing Fortran file: $f\n" if $opt_v;
    ($t = $f) =~ s/$suffix/.o/;
    if ( $opt_p and $opt_p =~ m/NONE/ ) { }
    elsif ( $opt_p and $opt_p =~ m/LIB/ and $lib ) { $t = $lib."(".$t.")"; }
    elsif ( $opt_p and $opt_p !~ m/LIB/ ) { $t = $opt_p."/".$t; }
    elsif ( $odir ) { $t = $odir."/".$t; }
    else { $t = $d."/".$t; }
    push( @{$tgts{$d}}, $t );
    print STDERR "\tTarget: $t\n" if $opt_v;
    $src_prqs{$t} = $d."/".$f;
    push( @{$inc_prqs{$t}}, &get_inc_prqs($f,$opt_v) );
    print STDERR "\t\tInclude Prereqs: @{$inc_prqs{$t}}\n" if $opt_v;
    foreach $ff ( @{$inc_prqs{$t}} ) {
        push( @{$mod_prqs{$t}}, &get_mod_prqs($ff,$opt_v) );
    }
    push( @{$mod_prqs{$t}}, &get_mod_prqs($f,$opt_v) );
    print STDERR "\t\tModule  Prereqs: @{$mod_prqs{$t}}\n" if $opt_v;
    print STDERR "\tSource: $t\n" if $opt_v;
    @mods = &get_mod_gens($f,$opt_v);
    print STDERR "\t\tModules Generated: @mods\n" if $opt_v;
    foreach $m ( @mods ) {
      if ( not exists($mod_sobj{$m}) ) {
        $mod_sobj{$m} = $t;
      }
      else {
        die "$sub: ERROR duplicate module source: $m $t $mod_sobj{$m}\n";
      }
    }
  }
#-process C/C++ sources
  if ( $opt_s ) {
    @files = &get_make_var( "Makefile", "SOURCEC", $opt_v );
  }
  else {
    my $file_csrc = qr/(?:\.(?:c|C)$)/; #suffixes for C/C++ sources
    opendir( DIR, $d) or die "$sub: ERROR opening $d for reading: $!";
    @files = grep( /$file_csrc/, readdir(DIR) );
    closedir( DIR );
  }
  foreach $f ( @files ) {
    print STDERR "Processing C/C++ file: $f\n" if $opt_v;
    ($t = $f) =~ s/$suffix/.o/;
    if ( $opt_p and $opt_p =~ m/NONE/ ) { }
    elsif ( $opt_p and $opt_p =~ m/LIB/ and $lib ) { $t = $lib."(".$t.")"; }
    elsif ( $opt_p and $opt_p !~ m/LIB/ ) { $t = $opt_p."/".$t; }
    elsif ( $odir ) { $t = $odir."/".$t; }
    else { $t = $d."/".$t; }
    push( @{$tgts{$d}}, $t );
    print STDERR "\tTarget: $t\n" if $opt_v;
    $src_prqs{$t} = $d."/".$f;
    push( @{$inc_prqs{$t}}, &get_inc_prqs($f,$opt_v) );
    print STDERR "\t\tInclude Prereqs: @{$inc_prqs{$t}}\n" if $opt_v;
  }
#-process includes
  if ( $opt_s ) {
    @files = &get_make_var( "Makefile", "SOURCEH", $opt_v );
  }
  else {
    my $file_isrc = qr/(?:\.(?:h|inc)$)/; #suffixes for include sources
    opendir( DIR, $d) or die "$sub: ERROR opening $d for reading: $!";
    @files = grep( /$file_isrc/, readdir(DIR) );
    closedir( DIR );
  }
  foreach $f ( @files ) {
    print STDERR "Processing include file: $f\n" if $opt_v;
    $t = $f;
    push( @{$tgts{$d}}, $t );
    print STDERR "\tTarget: $t\n" if $opt_v;
    push( @{$inc_prqs{$t}}, &get_inc_prqs($f,$opt_v) );
    print STDERR "\t\tINCLUDE Prereqs: @{$inc_prqs{$t}}\n" if $opt_v;
#    push( @{$mod_prqs{$t}}, &get_mod_prqs($f,$opt_v) );
#    print STDERR "\t\tModule  Prereqs: @{$mod_prqs{$t}}\n" if $opt_v;
  }
  print STDERR "Done processing in directory: $d\n" if $opt_v;
}
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# subroutine to reconcile module prerequisites for each object target
#-------------------------------------------------------------------------------
sub reconcile_mod_prqs {
  my $d; my $t; my $m; my %seen; my @objs;
  my $sub = $program.": ".(caller(0))[3];
  print STDERR "\n\nBegin reconciling module prerequisites\n" if $opt_v;
  foreach $t ( keys %mod_prqs ) {
    print STDERR "\tReconciling target: $t $#{$mod_prqs{$t}}\n" if $opt_v;
    #match module name with source object
    @objs=();
    foreach $m ( @{$mod_prqs{$t}} ) {
      if ( exists($mod_sobj{$m}) ) {
        print STDERR "\t\tmodule $m matched with $mod_sobj{$m}\n" if $opt_v;
        if ( $t ne $mod_sobj{$m} ) {
          push(@objs, $mod_sobj{$m});
        }
        else {
          print STDERR "\t\t\tcircular reference, prerequisite removed\n" if $opt_v;
        }
      }
      else {
        print STDERR "\t\tmodule $m not matched, prerequisite removed\n" if $opt_v;
      }
    }
    #remove duplicate prereqs
    %seen=(); @objs = grep {!$seen{$_}++} @objs;
    #set prereqs or delete empty target
    if ( @objs ) {
      print STDERR "\t\tNew module prereqs: @objs\n" if $opt_v;
      @{$mod_prqs{$t}}=@objs;
    }
    else {
      print STDERR "\t*** target $t removed from module prereqs ***\n" if $opt_v;
      delete($mod_prqs{$t});
    }
  }
  print STDERR "Done reconciling module prerequisites\n" if $opt_v;
}
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# subroutine to reconcile include prerequisites for each object target
#-------------------------------------------------------------------------------
sub reconcile_inc_prqs {
  my $sub = $program.": ".(caller(0))[3];
  my $file_objs = qr/\.(?:o|o\))$/;
  my $t; my $p; my @ps; my %seen; my %incv;
  print STDERR "\n\nBegin reconciling include prerequisites\n" if $opt_v;
  foreach $t ( grep( /$file_objs/, keys %inc_prqs ) ) {
    print STDERR "\tReconciling target: $t $#{$inc_prqs{$t}}\n" if $opt_v;
    %incv = (); #used to mark visited includes
    foreach $p ( @{$inc_prqs{$t}} ) {
      print STDERR "\t\tCall spanning_inc_prqs_tree for $p\n" if $opt_v;
      push(@{$inc_prqs{$t}}, spanning_inc_prqs_tree($p, \%incv));
    }
    #remove duplicate prereqs
    %seen=(); @{$inc_prqs{$t}} = grep {!$seen{$_}++} @{$inc_prqs{$t}};
    #remove include prerequisites that are not listed as a target
    @ps=();
    print STDERR "\t\tOld include prereqs: @{$inc_prqs{$t}}\n" if $opt_v;
    foreach $p ( @{$inc_prqs{$t}} ) {
      if ( exists($inc_prqs{$p}) ) { push(@ps, $p); }
    }
    $inc_prqs{$t} = [@ps];
    print STDERR "\t\tNew include prereqs: @{$inc_prqs{$t}}\n" if $opt_v;
    #remove targets that do not have include prerequisites
    if ( not @{$inc_prqs{$t}} ) {
      print STDERR "\t*** target $t removed from include prereqs ***\n" if $opt_v;
      delete($inc_prqs{$t});
    }
  }
  print STDERR "Done reconciling include prerequisites\n" if $opt_v;
}
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# recursive subroutine to traverse the spanning tree of include prerequisites
#-------------------------------------------------------------------------------
sub spanning_inc_prqs_tree {
  my ($p, $incv) = @_;
  my $sub = $program.": ".(caller(0))[3];
  my @incs = ();
  my $i;
  if ( exists($$incv{$p}) ) { return @incs; }
  $$incv{$p}=1;
  if ( exists($inc_prqs{$p}) ) {
    print STDERR "\t\t\tReconciling include: $p $#{$inc_prqs{$p}}\n" if $opt_v;
    foreach $i ( @{$inc_prqs{$p}} ) {
      push(@incs, $i);
      print STDERR "\t\t\tCall spanning_inc_prqs_tree for $i\n" if $opt_v;
      push(@incs, spanning_inc_prqs_tree($i, \%$incv));
    }
  }
  return @incs;
}
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# subroutine to output the prerequisite list for each object target
#-------------------------------------------------------------------------------
sub output_prqs {
  my ($file) = @_;
  my $sub = $program.": ".(caller(0))[3];
  my $file_objs = qr/\.(?:o|o\))$/;
  my $curdir = qr/\$\(CURDIR\)/;
  my $d; my $o; my $line; my $odir;
#-open file
  open FILO, ">$file" or die "$sub: error opening file $file for output: $!\n";
#-output target prerequisite list
  print STDERR "\n\nBegin outputting target preq list\n" if $opt_v;
  foreach $d ( keys %tgts ) {
    print STDERR "Directory: $d\n" if $opt_v;
    print FILO "#\n# $d\n#\n";
    foreach $o ( grep( /$file_objs/, @{$tgts{$d}} ) ) {

#      print STDERR "\t$o\n"
#print STDERR "\t$#mod_prqs{$o} \n"
#print STDERR "\t$#inc_prqs{$o} \n"

      if ( exists($mod_prqs{$o}) or exists($inc_prqs{$o}) ) {
        print STDERR "\tOutput target: $o $#{$inc_prqs{$o}} $#{$mod_prqs{$o}}\n" if $opt_v;
       #$line = $o." : ".$src_prqs{$o};
        $line = $o." :";
        if ( exists($mod_prqs{$o}) ) {
          $line = join(" ",$line,@{$mod_prqs{$o}});
        }
        if ( exists($inc_prqs{$o}) ) {
          $line = join(" ",$line,@{$inc_prqs{$o}});
        }
        print FILO "$line\n";
      }
      else {
        print STDERR "\tSkip target: $o $#{$inc_prqs{$o}} $#{$mod_prqs{$o}}\n" if $opt_v;
      }
    }
  }
  print STDERR "Done outputting target preq list\n" if $opt_v;
#-output search path (VPATH) settings
  if ( $opt_l ) {
    print STDERR "\n\nBegin outputting VPATH settings\n" if $opt_v;
    print FILO "#\n# VPATH settings\n#\n";
    foreach $d ( keys %tgts ) {
      print STDERR "\tDirectory: $d\n" if $opt_v;
      print FILO "VPATH += $d\n";
    }
    print STDERR "Done outputting VPATH settings\n" if $opt_v;
  }
#-close file
  close FILO;
}
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
# PERL Utility Functions -- these can be placed in a separate file
#-------------------------------------------------------------------------------

use strict;
use Cwd 'getcwd', 'realpath';

#this is needed when these functions are placed in a separate file (???)
my %delim = (qr/'/=>qr/'/, qr/"/=>qr/"/, qr/</=>qr/>/);


#-------------------------------------------------------------------------------
# subroutine to traverse directory trees and execute provided function
#-------------------------------------------------------------------------------
sub traverse_dir_trees {
  my ($func, $tdirs, $sdirs, $opt_d, $opt_v) = @_;
  my $sub = $program.": ".(caller(0))[3];
  my $d; my $dir; my @dirs=(); my $cdir;
  foreach $d ( @$tdirs ) {
    $dir = realpath($d);
    $cdir = getcwd();
    chdir $dir;
    &$func( $dir );
    if ( $opt_d and -e "Makefile" ) {
      @dirs = &get_make_var( "Makefile", "DIRS", $opt_v );
    }
    else {
      opendir( DIR, $dir) or die "$sub: ERROR opening $dir for reading: $!";
      @dirs = grep( -d && !/$sdirs/, readdir(DIR) );
      closedir( DIR );
    }
    if ( @dirs ) {
      &traverse_dir_trees( \&$func, \@dirs, $sdirs, $opt_d, $opt_v );
    }
    chdir $cdir;
  }
}
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# subroutine to scan a Makefile and obtain the list of values associated
# with a make variable
#-------------------------------------------------------------------------------
sub get_make_var {
  my ($file, $var, $opt_v) = @_;
  my $sub = $program.": ".(caller(0))[3];
  my $split_on_equal = qr/\s*[\+\:]?=\s*/;
  my $split_on_space = qr/\s+/;
  my $cont_line = qr/\\$/;
  my @list=(); my $line; my $t; my $s;
  print STDERR "Begin scanning Makefile for $var\n" if $opt_v;
  if ( -e $file ) {
  open FILE, "<$file" or die "$sub: ERROR opening file $file for input: $!\n";
  LINE: while ( defined( $line = <FILE> ) ) {
    chomp $line;
    #collect continued lines into single line
    if ( $line =~ s/$cont_line//) {
      $line .= <FILE>;
      redo unless eof(FILE);
    }
    #match line with variable name
    if ( $line =~ m/^\s*$var\s*([\+\:]?)=/ ) {
      print STDERR "\tline match: $line\n" if $opt_v;
      if ( $1 =~ m/\+/ or $line =~ m/\$\($var\)/ ) { #append
        $line =~ s/\$\($var\)//; #remove variable references
        ($t, $s) = split(/$split_on_equal/, $line, 2);
        print STDERR "\t\tappend: $s\n" if $opt_v;
        push( @list, split(/$split_on_space/, $s) );
      }
      else { #reset
        ($t, $s) = split(/$split_on_equal/, $line, 2);
        print STDERR "\t\treset: $s\n" if $opt_v;
        @list = split(/$split_on_space/, $s);
      }
    }
  }
  close FILE;
  print STDERR "\t$var = @list\n" if $opt_v;
  print STDERR "Done scanning Makefile for $var\n" if $opt_v;
  } else {
  print STDERR "Makefile does not exist. Return empty for $var\n" if $opt_v;
  }
  if ( wantarray() ) { return @list; }
  else { my $scalar = join(" ",@list); return $scalar; }
}
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# subroutine to scan a file and build a list of include prerequisites
#-------------------------------------------------------------------------------
sub get_inc_prqs {
  my ($file, $opt_v) = @_;
  my $sub = $program.": ".(caller(0))[3];
  my @list=(); my %cnt; my $line;
  my %delim = (qr/'/=>qr/'/, qr/"/=>qr/"/, qr/</=>qr/>/);
  open FILE, "<$file" or die "$sub: ERROR opening file $file for input: $!\n";
  foreach $line ( <FILE> ) {
    if ( $line =~ m/^[\#\s]*include\s*(['""'<])([\w\.\/]*)$delim{\1}/i ) {
      if ( $2 ) {
        print STDERR "\t\tinclude line match: $line" if $opt_v;
        $cnt{$2}++; #add to list
      }
    }
  }
  close FILE;
  @list = sort keys %cnt;
  return @list;
}
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# subroutine to scan a file and build a list of module prerequisites
#-------------------------------------------------------------------------------
sub get_mod_prqs {
  my ($file, $opt_v) = @_;
  my $sub = $program.": ".(caller(0))[3];
  my @list=(); my %cnt; my $line;
  open FILE, "<$file" or die "$sub: ERROR opening file $file for input: $!\n";
  foreach $line ( <FILE> ) {
    if ( $line =~ m/^\s*use\s+(\w*)/i ) {
      if ( $1 ) {
        print STDERR "\t\tuse line match: $line" if $opt_v;
        $cnt{lc($1)}++; #add to list (lowercase)
      }
    }
  }
  close FILE;
  @list = sort keys %cnt;
  return @list;
}
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# subroutine to scan a file and build a list of generated modules
#-------------------------------------------------------------------------------
sub get_mod_gens {
  my ($file, $opt_v) = @_;
  my $sub = $program.": ".(caller(0))[3];
  my @list=(); my $nm; my %cnt; my $line;
  open FILE, "<$file" or die "$sub: ERROR opening file $file for input: $!\n";
  foreach $line ( <FILE> ) {
    if ( $line =~ m/^\s*module\s+procedure/i ) {
        next
    }
    if ( $line =~ m/^\s*module\s+(\w*)/i ) {
      if ( $1 ) {
        print STDERR "\t\tmodule line match: $line" if $opt_v;
        $cnt{lc($1)}++; #add to list (lowercase)
      }
    }
  }
  close FILE;
  @list = sort keys %cnt;
  return @list;
}
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# END PERL Utility Functions
#-------------------------------------------------------------------------------
