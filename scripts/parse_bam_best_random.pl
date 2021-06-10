#!/usr/bin/env perl

# For mutiple hits, the default setting of Bowtie2 reports one random hit from alignments with equivalent of MAPQ score.
# However, alignments with equivalent MAPQ scores are not necessarily equal best sometimes. For instance, alignment with one mismatch might have similar MAPQ score to those alignments without mismatch.
# Therefore, the better measurement of "equivalent mapping" is "alignment score" (defined as AS tag) instead of MAPQ score.
# This script is used to randomly select one alignment for read aligned to mulitple locations with equivalent alignment score (recorded under tag AS).

use strict;
use warnings;
use IO::File;
use IO::Handle;
use List::Util 'sum';
use Math::Random qw(:all);

my $splitTag = $ENV{"SPLITTAG"} || "";
my $outALL = $ENV{"OUTALL"} || "No";

# input is sam file sorted by sequence name
# number of hits is recorded under tag XN
# matchlen of read is recorded under tag XL
# copy number of collapse reads tag XC

my $mismatchCutoff=$ENV{"MISMATCH"} || 1;
my $fh;
if (@ARGV==0) {
        $fh=IO::Handle->new();
        $fh->fdopen(fileno(STDIN),"r");
}
else {
        $fh=IO::File->new();
        $fh->open($ARGV[0],"r");
}

my $last_id;
my @reads_pool;
while(my $line=$fh->getline()){
        chomp($line);
        if($line=~ m/^@/sg){
                print $line."\n";
                next;
        }
        if($splitTag ne "" && $line eq $splitTag){
                next;
        }
        my @ary=split(/\t/,$line);
        if($ary[5]=~ m/I/){
                pos($ary[5])=0;
        }
        my @matches = ($ary[5]=~ m/(\d+)M/g);
        my $matchlen = sum(@matches);
        $line = $line."\t"."XL:i:$matchlen";
        if($ary[0]=~ m/.*-(.*)/){
                $line = $line."\t"."XC:i:$1";
        }
        my ($mismatch)=($line=~ m/NM:i:(\d+)\s/sg);
        pos($line)=0;
        if(!defined($mismatch)){
                print $line."\n";
        }
        next if $mismatch>$mismatchCutoff;

        if(defined($last_id) && $last_id ne $ary[0]){
                my $maxscore=-1000;
                grep {$maxscore = $_->{"score"} if $_->{"score"} > $maxscore} @reads_pool;
                @reads_pool = grep {$_->{"score"} == $maxscore} @reads_pool;
                my $numHits = scalar(@reads_pool);
                if($numHits >1){
                        if($outALL ne "Yes"){
                                my $strSeed = $last_id;
                                $strSeed=~ s/[^a-zA-Z0-9]//sg;
                                $strSeed = length($strSeed)>20 ? substr($strSeed,-20) : $strSeed;
                                random_set_seed_from_phrase($strSeed);
                                @reads_pool = ($reads_pool[random_uniform_integer(1, 1, $numHits)-1]);
                        }
                }
                foreach my $record (@reads_pool){
                        my $outline = $record->{"read"};

                        $outline =~ s/XN:i:(\d)/XN:i:$numHits/;
                        $outline =~ s/^(.*?)\t/$1--XN:$numHits\t/;
                        print $outline."\n";
                }
                @reads_pool=();
        }
        $last_id=$ary[0];
        my ($score)=($line=~ m/AS:i:(.*?)\s/sg);
        if(!defined($score)){
                 print "score here\n";
         }
        push(@reads_pool,{"read"=> $line,"score"=>$score});
}

my $maxscore=-1000;
grep {$maxscore = $_->{"score"} if $_->{"score"} > $maxscore} @reads_pool;
@reads_pool = grep {$_->{"score"} == $maxscore} @reads_pool;
my $numHits = scalar(@reads_pool);
if($numHits >1){
        my $strSeed = $last_id;
        $strSeed=~ s/[^a-zA-Z0-9]//sg;
        $strSeed = length($strSeed)>20 ? substr($strSeed,-20) : $strSeed;
        random_set_seed_from_phrase($strSeed);
        @reads_pool = ($reads_pool[random_uniform_integer(1, 1, $numHits)-1]);
}

if(scalar(@reads_pool)==1){
        my $outline = $reads_pool[0]->{"read"};
        $outline =~ s/XN:i:(\d)/XN:i:$numHits/;
        $outline =~ s/^(.*?)\t/$1--XN:$numHits\t/;
        print $outline."\n";
} else {
   foreach my $record (@reads_pool){
        my $outline = $record->{"read"};
        $outline =~ s/XN:i:(\d)/XN:i:$numHits/;
        $outline =~ s/^(.*?)\t/$1--XN:$numHits\t/;
        print $outline."\n";
    }
}
