#!/usr/bin/env perl
use strict;
use warnings;
use IO::Handle;
use IO::File;

my $fh;
if(@ARGV==0){
        $fh = IO::Handle->new();
        $fh->fdopen(fileno(STDIN),"r");
} else {
        $fh = IO::File->new();
        $fh->open($ARGV[0],"r");
}

my $id = $ARGV[1];

my $totalreads;
my $unmapped;
my $mapped;
my $uniquemapped;
while(my $line=$fh->getline()){
        if($line=~ m/(.*?) reads; of these/g){
                $totalreads=$1;
        }
        if($line=~ m/(\d+) .*aligned 0 times/g){
                $unmapped=$1;
        }
        if($line=~ m/(\d+) .*aligned exactly 1 time/g){
                $uniquemapped=$1;
        }
        if($line=~ m/(\d+) .*aligned >1 times/g){
                $mapped=$1;
        }
}
$mapped=$mapped+$uniquemapped;

print $id."\t".$totalreads."\t".$mapped."\t".$uniquemapped."\t".$unmapped."\n";
