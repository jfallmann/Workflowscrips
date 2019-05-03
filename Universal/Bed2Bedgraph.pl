#!/usr/bin/env perl
#Script Bed2Bedgraph.pl;
#Last changed Time-stamp: <2019-05-03 12:59:17 fall> by joerg

#### use things ###
use strict;
use warnings;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use autodie;
use PerlIO::gzip;
use Cwd;
use File::Path qw (make_path);
use Pod::Usage;
use Data::Dumper;
### use own modules
use FindBin::Real qw(Bin); # locate this script
use lib Bin() . "/../lib";
use Collection;

my $VERBOSE=0;
my ( $file, $odir, $species, $forward, $reverse, $type, $peak, $chromsize, $track, $conv, $anno, $scale, $u );


###############################################################################
# start reading command line options
###############################################################################

BEGIN{pod2usage(-verbose => 0)
		  unless GetOptions(
			  "file|f=s"	  => \$file,
			  "odir|o=s"	  => \$odir,
			  "forward|x=s"   => \$forward,
			  "reverse|y=s"   => \$reverse,
			  "type|t=s"	  => \$type,
			  "peak|p=s"	  => \$peak,
			  "chrom|c=s"	  => \$chromsize,
			  "track|a=s"	  => \$track,
			  "converted|v=s" => \$conv,
			  "anno|o=s"	  => \$anno,
			  "species|s=s"   => \$species,
			  "scale|l=s"	  => \$scale,
			  "ustart|u=s"    => \$u,
			  "help|h"		  => sub{pod2usage(-verbose => 1)},
			  "man|m"		  => sub{pod2usage(-verbose => 2)},
			  "verbose"		  => sub{ $VERBOSE++ }
		  );
}

##################
# Main
##################

#my $dir   = cwd();
#if(!$odir){$odir = "Bedgraph/$file";}
#my $outdir = $odir;
#
#if (!-d "$outdir"){
#	make_path("$outdir");
#}

die "You must provide filename and file with chrom-sizes\n" unless ($file && $chromsize);

###process call
my $pid = $$;
(my $job = `cat /proc/$pid/cmdline`)=~ s/\0/ /g;
print STDERR $job,"\n";

print STDERR "Analyzing Peak bedfile\n" if ($peak || $file =~ /bed_spk/);
print STDERR "Generating track\n" if ($track && $track eq "track");
print STDERR "Analyzing already genomic positions\n" if ($conv && $conv eq "on");
print STDERR "Adding annotation information\n" if ($anno && $anno eq "anno");
$anno ='' unless (defined $anno);
#print "FILE: $file\n";
$species = 'dummy' unless (defined $species);
my $spec;
if ($species eq "Human" or $species eq "human"){$spec				 = "Homo_sapiens";}
elsif ($species eq "Mouse" or $species eq "mouse"){$spec			 = "Mus_musculus";}
elsif ($species eq "Zebrafish" or $species eq "zebrafish"){$spec   = "Danio_rerio";}
elsif ($species eq "Drosophila" or $species eq "drosophila"){$spec = "Drosophila_melanogaster";}
elsif ($species eq "Worm" or $species eq "worm"){$spec			 = "Caenorhabditis_elegans";}
else{$spec														 = $species;}

my ($sizes, $tag)=Collection::fetch_chrom_sizes("$spec","$chromsize");

my ($covplus, $covminus, $annop, $annom, $totalreads) = Collection::bed_to_coverage($file, $anno, $sizes, $peak, $conv);

if (defined $scale){
	$scale /= $totalreads;
}
else{
	$scale = 1;
}

if (%{$covplus}){
	foreach my $key (sort {$covplus->{$a} <=> $covplus->{$b}} keys %{$covplus}){
		my $chr = $key;
		$chr = 'chr'.$chr if ($tag);
		foreach my $pos (sort {$covplus->{$key}->{$a} cmp $covplus->{$key}->{$b}} keys %{$covplus->{$key}}){
			if ($track && $track eq "track"){
				my $cov=$covplus->{$key}->{$pos}*$scale;
				my $annotation=$annop->{$key}->{$pos} if ($anno && $anno eq "anno");
				open (my $OUT, ">>:gzip", $forward);#"chr".$chrom."\.$type\.fw.track.gz");
				if ($anno && $anno eq "anno"){
					print $OUT "$chr\t$pos\t".($pos+1)."\t$cov\t$annotation\n";
				}
				else{
					print $OUT "$chr\t$pos\t".($pos+1)."\t$cov\n";
				}
				close ($OUT);
			}
			else{
				my $cov=$covplus->{$key}->{$pos}*$scale;
				my $annotation=$annop->{$key}->{$pos} if ($anno && $anno eq "anno");
				open (my $OUT, ">>:gzip", $forward);# "chr".$chrom."\.$type\.fw.gz");
				if ($anno && $anno eq "anno"){
					print $OUT "$chr\t$pos\t".($pos+1)."\t$cov\t$annotation\n";
				}
				else{
					print $OUT "$chr\t$pos\t".($pos+1)."\t$cov\n";
				}
				close ($OUT);
			}
		}
	}
}
else{
	open (my $OUT, ">>:gzip", $forward);# "chr".$chrom."\.$type\.fw.gz");
	foreach my $chr (keys %{$sizes}){
		$chr = 'chr'.$chr if ($tag);
		print $OUT join("\t",$chr,0,0,0)."\n";
	}
	close($OUT);
}

if(%{$covminus}){
	foreach my $key (sort {$covminus->{$a} <=> $covminus->{$b}} keys %{$covminus}){
		my $chr = $key;
		$chr = 'chr'.$chr if ($tag);
		foreach my $pos (sort {$covminus->{$key}->{$a} <=> $covminus->{$key}->{$b}} keys %{$covminus->{$key}}){
			if ($track && $track eq "track"){
				my $cov=$covminus->{$key}->{$pos}*$scale;
				my $annotation=$annom->{$key}->{$pos} if ($anno && $anno eq "anno");
				open (my $OUT, ">>:gzip", $reverse);# "chr".$chrom."\.$type\.re.track.gz");
				if ($anno && $anno eq "anno"){
					print $OUT "$chr\t$pos\t".($pos+1)."\t-$cov\t$annotation\n";
				}
				else{
					print $OUT "$chr\t$pos\t".($pos+1)."\t-$cov\n";
				}
				close ($OUT);
			}
			else{
				my $cov=$covminus->{$key}->{$pos}*$scale; # at least one occurence
				my $annotation=$annom->{$key}->{$pos} if ($anno && $anno eq "anno");
				my @tmp=split(/\t/,$key);
				my $chrom=$tmp[0];
				open (my $OUT, ">>:gzip", $reverse);#"chr".$chrom."\.$type\.re.gz");
				if ($anno && $anno eq "anno"){
					print $OUT "$chr\t$pos\t".($pos+1)."\t$cov\t$annotation\n";
				}
				else{
					print $OUT "$chr\t$pos\t".($pos+1)."\t$cov\n";
				}
				close ($OUT);
			}
		}
	}
}
else{
	open (my $OUT, ">>:gzip", $reverse);
	foreach my $chr (keys %{$sizes}){
		$chr = 'chr'.$chr if ($tag);
		print $OUT join("\t",$chr,0,0,0)."\n";
	}
	close($OUT);
}
