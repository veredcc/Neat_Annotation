#Author: Vered Chalifa-Caspi
#Bioinformatics Core Facility, Ben-Gurion University
#veredcc@bgu.ac.il

#This script extracts annotation and mapping between genes-RNAs-proteins from NCBI's GTF file
#The GTF file can be downloaded from an annotated reference genome from NCBI's "Genomes" database.
#The script creates three files: per genes, per transcripts, per proteins.
#A subsequent script, Construct_integrated_annot_file_from_RefSeq_gtf.pl, integrates them to one annotation file 
#for use in RNA-Seq analysis

use strict;
use warnings;

#input files
my $ref_genome_gtf = "Input/genomic.gtf";

#output dir
my $out_dir = "Results/01.Annotation";

#output files
my $genes_file              = $out_dir . "/" . "genes.txt";
my $transcripts_file        = $out_dir . "/" . "transcripts.txt";
my $proteins_file           = $out_dir . "/" . "proteins.txt";

#data structures

my %cdss;  #cdss appear per exon with redundant information. this removes the redundancy

#go over gtf file

open (GTF,    $ref_genome_gtf)      or die $!;
open (GENES,  ">$genes_file")       or die $!;
open (TRANSC, ">$transcripts_file") or die $!;
open (PROT,   ">$proteins_file")    or die $!;

#print header line in result files
print GENES  join ("\t", "gene_id", "gene_biotype", "description"), "\n";
print TRANSC join ("\t", "transcript_id", "transcript_biotype", "gene_id", "product"), "\n";
print PROT   join ("\t", "transcript_id", "protein_id", "gene_id", "product"), "\n";

while (<GTF>) {
	chomp;
	next if /^#/;
	my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split(/\t/);
	my %att = parse_attributes($attribute);
	if ($feature eq "gene") {
		$att{"description"} = "" unless $att{"description"};
		print GENES join ("\t", $att{"gene_id"}, $att{"gene_biotype"}, $att{"description"}), "\n";
	} elsif ($feature eq "transcript") {
		$att{"product"} = "" unless $att{"product"};
		print TRANSC join ("\t", $att{"transcript_id"}, $att{"transcript_biotype"}, $att{"gene_id"}, $att{"product"}), "\n";
	} elsif ($feature eq "CDS") {
		$att{"product"} = "" unless $att{"product"};
		my $data = join ("\t", $att{"transcript_id"}, $att{"protein_id"}, $att{"gene_id"}, $att{"product"});
		unless ($cdss{$data}) {
			$cdss{$data} = 1;
			print PROT $data, "\n";
		}
	} 
}

close (GENES);
close (TRANSC);
close (PROT);
close (GTF);

print "Done\n";


#subroutines
############

sub parse_attributes {
    my ($attribute) = @_;
    #print "attribute: $attribute\n\n";

    $attribute =~ s/; $//;  # Remove the trailing '; ' if present
    
    my %struct;

    # Use a regex to match key-value pairs that respect quoted terms
    while ($attribute =~ /(\w+)\s+"(.*?)"(?:;|$)/g) {
        my $key = $1;
        my $value = $2;
        #print "key: $key\n";
        #print "value: $value\n";
        $struct{$key} = $value;
    }

    #print "\n";
    return %struct;
}
