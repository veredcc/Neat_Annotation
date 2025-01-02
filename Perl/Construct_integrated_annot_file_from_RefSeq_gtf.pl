#Author: Vered Chalifa-Caspi
#Bioinformatics Core Facility, Ben-Gurion University
#veredcc@bgu.ac.il

#This script takes the output files of Extract_annotation_and_mapping_from_RefSeq_gtf.pl
#and unifies the results to one annotation file suitable for RNA-Seq analysis

use strict;
use warnings;

#output dir
my $out_dir = "Results/01.Annotation";

#input_files
my $genes_file              = $out_dir . "/" . "genes.txt";
my $transcripts_file        = $out_dir . "/" . "transcripts.txt";
my $proteins_file           = $out_dir . "/" . "proteins.txt";

#result_files
my $genes_annot_file        = $out_dir . "/" . "gene_annotation.txt";
my $log_file                = $out_dir . "/" . "log.txt";

#data structures
my $transcripts_info;
my $proteins_info;
my $genes_info;
my $gene2transcript;
my $trascript2protein;



open (LOG, ">$log_file") or die $!;

#read transcripts file
open (TR, $transcripts_file) or die $!;

my $headings1 = <TR>;
while (<TR>) {
	chomp;
	my ($transcript_id, $transcript_biotype, $gene_id, $product) = split (/\t/);
	$transcripts_info->{$transcript_id}->{"transcript_biotype"} = $transcript_biotype;
	$transcripts_info->{$transcript_id}->{"product"} = $product;
	$gene2transcript->{$gene_id}->{$transcript_id} = 1;
}

close (TR);

#read proteins file
open (PROT, $proteins_file) or die $!;

my $headings2 = <PROT>;
while (<PROT>) {
	chomp;
	my ($transcript_id, $protein_id, $gene_id, $product) = split (/\t/);
	$proteins_info->{$protein_id}->{"protein_product"} = $product;
	$trascript2protein->{$transcript_id}->{$protein_id} = 1;	
}

close (PROT);

#read genes file
open (GENE, $genes_file) or die $!;

my $headings3 = <GENE>;
while (<GENE>) {
	chomp;
	my ($gene_id, $gene_biotype, $description) = split (/\t/);
	$genes_info->{$gene_id}->{"gene_biotype"} = $gene_biotype;
	$genes_info->{$gene_id}->{"description"} = $description;
}
close (GENE);

#create integrated annotation file
open (OUT, ">$genes_annot_file") or die $!;
print OUT join ("\t", 'Gene_ID', 'Transcript_IDs', 'Protein_IDs', 'Gene_Type', 'Gene_Desctription', 'Transcript_Products', 'Protein_Products', ), "\n";

my $gene_id;
foreach $gene_id (sort keys %{$genes_info}) {
	my $gene_biotype = $genes_info->{$gene_id}->{"gene_biotype"};
	my $gene_description = $genes_info->{$gene_id}->{"description"};
	my $transcript_id;
	my @transcript_ids;
	my @transcript_biotypes;
	my @transcript_products;
	my @protein_ids;
	my @protein_products;
	foreach $transcript_id (sort keys %{$gene2transcript->{$gene_id}}) {
		push (@transcript_ids, $transcript_id);
		push (@transcript_biotypes, $transcripts_info->{$transcript_id}->{"transcript_biotype"});
		push (@transcript_products, $transcripts_info->{$transcript_id}->{"product"});
		my $protein_id;
		foreach $protein_id (sort keys %{$trascript2protein->{$transcript_id}}) {
			push (@protein_ids, $protein_id);
			push (@protein_products, $proteins_info->{$protein_id}->{"protein_product"});
		}
	}
	my $transcript_ids = join (";", @transcript_ids);
	my $transcript_biotypes = join (";", @transcript_biotypes);
	my $transcript_products = join (";", @transcript_products);
	my $protein_ids = join (";", @protein_ids);
	my $protein_products = join (";", @protein_products);
	
	print OUT join ("\t", $gene_id, $transcript_ids, $protein_ids, $gene_biotype, $gene_description, $transcript_products, $protein_products), "\n";	
}

close (OUT);
close (LOG);

print "Done\n";