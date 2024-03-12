#Author: Vered Caspi
#Date:   13.12.2023

#This script prepares KEGG annotation files for pathway enrichment analysis in RNA-Seq or proteomics projects.

#The script maps between Uniprot IDs and KEGG IDs directly,
#using a conversion file Uniprot-KEGG from KEGG REST API

#Rules are applied as follows:
#1. Start from all KEGG gene IDs from KEGG info file (per organism) retrieved from KEGG REST API
#2. If there is a direct conversion to uniprot (from A above) take all uniprot IDs.
#   If among them are uniprot IDs with the same gene symbol as the KEGG ID, take only them.
#   If none of them have the same gene symbl as the KEGG ID, take all of them.

#The script produces two result files:
#1. A detailed file for QA purposes
#2. A 2-column conversion file with KEGG ID and uniprot ID
#In addition, it prepares the necessary files for KEGG enrichment analysis in R

#The conversion file is used for:
#1. preparing files for upload to KEGG mapping-color tool
#2. counting the no. of genes per cluster (or gene list) per pathway
#3. enrichment testing in R 

#Instructions on how to retrieve the script's input files from KEGG REST API:
#(adjust to the relevant organism)
#Get gene to pathway in mouse:
#https://rest.kegg.jp/link/pathway/mmu
#in human:
#https://rest.kegg.jp/link/pathway/hsa
#
#Get pathway to name in mouse:
#https://rest.kegg.jp/list/pathway/mmu
#in human:
#https://rest.kegg.jp/list/pathway/hsa
#
#Get all mouse genes with their type, name and chromosomal location
#https://rest.kegg.jp/list/mmu
#in human:
#https://rest.kegg.jp/list/hsa
#
#Get all KEGG mouse genes and their uniprot
#https://rest.kegg.jp/conv/uniprot/mmu
#in human:
#https://rest.kegg.jp/conv/uniprot/hsa
#
#KEGG Rest instructions:
#https://www.genome.jp/linkdb/linkdb_api.html
#or use the web page:
#https://www.genome.jp/linkdb/
#
#Get gene info from uniprot - instructions:
#uniprot.org: in the search box, on the right, press on "Advanced", choose your Organism[OS], press on "search"
#press on Customize columns, choose to show: Reviewed, Entry Name, Gene Names (primary), Protein names
#press "save", press on "Download", choose format - tsv, press on "Download"

use strict;
use warnings;

my $organism = "Homo sapiens"; #KEGG organism name, e.g. Homo sapiens, Mus musculus

#input files
my $gene2pathway_file              = "From_KEGGREST_API/gene2pathway.txt";
my $pathway2name_file              = "From_KEGGREST_API/pathway2name.txt";
my $keggGene_info_file             = "From_KEGGREST_API/genes_info.txt";
#my $kegg2ensembl_file_from_KEGG    = "From_KEGGREST_API/Mouse_gene2ensembl.txt";
#my $ensembl_from_biomart_file      = "From_BioMart/mart_export_21.5.2023.txt";
#my $enst2ensg_file_from_ensmart    = "From_BioMart/mart_export_ensg2enst_28.5.2023.txt";
#my $uniprot2enst_from_uniprot_file = "From_Uniprot/uniprot2ensembl_from_uniprot_28.5.2023.tsv";
my $keggGene2up_file               = "From_KEGGREST_API/genes2uniprot.txt";
#my $uniprot2kegg_file_from_kegg    = "From_KEGGREST_API/genes2uniprot.txt";
my $uniprot_info_file              = "From_Uniprot/uniprotkb_organism_id_9606_2023_12_18.tsv";

#output dir
my $output_dir = "Functional_annotation_results/";

#result files

my $all_info_for_QA_file      = $output_dir . "KEGG_gene2all_info.txt";     #for QA
my $kegg_gene_to_uniprot_file = $output_dir . "KEGG_gene2uniprot_gene.txt"; #for perl script to prepare files for KEGG map and color
my $kegg_path_to_uniprot_file = $output_dir . "KEGG_pathway2gene.tab";   #for R enrichment analysis
my $kegg_path_to_name_file    = $output_dir . "KEGG_pathway2name.tab";   #for R enrichment analysis

#general info
#############

#open general files
open (GENE2PATH,     $gene2pathway_file) or die $!;
open (PATH2NAME,     $pathway2name_file) or die $!;
open (GENE2NAME,     $keggGene_info_file) or die $!;
open (KEGG2UP,       $keggGene2up_file) or die $!;
open (UP2INFO,       $uniprot_info_file) or die $!;

#data structures
my $gene2path;
my $path2genes;
my $path2name;
my $KeggID2info;
my $kegg2up;
my $up2kegg;
my $up2symbols;
my $up2title;

my $kegg2up_combined;
my $path2up;

#read and store general info
############################

while (<GENE2PATH>) {
	chomp;
	my ($kegg_gene_id, $kegg_path) = split (/\t/);
	$kegg_path =~ s/^path://;
	$gene2path->{$kegg_gene_id}->{$kegg_path} = 1;
	$path2genes->{$kegg_path}->{$kegg_gene_id} = 1;
}
while (<PATH2NAME>) {
	chomp;
	my ($kegg_path, $path_description) = split (/\t/);
	$path_description =~ s/ - $organism .*$//;
	$path2name->{$kegg_path} = $path_description;
}


while (<GENE2NAME>) {
	chomp;
	my ($kegg_gene_id, $gene_type, $gene_loc, $gene_desc) = split (/\t/);
	my ($gene_symbols, $gene_title);
	if ($gene_desc =~ /;/) {
		($gene_symbols, $gene_title) = split (/; /, $gene_desc);
	} else {
		$gene_symbols = "";
		$gene_title = $gene_desc;
	}
	$KeggID2info->{$kegg_gene_id}->{"symbols"} = $gene_symbols;
	$KeggID2info->{$kegg_gene_id}->{"title"} = $gene_title;	
}

while (<KEGG2UP>) {
	chomp;
	my ($kegg, $up) = split (/\t/);
	$up =~ s/^up://;
	$kegg2up->{$kegg}->{$up} = 1;
	$up2kegg->{$up}->{$kegg} = 1;
}

my $header = <UP2INFO>;
while (<UP2INFO>) {
	chomp;
	my ($up_id, $review_status, $up_acc, $symbols, $title) = split (/\t/);
	if ($symbols) {
		my @symbols = split("; ", $symbols);
		my $symbol;
		foreach $symbol (@symbols) {
			$up2symbols->{$up_id}->{$symbol} = 1;
		}
	}
	$up2title->{$up_id} = $title;
}

close (GENE2PATH);
close (PATH2NAME);
close (GENE2NAME);
close (KEGG2UP);
close (UP2INFO);

#Main combined
##############

open (OUT, ">$all_info_for_QA_file") or die $!;
print OUT join ("\t", "KEGG ID", "Orig No. Uniprot IDs", "Final No. Uniprot IDs", "Source", "KEGG Symbol", "KEGG Names",  "KEGG Title", "Uniprot ID", "Uniprot Symbol", "Uniprot Title"), "\n";

open (OUT_SHORT, ">$kegg_gene_to_uniprot_file") or die $!;

my $keggid;
foreach $keggid (sort keys %{$KeggID2info}) {
	my $symbols = $KeggID2info->{$keggid}->{"symbols"};
	my $kegg_symbol;
	if ($symbols =~ /,/) {
		($kegg_symbol) = $symbols =~ /^(.*?),/;
	} else {
		$kegg_symbol = $symbols;
	}
	my $keggid_title = "";
	$keggid_title = $KeggID2info->{$keggid}->{"title"} if ($KeggID2info->{$keggid}->{"title"});
	
	my $source;
	my $nr_up_ids = scalar (keys %{$kegg2up->{$keggid}});
	if ($nr_up_ids > 0) {  #if there is a uniprot (from KEGG), use it
		$source = "direct";
		my @all_upids = sort keys %{$kegg2up->{$keggid}};
		my @filtered_upids = ();
		my $upid;
		foreach $upid (@all_upids) {
			if ($up2symbols->{$upid}->{$kegg_symbol}) {
				push (@filtered_upids, $upid);
			}			
		}
		my @final_upids;
		if (@filtered_upids) {
			@final_upids = @filtered_upids;
		} else {
			@final_upids = @all_upids;
		}
		my $final_no_upids1 = scalar (@final_upids);
		foreach $upid (@final_upids) {
			my $up_symbols = join (";", sort keys (%{$up2symbols->{$upid}}));
			my $up_title = $up2title->{$upid};
			print OUT join ("\t", $keggid, $nr_up_ids, $final_no_upids1, $source, $kegg_symbol, $symbols, $keggid_title, $upid, $up_symbols, $up_title), "\n";
			print OUT_SHORT join ("\t", $keggid, $upid), "\n";
			$kegg2up_combined->{$keggid}->{$upid} = 1;
		}	
	}				
}

close (OUT);
close (OUT_SHORT);

#create files for R enrichment analysis
#######################################

#create kegg path 2 ensembl genes

open (OUT, ">$kegg_path_to_uniprot_file") or die $!;
print OUT join ("\t", "v1", "index"), "\n";

my $path;
foreach $path (sort keys %{$path2genes}) {
	my $kegg;
	foreach $kegg (sort keys %{$path2genes->{$path}}) {
		my $up;
		foreach $up (sort keys %{$kegg2up_combined->{$kegg}}) {  #here is the bug
			$path2up->{$path}->{$up} = 1;
		}
	}
	my $up;
	foreach $up (sort keys %{$path2up->{$path}}) {
		print OUT join ("\t", "path:" . $path, $up), "\n";
	}		
}

close (OUT);

#create path 2 name
###################

open (OUT, ">$kegg_path_to_name_file") or die $!;
print OUT join ("\t", "pathway", "info"), "\n";

my $pathid;
foreach $pathid (sort keys %{$path2name}) {
	my $name = $path2name->{$pathid};
	my $path = "path:" . $pathid;
	print OUT join ("\t", $path, $name), "\n";
}
close (OUT);

print "Done\n";
