#!/usr/bin/perl

use Bio::SeqIO;

my $seqio = Bio::SeqIO->new(-file => $ARGV[0], '-format' => 'Fasta');
open STDOUT, '>', "$ARGV[1]" or die;
while(my $seq = $seqio->next_seq) {
	$header = $seq->display_id;
	$sequence = $seq->seq;
	$nucl_index = 0;
	$result = '';
        @splitted_sequence = split('',"$sequence");
	while (@splitted_sequence) {
		$nucleotide = shift(@splitted_sequence);
		if ($nucl_index % 60 == 0){
			print("\n");
			$number_str = $nucl_index+1;
			$number_len = length("$nucl_index+1");
			$number_add = " " x (11 - $number_len);
			print($number_add);
			print($number_str);
			print(' ');
		}
		else {
			if ($nucl_index % 10 == 0){
				print(" ");
			};
		};
		if (index('ATCG',$nucleotide) != -1) {
			print($nucleotide);
		}
		else {
			print('N');
		};
		$nucl_index++;
	};
};
