#!/usr/bin/perl
use strict;
use warnings;

=pod

=head1 Analyse sequence

The program is intended to use to analyse if a sequence is indeed coding or non-coding.
It should be run as:

    $perl deliver_2.pl intronic.fasta coding.fasta unkown.fasta
=cut

my $ref_intronic = $ARGV[0];
my $ref_exonic = $ARGV[1];
my $unkown_seq = $ARGV[2];

print STDERR &log_time(), "Starting analysis\n";
&main(\$ref_intronic, \$ref_exonic, $unkown_seq);
print STDERR &log_time(), "Ending analysis\n";

exit(0);

##### Functions ######

=over

=back

=head2 read_fasta

Given a file path reads the file if it is in fasta format

    my $fasta_ref = read_fasta("intronic.fasta");

Returns a reference to a hash with the identifier as a key and the sequence as
a value. If the any sequence contains a non IUPAC symbol for DNA sequence it
throws an error.

=cut

sub read_fasta{
    # Returns a hash with the key being the id and the value the sequence
    local $/ = ">";
    my $file_fasta = shift;
    my %seqs;
    my $count = 0; # record counter
    open(FASTA, "$file_fasta") or die("Unable to open file $file_fasta\n");
    <FASTA>; # The first element is empty! remove it silently

    print STDERR &log_time(), "Reading $file_fasta fasta file\n";
    while (<FASTA>){
        chomp;
        my ($id, @seq) = split /\n/;
        my $id_name = &id_name($id);

        print STDERR &log_time(), "Reading sequence $id_name\n";
        $count++;

        my $seq = join("", @seq);
        # Prevents to incorporate any sequence not using the IUPAC code for DNA
        # Case insensitive! some sequences have undercase letters,
        # Prevents sequences with less than 6 nucleotids
        if ($seq =~ /[^ACGTRYSWKMBDHVN.]/ig){
            my $s = "does not use the IUPAC code for DNA\n";
            print STDERR &log_time(), $id_name, $s;
            next;
        } elsif (length $seq < 6){
            print STDERR &log_time(), $id_name, "is shorter than an hexamer\n";
            next;
        }
        $seqs{$id} = uc($seq);
    };
    close(FASTA);
    print STDERR &log_time(), "Just read $count FASTA sequences\n";
    return \%seqs;
};

=head2 hexamer_freq

Calculates the frequency of hexamers of the values of the hash

Requires a reference to a hash with the fasta file.

    my %hexamer_freq = &hexamer_freq($fasta_ref);

Returns a reference of a hash with the hexamer frequencies of the hash.

=cut

sub hexamer_freq{
    # Given a sequence find the number of hexamers of each type
    # First argument a sequence, second a hash of hexamers
    my (%counter, %freq);
    my $hash_ref = shift;
    my $total = 0;
    while(my ($id, $seq) = (each %{$hash_ref})){
        my $id_name = &id_name($id);


        # hexamers overlapping 3 nucleotides (1 codon position)
        for(my $i = 0; $i < length $seq; $i += 3){
            my $hexamer = substr($seq, $i, 6);
            # Avoid calculating hexamer for sequences shorter than a hexamer.
            if (length($hexamer) < 6){
                last;
            };
            $total++;
            if (defined $counter{$hexamer}){
                $counter{$hexamer}++;
            } else{
                $counter{$hexamer} = 1;
            }
        }
        print STDERR &log_time(),  $id_name, "'s hexamer counted\n";
    }

    # Calculates the frequency
    foreach my $key (keys %counter){
        $freq{$key} = $counter{$key}/$total;
    }
    return \%freq;
};

=head2 main

Calculates if a given sequence is intronic or coding.

First argument the reference intronic fasta file.
Second argument the reference coding fasta file.
Third argument the fasta file to analyse.

    main(\$intronic.fasta, \$coding.fasta, $unkown.fasta);

Prints if a sequence is intronic or coding
=cut

sub main{
    # Arguments
    my ($intronic, $coding, $unknown) = @_;

    # We read only once the sequence of interes and pass the referencee
    my $unknown_ref = &read_fasta($unknown);
    my $freq_intronic = &hexamer_freq( &read_fasta($$intronic));
    my $freq_coding = &hexamer_freq( &read_fasta($$coding));

    while( my ($id,  $sequen) = each %{$unknown_ref}){
        my $n_hexamers = 0;
        my $score = 0;
        for(my $i = 0; $i < length $sequen; $i += 3){
            my $hex = substr($sequen, $i, 6);
            # Stop if the length is shorter than 6 nucleotides
            if (length($hex) < 6){
                last;
            };
            $n_hexamers++;

            if (defined ${$freq_coding}{$hex} and defined ${$freq_intronic}{$hex}){
                $score += log2(${$freq_coding}{$hex} / ${$freq_intronic}{$hex});
            } else {
                print STDERR &log_time, "Unknown hexamer $hex. Use a better trainer\n"
            }
        }

        # Decide if a sequence is coding or intronic
        my $hexamer_score = $score/$n_hexamers;
        #print "$hexamer_score\n";
        if ($hexamer_score < 0){
            print "Intronic $id\n";
        }elsif($hexamer_score > 0){
            print "Coding $id\n";
        }elsif($hexamer_score == 0){
            my $msg = "Unable to decide if $id is or not a coding sequence\n";
            print $msg;
        }else {
            print STDERR &log_time, "Error in $id\n";
        };
    };
};

=head2 log_time

Prepares the localtime to be used for the logger

    log_time()
    2015/11/23 19:40:44
=cut

sub log_time{
    # Obtain the time to use in the logger
    my ($sec,$min,$hour,$mday,$mon,$year, @others)= localtime(time);
    my $nice_timestamp = sprintf ( "%04d/%02d/%02d %02d:%02d:%02d ",
                                   $year+1900,$mon+1,$mday,$hour,$min,$sec);
    return $nice_timestamp;
};

=head2 id_name

Splits the name by spaces and returns the first one.

Used to reduce space when spliting

    print($id);
    NM_001103386.01.e12_cds11 chrX 23878 11577 11716 FWD(+) 140bp  frame: 1
    my $name = id_name($id);
    print($name);
    NM_001103386.01.e12_cds11;

=cut

sub id_name{
    # Cut the name to something easiler to remember.
    my ($id )= @_;
    my @id_seq = split(/\s+/, $id);
    # There are some sequence of the same RefSeq identifier
    #my @id_name = split(/\./, $id_seq[0]);
    return $id_seq[0];
}

=head2 log2

Calculates the logarithm in base 2 of a given number

    print 1 == log2(2)
    1

=cut

sub log2 {
    # Calculates the logarithm of 2 of a number
    my $n = shift;
    return log($n)/log(2);
};
