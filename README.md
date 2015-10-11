# jellyfish-join

## Purpose

Given a [jellyfish](https://github.com/gmarcais/Jellyfish) database and a fasta or fastq file, perform a join of the kmers found in each sequence with their corresponding counts in the database.

This is very similar to the [query_per_sequnce](https://github.com/gmarcais/Jellyfish/tree/master/examples/query_per_sequence) example in jellyfish except that this will load the jellyfish database into memory (the same hash table it uses when counting). In my experience, this is much faster than the hash-aided binary search for a large number of queries.

## Requirements
  - Boost 1.54+
  - Jellyfish 2+
  - c++11

## Notes
Fasta files should be in a 2 line format (ie. no line wrapping)
