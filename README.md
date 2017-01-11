## Synopsis

A collection of functions useful for reading FASTA specified files,
parsing them, and reconstructing the fragments contained within
into the shortest possible super string

Includes test cases in tests.py and helper methods in utils.py

## Code Example

To reconstruct a fragment from a FASTA specified file, simply call `reconstruct_from_file(filename)`

Example
```
reconstruct_from_file("./data.txt")
```
## Approach

The program begins by reading in a FASTA specified from disk, using a simple state machine to parse the file. Each subfragment is stored into a dictionary, where the key is the subfragment name and the value is the subfragment itself.

Next, we generate a list of all possible left and right stubs. A left stub is a substring starting with the first character. These represent possible overlap regions on the lefthand side. For example, GATC could match G, GA, GAT, and GATC on the front end. Similarly, a right stub is a substring ending with the last character. For example, GATC could match C, TC, ATC, and GATC on the back end. Stub generation can be limited by the `min_length` parameter. We know fragments must overlap by at least half their length, so `min_length` can be set to half the length of the smallest subfragment

Given the set of all possible stubs, we create a map of each possible stub to the names of the subfragments that contain those stubs. For example, for the subfragment named "foo", consiting of "GATC",
- `right_stubs` would contain `{"GATC" : ["foo"],  "ATC" : ["foo"], "TC" : ["foo"], "C" : ["foo"]}`
- and `left_stubs` would contain `{"GATC" : ["foo"],  "GAT" : ["foo"], "GA" : ["foo"], "G" : ["foo"]}.`

Right and left stubs from subsequent fragments would be added to this list

After that, we have everything we need to reconstruct the fragment. We know the fragment order is guaranteed to be unique, so we pick a random subfragment to start. First, we try to add fragments on to the right side of the fragment. We sample successively shorter subfragments (always containing the last character). We're trying to match against the front half of the next fragment, so we look up our subfragments in our left_stubs map. We'll need to be sure to respect the overlap requirements, which specify that fragments must mutually overlap by at least half their length. After we have a match, we add it to our in-order list of assembled fragments, remove it from the list of unused fragments, and repeat the right-matching process with the fragment we just found. This continues until we run out of fragments

It's likely the random fragment we picked was not the first fragment in the superstring, so we need to repeat the process looking leftwards as well. We'll take successively shorter subfragments that always contain the first character and attempt to match them against the left half of the remaining fragments, contained in our right_stubs map.

At the end, we'll have a list of fragments in order, but with all of their overlapping sequences remaining. We'll zip them together, removing the duplicated nucleotides, and return the superstring.


## Installation

No installation should be necessary. The code relies only on standard python packages
- sys
- argparse
- random
- os
- textwrap

## Tests

The tests are designed to be platform agnostic, and include a command line utility to run them.

```
peter@ubuntu ~/rosalind % python tests.py -h
usage: tests.py [-h] [-m MODE] [-n NUM_ITERATIONS] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -m MODE, --mode MODE  The test mode to run. One of 'unit', 'stress' 'e2e',
                        'random' or 'all'.
  -n NUM_ITERATIONS, --num_iterations NUM_ITERATIONS
                        For stress and random tests, the number ofiterations
                        to run
  -v, --verbose         For stress and random tests, verbose output

```
Test modes
- **Unit** -- runs unit tests
- **Stress** -- runs a repeated number of large fragments
- **e2e** -- generates a random fragment, chops it up, writes it to disk, reads it back in, and reconstructs it
- **Random** -- repeatedly generates a random fragment, chops it up, and reconstructs it
- **All** -- runs all tests
