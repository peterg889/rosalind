"""
    A collection of tests for driver.py, along with a command line utility
    for help. Test framework agnostic -- built to be used directly by
    build agents
"""

from driver import (parse_file, generate_right_overlaps,
                   generate_left_overlaps, generate_overlap_map,
                   is_in_list, get_matching_item, match_right,
                   match_left, merge_strings, assemble_fragments,
                   names_to_fragments, reconstruct, reconstruct_from_file)
from utils import (generate_random_fragment, generate_test_subfragments,
                   write_fragments_to_disk)
import argparse
import random
import os

def test_parse_file():
    TEST_FILENAME = "./test1.txt"
    fragment_name_map, min_len = parse_file(TEST_FILENAME)

    assert fragment_name_map["foo"] == "GCG"

    assert fragment_name_map["bar"] == \
                  "TAGTGTGGGAAGCTTAACAAGCTGCGAACTGCAACTCAAACTTAAGACACCGCTCCGAAC"

    assert fragment_name_map["xyz"] == \
                  "TAGTGTGGGAAGCTTAACAATAGTGTGGGAAGCTTAACAA"

    assert min_len == 3

def test_generate_right_overlaps():
    TEST_FRAGMENT = 'GATCGATC'

    # min_length is greater than fragment length
    fragments = generate_right_overlaps(TEST_FRAGMENT, 9)
    assert fragments == []

    # min_length equals fragment length
    fragments = generate_right_overlaps(TEST_FRAGMENT, 8)
    assert fragments == ['GATCGATC']

    fragments = generate_right_overlaps(TEST_FRAGMENT, 6)
    assert fragments == ['GATCGATC', 'ATCGATC', 'TCGATC']

    fragments = generate_right_overlaps(TEST_FRAGMENT)
    assert fragments == ['GATCGATC', 'ATCGATC', 'TCGATC', \
                         'CGATC', 'GATC', 'ATC', 'TC', 'C']

def test_generate_left_overlaps():
    TEST_FRAGMENT = 'GATCGATC'

    # min_length is greater than fragment length
    fragments = generate_left_overlaps(TEST_FRAGMENT, 9)
    assert fragments == []

    # min_length equals fragment length
    fragments = generate_left_overlaps(TEST_FRAGMENT, 8)
    assert fragments == ['GATCGATC']

    fragments = generate_left_overlaps(TEST_FRAGMENT, 6)
    assert fragments == ['GATCGA', 'GATCGAT', 'GATCGATC']

    fragments = generate_left_overlaps(TEST_FRAGMENT)
    assert fragments == ['G', 'GA', 'GAT', 'GATC', 'GATCG', \
                         'GATCGA', 'GATCGAT', 'GATCGATC']

def test_generate_overlap_map():

    # base case
    fragments = {'foo': 'GATCG'}
    right_stubs, left_stubs = generate_overlap_map(fragments)

    assert right_stubs == {'ATCG': ['foo'], 'GATCG': ['foo'], 'CG': ['foo'], \
                           'G': ['foo'], 'TCG': ['foo']}
    assert left_stubs == {'GATCG': ['foo'], 'GAT': ['foo'], 'GA': ['foo'], \
                          'G': ['foo'], 'GATC': ['foo']}

    # test min_length
    right_stubs, left_stubs = generate_overlap_map(fragments, 4)

    assert right_stubs == {'ATCG': ['foo'], 'GATCG': ['foo']}
    assert left_stubs == {'GATCG': ['foo'], 'GATC': ['foo']}

    # test shared overlaps
    fragments = {'foo': 'GATCG', 'bar': 'AAACG'}
    right_stubs, left_stubs = generate_overlap_map(fragments)

    assert right_stubs == {'AAACG': ['bar'], 'ACG': ['bar'], \
                           'G': ['foo', 'bar'], 'AACG': ['bar'], \
                           'CG': ['foo', 'bar'], 'ATCG': ['foo'], \
                           'GATCG': ['foo'], 'TCG': ['foo']}

    assert left_stubs == {'A': ['bar'], 'AA': ['bar'], 'AAACG': ['bar'], \
                          'AAA': ['bar'], 'G': ['foo'], 'GAT': ['foo'], \
                          'AAAC': ['bar'], 'GATCG': ['foo'], 'GA': ['foo'], \
                          'GATC': ['foo']}


    # test multiple fragments
    fragments = {'foo': 'GATCG', 'bar': 'AAACG', 'xyz': "GACG"}
    right_stubs, left_stubs = generate_overlap_map(fragments, 2)

    assert right_stubs == {'AAACG': ['bar'], 'ACG': ['xyz', 'bar'], \
                           'AACG': ['bar'], 'CG': ['xyz', 'foo', 'bar'], \
                           'ATCG': ['foo'], 'GATCG': ['foo'], 'TCG': ['foo'], \
                           'GACG': ['xyz']}
    assert left_stubs == {'AA': ['bar'], 'AAACG': ['bar'], 'AAA': ['bar'], \
                          'GAT': ['foo'], 'AAAC': ['bar'], 'GATCG': ['foo'], \
                          'GA': ['xyz', 'foo'], 'GAC': ['xyz'], \
                          'GACG': ['xyz'], 'GATC': ['foo']}

def test_is_in_list():
    little_list = []
    big_list = ['a']
    assert is_in_list(big_list, little_list) == False

    little_list = ['c']
    big_list = ['a', 'b', 'c']
    assert is_in_list(big_list, little_list) == True

    little_list = ['z']
    big_list = ['a', 'b', 'c']
    assert is_in_list(big_list, little_list) == False

    little_list = ['z', 'a']
    big_list = ['a', 'b', 'c']
    assert is_in_list(big_list, little_list) == True

def test_get_matching_item():
    little_list = []
    big_list = ['a']
    assert get_matching_item(big_list, little_list) == None

    little_list = ['c']
    big_list = ['a', 'b', 'c']
    assert get_matching_item(big_list, little_list) == 'c'

    little_list = ['z']
    big_list = ['a', 'b', 'c']
    assert get_matching_item(big_list, little_list) == None

    little_list = ['z', 'a']
    big_list = ['a', 'b', 'c']
    assert get_matching_item(big_list, little_list) == 'a'

def test_match_right():
    fragments = {'foo': 'GATCGT', 'bar': 'TCGTAA'}
    fragment_names = ['bar']
    _, left_stubs = generate_overlap_map(fragments)
    match = match_right(fragments['foo'], left_stubs, fragment_names)
    assert match == 'bar'

    fragments = {'foo': 'GATCGT', 'bar': 'TCGTAA', 'xyz' : 'GTAATC'}
    fragment_names = ['bar', 'xyz']
    _, left_stubs = generate_overlap_map(fragments)
    match = match_right(fragments['foo'], left_stubs, fragment_names)
    assert match == 'bar'

    fragments = {'foo': 'GATCGT', 'bar': 'TCGTAA', 'xyz' : 'GTAATC'}
    fragment_names = ['foo', 'xyz']
    _, left_stubs = generate_overlap_map(fragments)
    match = match_right(fragments['bar'], left_stubs, fragment_names)
    assert match == 'xyz'

    # make sure the longest match is selected
    fragments = {'foo': 'GTAACGG', 'bar': 'TCGTAAA', 'xyz' : 'GTAAATC'}
    fragment_names = ['foo', 'xyz']
    _, left_stubs = generate_overlap_map(fragments)
    match = match_right(fragments['bar'], left_stubs, fragment_names)
    assert match == 'xyz'

    fragments = {'foo': 'GATCGT', 'bar': 'TCGTAA', 'xyz' : 'GTAATC'}
    fragment_names = ['foo', 'bar']
    _, left_stubs = generate_overlap_map(fragments)
    match = match_right(fragments['xyz'], left_stubs, fragment_names)
    assert match == None

    # check fragments won't match that overlap by exactly half
    fragments = {'foo': 'GATCGT', 'bar': 'TCGTAA', 'abc': 'TAAGGG'}
    fragment_names = ['foo', 'abc']
    _, left_stubs = generate_overlap_map(fragments)
    match = match_right(fragments['bar'], left_stubs, fragment_names)
    assert match == None

    # check greater than half length condition holds for odd lengths too
    fragments = {'foo': 'GATCGT', 'def': 'TCGATAA', 'hij': 'TAAAGGG'}
    fragment_names = ['foo', 'hij']
    _, left_stubs = generate_overlap_map(fragments)
    match = match_right(fragments['def'], left_stubs, fragment_names)
    assert match == None

def test_match_left():
    fragments = {'foo': 'GATCGT', 'bar': 'TCGATC'}
    fragment_names = ['bar']
    right_stubs, _ = generate_overlap_map(fragments)
    match = match_left(fragments['foo'], right_stubs, fragment_names)
    assert match == 'bar'

    fragments = {'foo': 'GATCGT', 'bar': 'TCGATC'}
    fragment_names = ['foo']
    right_stubs, _ = generate_overlap_map(fragments)
    match = match_left(fragments['bar'], right_stubs, fragment_names)
    assert match == None

    # make sure the longest match is selected
    fragments = {'foo': 'AGATCGT', 'bar': 'TCAGATC', 'xyz': 'GATGATC'}
    fragment_names = ['bar', 'xyz']
    right_stubs, _ = generate_overlap_map(fragments)
    match = match_left(fragments['foo'], right_stubs, fragment_names)
    assert match == 'bar'

    fragments = {'foo': 'GATCGT', 'bar': 'TCGATC'}
    fragment_names = ['bar']
    right_stubs, _ = generate_overlap_map(fragments)
    match = match_left(fragments['foo'], right_stubs, fragment_names)
    assert match == 'bar'

def test_merge_strings():

    # empty string cases
    first = ""
    second = "bbbcc"
    merge = merge_strings(first, second)
    assert merge == "bbbcc"

    first = "aaabb"
    second = ""
    merge = merge_strings(first, second)
    assert merge == "aaabb"

    first = "aaabb"
    second = "bbbcc"
    merge = merge_strings(first, second)
    assert merge == "aaabbbcc"

    first = "aaabb"
    second = "ccdd"
    merge = merge_strings(first, second)
    assert merge == "aaabbccdd"

    first = "AAAA"
    second = "ABB"
    merge = merge_strings(first, second)
    assert merge == "AAAABB"

    first = "AA"
    second = "AAAABB"
    merge = merge_strings(first, second)
    assert merge == "AAAABB"

    first = "GATCGATC"
    second = "GATCAAAA"
    merge = merge_strings(first, second)
    assert merge == "GATCGATCAAAA"

    first = "GATATATATACCC"
    second = "ATATATATACCCG"
    merge = merge_strings(first, second)
    assert merge == "GATATATATACCCG"

    first = "GATTTTTAAACCCAAA"
    second = "TTTTAAACCCAAAGAT"
    merge = merge_strings(first, second)
    assert merge == "GATTTTTAAACCCAAAGAT"

    first = "A"
    second = "AAAAAAAA"
    merge = merge_strings(first, second)
    assert merge == "AAAAAAAA"

    first = "A"
    second = "AAAABBBBB"
    merge = merge_strings(first, second)
    assert merge == "AAAABBBBB"

    first = "AAABAAA"
    second = "AAABBBBBBBB"
    merge = merge_strings(first, second)
    assert merge == "AAABAAABBBBBBBB"

def test_assemble_fragments():
    fragment_list = ["AAA", "ABB"]
    combined = assemble_fragments(fragment_list)
    assert combined == "AAABB"

    fragment_list = ["AAA", "ABB", "BBBBCC"]
    combined = assemble_fragments(fragment_list)
    assert combined == "AAABBBBCC"

    fragment_list = ["AAA", "XYZ"]
    combined = assemble_fragments(fragment_list)
    assert combined == "AAAXYZ"

    fragment_list = ["AAA", "AAAABBB"]
    combined = assemble_fragments(fragment_list)
    assert combined == "AAAABBB"

    fragment_list = ["AAAAA", "AAA"]
    combined = assemble_fragments(fragment_list)
    assert combined == "AAAAA"

    fragment_list = ["AAAAA", ""]
    combined = assemble_fragments(fragment_list)
    assert combined == "AAAAA"

    fragment_list = ["", "AAA"]
    combined = assemble_fragments(fragment_list)
    assert combined == "AAA"

def test_names_to_fragments():
    names = ["foo", "bar"]
    fragment_dict = {"foo" : "GAT", "bar" : "CTA"}

    fragment_list = names_to_fragments(names, fragment_dict)
    assert fragment_list == ['GAT', 'CTA']

def run_unit_tests():
    test_parse_file()
    test_generate_right_overlaps()
    test_generate_left_overlaps()
    test_generate_overlap_map()
    test_is_in_list()
    test_get_matching_item()
    test_match_right()
    test_match_left()
    test_merge_strings()
    test_assemble_fragments()
    test_names_to_fragments()

def stress_test_driver():
    """
        Generates a fragment of length TEST_LENGTH, then creates subfragments of
        length CHOP_LENGTH, each offset by STEP_SIZE. Tests to ensure the
        reconstructed fragment is the same as the initial fragment
    """
    TEST_LENGTH = 20000
    CHOP_LENGTH = 1000
    STEP_SIZE = 400


    test_fragment = generate_random_fragment(TEST_LENGTH)

    fragment_dict = generate_test_subfragments(test_fragment,
                                               CHOP_LENGTH, STEP_SIZE)

    reconstructed = reconstruct(fragment_dict, CHOP_LENGTH / 2)

    assert test_fragment == reconstructed

def random_test_driver():
    """
        Similar to stress_test_driver, but with random length, chop_length, and
        step_size. Boundaries are put in place to make sure generated test cases
        are sane, though may randomly generate some invalid test cases when
        there are identical subfragments generated
    """

    length = random.randrange(100, 2000)
    chop_length = random.randrange(length/10, length/2)
    step_size = random.randrange(5, max(chop_length/2 - 1, 6))

    test_fragment = generate_random_fragment(length)

    fragment_dict = generate_test_subfragments(test_fragment, chop_length,
                                               step_size, True)

    reconstructed = reconstruct(fragment_dict, chop_length / 2)

    if test_fragment != reconstructed:
        print test_fragment
        print reconstructed
    assert test_fragment == reconstructed

    return length, chop_length, step_size

def run_stress_test(n, verbose=False):
    for i in range(n):
        if verbose:
            print i
        stress_test_driver()

def run_random_test(n, verbose=False):
    for i in range(n):
        if verbose:
            print i
        random_test_driver()

def run_e2e_test(cleanup=True):

    """
        Generates a fragment of length TEST_LENGTH, then creates subfragments of
        length CHOP_LENGTH, each offset by STEP_SIZE. Writes file out to disk
        with name TEST_FILE. Reads TEST_FILE, then tests to ensure the
        reconstructed fragment is the same as the initial fragment

        Args:
            cleanup (bool) : if True, cleans up the file it created
                             if False, leaves the file
                             defaults to True
    """
    TEST_LENGTH = 20000
    CHOP_LENGTH = 1000
    STEP_SIZE = 400

    TEST_FILE = "test.txt"

    initial_test_fragment = generate_random_fragment(TEST_LENGTH)

    fragment_dict = generate_test_subfragments(initial_test_fragment,
                                               CHOP_LENGTH, STEP_SIZE)

    write_fragments_to_disk(fragment_dict, TEST_FILE)

    reconstructed = reconstruct_from_file(TEST_FILE)

    assert initial_test_fragment == reconstructed

    if cleanup:
        os.remove(TEST_FILE)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--mode', default='unit', type=str,
                        help='The test mode to run. One of \'unit\',' +
                        ' \'stress\' \'e2e\', \'random\' or \'all\'.')
    parser.add_argument('-n', '--num_iterations', default=100, type=int,
                        help="For stress and random tests, the number of" +
                        "iterations to run")
    parser.add_argument('-v', '--verbose', dest="verbose", action='store_true',
                        help="For stress and random tests, verbose output")

    parser.set_defaults(verbose=False)

    args = parser.parse_args()

    n = args.num_iterations
    verbose = args.verbose
    mode = args.mode.lower()

    num_tests = 0

    if mode == "unit":
        run_unit_tests()
        num_tests = 11
    elif mode == "stress":
        run_stress_test(n, verbose)
        num_tests = n
    elif mode == "random":
        run_random_test(n, verbose)
        num_tests = n
    elif mode == "e2e":
        run_e2e_test()
        num_tests = 1
    elif mode == "all":
        if verbose:
            print "running unit tests"
        run_unit_tests()
        if verbose:
            print "running stess tests"
        run_stress_test(n, verbose)
        if verbose:
            print "running random tests"
        run_random_test(n, verbose)
        if verbose:
            print "running e2e tests"
        run_e2e_test()

        num_tests = 11 + 2 * n + 1
    else:
        print "unrecognized test mode " + mode
        exit()

    test_str = ""
    if num_tests > 1:
        test_str = "tests"
    else:
        test_str = "test"

    print mode + " " + test_str + " completed successfully!"
    print str(num_tests) + " " + test_str + " run"

if __name__ == "__main__":
    main()


