"""
    A collection of functions useful for reading FASTA specified files,
    parsing them, and reconstructing the fragments contained within
    into the shortest possible super string
"""
import sys

def parse_file(filename):
    """
        Given a filename in FASTA format
        (https://en.wikipedia.org/wiki/FASTA_format), read and parse into a dict
        mapping fragment name to fragment content. Removes any new lines and
        white space as necessary

        Args:
            filename (str) : the file to read in

        Returns:
            fragment_name_map (dict) : dictionary mapping fragment
                                       name to fragment content

            min_len (int)    : the smallest length fragment encountered, useful
                               for limiting later subfragment generation
    """
    try:
        fragment_file = open(filename)
    except IOError as e:
        print "I/O error({0}): {1}".format(e.errno, e.strerror)

    fragment_name_map = {}
    name = None
    fragment = ''

    min_len = sys.maxint

    for line in fragment_file:

        # we're on a name line
        if line.startswith(">"):
            if name is not None:
                fragment_name_map[name] = fragment
                if len(fragment) < min_len:
                    min_len = len(fragment)

            name = line[1:].rstrip()
            fragment = ''
        # we're not on a name line
        else:
            fragment = fragment + line.rstrip()

    # don't forget the last one
    # assumes last line was not a name lines
    fragment_name_map[name] = fragment

    return fragment_name_map, min_len

def generate_right_overlaps(fragment, min_len=1):
    """ Generates a list of possible subfragments that could match this
        fragment, taken from the right-most side. min_len restricts the minimum
        subfragment length, inclusive. Fragments must always contain the
        last character

        Example:

        GATC

        could match

        C, TC, ATC, and GATC

        from right hand side

        Args:
            fragment (str) : the fragment to generate subfragments from
            min_len (int)  : the minimum subfragment length to generate,
                             inclusive. defaults to 1

        Returns:
            right_overlaps (list) : a list of all possible subfragments,
                                    taken from the right hand side
    """
    right_overlaps = []
    i = 0
    while i < len(fragment) - min_len + 1:
        overlap = fragment[i:]
        right_overlaps.append(overlap)
        i = i+1

    return right_overlaps

def generate_left_overlaps(fragment, min_len=1):
    """ Generates a list of possible subfragments that could match this
        fragment, starting with the left-most side. min_len restricts the
        minimum subfragment length, inclusive. Fragments must always contain the
        first character

        Example:

        GATC

        could match

        G, GA, GAT, GATC

        from left hand side

        Args:
            fragment (str) : the fragment to generate subfragments from
            min_len (int)  : the minimum subfragment length to generate,
                             inclusive. defaults to 1

        Returns:
            left_overlaps (list) : a list of all possible subfragments, taken
                                   from the left hand side
    """
    left_overlaps = []
    i = 1
    while i < len(fragment) + 1:
        overlap = fragment[:i]
        if len(overlap) >= min_len:
            left_overlaps.append(overlap)
        i = i+1

    return left_overlaps

def generate_overlap_map(fragments, min_len=1):
    """
        Given a dictionary that contains a mapping of fragment names to
        fragment contents, generates two dictionaries, useful for piecing
        overlapping fragments together.

        For each fragment, first generate right subfragments with length
        greater than min_len. For each possible sub fragment, keep a list of
        fragments that contain that subfragment, stored in right_stubs. Repeat
        process for all left subfragments, storing in left_stubs

        For example, if you called it with fragments = {"foo", "GATC"} and
        min_len = 2, right_stubs would contain
        {"GATC" : ["foo"],  "ATC" : ["foo"], "TC" : ["foo"]}

        and left_stubs would contain
        {"GATC" : ["foo"],  "GAT" : ["foo"], "GA" : ["foo"]}

        Args:
            fragments (dict) : dictionary mapping fragment names to
                               fragment contents
            min_len (int)    : optional parameter that limits the minimum length
                               subfragments to be generated, inclusive

        Returns:
            right_stubs (dict) : Dictionary mapping possible right subfragments
                                 to a list of fragment names that contain those
                                 subfragments
            left_stubs (dict)  : Dictionary mapping possible left subfragments
                                 to a list of fragment names that contain those
                                 subfragments
    """
    right_stubs = {}
    left_stubs = {}
    for fragment_name in fragments:

        # generate right subfragments
        right_overlaps = generate_right_overlaps(fragments[fragment_name],
                                                 min_len)

        for overlap in right_overlaps:
            if overlap in right_stubs:
                right_stubs[overlap].append(fragment_name)
            else:
                right_stubs[overlap] = [fragment_name,]

        # generate left subfragments
        left_overlaps = generate_left_overlaps(fragments[fragment_name],
                                               min_len)

        for overlap in left_overlaps:
            if overlap in left_stubs:
                left_stubs[overlap].append(fragment_name)
            else:
                left_stubs[overlap] = [fragment_name,]

    return right_stubs, left_stubs

def is_in_list(big_list, little_list):
    """
        Helper method. Returns True if there exists an element in little_list
        that's also in big_list and false otherwise

        Args:
            big_list (list)    : the list to search in
            little_list (list) : the set of things to look for

        Returns:
            True if there exists an element in little_list that's also in
            big_list False otherwise
    """
    for item in little_list:
        if item in big_list:
            return True

    return False

def get_matching_item(big_list, little_list):
    """
        Helper method. Returns the first item in little_list that's contained
        within big_list. Returns None if there is no item in little_list that's
        also in big_list

        Args:
            big_list (list)    : the list to search in
            little_list (list) : the set of things to look for

        Returns:
            item (str) : the first item in little_list that's also in big_list
            None       : if there are no items in little_list that are also
                         in big_list
    """
    for item in little_list:
        if item in big_list:
            return item
    return None

def match_right(right_fragment, left_stubs, remaining_fragments):
    """
        Given a fragment, a dictionary of left stubs (generated by
        generate_overlap_map), and a list of fragments that haven't
        been used yet, pick the remaining fragment with the longest possible
        overlap. Respect minimum overlap length requirements as specified in
        the instructions -- fragments must each overlap by a minimum of half
        their length

        Tries to match the largest possible fragment first. For example, with
        right_fragment = "GATC", it would try to first match "GATC", then "ATC".
        "TC" and "C" would be considered but skipped because they don't
        satisfy the minimum length requirements. Returns None if there are no
        matching fragments in remaining_fragments

        Args :
            right_fragment (str) : the fragment to use when looking for matches
            left_stubs (dict)    : a dictionary matching possible left_stubs to
                                   the fragments that contain them, generated by
                                   generate_overlap_map

            remaining_fragments (list)  : A list of strings representing the
                                          names of fragments that haven't been
                                          assembled yet

        Returns:
            match (str) : The name of a fragment that is the longest possible
                          left match for right_fragment among the fragments
                          we haven't used yet
                          Returns None if no such fragment exists
    """
    i = 0
    while i < len(right_fragment) + 1:
        subfragment = right_fragment[i:]
        match = None

        # if the subfragment exists in the map, and there's a fragment
        # in the list of fragments that contain that sub fragment that we
        # haven't used yet, and both length requirements are satisfied,
        # we have a match
        if subfragment in left_stubs \
        and is_in_list(remaining_fragments, left_stubs[subfragment]) \
        and len(subfragment) > len(right_fragment)/2 \
        and len(subfragment) > len(left_stubs[subfragment]) / 2:
            match = get_matching_item(remaining_fragments,
                                      left_stubs[subfragment])
            break

        i = i + 1

    return match

def match_left(left_fragment, right_stubs, remaining_fragments):
    """
        Given a fragment, a dictionary of right_stubs (generated by
        generate_overlap_map), and a list of fragments that haven't
        been used yet, pick the remaining fragment with the longest possible
        overlap. Respect minimum overlap length requirements as specified in
        the instructions -- fragments must each overlap by a minimum of half
        their length

        Tries to match the largest possible fragment first. For example, with
        left_fragment = "GATC", it would try to first match "GATC", then "GAT".
        "GA" and "G" would be considered but skipped because they don't
        satisfy the minimum length requirements. Returns None if there are no
        matching fragments in remaining_fragments

        Args :
            right_fragment (str) : the fragment to use when looking for matches
            left_stubs (dict)    : a dictionary matching possible left_stubs to
                                   the fragments that contain them, generated by
                                   generate_overlap_map

            remaining_fragments (list)  : A list of strings representing the
                                          names of fragments that haven't been
                                          assembled yet

        Returns:
            match (str) : The name of a fragment that is the longest possible
                          right match for right_fragment among the fragments
                          we haven't used yet
                          Returns None if no such fragment exists
    """

    i = 0
    match = None
    i = len(left_fragment)

    while i > 0:
        subfragment = left_fragment[:i]

        # if the subfragment exists in the map, and there's a fragment
        # in the list of fragments that contain that sub fragment that we
        # haven't used yet, and both length requirements are satisfied,
        # we have a match
        if subfragment in right_stubs \
        and is_in_list(remaining_fragments, right_stubs[subfragment]) \
        and len(subfragment) > len(left_fragment)/2 \
        and len(subfragment) > len(right_stubs[subfragment]) / 2:

            match = get_matching_item(remaining_fragments,
                                        right_stubs[subfragment])
            break
        i = i - 1

    return match

def merge_strings(a, b):
    """
        Given two strings a and b, returns the shortest possible superstring
        while preserving ordering

        The algorithm attempts to match the longest possible substring first,
        then tries progressively shorter substrings. For instance, given "AAAB"
        and "ABBA", the algorithm would try

        AAAB
        ABBA  -- no match

        AAAB
         ABBA -- no match

        AAAB
          ABBA -- match!

        It constructs the shortest possible superstring by adding the
        non overlapping regions from each string to the shared overlapping area

        Args:
            a (str) : first string
            b (str) : second string

        Returns:
            ret (str) : shortest possible superstring of a and b, preserving
                        order
    """
    if not a:
        return b
    if not b:
        return a

    alen = len(a)
    blen = len(b)

    overlap_len = min(alen, blen)

    aoverlap = a[-overlap_len:]
    boverlap = b[:overlap_len]

    while overlap_len > 0 and aoverlap != boverlap:
        overlap_len = overlap_len - 1
        aoverlap = a[-overlap_len:]
        boverlap = b[:overlap_len]

    ret = a[:-overlap_len] + a[-overlap_len:] + b[overlap_len:]

    return ret

def assemble_fragments(fragment_list):
    """
        Given a list of fragment strings are in the correct order, returns the
        shortest possible superstring. For example, ["GATC", "ATCG", "TCGA"]

        would return

        GATC   +
         ATCG  +
          TCGA +
        --------
        GATCGA

        Args:
            fragment_list (list) : a list of fragment strings in the correct
                                   order

        Returns:
            super_fragment (str) : the shortest possible superstring,
                                   preserving order

    """
    super_fragment = ""

    for fragment in fragment_list:
        super_fragment = merge_strings(super_fragment, fragment)

    return super_fragment

def names_to_fragments(fragment_names, fragment_dict):
    """
        Helper function. Given a list of fragment names and a dictionary
        mapping fragment names to fragment contents, return a list of fragment
        contents in the same order

        Args:
            fragment_names (list) : list of fragment name strings
            fragment_dict (dict)  : dict mapping fragment names to fragment
                                    contents

        Returns:
            ret_list (list) : list of fragment contents, in the same order as
                              fragment_names

    """

    ret_list = []
    for fragment in fragment_names:
        ret_list.append(fragment_dict[fragment])
    return ret_list

def reconstruct(fragment_dict, min_len=1):
    """
        Given a dictionary mapping fragment names to fragment contents, return
        the shortest possible superstring

        Starts by picking a random fragment, then greedily trying to match
        the fragment with the longest overlap to the right, and then to the
        left

        Args:
            fragment_dict (dict) : a dictionary mapping fragment names to
                                   fragment contents
            min_len (int)        : the minimum length subfragment to consider
                                   optional, defaults to 1

        Returns:
            super_fragment (str) : the shortest possible superstring containing
                                   all of the fragments
    """

    # for the set of fragments as a whole, we know the minimum subfragment
    # length we need to generate is half the length of the smallest fragment
    if min_len > 2:
        min_len = min_len / 2

    right_stubs, left_stubs = generate_overlap_map(fragment_dict, min_len)

    remaining_fragments = fragment_dict.keys()

    # pick a fragment to build off of. it doesn't matter which fragment
    # is picked -- the potential ordering is unique, so by adding on to the
    # left and right sides we can reconstruct the whole fragment
    assembled_fragment_list = [remaining_fragments[0],]

    # to start, both the left and right edges are the same as the initial
    # fragment
    left_edge = remaining_fragments[0]
    right_edge = remaining_fragments[0]

    # remove the seed fragment from the list of remaining
    remaining_fragments.remove(left_edge)


    # starting with the right edge, match and add fragments to the right
    # side of the fragment until you reach the end
    while right_edge != None:

        right_fragment = fragment_dict[right_edge]

        next_fragment_name = match_right(right_fragment, left_stubs,
                                        remaining_fragments)

        assembled_fragment_list.append(next_fragment_name)
        if next_fragment_name is not None:
            remaining_fragments.remove(next_fragment_name)
        right_edge = next_fragment_name

    # do the same thing for the left edge. Match and  add fragments to the left
    # side of the fragment until you reach the end
    while left_edge != None:
        left_fragment = fragment_dict[left_edge]

        next_fragment_name = match_left(left_fragment, right_stubs,
                                        remaining_fragments)

        # prepend to list
        assembled_fragment_list.insert(0, next_fragment_name)

        if next_fragment_name is not None:
            remaining_fragments.remove(next_fragment_name)

        left_edge = next_fragment_name

    # remove leading and trailing 'None's
    clean_list = [x for x in assembled_fragment_list if x is not None]

    # translate in order list of names to in order list of fragments
    in_order_fragments = names_to_fragments(clean_list, fragment_dict)

    # assemble list of fragments into shortest possible superstring
    super_fragment = assemble_fragments(in_order_fragments)

    return super_fragment

def reconstruct_from_file(filename):
    """
        Given a filename that contains fragments in the FASTA format
        ((https://en.wikipedia.org/wiki/FASTA_format), parse the file
        and return the shortest possible superstring that contains all the
        fragments

        Args:
            filename (str) : location of file on disk

        Returns:
            superstring (str) : shortest possible superstring containing all
                                fragments
    """

    fragments, min_len = parse_file(filename)

    superstring = reconstruct(fragments, min_len/2)

    return superstring

def main():
    """
        Main method. Emailed data file is stored in data.txt
    """
    superstring = reconstruct_from_file("./data.txt")

    print superstring

if __name__ == "__main__":

    main()
