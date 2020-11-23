#!/usr/bin/env python3

import os
import re
import unittest
import pytest
import glob
import sys
import logging
import argparse

import tqdm

from make_kreport import main as kreport_main

##
# Test suite for kraken_grouper.py
##


class GrouperTester(unittest.TestCase):
    def setUp(self) -> None:
        # Primitive variables
        self.ex_kraken_line = "C\ttest_id\t562\t52\t562:13 561:4 A:31 0:1 562:3\n"
        self.unclassified_kraken = "U\tc9f4363c-0b0e-4c90-9f56-e153c898f697_3\t0\t500\t0:409 255045:1 0:60"
        self.rm_pattern = "_\d+$"

        # File paths
        self.kraken_tab = "test_data/test.kraken"
        # self.kraken_tab = "/home/connor/Bioinformatics/Koonkie/Calysta/data/barcode01.kraken2"
        self.kreport = "test_data/test.kreport"
        self.kraken_one = "test_data/levelled_group.kraken"
        self.taxonomy = "test_data/ktaxonomy.tsv"
        self.output_prefix = "test_grouped"

        # Class instances
        self.k_match = KrakenMatch()
        self.kraken_matches = []
        test_lines = ["C\tc9f4363c-0b0e-4c90-9f56-e153c898f697_0\t9999998\t500\t0:45 9999998:5 0:155 2420306:1 0:18 1760:5 0:167 9999998:10 0:64",
                      "C\tc9f4363c-0b0e-4c90-9f56-e153c898f697_1\t131567\t500\t0:41 2144175:3 0:252 330879:5 0:72 1705566:5 0:92",
                      "C\tc9f4363c-0b0e-4c90-9f56-e153c898f697_2\t1883\t500\t0:238 1883:5 0:8 2:5 0:188 1330547:1 0:25",
                      "U\tc9f4363c-0b0e-4c90-9f56-e153c898f697_3\t0\t500\t0:409 255045:1 0:60",
                      "C\tc9f4363c-0b0e-4c90-9f56-e153c898f697_4\t1760\t500\t0:159 12916:2 0:21 384:1 0:18 1760:4 147645:1 2:5 0:22 1031538:1 0:143 1620421:3 0:75 1752398:2 0:1 2508168:3 0:9"]
        for kraken_line in test_lines:  # type: str
            km = KrakenMatch()
            km.load_match(kraken_line)
            self.kraken_matches.append(km)
        return

    def tearDown(self) -> None:
        for file_path in glob.glob(self.output_prefix + "*"):
            if os.path.isfile(file_path):
                os.remove(file_path)
        return

    def test_read_kraken_table(self):
        k_matches = read_kraken_table(self.kraken_one)
        self.assertEqual(12, len(k_matches))
        return

    def test_get_taxid_kmer_map(self):
        self.k_match.load_match(self.ex_kraken_line)
        taxid_counts = self.k_match.get_taxid_kmer_map()
        self.assertEqual(4, len(taxid_counts))
        self.assertEqual(16, taxid_counts["562"])
        return

    def test_group_kraken_matches(self):
        grouped_krakens = group_kraken_matches(kraken_matches=read_kraken_table(self.kraken_one),
                                               method="majority", remove_pattern=self.rm_pattern)
        self.assertEqual(1, len(grouped_krakens))
        k_match = grouped_krakens.pop(0)
        self.assertIsInstance(k_match, KrakenMatch)
        self.assertEqual("1079d4fe-4f91-4345-a1b7-9c6c0411bf86", k_match.seq_id)
        self.assertTrue(k_match.classified)
        self.assertEqual("1386", k_match.taxid)
        return

    def test_clean_taxid_counts(self):
        # Ensure '0' was removed
        taxid_counts = {'0': 49, '1228': 51, 'A': 3}
        clean_taxid_counts(taxid_counts)
        self.assertTrue('0' not in taxid_counts)
        self.assertTrue('A' not in taxid_counts)

        # Ensure '0' was retained
        taxid_counts = {'0': 180, '1228': 20}
        clean_taxid_counts(taxid_counts, min_zero_proportion=0.51)
        self.assertTrue('0' in taxid_counts)

        # Ensure nothing is changed
        taxid_counts = {'1': 180, '1228': 20, '454': 1208, '55': 666}
        clean_taxid_counts(taxid_counts)
        self.assertEqual(4, len(taxid_counts))
        return

    def test_find_majority_taxid(self):
        rep = find_majority_taxid({'561': 4, '560': 12, '562': 16, '564': 8, '1': 9})
        self.assertEqual('562', rep)
        return

    def test_find_lca_taxid(self):
        return

    def test_merge_kraken_matches(self):
        base_kraken_match, taxid_counts = merge_kraken_matches(self.kraken_matches)
        self.assertEqual(2282, taxid_counts['0'])
        self.assertEqual(2500, base_kraken_match.seq_len)
        self.assertTrue('' == base_kraken_match.seq_id)
        return

    def test_find_consensus_taxid_from_kraken_group(self):
        # Exit when the list of KrakenMatch instances is empty
        with pytest.raises(SystemExit):
            find_consensus_taxid_from_kraken_group(seq_id="test-seq.id", kraken_matches=[])

        # Ensure the k-mer match counts are being summed properly
        final_km = find_consensus_taxid_from_kraken_group(seq_id="test-seq.id", kraken_matches=self.kraken_matches)
        self.assertEqual('9999998', final_km.taxid)
        self.assertTrue(final_km.classified)
        # Ensure the final taxonomic classification is as expected
        return

    # Integrative test
    def test_kraken_grouper(self):
        kraken_grouper(["--kraken_output", self.kraken_tab,
                        "--taxonomy", self.taxonomy,
                        "--output_prefix", self.output_prefix])
        self.assertTrue(os.path.isfile(self.output_prefix + ".kreport"))
        return

##
# Ends tests
##


class MyFormatter(logging.Formatter):

    error_fmt = "%(levelname)s - %(module)s, line %(lineno)d:\n%(message)s"
    warning_fmt = "%(levelname)s:\n%(message)s"
    debug_fmt = "%(asctime)s\n%(message)s"
    info_fmt = "%(message)s"

    def __init__(self):
        super().__init__(fmt="%(levelname)s: %(message)s",
                         datefmt="%d/%m %H:%M:%S")

    def format(self, record):

        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._style._fmt

        # Replace the original format with one customized by logging level
        if record.levelno == logging.DEBUG:
            self._style._fmt = MyFormatter.debug_fmt

        elif record.levelno == logging.INFO:
            self._style._fmt = MyFormatter.info_fmt

        elif record.levelno == logging.ERROR:
            self._style._fmt = MyFormatter.error_fmt

        elif record.levelno == logging.WARNING:
            self._style._fmt = MyFormatter.warning_fmt

        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)

        # Restore the original format configured by the user
        self._style._fmt = format_orig

        return result


class KrakenMatch:
    def __init__(self):
        self.classified = True
        self.seq_id = ""
        self.taxid = '0'
        self.seq_len = 0
        self.kmer_matches = ""
        return

    def load_match(self, line) -> None:
        """
        Loads the KrakenMatch attributes based off the five tab-delimited fields of a Kraken classification table
        """
        fields = line.strip().split("\t")
        if fields[0] == "U":
            self.classified = False
        self.seq_id, self.taxid, self.seq_len, self.kmer_matches = fields[1:]
        self.seq_len = int(self.seq_len)
        return

    def get_taxid_kmer_map(self) -> dict:
        """
        Split the fifth field in the Kraken classification table into a dictionary,
         summing the number of k-mers attributed to each NCBI taxid

        Description:
        A space-delimited list indicating the LCA mapping of each $k$-mer in the sequence(s).
        For example, "562:13 561:4 A:31 0:1 562:3" would indicate that:

the first 13 $k$-mers mapped to taxonomy ID #562
the next 4 $k$-mers mapped to taxonomy ID #561
the next 31 $k$-mers contained an ambiguous nucleotide
the next $k$-mer was not in the database
the last 3 $k$-mers mapped to taxonomy ID #562
Note that paired read data will contain a "|:|" token in this list to indicate the end of one read and the beginning of another.

        When Kraken 2 is run against a protein database (see [Translated Search]),
         the LCA hitlist will contain the results of querying all six frames of each sequence.
        Reading frame data is separated by a "-:-" token.
        """
        taxid_kmer_map = {}
        for taxid_k_mers in self.kmer_matches.split(' '):
            taxid, count = taxid_k_mers.split(':')
            try:
                taxid_kmer_map[taxid] += int(count)
            except KeyError:
                taxid_kmer_map[taxid] = int(count)
        return taxid_kmer_map

    def tbl_format(self) -> str:
        buffer = []
        if self.classified:
            buffer.append('C')
        else:
            buffer.append('U')
        buffer += [self.seq_id, self.taxid, str(self.seq_len), self.kmer_matches]
        return "\t".join(buffer)


def get_options(sys_args):
    parser = argparse.ArgumentParser(description="A script for finding the"
                                                 " consensus Kraken taxonomic classification across levelled reads.",
                                     add_help=False)

    req_args = parser.add_argument_group("Required arguments")
    opt_args = parser.add_argument_group("Optional arguments")
    mis_args = parser.add_argument_group("Miscellaneous arguments")

    req_args.add_argument("-k", "--kraken_output", dest="k_file", required=True,
                          help="Path to a Kraken classification table, where each read's classification is on a row")
    req_args.add_argument('-t', '--taxonomy', dest='tax_file', required=True,
                          help='Output taxonomy file from make_ktaxonomy.py')

    opt_args.add_argument('-o', '--output_prefix', dest="out_prefix", default='grouped', required=False,
                          help="The path prefix to use, including the directories, "
                               "for writing the consensus Kraken classification files [ DEFAULT = './grouped' ]")
    opt_args.add_argument('-r', "--strip_regex", dest="re_pattern", required=False, default="_\d+$", type=str,
                          help="The regular expression pattern to use when converting split reads to the original"
                               " sequence name. [ DEFAULT = '_\d+$' ]")
    # opt_args.add_argument('--use-read-len', dest='use_read_len',
    #                       action='store_true', default=False, required=False,
    #                       help='Make report file using sum of read lengths [default: read counts]')

    mis_args.add_argument('--overwrite', action='store_true', default=False,
                          help='overwrites previously processed output folders')
    mis_args.add_argument('-v', '--verbose', action='store_true', default=False,
                          help='prints a more verbose runtime log')
    mis_args.add_argument("-h", "--help",
                          action="help", help="show this help message and exit")

    args = parser.parse_args(sys_args)
    return args


def prep_logging(log_file_name=None, verbosity=False) -> None:
    """
    Allows for multiple file handlers to be added to the root logger, but only a single stream handler.
    The new file handlers must be removed outside of this function explicitly

    :param log_file_name:
    :param verbosity:
    :return: None
    """
    if verbosity:
        logging_level = logging.DEBUG
    else:
        logging_level = logging.INFO

    # Detect whether a handlers are already present and return if true
    logger = logging.getLogger()
    if len(logger.handlers):
        return

    formatter = MyFormatter()
    # Set the console handler normally writing to stdout/stderr
    ch = logging.StreamHandler()
    ch.setLevel(logging_level)
    ch.terminator = ''
    ch.setFormatter(formatter)

    if log_file_name:
        output_dir = os.path.dirname(log_file_name)
        try:
            if output_dir and not os.path.isdir(output_dir):
                os.makedirs(output_dir)
        except (IOError, OSError):
            sys.stderr.write("ERROR: Unable to make directory '" + output_dir + "'.\n")
            sys.exit(3)
        logging.basicConfig(level=logging.DEBUG,
                            filename=log_file_name,
                            filemode='w',
                            datefmt="%d/%m %H:%M:%S",
                            format="%(asctime)s %(levelname)s:\n%(message)s")
        logging.getLogger('').addHandler(ch)
        logging.getLogger('').propagate = False
    else:
        logging.basicConfig(level=logging_level,
                            datefmt="%d/%m %H:%M:%S",
                            format="%(asctime)s %(levelname)s:\n%(message)s")
    return


def clean_taxid_counts(taxid_counts, min_zero_proportion=0.99) -> None:
    """ Remove ambiguous and missing k-mers from the dictionary, unless it's the majority """
    try:
        taxid_counts.pop('A')  # Remove ambiguous k-mers
    except KeyError:
        pass
    try:
        if taxid_counts['0']/sum(taxid_counts.values()) <= min_zero_proportion:
            taxid_counts.pop('0')  # Don't include the any k-mers that were not in the database
    except KeyError:
        pass
    return


def find_lca_taxid(taxid_counts: dict, alpha=0.51) -> str:
    # TODO: Implement majority-rule LCA from the component KrakenMatch k-mer matches
    return ""


def find_majority_taxid(taxid_counts: dict) -> str:
    return max(taxid_counts, key=lambda k: taxid_counts[k])


def merge_kraken_matches(kraken_matches: list) -> (KrakenMatch, dict):
    base_kraken_match = KrakenMatch()
    taxid_counts = {}
    while kraken_matches:
        k_match = kraken_matches.pop()  # type: KrakenMatch
        # Add the seq_length
        base_kraken_match.seq_len += int(k_match.seq_len)
        # Add the k-mer matches to the string
        base_kraken_match.kmer_matches += ' ' + k_match.kmer_matches
        for k, v in k_match.get_taxid_kmer_map().items():  # type: (str, int)
            try:
                taxid_counts[k] += v
            except KeyError:
                taxid_counts[k] = v
    return base_kraken_match, taxid_counts


def find_consensus_taxid_from_kraken_group(seq_id: str, kraken_matches: list, method="majority") -> KrakenMatch:
    """
    When reads are split into subsequences and classified individually by Kraken, they should be reconstituted into the
    original sequence and their individual k-mer classification should be consolidated.

    This is the main function for regrouping/consolidating the individual k-mer matches for subsequences.
    It operates on a list of KrakenMatch instances, all from the same read/supersequence, and groups the k-mer matches
    into a single dictionary mapping NCBI taxids to counts.
    This is then used to find the most likely taxonomic assignment, either through a simple majority-rule or by an LCA
    voting method.

    :param seq_id: The common supersequence name of all subsequences
    :param kraken_matches: A list of KrakenMatch instances derived from the same supersequence
    :param method: The name of the method to use to group the sub-KrakenMatch k-mer classifications,
     either 'majority' or 'lca'
    :return: A KrakenMatch instance that represents the whole sequence's Kraken classification
    """
    if len(kraken_matches) == 0:
        logging.error("No subsequences were found for read {}.\n".format(seq_id))
        sys.exit(7)

    base_kraken_match, taxid_counts = merge_kraken_matches(kraken_matches)
    base_kraken_match.seq_id = seq_id

    # Filter out irrelevant taxid codes
    clean_taxid_counts(taxid_counts)

    if len(taxid_counts) == 0:
        logging.debug("No taxid counts remained in read '{}' after removing ambiguous and missing taxid's.\n"
                      "".format(seq_id))
        base_kraken_match.taxid = '0'
        base_kraken_match.classified = False
        return base_kraken_match

    # Determine the final taxonomic classification of the whole sequence
    if method == "majority":
        base_kraken_match.taxid = find_majority_taxid(taxid_counts)
    elif method == "lca":
        base_kraken_match.taxid = find_lca_taxid(taxid_counts)
    else:
        logging.error("Unknown method ('{}') requested for grouping kraken matches of the same molecule.\n"
                      "".format(method))
        sys.exit(3)

    # Check whether the sequence was classified or not
    if base_kraken_match.taxid == '0':
        base_kraken_match.classified = False

    return base_kraken_match


def group_kraken_matches(kraken_matches: list, remove_pattern: str, method="majority") -> list:
    logging.info("Grouping kraken read classifications by sequence name...")

    # Set up the progress bar
    pbar = tqdm.tqdm(ncols=100)
    pbar.total = len(kraken_matches)

    # Ensure the method for merging/grouping the KrakenMatch instances is recognized
    if method not in ["majority", "lca"]:
        logging.error("Unknown method ('{}') requested for grouping kraken matches of the same molecule.\n"
                      "".format(method))
        sys.exit(3)

    grouped_matches = []
    group = []
    previous_seq_id = ""
    strip_re = re.compile(r"({})".format(remove_pattern))
    kraken_matches = sorted(kraken_matches, key=lambda x: x.seq_id)
    while kraken_matches:
        k_match = kraken_matches.pop(0)  # type: KrakenMatch
        if not strip_re.search(k_match.seq_id):
            logging.error("Sequence name '{}' did not match the provided regular expression pattern.\n")
            sys.exit(5)

        group_seq_id = strip_re.sub('', k_match.seq_id)
        if previous_seq_id and group_seq_id != previous_seq_id:
            grouped_matches.append(find_consensus_taxid_from_kraken_group(previous_seq_id, group, method))
            group = [k_match]
        else:
            group.append(k_match)
        previous_seq_id = group_seq_id
        pbar.update(1)

    # Group the last of the KrakenMatch instances
    grouped_matches.append(find_consensus_taxid_from_kraken_group(previous_seq_id, group, method))

    pbar.close()

    return grouped_matches


def read_kraken_table(k_path: str) -> list:
    logging.info("Reading kraken read-level classification table %s...\n" % k_path)
    sys.stdout.flush()

    # Save counts per taxid
    read_count = 0
    kraken_matches = []
    sys.stdout.write("\t%i million reads processed" % read_count)
    k_handler = open(k_path, 'r')
    for line in k_handler:
        k_match = KrakenMatch()
        k_match.load_match(line)
        kraken_matches.append(k_match)

        read_count += 1
        if read_count % 1000 == 0:
            sys.stdout.write('\r\t%0.3f million reads processed' % float(read_count / 1000000.))
            sys.stdout.flush()
    k_handler.close()
    sys.stdout.write('\r\t%0.3f million reads processed\n' % float(read_count / 1000000.))
    sys.stdout.flush()

    return kraken_matches


def write_kraken_table(kraken_matches: list, output_file: str) -> None:
    try:
        out_handler = open(output_file, 'w')
    except IOError:
        logging.error("Unable to write output file '{}'. Does the output directory exist?\n".format(output_file))
        sys.exit(7)

    buffer = []
    for k_match in kraken_matches:  # type: KrakenMatch
        buffer.append(k_match.tbl_format())
        if len(buffer) >= 1E4:
            out_handler.write("\n".join(buffer) + "\n")
            buffer.clear()

    out_handler.write("\n".join(buffer) + "\n")
    out_handler.close()
    return


def kraken_grouper(sys_args):
    args = get_options(sys_args)
    prep_logging()

    kraken_match_insts = read_kraken_table(k_path=args.k_file)
    grouped_matches = group_kraken_matches(kraken_matches=kraken_match_insts,
                                           remove_pattern=args.re_pattern)

    # Generating the output files
    kraken_file = args.out_prefix + ".kraken"
    kreport_file = args.out_prefix + ".kreport"
    write_kraken_table(grouped_matches, kraken_file)

    # Call main() in make_kreport.py
    kreport_main(["--kraken", kraken_file,
                  "--output", kreport_file,
                  '--taxonomy', args.tax_file,
                  "--use-read-len"])

    logging.info("Kraken-grouper completed successfully.\n")

    return


if __name__ == "__main__":
    kraken_grouper(sys.argv[1:])
