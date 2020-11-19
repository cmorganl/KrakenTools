#!/usr/bin/env python3

import os
import re
import unittest
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
        self.ex_kraken_line = "U\ttest_id\t562\t52\t562:13 561:4 A:31 0:1 562:3\n"
        self.rm_pattern = "_\d+$"

        # File paths
        self.kraken_tab = "test_data/test.kraken"
        self.kreport = "test_data/test.kreport"
        self.kraken_one = "test_data/levelled_group.kraken"
        self.taxonomy = "test_data/ktaxonomy.tsv"

        # Class instances
        self.k_match = KrakenMatch()
        return

    def tearDown(self) -> None:
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

    def test_find_lca_taxid(self):
        return

    def test_kraken_grouper(self):
        kraken_grouper(["--kraken_output", self.kraken_tab,
                        "--taxonomy", self.taxonomy,
                        "--output_prefix", "test_grouped"])
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
        self.taxid = 0
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


def find_lca_taxid(taxid_counts: dict, alpha=0.51) -> str:
    # TODO: Implement majority-rule LCA from the component KrakenMatch k-mer matches
    return ""


def find_majority_taxid(taxid_counts: dict) -> str:
    return max(taxid_counts, key=lambda k: taxid_counts[k])


def find_consensus_taxid_from_kraken_group(seq_id: str, kraken_matches: list, method="majority") -> KrakenMatch:
    base_kraken_match = KrakenMatch()
    base_kraken_match.seq_id = seq_id
    taxid_counts = {}
    classified = False
    while kraken_matches:
        k_match = kraken_matches.pop()  # type: KrakenMatch
        # Check that the subsequence was classified
        if k_match.classified:
            classified = True
        # Add the seq_length
        base_kraken_match.seq_len += int(k_match.seq_len)
        # Add the k-mer matches to the string
        base_kraken_match.kmer_matches += ' ' + k_match.kmer_matches
        taxid_counts.update(k_match.get_taxid_kmer_map())

    base_kraken_match.classified = classified
    # Filter out irrelevant taxid codes
    taxid_counts.pop('0')  # Don't include the any k-mers that were not in the database
    try:
        taxid_counts.pop('A')  # Remove ambiguous k-mers
    except KeyError:
        pass
    if method == "majority":
        base_kraken_match.taxid = find_majority_taxid(taxid_counts)
    elif method == "lca":
        base_kraken_match.taxid = find_lca_taxid(taxid_counts)
    else:
        logging.error("Unknown method ('{}') requested for grouping kraken matches of the same molecule.\n"
                      "".format(method))
        sys.exit(3)

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
