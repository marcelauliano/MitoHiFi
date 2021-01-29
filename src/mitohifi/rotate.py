#!/usr/bin/env python

import argparse
import gzip
import os
import sys

from Bio import SeqIO
from Bio.Seq import Seq


def parseArgs():
    parser = argparse.ArgumentParser(
        description="rotates contigs (or " "scaffolds) individually", add_help=False
    )
    req = parser.add_argument_group("Required")
    req.add_argument(
        "-i",
        "--infile",
        metavar="FILE",
        required=True,
        help="input FastA file, optionally gunzip compressed",
    )
    opt = parser.add_argument_group("Optional")
    opt.add_argument(
        "-h", "--help", action="help", help="show this help message and exit"
    )
    opt.add_argument(
        "-o",
        "--outfile",
        default=None,
        metavar="FILE",
        help="output rotated FastA file [stdout]",
    )
    opt.add_argument(
        "-r",
        "--rotate",
        type=int,
        default=50000,
        metavar="INT",
        help="rotation size (in bp); negative values turn counterclockwise "
        "from 3' to 5' [50,000]",
    )
    opt.add_argument(
        "--force",
        default=False,
        action="store_true",
        help="force rotations for all sequence records",
    )
    opt.add_argument(
        "--skip-failures",
        default=False,
        action="store_true",
        help="skip records that are not longer than the rotation size; "
        "output all records that can be rotated",
    )
    return parser.parse_args()


def main():
    args = parseArgs()
    infile = os.path.abspath(os.path.expanduser(args.infile))
    if infile.endswith(".gz"):
        infile = gzip.open(infile)
    rotation_size = args.rotate

    # Identify rotation direction and make sure it's not zero bp
    if rotation_size > 0:
        sys.stderr.write(
            "INFO: rotating each sequence record {} bp "
            "clockwise...\n\n".format(rotation_size)
        )
    elif rotation_size < 0:
        sys.stderr.write(
            "INFO: rotating each sequence record {} bp "
            "counterclockwise...\n\n".format(abs(rotation_size))
        )
    else:
        sys.stderr.write(
            "ERROR: the number of base pairs to rotate each "
            "contig must be an integer greater than or less than 0.\n"
        )
        sys.exit(1)

    # Verify input sequences exceed rotation size
    records_lengths = []
    for record in SeqIO.parse(infile, "fasta"):
        records_lengths.append(len(record.seq))
    if any(x <= abs(rotation_size) for x in records_lengths):
        if args.force or args.skip_failures:
            sys.stderr.write(
                "INFO: at least one input record is not longer "
                f"than {abs(rotation_size)} bp rotation size.\n\n"
            )
        else:
            sys.stderr.write(
                "ERROR: at least one input record is not longer "
                f"than {abs(rotation_size)} bp. To forcibly rotate anyways, re-run with --force "
                "or specify a smaller rotation size.\n\n"
            )
            sys.exit(1)

    # Rotate each sequence
    rotated_records = []
    for record in SeqIO.parse(infile, "fasta"):
        init_seq, seq_length = str(record.seq), len(record.seq)
        if seq_length > abs(rotation_size):
            rotd_seq = init_seq[rotation_size:] + init_seq[:rotation_size]
            record.seq = Seq(rotd_seq)
            rotated_records.append(record)
            sys.stderr.write(
                "INFO: {} sequence ({} bp) was rotated {} bp.\n".format(
                    record.id, seq_length, abs(rotation_size)
                )
            )
        else:
            if args.force:
                spin_size = rotation_size % len(record.seq)
                rotated_seq = init_seq[spin_size:] + init_seq[:spin_size]
                record.seq = Seq(rotated_seq)
                rotated_records.append(record)
                sys.stderr.write(
                    "INFO: {} sequence ({} bp) is not longer "
                    "than the specified {} bp rotation size. It was rotated "
                    "{} bp.\n".format(
                        record.id, seq_length, abs(rotation_size), spin_size
                    )
                )
            elif args.skip_failures:
                continue
            else:
                sys.stderr.write(
                    f"INFO: {record.id} sequence ({seq_length} bp) is not longer "
                    f"than the specified {abs(rotation_size)} bp rotation size. To forcibly "
                    "rotate anyways, re-run with --force or specify a "
                    "smaller rotation size.\n"
                )
                sys.exit(1)

    # Output records
    if args.outfile is None:
        ofh = sys.stdout
    else:
        ofh = os.path.abspath(os.path.expanduser(args.outfile))
    SeqIO.write(rotated_records, ofh, "fasta")


if __name__ == "__main__":
    main()
