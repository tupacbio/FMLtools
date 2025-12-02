#!/usr/bin/env python3

# offset = how many bases in the 3' direction the first base of the read should be from the 5' end of the motif
# e.g. for LpnPI and FspEI:
# forward strand site, forward strand read: site is external to the fragment, so the read should start at +14 from the beginning of the motif
# forward strand site, reverse strand read: site is internal to the fragment, so the read should start at +17 from the from the beginning of the motif
# reverse strand: opposite
# MspJI, raw CpG: +13, +16

import collections, os, click, pysam
import genomerator as gen

DEFAULT_SAME_STRAND_OFFSET = 13  # for MspJI and generic CpG
DEFAULT_OPPOSITE_STRAND_OFFSET = 16  # for MspJI and generic CpG
DEFAULT_WOBBLE = 1  # documented for MspJI


@click.command()
@click.option(
    "-s",
    "--same_strand_offset",
    type=int,
    default=DEFAULT_SAME_STRAND_OFFSET,
    help="how many bases in the 3' direction the first base of the read should be from the 5' end of the motif, if the read is on the same strand as the target region (default %i; use 14 for LpnPI/FspEI motifs)"
    % DEFAULT_SAME_STRAND_OFFSET,
)
@click.option(
    "-o",
    "--opposite_strand_offset",
    type=int,
    default=DEFAULT_OPPOSITE_STRAND_OFFSET,
    help="3'-oriented offset if read and target region are on opposite strands (default %i; use 16 for LpnPI/FspEI motifs)"
    % DEFAULT_OPPOSITE_STRAND_OFFSET,
)
@click.option(
    "-w",
    "--wobble",
    type=int,
    default=DEFAULT_WOBBLE,
    help="how far downstream from the expected position the cut site may 'wobble', opposite strand only (default %i)"
    % DEFAULT_WOBBLE,
)
@click.option(
    "-z", "--zeroes", is_flag=True, help="output all motif counts including zeroes"
)
@click.option("-q", "--quiet", is_flag=True, help="don't show progress bar")
@click.argument("region_bed", type=click.File("r"))
@click.argument("alignment_bam", required=False, default="-")
@click.argument(
    "out_bed", required=False, type=click.File("w"), default=click.open_file("-", "w")
)
def count_methyl_hits(
    region_bed,
    out_bed,
    alignment_bam,
    same_strand_offset,
    opposite_strand_offset,
    wobble,
    zeroes,
    quiet,
):
    """
    Given a BAM file of alignments and a BED file of target regions (restriction motif sites), count the number of alignments that start a specified distance from target regions (correct digestions).
    """

    sam_file = pysam.Samfile(alignment_bam, "rb")
    alignment_stream = gen.SamStream(
        filter(
            lambda alignment: not alignment.is_duplicate
            and not alignment.is_secondary
            and not alignment.is_supplementary,
            sam_file,
        ),
        references=sam_file.references,
        assert_sorted=True,
    )

    motif_stream = (
        gen.GenomeFeature.from_genomefeature(x, data=0)
        for x in gen.BedStream(source=region_bed, references=sam_file.references)
    )

    motif_counter = gen.OperationGenerator(
        a=motif_stream,
        b=alignment_stream,
        match=lambda motif, alignment: (
            motif.is_reverse == alignment.is_reverse
            and motif.start_offset(alignment) == same_strand_offset
        )
        or (
            motif.is_reverse != alignment.is_reverse
            and opposite_strand_offset
            <= motif.start_offset(alignment)
            <= opposite_strand_offset + wobble
        ),
        a_is_passed=lambda motif, alignment: alignment
        > motif.right + max(same_strand_offset, opposite_strand_offset),
        b_is_passed=lambda motif, alignment: motif
        > alignment.right + max(same_strand_offset, opposite_strand_offset),
    )

    with click.progressbar(
        length=1, file=(open(os.devnull, "w") if quiet else None)
    ) as bar:
        last_progress = 0
        for motif_count in motif_counter:
            if zeroes or motif_count.data > 0:
                bed_line = motif_count.bed(
                    reference_names=sam_file.references, name=motif_count.data
                )
                fields = bed_line.strip().split("\t")

                current_start = int(fields[1])

                if motif_count.is_reverse:
                    # Reverse strand (-): Trim first 2 bases by increasing start by 2
                    fields[1] = str(current_start + 2)

                else:
                    # Forward strand (+): Set End to be Start + 2
                    fields[2] = str(current_start + 2)

                # Write the modified line
                out_bed.write("\t".join(fields) + "\n")

            if not quiet:
                new_progress = motif_count.progress(sam_file.lengths)
                bar.update(new_progress - last_progress)
                last_progress = new_progress
        if not quiet:
            bar.update(1)

    click.echo("%i total reads" % alignment_stream.count, err=True)
    click.echo("%i aligned to motifs" % motif_counter.count_b_hits, err=True)


if __name__ == "__main__":
    count_methyl_hits()

