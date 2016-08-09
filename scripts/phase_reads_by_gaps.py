import argparse
from collections import defaultdict
import csv
import networkx as nx
import operator
import pprint
import pybedtools
import pysam
import sys

#import matplotlib as mpl
# mpl.use('Agg')
# import matplotlib.pyplot as plt
# from matplotlib import pylab

import logging

# create logger
logger = logging.getLogger("read_phaser")
logger.setLevel(logging.DEBUG)

# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)

# Constants related to zero-based indexing of fields from the SV call format.
QUERY_START = 10
QUERY_END = 11
QUERY_LENGTH = 13


def get_node_size(node):
    return int(node[4])


def find_consensus_calls(graph, name_column, read_graph):
    """
    Find all connected subgraphs such that calls with only self-self overlaps
    will exist as singleton graphs while all nodes that overlap each other
    directly or transitively will be clustered in the same graph.
    """
    for subgraph in nx.connected_component_subgraphs(graph):
        # Collect all nodes in this group by their start positions.
        number_of_nodes = len(subgraph.nodes())
        if number_of_nodes > 1:
            # Add an edge between all read pairs in this subgraph.
            for i in xrange(number_of_nodes):
                read_name_for_i = subgraph.nodes()[i][name_column]
                for j in xrange(number_of_nodes):
                    read_name_for_j = subgraph.nodes()[j][name_column]

                    if i != j:
                        # If an edge already exists between two reads, increment
                        # the edge weight by one. Otherwise, add the edge.
                        if read_graph.has_edge(read_name_for_i, read_name_for_j):
                            read_graph[read_name_for_i][read_name_for_j]["weight"] += 1
                        else:
                            read_graph.add_edge(read_name_for_i, read_name_for_j, weight=1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("reads", help="BAM of reads to phase or FOFN file with list of BAM paths to use")
    parser.add_argument("gaps", help="tabixed BED file of gaps found in reads")
    parser.add_argument("region", help="region to phase by gaps")
    parser.add_argument("phased_prefix", help="prefix for one or more BAM files with phased reads")
    parser.add_argument("--min_proportion_of_phased_reads", type=float, default=0.25, help="minimum proportion of total reads in the largest phased subgraph")
    parser.add_argument("--reciprocal_overlap", type=float, default=0.5, help="proportion of reciprocal overlap required to consider two calls the same")
    parser.add_argument("--name_column", type=int, default=None, help="column to use as the name field")
    parser.add_argument("--type_column", type=int, default=None, help="column to use as the type of event")
    args = parser.parse_args()

    if args.reads.endswith(".fofn"):
        with open(args.reads, "r") as fh:
            bams = [pysam.AlignmentFile(filename.rstrip(), "rb") for filename in fh]
    else:
        bams = [pysam.AlignmentFile(args.reads, "rb")]

    gaps = pysam.TabixFile(args.gaps)
    read_graph = nx.Graph()

    # Load all calls from the given BED file. Loading directly into a BedTool
    # from the filename is not guaranteed to work if the input BED file has the
    # same number of columns as a standard bedN format (e.g., bed12) and any of
    # those columns contain values of a type that differs from the standard.
    region_reads = dict([(read.query_name, read)
                         for bam in bams
                         for read in bam.fetch(region=args.region)])
    logger.debug("Found %i reads in region %s", len(region_reads), args.region)
    calls = pybedtools.BedTool([record for record in gaps.fetch(region=args.region)])
    logger.debug("Found %i gaps in region %s", len(calls), args.region)

    if len(calls) > 0:
        columns_per_call = len(list(calls[0]))

        # Intersect the given calls with themselves. Self-self overlaps will be
        # reported along with any other overlaps matching the given reciprocal
        # overlap proportion.
        intersected_calls = calls.intersect(b=calls, f=args.reciprocal_overlap, r=True, wao=True, sorted=True)
        logger.debug("Found %i intersected calls in region %s", len(intersected_calls), args.region)

        # Create a graph connecting all calls that share a reciprocal overlap.
        current_contig = None
        for call in intersected_calls:
            if current_contig != call[0]:
                # If this isn't the first pass through the calls and we've found a
                # new contig, find the consensus calls and print them.
                if not current_contig is None:
                    find_consensus_calls(graph, args.name_column, read_graph)

                # If we've switched to a different contig, create a new graph for
                # that contig.
                current_contig = call[0]
                graph = nx.Graph()

            left = tuple(call[:columns_per_call])

            # Omit the final column that contains the total bases overlapping
            # between inputs.
            right = tuple(call[columns_per_call:-1])

            # Only keep edges if the overlap is of the same event type.
            if left[args.type_column] == right[args.type_column]:
                graph.add_edge(left, right)
        else:
            # If we've finished processing the last call, get consensus calls for
            # the final contig's graph.
            find_consensus_calls(graph, args.name_column, read_graph)

        # # Prune edges from the read graph with a weight corresponding to fewer
        # # than 2 events supporting a link between reads (weight < 4).
        # print("Pruning read graph")
        # MIN_WEIGHT = 2
        # pruned_edges = 0
        # for edge in read_graph.edges():
        #     if read_graph[edge[0]][edge[1]]["weight"] < MIN_WEIGHT:
        #         read_graph.remove_edge(*edge)
        #         pruned_edges += 1

        # print ("Pruned %i edges" % pruned_edges)

        print("Region reads: %i" % len(region_reads))
        max_subgraph_reads = 0
        max_subgraph = None
        for subgraph in nx.connected_component_subgraphs(read_graph):
            if len(subgraph.nodes()) > max_subgraph_reads:
                max_subgraph_reads = len(subgraph.nodes())
                max_subgraph = subgraph

            print("Subgraph with %i nodes" % len(subgraph.nodes()))

        print("Max subgraph reads: %i" % max_subgraph_reads)
        other_reads = set(region_reads) - set(max_subgraph.nodes())
        print("Other reads: %i" % len(other_reads))
        # nx.draw_spring(max_subgraph)
        # plt.savefig("read_graph.png")

        # If there are enough reads in the largest subgraph to warrant a
        # separate assembly, write out two separate sets of reads. Otherwise,
        # write out the original set of reads.
        if (len(max_subgraph.nodes()) / float(len(region_reads))) >= args.min_proportion_of_phased_reads:
            print("Writing out phased reads")

            phase_set1 = pysam.AlignmentFile("%s_1.bam" % args.phased_prefix, "wb", template=bams[0])
            phase_set2 = pysam.AlignmentFile("%s_2.bam" % args.phased_prefix, "wb", template=bams[0])

            # Write out phased reads from the largest subgraph.
            for read in max_subgraph.nodes():
                phase_set1.write(region_reads[read])

            # Write out all other reads which belong to the second group.
            for read in other_reads:
                phase_set2.write(region_reads[read])

            phase_set1.close()
            phase_set2.close()
        else:
            # Write out all reads to one file.
            print("Writing out original reads")
            phase_set1 = pysam.AlignmentFile("%s_1.bam" % args.phased_prefix, "wb", template=bams[0])
            for read in region_reads:
                phase_set1.write(region_reads[read])

            phase_set1.close()
