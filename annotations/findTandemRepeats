#!/usr/bin/env python
"""Find tandem repeats in a genome"""

from __future__ import division

import argparse
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)
signal.signal(signal.SIGINT, signal.SIG_DFL)
collections = None
heapq=None
pysam=None
args=None

"""Usage"""
argparser = argparse.ArgumentParser("""Find tandem repeats in a genome""",
epilog="""
OVERVIEW
    Find tandem repeats in a genome by local kmer matching.

    The kmer of length `--k` at a position X in the genome is matched against
    all kmers from X+1 to X+`--maxPeriod`.  A kmer match at X+P is considered
    support for a tandem repeat run of period P from X to X+P.

    Nearby kmer matches at similar period (within 10% of P) are merged into
    the tandem repeat run.  When no kmer match is seen for 100 bp, the tandem
    repeat run is ended.

    Runs with at least `--minMatchBp` basepairs, `--minMatchFraction` fraction
    of possibly matching basepairs in kmer matches, that span at least
    `--minLength` basepairs of the genome, and that involve at least `--minCopies`
    copies of the repeat unit are output, in sorted order.
""", formatter_class=argparse.RawDescriptionHelpFormatter)


argparser.add_argument("reffa", metavar="ref.fa")
argparser.add_argument("outbed", metavar="out.bed")
argparser.add_argument("--chrom", metavar="chr", help="Limit analysis to this chromosome")
argparser.add_argument("--k", metavar="N", type=int,
    help="kmer size to use for repeat matching (default: %(default)s)", default=8)
argparser.add_argument("--maxPeriod", metavar="N", type=int,
    help="maximum repeat period size to consider, bp (default: %(default)s)", default=2000)
argparser.add_argument("--minCopies", metavar="F", type=float,
    help="minimum copies of tandem repeat (default: %(default)s)", default=1.5)
argparser.add_argument("--minLength", metavar="N", type=int,
    help="minimum length of a tandem repeat run, bp (default: %(default)s)", default=50)
argparser.add_argument("--minMatchBp", metavar="N", type=int,
    help="minimum basepairs in kmer matches in the tandem repeat run (default: %(default)s)", default=30)
argparser.add_argument("--minMatchFraction", metavar="F", type=float,
    help="minimum fraction of basepairs in kmer matches in the tandem repeat run (default: %(default)s)", default=0.25)
argparser.add_argument("--merge", action="store_true", help="merge overlapping tandem repeats")

class TrBedRecord:
    """A "tandem repeat run" from chrom:chromStart-chromEnd with a specified period.
       The score is the number of kmer matches that support the run."""
    def __init__(self, chrom, chromStart, chromEnd, period):
        self.chrom = chrom
        self.chromStart = chromStart
        self.chromEnd = chromEnd
        self.period = period
        self.score = args.k

    def __lt__(self, rhs):
        return (self.chrom, self.chromStart, self.chromEnd, self.period) < (
            rhs.chrom, rhs.chromStart, rhs.chromEnd, rhs.period)


class BedMergingWriter():
    """A writer that merges overlapping BED records.
    The records must be provided in sorted order.
    """
    def __init__(self, writer):
        self.writer = writer
        self.activeChrom = None
        self.activeChromStart = None
        self.activeChromEnd = None

    def write(self, line):
        # write a new BED line
        chrom,chromStart,chromEnd = line.split()[:3]
        chromStart = int(chromStart)
        chromEnd = int(chromEnd)

        # If this line does not overlap the active record then
        # write the active record
        if chrom != self.activeChrom or chromStart > self.activeChromEnd:
            if self.activeChrom is not None:
                self.writer.write("%s\t%d\t%d\n" % (self.activeChrom, self.activeChromStart, self.activeChromEnd))
            self.activeChrom = chrom
            self.activeChromStart = chromStart
            self.activeChromEnd = chromEnd
        self.activeChromEnd = max(chromEnd, self.activeChromEnd) # extend the active record


    def close(self):
        # Write the final record
        if self.activeChrom is not None:
           self.writer.write("%s\t%d\t%d\n" % (self.activeChrom, self.activeChromStart, self.activeChromEnd))
        self.writer.close()


class SortedWriter():
    """A writer that turns almost-sorted data into fully-sorted data.
    Each line is assigned a key, and lines are written in order by key.
    """
    def __init__(self, writer):
        self.writer = writer
        self.active = [] # heap of active lines to write as tuples: (key, line)
        self.ix = 0      # serial index to maintain a stable sort


    def write(self, key, line):
        # add line to the active heap
        heapq.heappush(self.active, (key, self.ix, line))
        self.ix = self.ix + 1


    def flush(self, flushkey):
        # Flush any records with key < flushkey.
        while self.active and self.active[0][0] < flushkey:
            key, ix, line = heapq.heappop(self.active)
            self.writer.write("%s" % (line))


    def close(self):
        # Write remaining lines.
        while self.active:
            key, ix, line = heapq.heappop(self.active)
            self.writer.write("%s" % (line))
        self.writer.close()


class TrBedMerger:
    """Merge tandem repeat (evidence) with similar periods into tandem repeat calls.
    Write tandem repeats to the provided file, in sorted order."""
    def __init__(self, fn):
        writer = open(fn, "w")
        if args.merge:
            writer = BedMergingWriter(writer)
        self._sortedwriter = SortedWriter(writer)
        self._activeTrs = list()
        self._chromIx = dict() # map from chrom name to serial index
        self._prevTrBedKey = (-1,-1)


    def _filterorwrite(self, trBed):
        """Write a tandem repeat record, if it passes filter criteria.
        The tandem repeats are not written in perfect sorted order by chrom,chromStart; so,
        a sorted writer is used to provide sorted output."""
        L = trBed.chromEnd - trBed.chromStart # length of tandem repeat run
        P = trBed.period # period of repeat
        if (L >= args.minLength and (L/P) >= args.minCopies and
           trBed.score >= args.minMatchBp and (trBed.score/max(P,(L-P)) >= args.minMatchFraction)):
            self._sortedwriter.write((self._chromIx[trBed.chrom], trBed.chromStart),
                "%s\t%d\t%d\tcopies=%0.1f,period=%d,score=%d\n" % (trBed.chrom, trBed.chromStart, trBed.chromEnd, L/P, trBed.period, trBed.score))


    def _kmergap(self, bed):
        """Acceptable gap between adjacent matching kmers in a repeat."""
        return 100


    def addTrBed(self, trBed):
        """Merge a new (putative) tandem repeat.
        The trBed objects are provided in order by chrom,chromStart."""
        # Assign the chromosome an index if it has not been seen before.
        if trBed.chrom not in self._chromIx:
            self._chromIx[trBed.chrom] = len(self._chromIx)

        trBedKey = (self._chromIx[trBed.chrom], trBed.chromStart)
        assert(self._prevTrBedKey <= trBedKey)
        self._prevTrBedKey = trBedKey

        # Flush self._activeTrs of intervals that can not be extended by a later kmer match.
        activeheap = []
        for activeTr in sorted(self._activeTrs):
            if activeTr.chrom != trBed.chrom:
                self._filterorwrite(activeTr)
            elif (activeTr.chromEnd + self._kmergap(activeTr)) < (trBed.chromStart + activeTr.period):
                self._filterorwrite(activeTr)
            else: # retain activeTr as active
                heapq.heappush(activeheap, activeTr)
        self._activeTrs = activeheap # reset the active heap

        # Merge the new tandem repeat interval with any active TR with a similar period (within 10%),
        # or add it as a new tandem repeat.
        merged = False
        for activeTr in self._activeTrs:
            kmergap = self._kmergap(activeTr)
            if ((trBed.chromEnd - activeTr.chromEnd) <= kmergap) and (
                max(activeTr.period, trBed.period) / min(activeTr.period, trBed.period) < 1.1):
                if (trBed.chromEnd > activeTr.chromEnd):
                    activeTr.score += min(args.k, trBed.chromEnd - activeTr.chromEnd) # bases in kmer matches
                    activeTr.chromEnd = trBed.chromEnd
                merged = True
        if not merged:
            heapq.heappush(self._activeTrs, trBed)

        # No future TR will start before the head of the heap.
        self._sortedwriter.flush((self._chromIx[self._activeTrs[0].chrom], self._activeTrs[0].chromStart))


    def close(self):
        # Flush remaining items
        for activeTr in sorted(self._activeTrs):
            self._filterorwrite(activeTr)
        # Clear the active list.
        self._activeTrs = []
        # Close the sorted writer to flush it.
        self._sortedwriter.close()


def main():
    global args; args = argparser.parse_args()
    global collections; import collections
    global heapq; import heapq
    global pysam; import pysam

    # Open input files
    reffa = pysam.Fastafile(args.reffa)

    # Open output file.  Use the "TrBedMerger" object
    # to merge individual kmer links at a specified
    # period into tandem repeat runs.
    trBedMerger = TrBedMerger(args.outbed)
    k = args.k

    # Iterate over chromosomes in order.
    for chrom in reffa.references:
        if not args.chrom or chrom == args.chrom:
            seq = reffa.fetch(chrom).upper()
            chromSize = len(seq)
            if chromSize < k:
                continue

            # Map from kmer to a list of positions where the kmer occurs
            kmerToPos = collections.defaultdict(collections.deque)

            # Prepare the first window on the chromosome
            for i in range(0,min(args.maxPeriod, chromSize)):
                kmerToPos[seq[i:i+k]].append(i)

            # Look for tandem repeats that start at each position on the chromosome.
            # Find "tandem repeats" as kmer matches offset from [1,args.maxPeriod].
            for i in range(0,chromSize-k):
                chromStart = i
                chromStartKmer = seq[chromStart:chromStart+k]
                kmerToPos[chromStartKmer].popleft() # do not count a kmer matching itself
                if "N" not in chromStartKmer and kmerToPos[chromStartKmer]:
                    for ix,pos in enumerate(kmerToPos[chromStartKmer]):
                        # Consider only the first three period sizes to avoid duplicate
                        # processing of all multiples of fundamental period.
                        if ix == 3:
                            break
                        trBedMerger.addTrBed(TrBedRecord(chrom, chromStart, k+pos, pos-chromStart))

                # record the next kmer that will be in the window for the next position
                nextKmerStart = i + args.maxPeriod
                if nextKmerStart + k < chromSize:
                    kmerToPos[seq[nextKmerStart:nextKmerStart+k]].append(nextKmerStart)

    # Close files
    reffa.close()
    trBedMerger.close()


if __name__ == "__main__":
    main()
