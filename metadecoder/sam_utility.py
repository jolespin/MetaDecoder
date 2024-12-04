import os
import gzip
import struct
from math import ceil

def read_sam_header(input_sam):
    sorted_sam = False
    open_file = open(input_sam, 'r')
    for line in open_file:
        if line.startswith('@'):
            if line.startswith('@SQ'):
                for tag_value in line.rstrip('\n').split('\t'):
                    if tag_value.startswith('SN'):
                        sequence_id = tag_value[3 : ]
                    elif tag_value.startswith('LN'):
                        sequence_length = int(tag_value[3 : ])
                yield (sequence_id, sequence_length)
            elif line.startswith('@HD') and ('SO:coordinate' in line):
                sorted_sam = True
        else:
            break
    open_file.close()
    # assert sorted_sam, f"The sam file \"{input_sam}\" must be sorted by coordinate." #
    return None


def generate_block(input_alignment, blocks):
    '''
    Split a sam file into some parts.
    Parameters:
        input_sam: the input sam file (with @SQ headers).
        blocks: the number of blocks.
    Return:
        a generator of start and end positions of a block of a sam file.
    '''
    file_size = os.path.getsize(input_alignment)

    if input_alignment.endswith(".bam"):
        with open(input_alignment, "rb") as open_file:
            # Skip the header by moving to the first alignment
            magic = open_file.read(4)
            if magic != b"BAM\1":
                raise ValueError("Not a valid BAM file")

            l_text = struct.unpack("<i", open_file.read(4))[0]
            open_file.seek(l_text, 1)  # Skip header text
            n_ref = struct.unpack("<i", open_file.read(4))[0]

            # Skip reference sequence dictionary
            for _ in range(n_ref):
                l_name = struct.unpack("<i", open_file.read(4))[0]
                open_file.seek(l_name + 4, 1)

            block_start = open_file.tell()
            block_size = ceil((file_size - block_start) / blocks)

            while block_start < file_size:
                open_file.seek(block_start + block_size, 0)
                open_file.readline()  # Align to the next record
                block_end = open_file.tell()
                yield block_start, min(block_end, file_size)
                block_start = block_end
    else:

        block_start = 0

        open_file = open(input_alignment, 'rb')
        for line in open_file:
            if not line.startswith(b'@'):
                block_start = open_file.tell() - len(line)
                break
        # block_start has been defined #

        block_size = ceil((file_size - block_start) / blocks)
        while block_start < file_size:
            open_file.seek(block_size, 1)
            open_file.readline()
            block_end = open_file.tell()
            yield(block_start, min(block_end, file_size))
            block_start = block_end
    open_file.close()
    return None


def read_sam_file(input_sam, block_start, block_end, mapq = 0):
    '''
    Parameters:
        input_sam: the input sam file (with @SQ headers).
        block_start: the start position of a block.
        block_end: the end position of a block.
        mapq: MAPQ.
    Return:
        (read_id, ref_id, pos, cigar)
    '''

    open_file = open(input_sam, 'rb')
    open_file.seek(block_start, 0)
    for line in open_file:
        lines = line.rstrip(b'\n').split(b'\t')
        if (~ int(lines[1]) & 4) and (int(lines[4]) >= mapq):
            # read id, reference sequence, position, cigar #
            yield (lines[0].decode('ascii'), lines[2].decode('ascii'), int(lines[3]), lines[5].decode('ascii'))
        block_start += len(line)
        if block_start >= block_end:
            break
    open_file.close()
    return None



def read_bam_header(input_bam):
    """
    Reads the BAM header and sequence dictionary.
    Parameters:
        input_bam: the input BAM file path.
    Yields:
        sequence_id and sequence_length for each reference sequence in the header.
    """
    with gzip.open(input_bam, "rb") as open_file:
        magic = open_file.read(4)
        if magic != b"BAM\1":
            raise ValueError("Not a valid BAM file")

        # Read header text length
        l_text = struct.unpack("<i", open_file.read(4))[0]
        header_text = open_file.read(l_text).decode("utf-8")
        print("Header Text:", header_text)

        # Number of reference sequences
        n_ref = struct.unpack("<i", open_file.read(4))[0]

        # Parse each reference sequence
        for _ in range(n_ref):
            l_name = struct.unpack("<i", open_file.read(4))[0]
            name = open_file.read(l_name).decode("utf-8").rstrip("\x00")
            l_ref = struct.unpack("<i", open_file.read(4))[0]
            yield name, l_ref



def read_bam_file(input_bam, block_start, block_end, mapq=0):
    """
    Reads alignment data from a BAM block.
    Parameters:
        input_bam: the input BAM file (binary).
        block_start: the start position of a block.
        block_end: the end position of a block.
        mapq: minimum MAPQ value to filter.
    Yields:
        (read_id, ref_id, pos, cigar) for each alignment.
    """
    with open(input_bam, "rb") as open_file:
        open_file.seek(block_start, 0)

        while open_file.tell() < block_end:
            block_size_data = open_file.read(4)
            if not block_size_data:
                break
            block_size = struct.unpack("<i", block_size_data)[0]

            # Read alignment block
            alignment_data = open_file.read(block_size)
            ref_id, pos, _, mapq_val, _, l_read_name = struct.unpack("<iiBBBi", alignment_data[:16])
            read_name = alignment_data[16:16 + l_read_name - 1].decode("ascii")
            cigar = alignment_data[16 + l_read_name:].decode("ascii")

            if mapq_val >= mapq:
                yield read_name, ref_id, pos, cigar

