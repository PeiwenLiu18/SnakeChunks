def fastq_stats(fastq_file):
    read_count = 0
    total_bases = 0
    gc_count    = 0
    line_number = 0
    with open(fastq_file) as f:
        for line in f:
            line = line.rstrip()
            line_number += 1

            if line_number % 4 == 1:
                read_count += 1
            elif line_number % 4 == 2:
                total_bases += len(line)
                gc_count += line.count('G')
                gc_count += line.count('C')

    avg_read_length = total_bases / read_count
    gc_rate = (gc_count/total_bases)*100

    return (read_count, avg_read_length, gc_rate)

##def peaks_bed_stats(bed_file):


