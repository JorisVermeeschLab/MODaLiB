def find_microsatellite(seq, min_motif_length=2, max_motif_length=5):
    n = len(seq)

    counts = list(range(min_motif_length, max_motif_length + 1))
    best_i = 0
    best_n_repeats = 1
    best_motif_length = 1

    for i in range(n):
        for k in range(max_motif_length - min_motif_length + 1):
            motif_length = k + min_motif_length
            if seq[i] == seq[i - motif_length]:
                counts[k] += 1
            else:
                counts[k] = motif_length

            n_repeats = int(counts[k] / motif_length)
            if n_repeats >= 3:
                if n_repeats * motif_length > best_n_repeats * best_motif_length:
                    best_n_repeats = n_repeats
                    best_motif_length = motif_length
                    best_i = i
    
    start = best_i + 1 - best_motif_length * best_n_repeats
    end = best_i
    print('MOTIF: ', seq[best_i+1-best_motif_length:best_i+1])
    print('REPEATS: ', best_n_repeats)
    print('START: ', start)
    print('END: ', end)
    print('SEQ: ', seq[start:end+1])

    # TODO: extend left and right and account for mutations


input_sequence = "AATCTACTTATCTCTATCTATCTATCTATCTATCTCTCTTTTTTCTCCAAACTCTTACT"
print(input_sequence)
find_microsatellite(input_sequence)
