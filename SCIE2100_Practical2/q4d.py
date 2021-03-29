def findstartcodon(seq, reading_frame, direction=True):
    if direction is True:
        i = reading_frame - 1

        while i < len(seq):
            if seq[i:i + 3] == 'ATG':
                return i + 1
            else:
                i = i + 3
    else:
        back_seq = seq[::-1]
        i = reading_frame - 1
        while i < len(back_seq):
            if back_seq[i:i + 3] == 'ATG':
                return len(seq) - i
            else:
                i = i + 3
    return 0