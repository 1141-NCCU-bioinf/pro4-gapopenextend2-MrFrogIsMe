import pandas as pd

def load_pam_matrix(file_path):
    '''Load PAM scoring matrix from a file into a pandas DataFrame.'''
    pam_df = pd.read_csv(file_path, sep='\s+', comment='#', index_col=0)
    return pam_df

def load_fasta(file_path):
    labels, sequences = [], []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                labels.append(line.strip())
                sequences.append('')
            else:
                sequences[-1] += line.strip()
    return labels, sequences

def init_dp(seq1, seq2, aln, gap_open, gap_extend):
    NEG_INF = float('-inf')
    # 三個矩陣：M (match/mismatch), X (seq1 gap), Y (seq2 gap)
    dp_M = [[NEG_INF for _ in range(len(seq1) + 1)] for _ in range(len(seq2) + 1)]
    dp_X = [[NEG_INF for _ in range(len(seq1) + 1)] for _ in range(len(seq2) + 1)]
    dp_Y = [[NEG_INF for _ in range(len(seq1) + 1)] for _ in range(len(seq2) + 1)]

    # 三個方向矩陣
    dir_M = [[None for _ in range(len(seq1) + 1)] for _ in range(len(seq2) + 1)]
    dir_X = [[None for _ in range(len(seq1) + 1)] for _ in range(len(seq2) + 1)]
    dir_Y = [[None for _ in range(len(seq1) + 1)] for _ in range(len(seq2) + 1)]

    if aln == 'local':
        for i in range(len(seq2) + 1):
            dp_M[i][0] = dp_X[i][0] = dp_Y[i][0] = 0
        for j in range(len(seq1) + 1):
            dp_M[0][j] = dp_X[0][j] = dp_Y[0][j] = 0
        return (dp_M, dp_X, dp_Y), (dir_M, dir_X, dir_Y)

    # Global alignment initialization
    dp_M[0][0] = 0
    # X: Horizontal gap (seq1 的字元對到 seq2 的 gap)
    for j in range(1, len(seq1) + 1):
        dp_X[0][j] = gap_open + (j - 1) * gap_extend
        dir_X[0][j] = ('X', 0, j - 1)
    # Y: Vertical gap (seq2 的字元對到 seq1 的 gap)
    for i in range(1, len(seq2) + 1):
        dp_Y[i][0] = gap_open + (i - 1) * gap_extend
        dir_Y[i][0] = ('Y', i - 1, 0)

    return (dp_M, dp_X, dp_Y), (dir_M, dir_X, dir_Y)

def alignment(input_path, score_path, output_path, aln, gap_open, gap_extend):
    (label1, label2), (seq1, seq2) = load_fasta(input_path)
    scores = load_pam_matrix(score_path)
    (dp_M, dp_X, dp_Y), (dir_M, dir_X, dir_Y) = init_dp(seq1, seq2, aln, gap_open, gap_extend)

    NEG_INF = float('-inf')
    max_score = NEG_INF
    max_poses = []

    for i in range(1, len(seq2) + 1):
        for j in range(1, len(seq1) + 1):
            c1 = seq1[j-1]
            c2 = seq2[i-1]
            match_score = int(scores.loc[c2, c1])

            # M[i][j]: match/mismatch，從對角線來
            m_candidates = {
                ('M', i-1, j-1): dp_M[i-1][j-1] + match_score,
                ('X', i-1, j-1): dp_X[i-1][j-1] + match_score,
                ('Y', i-1, j-1): dp_Y[i-1][j-1] + match_score,
            }
            if aln == 'local':
                m_candidates[(None, None, None)] = 0
            dp_M[i][j] = max(m_candidates.values())
            dir_M[i][j] = max(m_candidates, key=m_candidates.get)
            if dir_M[i][j] == (None, None, None):
                dir_M[i][j] = None

            # X[i][j]: seq1 gap (水平移動，j-1 -> j)
            x_candidates = {
                ('M', i, j-1): dp_M[i][j-1] + gap_open,
                ('X', i, j-1): dp_X[i][j-1] + gap_extend,
                ('Y', i, j-1): dp_Y[i][j-1] + gap_open,
            }
            if aln == 'local':
                x_candidates[(None, None, None)] = 0
            dp_X[i][j] = max(x_candidates.values())
            dir_X[i][j] = max(x_candidates, key=x_candidates.get)
            if dir_X[i][j] == (None, None, None):
                dir_X[i][j] = None

            # Y[i][j]: seq2 gap (垂直移動，i-1 -> i)
            y_candidates = {
                ('M', i-1, j): dp_M[i-1][j] + gap_open,
                ('X', i-1, j): dp_X[i-1][j] + gap_open,
                ('Y', i-1, j): dp_Y[i-1][j] + gap_extend,
            }
            if aln == 'local':
                y_candidates[(None, None, None)] = 0
            dp_Y[i][j] = max(y_candidates.values())
            dir_Y[i][j] = max(y_candidates, key=y_candidates.get)
            if dir_Y[i][j] == (None, None, None):
                dir_Y[i][j] = None

            # Track max for local alignment
            if aln == 'local':
                cell_max = max(dp_M[i][j], dp_X[i][j], dp_Y[i][j])
                if cell_max > max_score:
                    max_score = cell_max
                    max_poses = []
                if cell_max == max_score and cell_max > 0:
                    if dp_M[i][j] == cell_max:
                        max_poses.append(('M', i, j))
                    if dp_X[i][j] == cell_max:
                        max_poses.append(('X', i, j))
                    if dp_Y[i][j] == cell_max:
                        max_poses.append(('Y', i, j))

    # Trace-back
    def trace_back(aln1, aln2, matrix, i, j, res):
        if matrix == 'M':
            direction = dir_M[i][j]
        elif matrix == 'X':
            direction = dir_X[i][j]
        else:
            direction = dir_Y[i][j]

        if direction is None:
            if i == 0 and j == 0:
                res.append((aln1, aln2))
            elif aln == 'local':
                res.append((aln1, aln2))
            return

        prev_matrix, pi, pj = direction

        if pi == i - 1 and pj == j - 1:
            # Diagonal: match/mismatch
            trace_back(seq1[j-1] + aln1, seq2[i-1] + aln2, prev_matrix, pi, pj, res)
        elif pi == i and pj == j - 1:
            # Horizontal: gap in seq2
            trace_back(seq1[j-1] + aln1, '-' + aln2, prev_matrix, pi, pj, res)
        elif pi == i - 1 and pj == j:
            # Vertical: gap in seq1
            trace_back('-' + aln1, seq2[i-1] + aln2, prev_matrix, pi, pj, res)

    res = []
    if aln == 'global':
        i, j = len(seq2), len(seq1)
        # 從三個矩陣中選最佳的結尾
        final_scores = [('M', dp_M[i][j]), ('X', dp_X[i][j]), ('Y', dp_Y[i][j])]
        best_matrix = max(final_scores, key=lambda x: x[1])[0]
        trace_back('', '', best_matrix, i, j, res)
    elif aln == 'local':
        for matrix, i, j in max_poses:
            trace_back('', '', matrix, i, j, res)
        if res:
            max_len = max(len(a[0]) for a in res)
            res = [a for a in res if len(a[0]) == max_len]

    res.sort(key=lambda x: (x[0], x[1]))
    with open(output_path, 'w') as f:
        for aln1, aln2 in res:
            f.write(label1 + '\n')
            f.write(aln1 + '\n')
            f.write(label2 + '\n')
            f.write(aln2 + '\n')
            if aln == 'global':
                break

if __name__ == "__main__":
    alignment("examples/test.fasta", "examples/pam250.txt", "examples/my_result_global.fasta", "global", -10, -2)
    alignment("examples/test.fasta", "examples/pam250.txt", "examples/my_result_local.fasta", "local", -10, -2)
