import itertools
import numpy as np
import pandas as pd


def generate_csv(input_data, labels, output_file):
    with open(input_data) as f:
        lines = f.readlines()
    lines = map(lambda line: line.rstrip('\n'), lines)
    lines = map(lambda line: line.split(' '), lines)
    lines = map(lambda line: list(filter(lambda chars: len(chars) != 0, line)), lines)
    lines = map(lambda line: list(map(lambda chars: tuple(map(lambda char: int(char), chars.split(':'))), line)), lines)
    flat_lines = list(itertools.chain.from_iterable(lines))

    cols = len(lines)
    rows = max(flat_lines)[0]+1

    arr = np.zeros((rows, cols))
    for col, line in enumerate(lines):
        for row, val in line:
            arr[row, col] = val

    labels = pd.read_csv(labels, header=None)

    label = []
    for idx, val in labels.iterrows():
        if (val == 1).bool():
            label.append('ca.%d' % (idx + 1))
        else:
            label.append('non_ca.%d' % (idx + 1))

    df = pd.DataFrame(arr)
    df.columns = label

    df.to_csv(output_file, index=False)


generate_csv('dexter_train.data', 'dexter_train.labels', 'dexter_train.csv')
generate_csv('dexter_valid.data', 'dexter_valid.labels', 'dexter_valid.csv')
