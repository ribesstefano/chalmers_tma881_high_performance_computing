import argparse
import random

def get_rand_coord():
    rand_num = random.uniform(-10, 10)
    pos = True if rand_num > 0 else False
    rand_num_repr = f'{abs(rand_num):05.3f}'
    if rand_num < 10:
        rand_num_repr = '0' + rand_num_repr
    if pos:
        rand_num_repr = '+' + rand_num_repr
    else:
        rand_num_repr = '-' + rand_num_repr
    return rand_num_repr[:7]

def main():
    parser = argparse.ArgumentParser(description='Generate random cells.')
    parser.add_argument('-n', '--num_cells', metavar='N', type=int,
                        dest='num_cells', default=1000,
                        help='Number of cells to generate. Default: 1000')
    parser.add_argument('-f', '--filename', dest='filename',
                        default='cells.txt',
                        help='Output filename. Default: cells.txt')
    args = parser.parse_args()

    with open(args.filename, 'w') as f:
        for _ in range(args.num_cells):
            line = f'{get_rand_coord()} {get_rand_coord()} {get_rand_coord()}\n'
            f.write(line)


if __name__ == '__main__':
    main()