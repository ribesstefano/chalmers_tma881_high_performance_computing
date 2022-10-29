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
                        default='cells',
                        help='Output filename. Default: cells')
    parser.add_argument('-d', '--default', dest='default',
                        action='store_true',
                        default=False,
                        help='Generate default cells text. Default: False')
    args = parser.parse_args()

    with open(args.filename, 'w') as f:
        if args.default:
            default_cells = '''+01.330 -09.035 +03.489
-03.718 +02.517 -05.995
+09.568 -03.464 +02.645
-09.620 +09.279 +08.828
+07.630 -02.290 +00.679
+04.113 -03.399 +05.299
-00.994 +07.313 -06.523
+03.376 -03.614 -06.657
+01.304 +09.381 -01.559
-04.238 -07.514 +08.942

'''
            f.write(default_cells)
        else:
            for _ in range(args.num_cells):
                line = f'{get_rand_coord()} {get_rand_coord()} {get_rand_coord()}\n'
                f.write(line)


if __name__ == '__main__':
    main()