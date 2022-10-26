import argparse, os, subprocess, shutil, json
from itertools import product

# scp benchmark_cell_distances.py hpcuser124@gantenbein.math.chalmers.se:chalmers_tma881_high_performance_computing/assignment_2/

# TEST_DATA_DIR = '/home/stefanor/Downloads/chalmers_tma881_high_performance_computing/reference_code/cell_distances/test_data/'
TEST_DATA_DIR = '/home/hpc2022/cell_distances/test_data'
devnull = subprocess.DEVNULL

def check_runtime(num_threads, cell_size, curr_path, time_bound_ms, repeat):
    bench = 'bench.csv'
    warmup = 1
    time_bound_total = 2 * time_bound_ms * (warmup + repeat)

    # Link test files from base to current directory
    if os.path.exists(curr_path + '/cells'):
        os.remove(curr_path + '/cells')
    os.symlink(TEST_DATA_DIR + f'/cell_e{cell_size}', curr_path + '/cells')
    # Run program
    cell_distances_cmd = f'./cell_distances -t{num_threads}'
    cmd = f'hyperfine --export-csv {bench} --time-unit millisecond '
    cmd += f'--warmup {warmup} --max-runs {repeat} \"{cell_distances_cmd}\"'
    if time_bound_total > 0:
        try:
            subprocess.Popen(cmd, shell=True, cwd=curr_path, stdout=devnull, stderr=devnull).wait(time_bound_total)
        except subprocess.TimeoutExpired:
            print(f'TOO SLOW FOR {cell_distances_cmd} ON cell_e{cell_size}')
            return
    else:
        subprocess.Popen(cmd, shell=True, cwd=curr_path, stdout=devnull, stderr=devnull).wait()
    # Reading file
    try:
        with open(os.path.join(curr_path, bench), 'r') as f:
            time = float(f.readlines()[-1].split(',')[-2])
            if time > time_bound_ms and time_bound_ms > 0:
                print(f'TOO SLOW ({time} ms) FOR {cell_distances_cmd} ON FILE cell_e{cell_size}')
                return
            else:
                return time
    except:
        print(f'RUNTIME ERROR FOR {cell_distances_cmd} ON FILE cell_e{cell_size}')

def check_runtimes(path):
    return (     check_runtime( 1, 4, path, 0.26, 10)
                 and check_runtime( 5, 5, path, 10.0,    5)
                 and check_runtime(10, 5, path, 5.3, 10)
                 and check_runtime(20, 5, path, 2.8, 10)
                 )

def check_correctness(distreffiles, distfile):
    distrefs = [dict() for f in distreffiles]
    dist = dict()
    for (distreffile,distref) in zip(distreffiles,distrefs):
        for l in open(distreffile, 'r').readlines():
            l = l.split(' ')
            distref[l[0].strip()] = int(l[1])
    for l in open(distfile, 'r').readlines():
        l = l.split(' ')
        dist[l[0].strip()] = int(l[1])
    invtolerance = 50
    abstolerance = 3
    first_wrong = [None for _ in distrefs]
    for (dx,distref) in enumerate(distrefs):
        for k in distref.keys():
            if ( ( k not in dist and distref[k] > abstolerance ) or 
                     ( k in dist and abs(dist[k] - distref[k]) * invtolerance > distref[k] and abs(dist[k] - distref[k]) > abstolerance ) ) :
                first_wrong[dx] = k
                break
    if any(k is None for k in first_wrong):
        return True
    print('NONE OF THE REFERENCE SOLUTIONS MATCHES:')
    for (k,distref,distreffile) in zip(first_wrong,distrefs,distreffiles):
        print('WHEN COMPARING WITH {}'.format(distreffile))
        print('WRONG VALUE AT {}'.format(k))
        print('DETECTED VALUE IS {}'.format(dist[k] if k in dist else 'absent'))
        print('SHOULD BE WITHIN 2%-RELATIVE OR 3%-ABSOLUTE DISTANCE FROM EXPECTED {}'.format(distref[k]))

    return False

def run_check_build(working_dir):
    # Check files format
    is_valid_file = lambda f: (
             f in ['makefile', 'Makefile'] or
             f.endswith('.cc') or f.endswith('.c') or
             f.endswith('.hh') or f.endswith('.h') )
    is_valid_folder = lambda folder: (
        all( all(map(is_valid_file, files))
                 for (root, _, files) in os.walk(folder) ) )

    print('checking for additional files after cleaning...')
    if not is_valid_folder(working_dir):
        print('ADDITIONAL FILES IN TAR')
        return False

    # build clean build
    print('bulding and cleaning...')
    subprocess.Popen(['make', 'cell_distances'], cwd=working_dir).wait()
    subprocess.Popen(['make', 'clean'], cwd=working_dir).wait()
    if not is_valid_folder(working_dir):
        print('ADDITIONAL FILES AFTER BUILD CLEAN')
        return False

    print('bulding...')
    subprocess.Popen(['make', 'cell_distances'], cwd=working_dir, stdout=devnull).wait()
    return True

def run_check_correctness():
    print('checking correctness...')
    if os.path.exists(working_dir + '/cells'):
        os.remove(working_dir + '/cells')
    os.symlink('/home/stefanor/Downloads/chalmers_tma881_high_performance_computing/reference_code/cell_distances/test_data/cell_e5', working_dir + '/cells')
    dist_file = 'distances/' + stem
    with open(dist_file, 'w') as f:
        subprocess.Popen(['./cell_distances', '-t8'], cwd=working_dir, stdout = f).wait(2000)
    if os.stat(dist_file).st_size == 0:
        print('PROGRAM DOES NOT WRITE TO STDOUT')
        return False
    else:
        return check_correctness(['/home/stefanor/Downloads/chalmers_tma881_high_performance_computing/reference_code/cell_distances/test_data/dist_e5',
                                                            '/home/stefanor/Downloads/chalmers_tma881_high_performance_computing/reference_code/cell_distances/test_data/dist_e5_alt'],
                                                            dist_file)

def build_program(block_size_x=2048, block_size_y=2048, working_dir='benchmarking'):
    subprocess.run(['mkdir', '-p', working_dir])
    # Build program
    subprocess.run(['cp', '-rf', 'Makefile', working_dir])
    subprocess.run(['cp', '-rf', 'cell_distances.c', working_dir])
    macros = f'BLOCK_DIMS=-DBLOCK_SIZE_X={block_size_x} -DBLOCK_SIZE_Y={block_size_y}'
    subprocess.Popen(['make', '-B', 'cell_distances', macros], cwd=working_dir, stdout=devnull).wait()

def main():
    # Create test directory
    working_dir = 'benchmarking'
    benchmark_filename = 'bench.json'
    # # Run default benchmarks
    # block_size_x = 2048
    # block_size_y = 2048
    # build_program(block_size_x, block_size_y, working_dir)
    # time_ms = check_runtime(num_threads=1, cell_size=4, curr_path=working_dir, time_bound_ms=0.26, repeat=10)
    # time_ms = check_runtime(num_threads=5, cell_size=5, curr_path=working_dir, time_bound_ms=10.0, repeat=5)
    # time_ms = check_runtime(num_threads=10, cell_size=5, curr_path=working_dir, time_bound_ms=5.3, repeat=10)
    # time_ms = check_runtime(num_threads=20, cell_size=5, curr_path=working_dir, time_bound_ms=2.8, repeat=10)

    # Run benchmarks
    benchmark_results = {}
    processed_combs = {}
    pow_2s = [int(2**x) for x in range(7, 11)]
    pow_2s = list(range(int(2**7), int(2**11), 32))

    for x, y in product(pow_2s, pow_2s):
        print(f'RUNNING ({x}, {y})...')
        build_program(block_size_x=x, block_size_y=y, working_dir=working_dir)
        time_ms = check_runtime(num_threads=20, cell_size=5, curr_path=working_dir, time_bound_ms=0.0, repeat=5)
        print(f'({x}, {y}) = {time_ms} ms')
        benchmark_results[f'{x},{y}'] = (x, y, time_ms)
        # if (x, y) not in processed_combs:
        #   run benchmarcks
        # processed_combs[(x, y)] = None
        # processed_combs[(y, x)] = None

    with open(os.path.join(working_dir, benchmark_filename), 'w') as fp:
        json.dump(benchmark_results, fp, indent=4)

if __name__ == '__main__':
    main()