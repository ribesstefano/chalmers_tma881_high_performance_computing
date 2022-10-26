'''
The analysis contained in this script has been carried out on the following
machine hardware:
# ==============================================================================
>>> lscpu
# ==============================================================================
Architecture:            x86_64
  CPU op-mode(s):        32-bit, 64-bit
  Address sizes:         46 bits physical, 48 bits virtual
  Byte Order:            Little Endian
CPU(s):                  112
  On-line CPU(s) list:   0-111
Vendor ID:               GenuineIntel
  Model name:            Intel(R) Xeon(R) Platinum 8180 CPU @ 2.50GHz
    CPU family:          6
    Model:               85
    Thread(s) per core:  2
    Core(s) per socket:  28
    Socket(s):           2
    Stepping:            4
    BogoMIPS:            5000.00
    Flags:               fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush dts acpi mmx fxsr sse sse2 ss ht tm pbe syscall nx pdpe1gb rdtscp lm constant_tsc art arch_perfmon pebs bt
                         s rep_good nopl xtopology nonstop_tsc cpuid aperfmperf pni pclmulqdq dtes64 monitor ds_cpl vmx smx est tm2 ssse3 sdbg fma cx16 xtpr pdcm pcid dca sse4_1 sse4_2 x2apic movbe popcnt tsc_de
                         adline_timer aes xsave avx f16c rdrand lahf_lm abm 3dnowprefetch cpuid_fault epb cat_l3 cdp_l3 invpcid_single pti intel_ppin ssbd mba ibrs ibpb stibp tpr_shadow vnmi flexpriority ept vpi
                         d ept_ad fsgsbase tsc_adjust bmi1 hle avx2 smep bmi2 erms invpcid rtm cqm mpx rdt_a avx512f avx512dq rdseed adx smap clflushopt clwb intel_pt avx512cd avx512bw avx512vl xsaveopt xsavec x
                         getbv1 xsaves cqm_llc cqm_occup_llc cqm_mbm_total cqm_mbm_local dtherm ida arat pln pts pku ospke md_clear flush_l1d arch_capabilities
Virtualization features: 
  Virtualization:        VT-x
Caches (sum of all):     
  L1d:                   1.8 MiB (56 instances) 1887436.8 B
  L1i:                   1.8 MiB (56 instances) 1887436.8 B
  L2:                    56 MiB (56 instances) 58720256 B
  L3:                    77 MiB (2 instances) 80740352 B
NUMA:                    
  NUMA node(s):          2
  NUMA node0 CPU(s):     0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110
  NUMA node1 CPU(s):     1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,75,77,79,81,83,85,87,89,91,93,95,97,99,101,103,105,107,109,111
Vulnerabilities:         
  Itlb multihit:         KVM: Mitigation: VMX disabled
  L1tf:                  Mitigation; PTE Inversion; VMX conditional cache flushes, SMT vulnerable
  Mds:                   Mitigation; Clear CPU buffers; SMT vulnerable
  Meltdown:              Mitigation; PTI
  Mmio stale data:       Mitigation; Clear CPU buffers; SMT vulnerable
  Retbleed:              Mitigation; IBRS
  Spec store bypass:     Mitigation; Speculative Store Bypass disabled via prctl
  Spectre v1:            Mitigation; usercopy/swapgs barriers and __user pointer sanitization
  Spectre v2:            Mitigation; IBRS, IBPB conditional, RSB filling, PBRSB-eIBRS Not affected
  Srbds:                 Not affected
  Tsx async abort:       Mitigation; Clear CPU buffers; SMT vulnerable
# ==============================================================================
>>> lscpu --caches
# ==============================================================================
NAME ONE-SIZE ALL-SIZE WAYS TYPE        LEVEL  SETS PHY-LINE COHERENCY-SIZE
L1d       32K     1.8M    8 Data            1    64        1             64
L1i       32K     1.8M    8 Instruction     1    64        1             64
L2         1M      56M   16 Unified         2  1024        1             64
L3      38.5M      77M   11 Unified         3 57344        1             64
'''
import json
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
from functools import reduce
from itertools import product

# NOTE: The tests have been carried out running 20 parallel threads
num_threads = 20
elem_byte_size = 2
# NOTE: Physical cores
num_cores = 56
num_sockets = 2
threads_per_core = 2
threads_core_utilization = threads_per_core / num_threads
# Individual cache sizes in Bytes
L1 = 1887436.8 / num_cores
L2 = 58720256 / num_cores
L3 = 80740352 / num_sockets
# ==============================================================================
# Generate table of results
# ==============================================================================
with open('bench.json', 'r') as fp:
    data = json.load(fp)
table = []
col_names = ['Block Size X',
             'Block Size Y',
             'Execution Time [s]',
             'Size [KB]',
             'L1 Utilization [%]',
             'L2 Utilization [%]',
             'L3 Utilization [%]',
             'Block Dimensions',
            ]
for k in data.keys():
    block_x = data[k][0]
    block_y = data[k][1]
    exec_t = data[k][2]
    block_bytes = block_x * block_y * elem_byte_size
    bytes_per_threads = block_bytes * threads_core_utilization
    # NOTE: We only round for reporting purposes
    L1_utilization = round(bytes_per_threads / L1 * 100, 2)
    L2_utilization = round(bytes_per_threads / L2 * 100, 2)
    L3_utilization = round(bytes_per_threads / L3 * 100, 2)
    entry = (block_x,
             block_y,
             exec_t,
             round(block_bytes / 1024, 4),
             L1_utilization,
             L2_utilization,
             L3_utilization,
             f'{block_x} x {block_y}',
            )
    table.append(entry)
df = pd.DataFrame(table, columns=col_names)
# df = pd.read_json('bench_df.json')
df.to_json('bench_df.json', indent=4)
# ==============================================================================
# Bar plots
# ==============================================================================
# ax = sns.barplot(data=df, x='Execution Time [s]', y='Block Dimensions', orient='h')
ax = sns.barplot(data=df, x='Execution Time [s]', y='L2 Utilization [%]', orient='h')
# TODO: The following is only available for matplotlib>=3.4.0
# for container in ax.containers:
#     ax.bar_label(container, padding=10)
plt.grid(axis='x')
plt.title(f'Execution time vs. Block Dimensions ({num_threads} Threads)')
plt.tight_layout()
plt.savefig(f'l2_utilization_vs_exec_time_{num_threads}threads.pdf')
plt.show()
# ==============================================================================
# Ideal block size
# ==============================================================================
ideal_block_size = int(L2 / threads_core_utilization)
print(f'Ideal block size: {ideal_block_size} [B]')

def factors(n):    
    return set(reduce(list.__add__, 
                ([i, n // i] for i in range(1, int(n**0.5) + 1) if n % i == 0)))

fct = sorted(list(factors(ideal_block_size)))

cnt = 1
for x, y in product(fct, fct):
    if x * y == ideal_block_size:
        print(f'Best block dimensions candidate n.{cnt:3d} for {num_threads} threads: {x:10d} x {y:10d}')
        cnt += 1
# ==============================================================================
# Scatter plots
# ==============================================================================
# ax = sns.scatterplot(data=df, x='Block Size X', y='Block Size Y', hue='Execution Time [ms]', size='Execution Time [ms]', palette='vlag')
# ax = sns.scatterplot(data=df, x='Block Size X', y='Block Size Y', hue='Size [B]', size='Size [B]', palette='vlag')

# fig, ax = plt.subplots()
# plt.scatter(x_val, y_val, s=t_val)

# # ax.legend(bbox_to_anchor=(1.04, 0.5), loc='center left', borderaxespad=0)
# # ax.set(xscale='log', yscale='log')
# ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=5)

# # Show colorbar
# # norm = plt.Normalize(df['Execution Time [ms]'].min(), df['Execution Time [ms]'].max())
# norm = plt.Normalize(df['Size [B]'].min(), df['Size [B]'].max())
# sm = plt.cm.ScalarMappable(cmap='vlag', norm=norm)
# sm.set_array([])
# ax.figure.colorbar(sm)
# plt.scatter(x, y, s=t, label='Execution time [us]')
# plt.legend()