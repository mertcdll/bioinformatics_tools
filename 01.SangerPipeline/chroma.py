import matplotlib.pyplot as plt
from cfutils.parser import parse_abi
from cfutils.show import plot_chromatograph

seq = parse_abi("/home/genwork2/Mert/ab1fastq/input/ALPARSLAN_ALIK_REZTEP_SMARCD1_EX12_13F_D09.ab1")
peaks = seq.annotations["peak positions"]

fig, axes = plt.subplots(2, 1, figsize=(12, 6), sharex=True)
plot_chromatograph(
    seq,
    region=(0, len(peaks)),
    ax=axes[0],
    show_bases=True,
    show_positions=True,
    color_map=dict(zip("ATGC", ["C0", "C2", "C1", "C4"])),
)

# Save the plot as a PNG file
plt.savefig('/home/genwork2/Mert/ab1fastq/input/your_plot.png', dpi=3500, bbox_inches='tight')

plt.show()


peaks