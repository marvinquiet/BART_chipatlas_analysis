import os

def sort_bed(input_dir, output_dir):
    for peak in os.listdir(input_dir):
        input_peak = os.path.join(input_dir, peak)
        output_peak = os.path.join(output_dir, 'sorted_'+peak)

        cmd = "sort -k1,1 -k2,2n {} > {}".format(input_peak, output_peak)
        os.system(cmd)

if __name__ == '__main__':
    input_dir = "data/TFs_peaks_per_srx"
    output_dir = "data/sorted_TFs_peaks_per_srx"
    os.makedirs(output_dir, exist_ok=True)

    sort_bed(input_dir, output_dir)
