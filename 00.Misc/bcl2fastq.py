import os
import argparse

class BCLConverter:

    def __init__(self, **kwargs):
        self.runfolder_dir = kwargs.get('runfolder_dir', None)
        self.input_dir = kwargs.get('input_dir', None) 
        self.intensities_dir = kwargs.get('intensities_dir', None)
        self.interop_dir = kwargs.get('interop_dir', None)
        self.stats_dir = kwargs.get('stats_dir', None)
        self.reports_dir = kwargs.get('reports_dir', None)
        self.output_dir = kwargs.get('output_dir', None)
        self.sample_sheet_path = kwargs.get('sample_sheet_path', None)
        self.loading_thread = kwargs.get('loading_thread', None)
        self.processing_thread = kwargs.get('processing_thread', None)
        self.writing_thread = kwargs.get('writing_thread', None)
        self.fastq_compression_level = kwargs.get('fastq_compression_level', None)
        self.no_lane_splitting = kwargs.get('no_lane_splitting', None)
        self.barcode_mismatches = kwargs.get('barcode_mismatches', None)
        self.use_bases_mask = kwargs.get('use_bases_mask', None)
        self.tiles = kwargs.get('tiles', None)
        self.min_log_level = kwargs.get('min_log_level', None)
        self.minimum_trimmed_read_length = kwargs.get('minimum_trimmed_read_length', None)
        self.mask_short_adapter_reads = kwargs.get('mask_short_adapter_reads', None)
        self.adapter_stringency = kwargs.get('adapter_stringency', None)
        self.ignore_missing_bcls = kwargs.get('ignore_missing_bcls', None)
        self.ignore_missing_filter = kwargs.get('ignore_missing_filter', None)
        self.ignore_missing_positions = kwargs.get('ignore_missing_positions', None)
        self.ignore_missing_controls = kwargs.get('ignore_missing_controls', None)
        self.write_fastq_reverse_complement = kwargs.get('write_fastq_reverse_complement', None)
        self.with_failed_reads = kwargs.get("with_failed_reads", None)
        self.create_fastq_for_index_reads = kwargs.get("create_fastq_for_index_reads", None)
        self.find_adapters_with_sliding_window = kwargs.get("find_adapters_with_sliding_window", None)
        self.no_bgzf_compression = kwargs.get("no_bgzf_compression", None)



    def converter(self):
        command = f"bcl2fastq -R {self.runfolder_dir} -o {self.output_dir} --sample-sheet {self.sample_sheet_path} -r {self.loading_thread} -w {self.writing_thread} -p {self.processing_thread} --fastq-compression-level {self.fastq_compression_level}"
        
        if self.input_dir:
            command += f" -i {self.input_dir}"
        
        if self.intensities_dir:
            command += f" --intensities_dir {self.intensities_dir}"

        if self.interop_dir:
            command += f" --interop-dir {self.interop_dir}"

        if self.stats_dir:
            command += f" --stats-dir {self.stats_dir}"

        if self.reports_dir:
            command += f" --reports-dir {self.reports_dir}"

        if self.no_lane_splitting:
            command += " --no-lane-splitting"
        
        if self.barcode_mismatches:
            command += f" --barcode-mismatches {self.barcode_mismatches}"
        
        if self.use_bases_mask:
            command += f" --use-bases-mask {self.use_bases_mask}"

        if self.min_log_level:
            command += f" -l {self.min_log_level}"
        
        if self.tiles:
            command += f" --tiles {self.tiles}"

        if self.minimum_trimmed_read_length: 
            command += f" --minimum-trimmed-read-length {self.minimum_trimmed_read_length}"
        
        if self.mask_short_adapter_reads:
            command += f" --mask-short-adapter-reads {self.mask_short_adapter_reads}"

        if self.adapter_stringency:
            command += f" --adapter-stringency {self.adapter_stringency}" 
        
        if self.ignore_missing_bcls:
            command += f" --ignore-missing-bcls"
        
        if self.ignore_missing_filter:
            command += f" --ignore-missing-filter"
        
        if self.ignore_missing_positions:
            command += f" --ignore-missing-positions"
        
        if self.ignore_missing_controls:
            command += f" --ignore-missing-controls"
        
        if self.write_fastq_reverse_complement:
            command += f" --write-fastq-reverse-complement"

        if self.with_failed_reads:
            command += f" --with-failed-reads"

        if self.create_fastq_for_index_reads:
            command += f" --create-fastq-for-index-reads"

        if self.find_adapters_with_sliding_window:
            command += f" --find-adapters-with-sliding-window"

        if self.no_bgzf_compression:
            command += f" --no-bgzf-compression"

        os.system(command)

        return self.output_dir

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Converts BCLs to Fastq files in a specified directory")

    parser.add_argument("--runfolder_dir", required=True, help="Input directory containing BCL data")
    parser.add_argument("--input_dir", help = "Path to the input directory (=<runfolder-dir>/Data/Intensities/BaseCalls/)")
    parser.add_argument("--intensities_dir", help = "Path to the intensities directory (=<input-dir>/../) ")
    parser.add_argument("--interop_dir", help = "path to demultiplexing statistics directory (=<runfolder-dir>/InterOp/)")
    parser.add_argument("--stats_dir", help = "path to human-readable demultiplexing statistics directory (=<output-dir>/Stats/)")
    parser.add_argument("--reports_dir", help = "path to reporting directory (=<output-dir>/Reports/)")
    parser.add_argument("--output_dir", required=True, help="Output directory for Fastq files")
    parser.add_argument("--sample_sheet", required=True, help="Sample sheet path")
    parser.add_argument("--loading_thread", type=int, default=4, required=True, help="Number of loading threads")
    parser.add_argument("--processing_thread", type=int, required=True, help="Number of processing threads")
    parser.add_argument("--writing_thread", type=int, default=4, required=True, help="Number of writing threads")
    parser.add_argument("--fastq_compression_level", type=int, default=4, required=True, help="Fastq compression level")
    parser.add_argument("--barcode_mismatches", type=int, default=1, help="Number of allowed barcode mismatches per index")
    parser.add_argument("--use_bases_mask", help="Bases mask for sequencing")
    parser.add_argument("--no_lane_splitting", action="store_true", help="Disable lane splitting")
    parser.add_argument("--min_log_level", type=str, default = "INFO", help="Default = INFO, Recognized values: NONE, FATAL, ERROR, WARNING, INFO, DEBUG, TRACE")
    parser.add_argument("--tiles", type=str, help="comma-separated list of regular expressions to select only a subset of the tiles available in the flow-cell. Multiple entries allowed, each applies to the corresponding base-calls.")
    parser.add_argument("--minimum_trimmed_read_length", type=int, default=35, help=" minimum read length after adapter trimming")
    parser.add_argument("--mask_short_adapter_reads", type=int, default=22, help="smallest number of remaining bases (after masking bases below the minimum trimmed read length) below which whole read is masked")
    parser.add_argument("--adapter_stringency", type=float, default=0.9, help="adapter stringency")
    parser.add_argument("--ignore_missing_bcls", action="store_true", help="assume 'N'/'#' for missing calls")
    parser.add_argument("--ignore_missing_filter", action="store_true", help="assume 'true' for missing filters")
    parser.add_argument("--ignore_missing_positions", action="store_true", help=" assume [0,i] for missing positions, where i is incremented starting from 0")
    parser.add_argument("--ignore_missing_controls", action="store_true", help="(deprecated) assume 0 for missing controls")
    parser.add_argument("--write_fastq_reverse_complement", action="store_true", help="generate FASTQs containing reverse complements of actual data")
    parser.add_argument("--with_failed_reads", action="store_true", help="include non-PF clusters")
    parser.add_argument("--create_fastq_for_index_reads", action="store_true", help="create FASTQ files also for index reads")
    parser.add_argument("--find_adapters_with_sliding_window", action="store_true", help="find adapters with simple sliding window algorithm")
    parser.add_argument("--no_bgzf_compression", action="store_true", help="turn off BGZF compression for FASTQ files")


    args = parser.parse_args()

    bclconverter = BCLConverter(**vars(args))
    bclconverter.converter()


