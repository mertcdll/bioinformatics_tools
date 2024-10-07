import numpy as np 
import pandas as pd
from io import StringIO
from reportlab.lib.pagesizes import letter, inch
from reportlab.platypus import SimpleDocTemplate, Paragraph, Table, TableStyle
from reportlab.platypus import Image as PLImage
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY
from plotnine import *
import argparse
import os
import io
import concurrent.futures as cf

class QCReportGenerator:
    
    def __init__(self, input_path, output_path, max_workers):
        self.input_path = input_path
        self.output_path = output_path
        self.max_workers = max_workers

    def fastqc_sections(self, fastqc_file):
        
        with open(fastqc_file, 'r') as file:
            fastqc_metrics = file.read()

        fqc_sections = fastqc_metrics.strip().split('>>END_MODULE')

        sec_dataframes_dict = {}

        for section in fqc_sections[:-1]:
            lines = section.split('\n')
            sec_names = lines[1].strip('>>').split('\t')[0]

            if sec_names == 'Sequence Duplication Levels':
                sec_headers = lines[3].strip('#').split('\t')
                sec_rows = '\n'.join(lines[4:])
            else:
                sec_headers = lines[2].strip('#').split('\t')
                sec_rows = '\n'.join(lines[3:])

            df = pd.read_csv(StringIO(sec_rows), delimiter='\t', names=sec_headers)
            sec_dataframes_dict[sec_names] = df
        

        return sec_dataframes_dict
    
    def generate_pdf_report(self, input_file, output_file):

        sec_dataframes_dict = self.fastqc_sections(input_file)

        doc = SimpleDocTemplate(output_file, pagesize=letter)

        story = []

        styles = getSampleStyleSheet()
        styles.add(ParagraphStyle(name='justify', parent=styles['Normal'], alignment=TA_JUSTIFY))

        story.append(Paragraph("Genoks QC Report", styles['Title']))
        story.append(Paragraph("<br/><br/>", styles['Normal']))

        
        def add_section(title, content):
            story.append(Paragraph(title, styles['Heading1']))
            story.append(Paragraph(content, styles['justify']))
            story.append(Paragraph("<br/><br/>", styles['Normal']))

        # Basic Statistics Section
        bs = sec_dataframes_dict.get('Basic Statistics')
        bs_desc = '''The Basic Statistics module generates simple composition statistics for the analyzed file.
        Filename: The original name of the file being analyzed.
        File type: Indicates whether the file contains actual base calls or colorspace data that has been converted to base calls.
        Encoding: Specifies the ASCII encoding of quality values found in the file.
        Total Sequences: The total number of sequences processed, reported as both an actual count and an estimated count (currently the same). Future versions may allow analysis of a subset of sequences to estimate the total count, but this feature is currently disabled due to uneven distribution of problematic sequences.
        Filtered Sequences: In Casava mode, sequences flagged for filtering are removed from all analyses. The number of filtered sequences is reported separately, and the total sequences count above excludes these filtered sequences for the rest of the analysis.
        Sequence Length: Provides the length of the shortest and longest sequence in the dataset. If all sequences have the same length, only one value is reported.
        %GC: Represents the overall percentage of GC content across all bases in all sequences.'''
        
        if bs is not None and not bs.empty:

            add_section("Basic Statistics", bs_desc)
            table_data = [list(bs.columns)] + bs.values.tolist()
            table = Table(table_data)
            table.setStyle(TableStyle([('BACKGROUND', (0, 0), (-1, 0), colors.grey),
                                    ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                                    ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                                    ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                                    ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
                                    ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
                                    ('GRID', (0, 0), (-1, -1), 1, colors.black)]))
            story.append(table)
        
        else:
            add_section("Basic Statistics", bs_desc)
            story.append(Paragraph("There is no basic statistics table for this sequencing run.", styles['Normal']))


        story.append(Paragraph("<br/><br/>", styles['Normal']))

        # Per Base Sequence Quality Section
        pbsq = sec_dataframes_dict.get('Per base sequence quality')

        pbsq_desc = '''This graph provides a representation of quality values observed for all bases at each position within the FastQ file.

    The key features of the plot are as follows:

        Central Black Line (Median): The black line represents the median value of the quality scores at each position. The median is a measure of the central tendency of the data.

        Yellow Box (Inter-quartile Range): The yellow box represents the inter-quartile range, which spans from the 25th percentile (lower quartile) to the 75th percentile (upper quartile) of the quality scores. It gives a measure of the spread or variability of the data.

        Upper and Lower Whiskers: The whiskers extend from the box and represent the range of the data beyond the inter-quartile range. The upper whisker goes up to the 90th percentile, and the lower whisker goes down to the 10th percentile.

        Blue Line with Red Dots (Mean Quality): The blue line with red dots represents the mean (average) quality score at each position. The mean provides an overall measure of the average quality.

    The y-axis of the plot represents the quality scores, where higher scores indicate better base calls. The background of the graph is divided into color-coded regions.

        Green: Very good quality calls.
        Yellow: Calls of reasonable quality.
        Red: Calls of poor quality.'''

        if pbsq is not None and not pbsq.empty:
            
            pbsq['Baseint'] = pbsq["Base"].apply(lambda x: sum(map(int, x.split("-"))) / 2 if '-' in x else int(x))
            pbsq['Base'] = pd.Categorical(pbsq['Base'], categories=sorted(pbsq['Base'].unique()))
            pbsq['Baseint'] = pbsq['Baseint'].astype(int)
            pbsq['Basecat'] = pd.Categorical(pbsq['Baseint'], categories=sorted(pbsq['Baseint'].unique()))
            
            pbsq_chart = (
                ggplot(pbsq, aes(x='Basecat',ymin='10th Percentile', lower='Lower Quartile', middle='Median', upper='Upper Quartile', ymax='90th Percentile')) +
                geom_boxplot(stat= 'identity', width=0.5, color='black', fill='yellow', alpha=1) +
                geom_point(aes(x='Basecat', y='Mean'), color='red', size=0.5) +
                geom_line(aes(x='Basecat', y='Mean', group=1), color='blue') +
                geom_rect(aes(xmin=-float('inf'), xmax=float('inf'), ymin=float('-inf'), ymax=20), fill='red', alpha=0.02) +
                geom_rect(aes(xmin=-float('inf'), xmax=float('inf'), ymin=20, ymax=30), fill='yellow', alpha=0.02) +
                geom_rect(aes(xmin=-float('inf'), xmax=float('inf'), ymin=30, ymax=40), fill='green', alpha=0.02) +
                scale_y_continuous(breaks=list(range(0, 41, 2)), expand=(0.01, 0)) +
                theme_minimal() +
                ylim(0, 40) +
                labs(x='Position in read (bp)', y='', title='Quality Scores Across All Bases Sanger/Illumina 1.9 Encoding') +
                theme(plot_title=element_text(size=12, hjust=0.5, face='bold'), axis_text_x=element_text(size=8, rotation=90))
            )
            

            add_section("Per Base Sequence Quality", pbsq_desc)
            buffer=io.BytesIO()
            pbsq_chart.save(buffer, format='png')
            buffer.seek(0)
            story.append(PLImage(buffer, 5.5 * inch, 3 * inch))
        
        else:
            add_section("Per Base Sequence Quality", pbsq_desc)
            story.append(Paragraph("There is no per base sequence quality graph for this sequencing run.", styles['Normal']))
        
        story.append(Paragraph("<br/><br/>", styles['Normal']))

        # Per Tile Sequence Quality Section
        
        ptsq = sec_dataframes_dict.get('Per tile sequence quality')
        ptsq_desc = '''This graph will be visible in analysis results when an Illumina library that preserves its original sequence identifiers, which contain information about the flowcell tile from which each read originated is utilized. The graph allows us to examine the quality scores for each tile across all bases, helping identify any potential loss in quality associated with specific parts of the flowcell.

    The plot illustrates the variation from the average quality for each tile, with colors ranging from cold to hot. Cold colors indicate positions where the quality equals or surpasses the average for that base in the run, while hotter colors suggest that a tile exhibited lower quality compared to other tiles for that specific base. An ideal plot would appear predominantly blue.'''
        
        if ptsq is not None and not ptsq.empty:
        
            ptsq['Baseint'] = ptsq["Base"].apply(lambda x: sum(map(int, x.split("-"))) / 2 if '-' in x else int(x))
            ptsq['Base'] = pd.Categorical(ptsq['Base'], ordered=True, categories=ptsq['Base'].unique())
            ptsq['Baseint'] = ptsq['Baseint'].astype(int)
            ptsq['Tilecat'] = pd.Categorical(ptsq['Tile'], categories=sorted(ptsq['Tile'].unique()))
            
            ptsq_chart = (
                ggplot(ptsq, aes(x='Base', y='Tilecat', fill='Mean'))
                + geom_tile()
                + scale_fill_gradient(low='red', high='blue')
                + scale_y_discrete(expand=(0.01, 0), breaks=lambda x: x[::50]) 
                + labs(x='Position in read (bp)', y='', title='Quality per tile') 
                + theme(plot_title=element_text(size=12, hjust=0.5, face='bold'), axis_text_x=element_text(size=8, rotation=90), panel_background=element_rect(fill='white'))
            )
            
            add_section("Per Tile Sequence Quality", ptsq_desc)
            buffer=io.BytesIO()
            ptsq_chart.save(buffer, format='png')
            buffer.seek(0)
            story.append(PLImage(buffer, 5.5 * inch, 3 * inch))
        
        else:
            add_section("Per Tile Sequence Quality", ptsq_desc)
            story.append(Paragraph("There is no per tile sequence quality graph for this sequencing run.", styles['Normal']))


        story.append(Paragraph("<br/><br/>", styles['Normal']))

        # Per Sequence Quality Scores Section
        psqs = sec_dataframes_dict.get('Per sequence quality scores')
        psqs_desc = '''The per-sequence quality score report enables us to identify whether a specific group of sequences exhibits consistently low quality values. Such occurrences are typically observed when certain sequences are poorly imaged (e.g., located on the edge of the field of view). However, these sequences should only make up a small fraction of the total dataset.'''

        if psqs is not None and not psqs.empty:

            psqs_chart = (
                ggplot(psqs, aes(x='Quality', y='Count'))
                + geom_line(color='red')
                + scale_y_continuous(breaks=list(np.linspace(min(psqs['Count']), max(psqs['Count']), 20)), expand=(0.01, 0)) 
                + scale_x_continuous(breaks=list(np.arange(min(psqs['Quality']), max(psqs['Quality']), 2)))
                + theme_minimal()
                + geom_text(x=max(psqs['Quality'])-6, y= 0.9*max(psqs['Count']), label="Average Quality per Read", color="red", size=7, va="bottom")   
                + labs(x='Mean Sequence Quality (Phred Score)', y='', title='Quality score distribution over all sequences')
                + theme(plot_title=element_text(size=12, hjust=0.5, face='bold'), axis_text_x=element_text(size=8))
            )
            add_section("Per Sequence Quality Scores", psqs_desc)
            buffer=io.BytesIO()
            psqs_chart.save(buffer, format='png')
            buffer.seek(0)
            story.append(PLImage(buffer, 5.5 * inch, 3 * inch))

        else:
            add_section("Per Sequence Quality Scores", psqs_desc)
            story.append(Paragraph("There is no per sequence quality scores graph for this sequencing run.", styles['Normal']))
        
        story.append(Paragraph("<br/><br/>", styles['Normal']))

        # Per Base Sequence Content Section
        pbsc = sec_dataframes_dict.get('Per base sequence content')
        
        pbsc_desc = '''The "Per Base Sequence Content" graph illustrates the distribution of the four normal DNA bases at each position in a file. In an unbiased library, one would anticipate minimal variation among the different bases throughout the sequence run, resulting in parallel lines on the plot. The proportion of each base should generally align with the overall abundance of these bases in your genome, and they should not exhibit significant imbalances relative to each other.'''

        if pbsc is not None and not pbsc.empty:

            pbsc_melted = pd.melt(pbsc, id_vars='Base', var_name='Nucleotide', value_name='Value')
            pbsc_melted['Baseint'] = pbsc_melted["Base"].apply(lambda x: sum(map(int, x.split("-"))) / 2 if '-' in x else int(x))
            pbsc_melted['Baseint'] = pbsc_melted['Baseint'].astype(int)
            pbsc_chart = (
                ggplot(pbsc_melted, aes(x='Baseint', y='Value', color='Nucleotide'))
                + geom_line()
                + scale_x_continuous(name='Position in read(bp)', breaks=np.arange(min(pbsc_melted['Baseint']), max(pbsc_melted['Baseint']) + 5, 5))
                + scale_y_continuous(name=None, breaks=np.arange(0, 101, 10), limits=[0, 100])
                + labs(title='Sequence Content Across All Bases(Percentage)', color='Nucleotide')
                + theme_minimal()
                + theme(plot_title=element_text(size=12, hjust=0.5, face='bold'), axis_text_x=element_text(size=8, rotation=90), axis_title_y=element_blank())
            )
            add_section("Per Base Sequence Content", pbsc_desc)
            buffer=io.BytesIO()
            pbsc_chart.save(buffer, format='png')
            buffer.seek(0)
            story.append(PLImage(buffer, 5.5 * inch, 3 * inch))
        else:
            add_section("Per Base Sequence Content", pbsc_desc)
            story.append(Paragraph("There is no per base sequence content graph for this sequencing run.", styles['Normal']))

        story.append(Paragraph("<br/><br/>", styles['Normal']))

        # Per Sequence GC Content Section 
        psgcc = sec_dataframes_dict.get('Per sequence GC content')
        psgcc_desc = '''The Per Sequence GC Content plot provides us with insights into the distribution of GC content across all sequences within your run. For most applications, this distribution should be relatively normal, with the average GC content aligning with the desired value for the library. For reduced complexity libraries, it is important to avoid biases caused by having either too high or too low of a GC content. Different library preparation methods can result in different average GC contents, making this plot a useful resource for verifying that the GC content distribution aligns with the desired values.'''

        if psgcc is not None and not psgcc.empty:

            psgcc_chart = (
                ggplot(psgcc, aes(x='GC Content', y='Count'))
                + geom_line(color='red')
                + labs(title='GC distribution over all sequences', x='Mean GC content (%)')
                + scale_x_continuous(name='Mean GC Content (%)', breaks=np.arange(min(psgcc['GC Content']), max(psgcc['GC Content']) + 3, 3))
                + scale_y_continuous(name=None, breaks=list(np.linspace(min(psgcc['Count']), max(psgcc['Count']), 15)), expand=(0.01, 0))
                + theme_minimal()
                + theme(plot_title=element_text(size=12, hjust=0.5, face='bold'), axis_text_x=element_text(size=8, rotation=90), axis_title_y=element_blank())
            )   
            add_section("Per Sequence GC Content", psgcc_desc)
            buffer=io.BytesIO()
            psgcc_chart.save(buffer, format='png')
            buffer.seek(0)
            story.append(PLImage(buffer, 5.5 * inch, 3 * inch))

        else:
            add_section("Per Sequence GC Content", psgcc_desc)
            story.append(Paragraph("There is no per sequence gc content graph for this sequencing run.", styles['Normal']))


        story.append(Paragraph("<br/><br/>", styles['Normal']))

        # Per Base N Content Section
        pbnc = sec_dataframes_dict.get('Per base N content')
        pbnc_desc = '''The Per Base N Content plot illustrates the distribution of Ns (undetermined nucleotides) at each position in the file. Although these Ns will typically occur at random within the file, specific positions with a higher N content may indicate regions of low complexity, which could pose issues for some applications (e.g., when mapping reads to a reference genome). As such, the presence of consistently elevated Ns can be useful for troubleshooting problematic sequencing libraries.'''

        if pbnc is not None and not pbnc.empty:

            pbnc['Baseint'] = pbnc["Base"].apply(lambda x: sum(map(int, x.split("-"))) / 2 if '-' in x else int(x))

            pbnc_chart = (
                ggplot(pbnc, aes(x='Baseint', y='N-Count'))
                + geom_line(color='red')
                + labs(title='N-Count Across All Bases', x='Position in read', y='')
                + scale_x_continuous(breaks=np.arange(0, 153, 3))
                + scale_y_continuous(breaks=np.arange(0,101,10), limits = [0,100])
                + theme_minimal()
                + theme(plot_title=element_text(size=12, hjust=0.5, face='bold'), axis_text_x=element_text(size=8, rotation=90))
            )   
            add_section("Per Base N Content", pbnc_desc)
            buffer=io.BytesIO()
            pbnc_chart.save(buffer, format='png')
            buffer.seek(0)
            story.append(PLImage(buffer, 5.5 * inch, 3 * inch))

        else:
            add_section("Per Base N Content", pbnc_desc)
            story.append(Paragraph("There is no per base N content graph for this sequencing run.", styles['Normal']))


        story.append(Paragraph("<br/><br/>", styles['Normal']))

        # Sequence Length Distribution Section 
        sld = sec_dataframes_dict.get('Sequence Length Distribution')
        sld_desc = '''The Sequence Length Distribution graph illustrates the number of sequences with specific read lengths. Ideally, one would anticipate a relatively tight distribution centered around a specific length that aligns with library preparation and sequencing objectives. Deviations from this expected length could indicate issues with sample quality or with the library preparation itself.'''


        if sld is not None and not sld.empty:

            if '-' in str(sld["Length"].iloc[0]):
                sld['Lengthint'] = sld["Length"].apply(lambda x: sum(map(int, x.split("-"))) / 2 if '-' in x else int(x))
            else:
                sld['Lengthint'] = sld["Length"].astype(int)

            if '-' in str(sld["Length"].iloc[0]):

                sld_chart = (
                    ggplot(sld, aes(x='Lengthint', y='Count'))
                    + geom_line(color='red')
                    + scale_y_continuous(breaks=list(np.linspace(min(sld['Count']), max(sld['Count']), 20)), expand=(0.01, 0)) 
                    + scale_x_continuous(breaks=list(np.arange(min(sld['Lengthint']), max(sld['Lengthint']),5)))
                    + theme_minimal()
                    + geom_text(x=max(sld['Lengthint'])-20, y= 0.9*max(sld['Count']), label="Sequence Length", color="red", size=7, va="bottom")   
                    + labs(x='Sequence Length (bp)', y= '', title = 'Distribution of Sequence Lengths over All Sequences')
                    + theme(plot_title=element_text(size=12, hjust=0.5, face='bold'), axis_text_x=element_text(size=8, rotation=90))
                )
                
                

            else:
                
                sld_chart = (
                    ggplot(sld, aes(x='Lengthint', y='Count'))
                    + geom_point(color='red')
                    + scale_y_continuous(breaks=list(np.linspace(0, int(sld['Count'][0]), 20).astype(int)), expand=(0.01, 0)) 
                    + scale_x_continuous(breaks=list(np.arange(0, sld['Lengthint'][0] + 1, 5)))
                    + theme_minimal()
                    + geom_text(x=sld['Lengthint'][0], y=0.9 * int(sld['Count'][0]), label="Sequence Length", color="red", size=7, va="bottom")   
                    + labs(x='Sequence Length (bp)', y='', title='Distribution of Sequence Lengths over All Sequences')
                    + theme(plot_title=element_text(size=12, hjust=0.5, face='bold'), axis_text_x=element_text(size=8, rotation=90))
                )

            add_section("Sequence Length Distribution", sld_desc)
            buffer=io.BytesIO()
            sld_chart.save(buffer, format='png')
            buffer.seek(0)
            story.append(PLImage(buffer, 5.5 * inch, 3 * inch))
        
        else:
            add_section("Sequence Length Distribution", sld_desc)
            story.append(Paragraph("There is no sequence length distribution graph for this sequencing run.", styles['Normal']))
           

        story.append(Paragraph("<br/><br/>", styles['Normal']))


        # Adapter Content Section
        ac = sec_dataframes_dict.get('Adapter Content')
        ac_desc = '''The "Adapter Content" plot displays the presence and abundance of adapter sequences at each position in the read. Adapters are often used during library preparation to facilitate sequencing of the target regions. However, it is essential to ensure that adapter sequences are effectively removed or trimmed before downstream analyses, as they may interfere with alignment and lead to misinterpretation of the data.'''


        if ac is not None and not ac.empty:
            
            ac_melted = pd.melt(ac, id_vars='Position', var_name='Adapter Type', value_name='Value')

            ac_melted['Baseint'] = ac_melted["Position"].apply(lambda x: sum(map(int, x.split("-"))) / 2 if '-' in x else int(x))

            ac_chart = (
                ggplot(ac_melted, aes(x='Baseint', y='Value', color='Adapter Type'))
                + geom_line()
                + scale_x_continuous(name='Position in read(bp)', breaks=np.arange(min(ac_melted['Baseint']), max(ac_melted['Baseint']) + 5, 5))
                + scale_y_continuous(name=None, breaks = np.arange(0,101,10), limits=[0, 100])
                + labs(title='Adapter Content', color='Adapter Type')
                + theme_minimal()
                + theme(plot_title=element_text(size=12, hjust=0.5, face='bold'), 
                    axis_text_x=element_text(size=8, rotation=90), 
                    axis_title_y=element_blank(), 
                    legend_position='top',
                    legend_text=element_text(size=4.4),
                    legend_title=element_text(size=7))
        
            )   


            add_section("Adapter Content", ac_desc)
            buffer=io.BytesIO()
            ac_chart.save(buffer, format='png')
            buffer.seek(0)
            story.append(PLImage(buffer, 5.5 * inch, 3 * inch))

        else: 
            add_section("Adapter Content", ac_desc)
            story.append(Paragraph("There is no adapter content graph for this sequencing run.", styles['Normal']))
           
        story.append(Paragraph("<br/><br/>", styles['Normal']))

        # Overrepresented Sequences Section
        orseqs = sec_dataframes_dict.get('Overrepresented sequences')
        
        orseqs_desc = '''A typical high-throughput library comprises a diverse set of sequences, with no individual sequence representing an insignificant fraction of the whole. However, discovering a single sequence that is excessively overrepresented can indicate its high biological significance, contamination in the library, or lower diversity than expected.
            This module provides a list of sequences that make up more than 0.1% of the total sequences. To conserve memory, only sequences within the first 100,000 entries are tracked until the end of the file. Consequently, an overrepresented sequence that appears later in the file may be missed.
            For each overrepresented sequence, the program searches for matches in a database of common contaminants and reports the best hit, requiring hits to be at least 20 base pairs long with no more than 1 mismatch. Finding a hit does not conclusively identify the source of contamination but may provide valuable clues. Additionally, similar adapter sequences might yield reported hits, which may not be technically accurate but have closely related sequences.
            Due to the requirement for exact sequence matches over the entire length of a sequence, reads longer than 75 base pairs are truncated to 50 base pairs for analysis. Longer reads are more susceptible to sequencing errors, artificially inflating observed diversity and potentially underrepresenting highly duplicated sequences.'''
        
        if orseqs is not None and not orseqs.empty:
            
            add_section("Overrepresented Sequences", orseqs_desc)
            orseqs_table_data = [list(orseqs.columns)] + orseqs.values.tolist()
            orseqs_table = Table(orseqs_table_data)
            orseqs_table.setStyle(TableStyle([('BACKGROUND', (0, 0), (-1, 0), colors.grey),
                                            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                                            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                                            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                                            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
                                            ('LEFTPADDING', (0, 0), (0, -1), 1),
                                            ('RIGHTPADDING', (0, 0), (0, -1), 1),
                                            ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
                                            ('GRID', (0, 0), (-1, -1), 1, colors.black),
                                            ('FONTSIZE',(0, 0), (-1, -1), 6),
                                            ('LEFTMARGIN', (0, 0), (0, -1), 0),
                                            ('RIGHTMARGIN', (0, 0), (0, -1), 0), 
                                            ('FONTSIZE', (0, 0), (-1, 0), 10)
                                            ]))
            
            story.append(orseqs_table)
        else:
            add_section("Overrepresented Sequences", orseqs_desc)
            story.append(Paragraph("There are no overrepresented sequences for this sequencing run.", styles['Normal']))
    
        story.append(Paragraph("<br/><br/>", styles['Normal']))
        doc.build(story)

    def generate_reports_in_directory(self):
        with cf.ProcessPoolExecutor(max_workers=self.max_workers) as executor:
            workers = []

            for root, dirs, _ in os.walk(self.input_path):
                if root == self.input_path:
                    for folder in dirs:
                        if folder.endswith('_fastqc'):
                            fastqc_folder_path = os.path.join(root, folder)
                            input_file = os.path.join(fastqc_folder_path, 'fastqc_data.txt')
                        
                            report_prefix = folder.replace('_fastqc', '')
                        
                            output_file = os.path.join(self.output_path, f"{report_prefix}_qcreport.pdf")
                            workers.append(executor.submit(self.generate_pdf_report, input_file, output_file))
            for ws in cf.as_completed(workers):
                ws.result()
                pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate PDF reports from FastQC data.')
    parser.add_argument('--input_path', type=str, help='Path to the directory containing FastQC data folders')
    parser.add_argument('--output_path', type=str, help='Path to the output directory to be set')
    parser.add_argument('--cores', type=int, default=1, help='Number of cores to use for parallel execution')
    args = parser.parse_args()

    input_directory_abs = os.path.abspath(args.input_path)
    output_directory_abs = os.path.abspath(args.output_path)

    report_generator = QCReportGenerator(args.input_path, args.output_path, args.cores)
    report_generator.generate_reports_in_directory()
    
