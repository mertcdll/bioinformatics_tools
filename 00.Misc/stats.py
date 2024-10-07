import json
import os
import argparse
import pandas as pd
from tabulate import tabulate
from reportlab.lib.pagesizes import letter
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Spacer, paragraph,PageTemplate,Frame


class Stats_to_pdf:
    def __init__(self, json_file, output):
        self.json_file = json_file
        self.output = output
        self.data = None
        self.conversion_results = []
        self.sample_aggregated_metrics = {}
        self.table_style = [
            ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 8.5),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
        ]
        self.load_data()
        self.aggregate_sample_metrics()
        self.generate_reports()
        self.final()

    def load_data(self):
        with open(self.json_file) as f:
            self.data = json.load(f)
            self.conversion_results = self.data['ConversionResults']

    def create_table(self, headers, table_data):
        table = Table([headers] + table_data)
        table.setStyle(TableStyle(self.table_style))
        return table

    def aggregate_sample_metrics(self):
        for lane in self.conversion_results:
            for sample in lane['DemuxResults']:
                sample_id = sample["SampleId"]
                if sample_id not in self.sample_aggregated_metrics:
                    self.sample_aggregated_metrics[sample_id] = {
                        "SampleId": sample_id,
                        "SampleName": sample["SampleName"],
                        "IndexMetrics": [],
                        "NumberReads": 0,
                        "Yield": 0,
                        "ReadMetrics": []
                    }


                self.sample_aggregated_metrics[sample_id]["IndexMetrics"].extend(sample["IndexMetrics"])
                
                self.sample_aggregated_metrics[sample_id]["NumberReads"] += sample["NumberReads"]
                self.sample_aggregated_metrics[sample_id]["Yield"] += sample["Yield"]

                for read_metric in sample["ReadMetrics"]:
                    read_number = read_metric["ReadNumber"]
                    if len(self.sample_aggregated_metrics[sample_id]["ReadMetrics"]) < read_number:
                        self.sample_aggregated_metrics[sample_id]["ReadMetrics"].append({
                            "ReadNumber": read_number,
                            "Yield": 0,
                            "YieldQ30": 0,
                            "QualityScoreSum": 0,
                            "TrimmedBases": 0
                        })

                    self.sample_aggregated_metrics[sample_id]["ReadMetrics"][read_number - 1]["Yield"] += read_metric["Yield"]
                    self.sample_aggregated_metrics[sample_id]["ReadMetrics"][read_number - 1]["YieldQ30"] += read_metric["YieldQ30"]
                    self.sample_aggregated_metrics[sample_id]["ReadMetrics"][read_number - 1]["QualityScoreSum"] += read_metric["QualityScoreSum"]
                    self.sample_aggregated_metrics[sample_id]["ReadMetrics"][read_number - 1]["TrimmedBases"] += read_metric["TrimmedBases"]

    def generate_reports(self):
        for sample_id, sample_metrics in self.sample_aggregated_metrics.items():

            elements = []

            table_data0 = [[sample_metrics["SampleId"], sample_metrics["SampleName"],
                            sample_metrics["NumberReads"], sample_metrics["Yield"]]]
            headers0 = ["Sample id", "Sample name", "Number of Reads", "Yield"]
            t0 = self.create_table(headers0, table_data0)
            elements.extend([t0, Spacer(0, 20)])

            table_data = []
            headers = ["Index Sequence", "Mismatch Counts_0", "Mismatch Counts_1"]
            for metrics in sample_metrics["IndexMetrics"]:
                index_sequence = metrics["IndexSequence"]
                mismatch_counts = metrics["MismatchCounts"]
                mismatch_counts_0 = mismatch_counts["0"]
                mismatch_counts_1 = mismatch_counts["1"]
                table_data.append([index_sequence, mismatch_counts_0, mismatch_counts_1])

            t = self.create_table(headers, table_data)
            elements.extend([t, Spacer(0, 20)])

            table_data1 = []
            headers1 = ["Read Number", "Yield", "Yield Q30", "Quality Score Sum", "Trimmed Bases"]
            for metric in sample_metrics["ReadMetrics"]:
                read_number = metric["ReadNumber"]
                yield_value = metric["Yield"]
                yield_q30 = metric["YieldQ30"]
                quality_score_sum = metric["QualityScoreSum"]
                trimmed_bases = metric["TrimmedBases"]
                table_data1.append([read_number, yield_value, yield_q30, quality_score_sum, trimmed_bases])

            t1 = self.create_table(headers1, table_data1)
            elements.extend([t1, Spacer(0, 20)])


                                
            pdf_filename = os.path.join(self.output, f"{sample_id}.pdf")
            doc = SimpleDocTemplate(pdf_filename, pagesize=letter)
            doc.build(elements)

    
    def final(self):
        csv_data = []
        cnv_data_each_read = []
        table_data = []
        for sample in self.data['ConversionResults']:
            for elm in sample["DemuxResults"]:
                sample_id = elm["SampleId"]
                number_reads = elm.get("NumberReads", "")
                yield_mb = elm.get("Yield", "")
                read_metrics = elm.get("ReadMetrics", [])
                csv_data.append([sample_id, number_reads, yield_mb])
                for i, read in enumerate(read_metrics, start=1):
                    read_number = i
                    read_yield = read.get("Yield", "")
                    yield_q30 = read.get("YieldQ30", "")
                    quality_score_sum = read.get("QualityScoreSum", "")
                    trimmed_bases = read.get("TrimmedBases", "")
                    cnv_data_each_read.append([sample_id, read_number, read_yield, yield_q30, quality_score_sum, trimmed_bases])
                    table_data.append([sample_id, number_reads, yield_mb, read_number, read_yield, yield_q30, quality_score_sum, trimmed_bases])
        
        df = pd.DataFrame(csv_data)
        df.columns = ["sample_id", "Total_reads_number", "Total_yield_base"]
        sum_df = df.groupby("sample_id").agg({"Total_reads_number":"sum", "Total_yield_base": "sum"})
        sum_df.reset_index(inplace=True)
        df1 = pd.DataFrame(cnv_data_each_read)
        df1.columns = ["sample_id", "read_number", "yield_base", "yield_base_q30", "quality_score_sum", "trimmed_bases"]
        sum_df1 = df1.groupby(["sample_id", "read_number"]).agg({"yield_base": "sum", "yield_base_q30": "sum", "quality_score_sum": "sum", "trimmed_bases": "sum"})
        sum_df1.reset_index(inplace=True)

        with pd.ExcelWriter(os.path.join(self.output, f"QC_metrics.xlsx")) as writer:
                sum_df.to_excel(writer, sheet_name='Total reads and bases', index=False, header=True)
                sum_df1.to_excel(writer, sheet_name='Read details', index=False, header=True)

        table_style = TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0),8.5),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
        ])

        pdf = SimpleDocTemplate(os.path.join(self.output, f"QC_metrics.pdf"), pagesize=letter)

        headers = ["sample_id", "Total_reads_number", "Total_yield_base", "read_number", "yield_base", "yield_base_q30", "quality_score_sum", "trimmed_bases"]
        num_columns = len(table_data[0])
        num_rows = len(table_data)

        table = Table([headers] + table_data, colWidths=[75]*num_columns, rowHeights=30, repeatRows=1)
        table.setStyle(table_style)

        elements = []
        elements.append(table)
        pdf.build(elements)

stats_file_path = "/home/genwork2/Mert/statsupdate/Stats_nolanesplitting.json"
output_path = "/home/genwork2/Mert/statsupdate/statsout"
instance = Stats_to_pdf(stats_file_path, output_path)
