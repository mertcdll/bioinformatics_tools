import pdfkit

class HTMLtoPDFConverter:
    def __init__(self, html_file_path, pdf_file_path, css_file_path=None):
        self.html_file_path = html_file_path
        self.pdf_file_path = pdf_file_path
        self.css_file_path = css_file_path

    def convert_to_pdf(self):
        try:
            pdfkit.from_file(self.html_file_path, self.pdf_file_path, css=self.css_file_path)
            print("Conversion successful!")
        except Exception as e:
            print(f"Error during conversion: {e}")




if __name__ == "__main__":
    html_to_pdf_converter = HTMLtoPDFConverter(
        "/home/dell/Documents/IlyomeReport/report.html",
        "/home/dell/Documents/IlyomeReport/out2.pdf",
        css="/home/dell/Documents/IlyomeReport/style.A4.css"
    )

    html_to_pdf_converter.convert_to_pdf()
