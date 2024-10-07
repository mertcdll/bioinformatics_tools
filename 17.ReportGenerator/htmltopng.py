import os
from html2image import Html2Image

class HTMLtoImageConverter:
    def __init__(self, directory_path):
        self.directory_path = directory_path

    def convert_to_image(self, html_file):
        html_path = os.path.join(self.directory_path, html_file)
        png_file = html_file.replace(".html", ".png")
        hti = Html2Image(custom_flags=['--virtual-time-budget=1000000000', '--hide-scrollbars'])
        hti.output_path = self.directory_path
        hti.size = (2000, 250)
        hti.screenshot(html_file=html_path, save_as=png_file)
        print(f"Conversion of {html_file} to {png_file} successful!")

if __name__ == "__main__":
    directory_path = "/home/dell/Documents/IlyomeReport/"
    
    html_files = [f for f in os.listdir(directory_path) if f.startswith("metric_") and f.endswith(".html")]

    html_to_image_converter = HTMLtoImageConverter(directory_path)

    for html_file in html_files:
        html_to_image_converter.convert_to_image(html_file)
