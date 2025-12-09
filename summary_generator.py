import os
from datetime import datetime

# Get current script directory
script_dir = os.path.dirname(os.path.abspath(__file__))

# Get list of PNG files
png_files = [f for f in os.listdir(script_dir) if f.lower().endswith('.png')]

# Create HTML content
html_lines = ['<html>', '<body>']
for img in png_files:
    html_lines.append(f'<div><img src="../{img}" style="max-width:100%; display:block; margin-bottom:20px;"></div>')
html_lines += ['</body>', '</html>']
html_content = '\n'.join(html_lines)

# Create output directory "2summary"
summary_dir = os.path.join(script_dir, 'summaries')
os.makedirs(summary_dir, exist_ok=True)

# Generate filename with timestamp
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
html_filename = f"{timestamp}_Summary.html"
html_path = os.path.join(summary_dir, html_filename)

# Write HTML file
with open(html_path, 'w', encoding='utf-8') as f:
    f.write(html_content)

print(f"HTML summary saved to: {html_path}")