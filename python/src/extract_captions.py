import re

with open('pmc_page.html', encoding='utf-8') as f:
    text = f.read()

for f_id in ['F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7']:
    pattern = rf'id="{f_id}".*?<figcaption>(.*?)</figcaption>'
    m = re.search(pattern, text, re.DOTALL)
    if m:
        clean = re.sub('<.*?>', '', m.group(1))
        # clean extra whitespace
        clean = re.sub(r'\s+', ' ', clean)
        print(f"=== FIGURE {f_id} ===")
        print(clean.strip())
        print()

