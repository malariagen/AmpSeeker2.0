import json, glob, os, re
import plotly.graph_objects as go
import plotly.io as pio

os.makedirs('results/plots', exist_ok=True)
for nb_path in glob.glob('*.ipynb'):
    nb = json.load(open(nb_path))
    fig_num = 0
    nb_name = nb_path.replace('.ipynb', '')
    for cell in nb.get('cells', []):
        for out in cell.get('outputs', []):
            data = out.get('data', {})
            fig = None

            # Approach 1: native Plotly JSON output
            if 'application/vnd.plotly.v1+json' in data:
                try:
                    fig = go.Figure(data['application/vnd.plotly.v1+json'])
                except Exception as e:
                    print(f'Warning (JSON): {e}')

            # Approach 2: plotly.offline HTML div output
            elif 'text/html' in data:
                html = ''.join(data['text/html']) if isinstance(data['text/html'], list) else data['text/html']
                match = re.search(r'var figure = ({.*?"data":.*?})\s*;?\s*</script>', html, re.DOTALL)
                if not match:
                    match = re.search(r'Plotly\.newPlot\([^,]+,\s*({.*?})\s*\)', html, re.DOTALL)
                if match:
                    try:
                        fig = go.Figure(json.loads(match.group(1)))
                    except Exception as e:
                        print(f'Warning (HTML): {e}')

            if fig is not None:
                fig_num += 1
                try:
                    fig.write_image(f'results/plots/{nb_name}_figure_{fig_num}.pdf', format='pdf')
                    fig.write_image(f'results/plots/{nb_name}_figure_{fig_num}.png', format='png', scale=2)
                    print(f'Saved: {nb_name}_figure_{fig_num}.pdf/.png')
                except Exception as e:
                    print(f'Warning (save): {e}')

print(f'PDF/PNG save complete')
