from pathlib import Path

import bokeh.embed
import bokeh.resources
import yaml

from jinja2 import Template

from hits.visualize.interactive import parallel_coordinates

def make_modal(data):
    if isinstance(data, (str, Path)):
        data = yaml.safe_load(Path(data).read_text())

    modal_template = Template(Path('modal.jinja').read_text())

    modal = modal_template.render(data)

    return modal

def save_bokeh_html(layout, html_fn, description_data=None, modal_data=None):
    template_fn = 'bokeh_template.jinja'
    template = Path(template_fn).read_text()

    navbar_fn = 'navbar.html'
    template_variables = {
        'navbar': Path(navbar_fn).read_text(),
    }

    if description_data is not None:
        template_variables['description_data'] = description_data
    
    if modal_data is not None:
        template_variables['modal'] = make_modal(modal_data)
    
    html = bokeh.embed.file_html(layout, 
                                 bokeh.resources.CDN,
                                 template=template,
                                 template_variables=template_variables,
                                )

    Path(html_fn).write_text(html)
    
def save_PC_html(df, fn, width=1180, modal_data=None, description_data=None, **plot_kwargs):
    html_before = '''\
<html lang="en">

<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">

  <link rel="icon" href="images/icon.png">

  <title>Repair-seq data browser</title>

  <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-1BmE4kWBq78iYhFldvKuhfTAU6auU8tT94WrHftjDbrCEXSU1oBoqyl2QvZ6jIW3" crossorigin="anonymous">
  <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js" integrity="sha384-ka7Sk0Gln4gmtz2MlQnikT1wXgYsOg+OMhuP+IlRH9sENBO0LRn5q+8nbTov4+1p" crossorigin="anonymous"></script>
'''
    html_between = '''\
    <link href="styles/style.css" rel="stylesheet">
</head>
<body>
''' + Path('navbar.html').read_text()
    
    if modal_data is not None:
        help_button = '''\
<a style="margin-left: auto" href="#Modal" data-bs-toggle="modal" data-bs-target="#Modal">
<img id="help-button" src="images/question-square.svg" alt="Bootstrap" width="40" height="40">
</a>
'''
        #help_button_description = 'Click the question mark button for a guide to interacting with the plot.'
        help_button_description = ''
        modal = make_modal(modal_data)
    else:
        help_button = ''
        help_button_description = ''
        modal = ''
    
    if description_data is not None:
        html_between += f'''\
<div class="container" style="width: {width}px; margin-left: 50px; margin-top: 10px;">
<div style="display: flex">
<h2>
{description_data["title"]}
</h2>
{help_button}
</div>
{description_data["details"]}
{help_button_description}
</div>
'''
    html_after = f'''\
{modal}
</body>
</html>
'''
    
    extra_injections = {
        'html_before': html_before,
        'html_between': html_between,
        'html_after': html_after,
    }
    
    parallel_coordinates(df,
                         extra_injections=extra_injections,
                         save_as=fn,
                         table_width=width,
                         width=1000,
                         height=600,
                         **plot_kwargs,
                        )