from plotly.graph_objs import *
from plotly.offline import plot as offpy
import plotly.graph_objs as go
import numpy as np
import math
import os




def bar_chart_plot(x,y):
    data = [go.Bar(
        x=x,
        y=y
    )]

    fig = go.Figure(data=data)

    offpy(fig, filename="dd", auto_open=True, show_link=False)


def visualize_evolution(psi, phi, words, num_topics):
    topic_words = []
    for i in range(num_topics):
        words_indices = np.argsort(phi[i, :])[:10]
        topic_words.append([words[j] for j in words_indices])

    xs = np.linspace(0, 1, num=1000)
    data = []
    for i in range(len(psi)):
        ys = [math.pow(1 - x, psi[i][0] - 1) * math.pow(x, psi[i][1] - 1) / scipy.special.beta(psi[i][0], psi[i][1]) for
              x in xs]
        trace = go.Scatter(x=xs, y=ys, name=', '.join(topic_words[i]))
        data.append(trace)

    layout = go.Layout(
        xaxis=dict(
            autorange=True,
            showgrid=True,
            zeroline=True,
            showline=False,
            autotick=True,
            ticks='',
            showticklabels=False
        ),
        yaxis=dict(
            autorange=True,
            showgrid=True,
            zeroline=True,
            showline=False,
            autotick=True,
            ticks='',
            showticklabels=False
        ),
    )
    fig = go.Figure(data=data, layout=layout)
    offpy(fig, filename="visualize_evolution.html", auto_open=True, show_link=False)


def visualize_associations(X, Y, Z):
    trace = go.Heatmap(z=Z,
                       x=X,
                       y=Y)
    data = [trace]
    fig = go.Figure(data=data)
    offpy(fig, filename="associations.html", auto_open=True, show_link=False)


def show_image(image_url):
    import webbrowser
    new = 2
    html = "<html><head></head><body><img src="+image_url+" height='750' width='1500'></body></html>"
    dir_path = os.path.dirname(os.path.realpath(__file__))
    with open(dir_path+"/show_image.html", 'w') as file_writer:
        file_writer.write(html)

    dir_path = os.path.dirname(os.path.realpath(__file__))
    out_file_url = "file://"+dir_path+"/show_image.html"
    webbrowser.open(out_file_url, new=new)

if __name__ == "__main__":
    show_image(1)

