from flask import Flask, render_template, redirect, request, url_for, send_file, Response
from flask_bootstrap import Bootstrap
from markupsafe import Markup
import pandas as pd
import functions_db_query as dbq

import io, random
import matplotlib.pyplot as plt 
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

plt.switch_backend('Agg')

# create a flask application object
app = Flask(__name__)

# flask_wft requires encryption key
app.config['SECRET_KEY'] = 'chromosome22'
Bootstrap(app)

# search form in home page
@app.route('/', methods=['GET','POST'])
def index():
	if request.method == "POST":
		search_type = request.form.get("search_type")
		search_value = request.form.get("search_value")
		return redirect(url_for('search', search_type=search_type, search_value=search_value))
	return render_template('index.html')

# return result data base on search
@app.route('/search/<search_type>/<search_value>', methods = ['GET'])
def search(search_type, search_value):
	data = dbq.search_db(search_type, search_value)
	return render_template('result.html', data=data, search_type=search_type, search_value=search_value)

# select statistics and populations form
@app.route('/stats_pop_select', methods=['POST'])
def stats_pop():
	search_type = request.form["search_type"]
	search_value = request.form["search_value"]
	data = dbq.search_db(search_type, search_value) # data from db in list of tuples
	pop_list = request.form.getlist("population")
	stats_list = request.form.getlist("summarystats")
	data_df = dbq.to_df(data) # put data into dataframe for stats
	stats_df = dbq.calc_stats(data_df, stats_list, pop_list)
	stats_df.to_csv('stats.txt', sep=',', index=False, header=True)
	fst_df = dbq.calcFst(data_df, pop_list)

	model_plot=''
	if "shannon" in stats_list:
		plot_url = dbq.get_shannon_plot(stats_df)
		model_plot = Markup('<img src="data:image/png;base64,{}" width: 360px; height: 288px>'.format(plot_url))
	
	
	return render_template('stats_pop.html', data=data, search_type=search_type, search_value=search_value,
							tables=[stats_df.to_html(classes='data')], fsts=[fst_df.to_html(classes='data')],
							model_plot=model_plot)

# download txt file
@app.route('/download')
def download_file():
	path = "stats.txt"
	return send_file(path, as_attachment=True)


# @app.route('/images/<cropzonekey>')
# def images(cropzonekey):
#     return render_template("images.html", title=cropzonekey)

# @app.route('/fig/<cropzonekey>')
# def fig(cropzonekey):
#     fig = draw_polygons(cropzonekey)
#     img = StringIO()
#     fig.savefig(img)
#     img.seek(0)
#     return send_file(img, mimetype='image/png')

# @app.route('/plot.png')
# def plot_png():
#     fig = create_figure()
#     output = io.BytesIO()
#     FigureCanvas(fig).print_png(output)
#     return Response(output.getvalue(), mimetype='image/png')

# def create_figure():
# 	fig = Figure()
# 	axis = fig.add_subplot(1, 1, 1)
# 	xs = range(100)
# 	ys = [random.randint(1, 50) for x in xs]
# 	axis.plot(xs, ys)
# 	return fig

# start the web server
if __name__ == '__main__':
	app.run(debug=True, port=5001)